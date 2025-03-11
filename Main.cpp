/*

888b     d888        d8888 888b     d888  .d88888b. Y88b   d88P  .d8888b.
8888b   d8888       d88888 8888b   d8888 d88P" "Y88b Y88b d88P  d88P  Y88b
88888b.d88888      d88P888 88888b.d88888 888     888  Y88o88P          888
888Y88888P888     d88P 888 888Y88888P888 888     888   Y888P         .d88P
888 Y888P 888    d88P  888 888 Y888P 888 888     888   d888b     .od888P"
888  Y8P  888   d88P   888 888  Y8P  888 888     888  d88888b   d88P"
888   "   888  d8888888888 888   "   888 Y88b. .d88P d88P Y88b  888".....
888       888 d88P     888 888       888  "Y88888P" d88P   Y88b 888888888

Authors:
	- Sirio Brunalti (sirio.brunialti@kaust.edu.sa)
	- Enrico Garavaglia
	- Tiziano Faravelli (tiziano.faravelli@polimi.it)
	- Eliseo Ranzi


HOW TO CITE
	If a publication is produced using this software the following references
	must be reported
	- E. Ranzi, A. Frassoldati, S. Granata, T. Faravelli, Wide-Range Kinetic Modeling Study of the Pyrolysis, Partial Oxidation, and Combustion of Heavy n-Alkanes, Ind. Eng. Chem. Res. 44 (2005) 5170-5183
	- E. Ranzi, T. Faravelli, P. Gaffuri, A. Sogaro, Low-Temperature Combustion - Automatic-Generation of Primary Oxidation Reactions and Lumping Procedures, Combustion and Flame 102 (1995) 179-192.
	- S. Brunialti, X. Zhang, T. Faravelli, A. Frassoldati, S.M. Sarathy, Automatically generated detailed and lumped reaction mechanisms for low- and high-temperature oxidation of alkanes, Proceedings of the Combustion Institute 39 (2023) 335-344.

WARNING
	This version of MAMOX2 and its source code are intended for educational purposes.
	Use for commercial purposes is not permitted. For any commercial issue please contact
	Sirio Brunialti (sirio.brunialti@kaust.edu.sa) and Tiziano Faravelli (tiziano.faravelli@polimi.it).

LIMITED WARRANTY
	The Software and related documentation are provided “AS IS” and without
	any warranty of any kind and Seller EXPRESSLY DISCLAIMS ALL WARRANTIES,
	EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
	OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
*/



#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#define NOMINMAX		// need to be defined before inclusion of windows.h because otherwise there is a problem with the functions std::min and std::max
#include <windows.h>
#include <algorithm>
#include <thread>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <vector>
#include <direct.h>
#include <iomanip>

#include  "CKMechReader.h"
#include "utilities.h"
#include "readConfig.h"
#include "reaction.h"
#include "kinox.h"
#include "chemkinOut.h"
#include "thermoOut.h"
//#include "cantera/base/ctml.h"
#include "Simulation.h"
#include "LumpedReaction.h"

#define MAX_SIMULATION_TIME 6400

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION DECLARATIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

/**
* Generate the list of all possible radicals that can be formed from the input fuel
*
* @param HC Alkane the child radicals are created for
* @param Rvec pointer to the vector that is going to be filled with the child radicals
* @param cToRMap pointer to the vector that is going to contain the map correlating the
*		 index of the carbon in the structure with the corresponding radical in Rvec
* @return void
*/
void generateR(Molecola HC, std::vector<Molecola>* Rvec, std::vector<int>* cToRMap);

/**
* Generate the list of all possible ROO that can be formed from the input fuel
*
* @param HC Alkane the child ROO are created for
* @param ROOvec pointer to the vector that is going to be filled with the child ROO
* @param cToROOMap pointer to the vector that is going to contain the map correlating the
*		 index of the carbon in the structure with the corresponding ROO in ROOvec
* @return void
*/
void generateROO(Molecola HC, std::vector<Molecola>* ROOvec, std::vector<int>* cToROOMap);

/**
* Generate the list of all possible QOOH that can be formed from the input fuel
*
* @param HC Alkane the child QOOH are created for
* @param QOOHvec pointer to the vector that is going to be filled with the child QOOH
* @param cToQOOHMap pointer to the vector that is going to contain the map correlating the
*		 indexes of OOH and radical (in the order [OOHpos,Rpos] in the structure
*		 with the corresponding QOOH in QOOHvec
* @return void
*/
void generateQOOH(Molecola HC, std::vector<Molecola>* QOOHvec,
	std::vector<std::vector<int>>* cToQOOHMap);

/**
* Generate the list of all possible OOQOOH that can be formed from the input fuel
*
* @param HC Alkane the child OOQOOH are created for
* @param OOQOOHvec pointer to the vector that is going to be filled with the child OOQOOH
* @param cToOOQOOHMap pointer to the vector that is going to contain the map correlating the
*		 indexes of OOH and OO* (in the order [OOHpos,OO*pos] in the structure
*		 with the corresponding OOQOOH in OOQOOHvec
* @return void
*/
void generateOOQOOH(Molecola HC, std::vector<Molecola>* OOQOOHvec,
	std::vector<std::vector<int>>* cToOOQOOHMap);

/**
* Generate all the initiation reactions for the fuel
*
* @param HC the fuel the reaction are generate for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> initiationReactions(Molecola HC, Kinox* k);

/**
* Generate all the H-abstraction reactions for the fuel and the specific abstracting molecule
*
* @param HC the fuel the reaction are generate for
* @param absR the molecule abstracting the hydrogen
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> hAbstractionReactions(Molecola HC, Molecola absR, Kinox* k);

/**
* Generate all the O2 + R -> ROO reactions for the list of radicals
*
* @param Rs the list of radicals to generate the reactions for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> O2AdditionToRReactions(std::vector<Molecola> Rs, Kinox* k);

/**
* Generate all the ROO -> R + O2 reactions for the list of ROO
*
* @param ROOs the list of ROO to generate the reactions for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> O2EliminationFromROOReactions(std::vector<Molecola> ROOs, Kinox* k);

/**
* Generate all the O2 + QOOH -> OOQOOH reactions for the list QOOHs
*
* @param QOOHs the list of QOOHs to generate the reactions for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> O2AdditionToQOOHReactions(std::vector<Molecola> QOOHs, Kinox* k);

/**
* Generate all the OOQOOH -> QOOH + O2 reactions for the list OOQOOHs
*
* @param OOQOOHs the list of OOQOOHs to generate the reactions for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> O2EliminationFromOOQOOHReactions(std::vector<Molecola> OOQOOHs, Kinox* k);

/**
* Generate all the isomerization reactions for all the radicals in the vector Rs
*
* @param Rs the vector of radicals the reactions are going to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> RIsomerizationReaction(std::vector<Molecola> Rs, Kinox* k);

/**
* Generate all the isomerization reactions for all the ROO in the vector ROOs
*
* @param ROOs the vector of ROO the reactions are going to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> ROOIsomerizationReactions(std::vector<Molecola> ROOs, Kinox* k);

/**
* Generate all the isomerization reactions for all the QOOH in the vector QOOHs
*
* @param QOOHs the vector of QOOH the reactions are going to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> QOOHIsomerizationReactions(std::vector<Molecola> QOOHs, Kinox* k);

/**
* Generate all the isomerization reactions for all the OOQOOH in the vector OOQOOHs
*
* @param OOQOOHs the vector of OOQOOH the reactions are going to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> OOQOOHIsomerizationReactions(std::vector<Molecola> OOQOOHs, Kinox* k);

/**
* Generate all the isomerization reactions for all the P(OOH)2 in the vector POOH2s
*
* @param POOH2s the vector of P(OOH)2 the reactions are going to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all the reactions
*/
std::vector<Reaction> POOH2IsomerizationReactions(std::vector<Molecola> POOH2s, Kinox* k);

/**
* Generate all the reactions of concerted elimination of HO2 from alkylperoxy radicals
* (ROO -> OLE + HO2) for the provided ROO
*
* @param ROOs vector of all the ROO the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> olefinsFromROOReactions(std::vector<Molecola> ROOs, Kinox* k);

/**
* Generate all the reactions of beta decomposition for the provided radicals
*
* @param Rs vector of all the radicals the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> RBetaDecompositioReactions(std::vector<Molecola> Rs, Kinox* k);

/**
* Generate all the reactions of radical decomposition to olefins from H abstraction
* from O2 (R + O2 -> OLE + HO2)
*
* @param Rs vector of all the radicals the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> olefinsFromRPlusO2Reactions(std::vector<Molecola> Rs, Kinox* k);

/**
* Generate all the reactions of beta-QOOH decomposition for the provided QOOHs, if a
* QOOH in the list is not a beta-QOOH it gets skipped.
*
* @param QOOHs vector of all the QOOH the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> betaQOOHDecompositionReaction(std::vector<Molecola> QOOHs, Kinox* k);

/**
* Generate all the reactions of gamma-QOOH decomposition for the provided QOOHs, if a
* QOOH in the list is not a gamma-QOOH it gets skipped.
*
* @param QOOHs vector of all the QOOH the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> gammaQOOHDecompositionReaction(std::vector<Molecola> QOOHs, Kinox* k);

/**
* Generate all the reactions of delta-QOOH decomposition for the provided QOOHs, if a
* QOOH in the list is not a delta-QOOH it gets skipped.
*
* @param QOOHs vector of all the QOOH the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> deltaQOOHDecompositionReaction(std::vector<Molecola> QOOHs, Kinox* k);

/**
* Generate all the reactions of QOOH decomposition to cyclic ethers (QOOH -> cEth + OH)
*
* @param QOOHs vector of all the QOOH the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> QOOHToCEthReactions(std::vector<Molecola> QOOHs, Kinox* k);

/**
* Generate all the reactions of P(OOH)2 decomposition to cyclic ethers (P(OOH)2 -> cEthOOH + OH)
*
* @param POOH2s vector of all the P(OOH)2 the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> POOH2ToCEthOOHReactions(std::vector<Molecola> POOH2s, Kinox* k);

/**
* Generate all the reactions of beta-P(OOH)2 decomposition for the provided P(OOH)2s, if a
* P(OOH)2 in the list is not a beta-P(OOH)2 it gets skipped.
*
* @param POOH2s vector of all the P(OOH2) the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> betaPOOH2DecompositionReaction(std::vector<Molecola> POOH2s, Kinox* k);

/**
* Generate all the reactions of gamma-P(OOH)2 decomposition for the provided P(OOH)2s, if a
* P(OOH)2 in the list is not a gamma-P(OOH)2 it gets skipped.
*
* @param POOH2s vector of all the P(OOH2) the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> gammaPOOH2DecompositionReaction(std::vector<Molecola> POOH2s, Kinox* k);

/**
* Generate all the reactions of delta-P(OOH)2 decomposition for the provided P(OOH)2s, if a
* P(OOH)2 in the list is not a delta-P(OOH)2 it gets skipped.
*
* @param POOH2s vector of all the P(OOH2) the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> deltaPOOH2DecompositionReaction(std::vector<Molecola> POOH2s, Kinox* k);

/**
* Generate all the cyclic ether-OOH decomposition reactions for the provided
* cyclic ethers-OOH
*
* @param cEthOOHs vector of all the cyclic ethers-OOH the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> cyclicEtherOOHDecompositionReactions(std::vector<Molecola> cEthOOHs,
	Kinox* k);

/**
* Generate all the reactions for the decomposition of OOQOOH species to OLE-OOH
*
* @param OOQOOHs vector of all the OOQOOH the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> OOQOOHToOLEOOHReactions(std::vector<Molecola>OOQOOHs, Kinox* k);

/**
* Generate all the reactions for the decomposition of OLE-OOH species
*
* @param OLEOOHs vector of all the OLEOOH the reactions have to be generated for
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> OLEOOHDecompositionReactions(std::vector<Molecola> OLEOOHs, Kinox* k);

/**
* Generate all the reactions for the formation of ketohydroperoxides
*
* @param OOQOOHs vector of all the OOQOOH to decompose to ketohydroperoxides
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> KHPFormationReactions(std::vector<Molecola> OOQOOHs, Kinox* k);

/**
* Generate all the ketohydroperoxides decomposition reactions
*
* @param KHPs vector of all the ketohydroperoxides to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> KHPDecompositionReactions(std::vector<Molecola> KHPs, Kinox* k);

/**
* Generate all the cyclic ethers decomposition reactions
*
* @param cEths vector of all the cyclic ethers to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> cyclicEthersDecompositionReactions(std::vector<Molecola> cEths,
	Kinox* k);

/**
* Generate all the reactions of olefins' conversion to allylic radicals
*
* @param OLEs vector of all the olefins to convert
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> allylicRadicalsFormationReactions(std::vector<Molecola> OLEs, Kinox* k);

/**
* Generate all the reactions of allylic radicals' conversion to alkenyl RO
*
* @param AllRs vector of all the allylic radicals' to convert
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> alkenylROFormationReactions(std::vector<Molecola> AllRs, Kinox* k);

/**
* Generate all the reactions of alkenyl RO decomposition
*
* @param AlkROs vector of all the alkenyl RO to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> alkenylRODecompositionReactions(std::vector<Molecola> AlkROs, Kinox* k);

/**
* Generate all the aldehydes decomposition reactions
*
* @param ALD vector of all the aldehydes to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> aldehydesDecompositionReactions(std::vector<Molecola> ALDs, Kinox* k);

/**
* Generate all the aldehyde olefins decomposition reactions
*
* @param ALDOLE vector of all the aldehyde olefins to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> aldehydeOlefinsDecompositionReactions(std::vector<Molecola> ALDOLEs, Kinox* k);

/**
* Generate all the ketones decomposition reactions
*
* @param KETOs vector of all the ketones to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> ketonesDecompositionReactions(std::vector<Molecola> KETOs, Kinox* k);

/**
* Generate all the ketones olefins decomposition reactions
*
* @param KETOOLEs vector of all the ketones olefins to decompose
* @param k pointer to the kinox object that provides the rate rules
* @return a vector containing all reactions
*/
std::vector<Reaction> ketonesOlefinsDecompositionReactions(std::vector<Molecola> KETOOLEs, Kinox* k);


/**
* Decompose the cyclic ether radical
* @param m the cyclic ether radical to decompose
* @return a vector containig all the decomposition products
*/
std::vector<Molecola> decomposeCEthR(Molecola m);

/**
* Decompose the linear ether RO
* @param m the linear ether RO to decompose
* @return a vector containig all the decomposition products
*/
std::vector<Molecola> decomposeLinEthRO(Molecola mol);

/**
* Decompose the RO species
* @param m the RO species to decompose
* @return a vector containig all the decomposition products of the RO in vec
*/
std::vector<Molecola> decomposeRO(Molecola m2);

/**
* Decompose the RO species and the eventually produced RO products
* @param m the RO species to decompose
* @return a vector containig all the decomposition products of the RO in vec
*/
std::vector<Molecola> fullyDecomposeRO(Molecola vec);

/**
* Decompose all the RO species in the vector and the eventually produced RO products
* @param m vector of species containing the RO to decompose
* @return a vector containig all the non RO species of vec and the decomposition products
*		  of the RO in vec
*/
std::vector<Molecola> fullyDecomposeRO(std::vector<Molecola> vec);

/**
* Get all the products that belong to a specific family of species from all the provided
* reactions. Duplicate species are skipped.
*
* @param reactions vector of reactions the products have to be taken from
* @param kind the kind of species to be taken
* @return a vector containing all the desired species
*/
std::vector<Molecola> getProducts(std::vector<Reaction> reactions, species kind);

/**
* Seek in all the products of the provided reactions for species belonging to the specified
* class of species and add to the molecule vector all the missing species.
*
* @param molVec pointer to the vector containing the list of species to be filled
* @param reactions the vector of reactions the products have to be taken from
* @param kind the kind of species to be taked
* @return integer equal to the number of new elements added to molVec
*/
int getAdditionalProducts(std::vector<Molecola>* molVec, std::vector<Reaction> reactions,
	species kind);

/**
* Convert the enumarator Carbonio to the corresponding number:	primary		-> 1
*																secondary	-> 2
*																tertiary	-> 3
*																quaternary	-> 4
* @param c the Carbonio enumerator to convert
* @return the integer
*/
int carbonioToInt(Carbonio c);

/**
* Process the vector of reactions and deals with duplicate reactions. If two reactions
* with same products and same reactants are found they are deemed to be duplicates. If
* they have also the same rate rule used (and same rate costants) they are merged into
* a single reaction and the reaction comment is modified to take into account the
* multiplication factor. If the reactions don't have the same rate rule they are flagged
* as duplicates.
*
* @param reacVec pointer to the vector to process
* @return int return the number of merged reactions
*/
int processDuplicateReactions(std::vector<Reaction>* reacVec);

/**
* Add all the new species in the reaction vector (both reactants and products) to the
* vector of molecules
*
* @param molVec pointer to the vector of molecules
* @param reacVec pointer to the vector of reactions
* @return void
*/
void addNewSpecies(std::vector<Molecola>* molVec, std::vector<Reaction>* reacVec);

/**
* Print the species in the provided file with the proper formatting
*
* @param outfile pointer to the ofstream in which the species have to be printed
* @param mols the vector of species to print
* @param label the string with the label to print as comment at the beginning
*		 of the list
* @param chemOut pointer to the ChemkinOut object providing the names
* @return void
*/
void printSpeciesInFile(std::ofstream* outfile, std::vector<Molecola> mols,
	std::string label, ChemkinOut* chemOut);

/**
* Print the species names in the provided file with the proper formatting
*
* @param outfile pointer to the ofstream in which the species have to be printed
* @param mols the vector of species' names to print
* @param label the string with the label to print as comment at the beginning
*		 of the list
* @return void
*/
void printSpeciesInFile(std::ofstream* outfile, std::vector<std::string> mols,
	std::string label);

/**
* Print the reaction in Chemkin format
*
* @param outfile pointer to the file the reaction is going to printed in
* @param raec the reaction to print
* @param chemOut pointer to the ChemkinOut object that names the molecules
* @return void
*/
void printReaction(std::ofstream* outfile, Reaction reac, ChemkinOut* chemOut);

/**
* * Print the reactions in Chemkin format
*
* @param outfile pointer to the file the reaction is going to printed in
* @param raecs the vector of reactions to print
* @param chemOut pointer to the ChemkinOut object that names the molecules
* @return void
*/
void printReactions(std::ofstream* outfile, std::vector<Reaction> reacs,
	ChemkinOut* chemOut);

/**
* Find if there is a reaction pathway consuming the species both in the base
* mechanism and in the list of generated reactions
*
* @param spec the molecule to find the consumption pathways for
* @param baseMechReacs the vector of baseReactions of the base mechanism
* @param totReactions the vector of reactions generated
* @param chemOut the pointer to the ChemkinOut object that provides the naming
* @return bool true if there is a consumption pathway otherwise false
*/
bool isThereDecompositionPath(Molecola spec,
	std::vector<baseReaction>* baseMechReacs,
	std::vector<Reaction>* totReactions, ChemkinOut* chemOut);

/**
* Return a string with the name of the kind of species
*
* @param kind the kind of species
* @return the string with the name of the kind of species
*/
std::string speciesToText(species kind);


/**
* Return true if the reaction is included in the list of base reactions.
* A reaction is deemed to be included if there is a reaction with the same
* products and reactants (it is checked also against reverse reactions).
* Kinetic parameters are not taken into account.
*
* @param reac the reaction to check if is included
* @param reacList the list of base reactions in which check if reac is present
* @param chemOut pointer to the ChemkinOut object used for naming the species
* @return bool true if the reaction is present otherwise false
*/
bool isIncluded(Reaction reac, std::vector<baseReaction>* reacList,
	ChemkinOut* chemOut);

/**
* Return a string with expanded environment variables (Windows)
*
* @param inString the string to expand
* @return the string with expanded environment variables
*/
std::string expandEnvVar(std::string inString);

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~ MAIN ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

int main(int argc, char* argv[])
{
	// Set window size and position
	HWND console = GetConsoleWindow();
	RECT rec;
	GetWindowRect(console, &rec);
	MoveWindow(console, rec.left, rec.top, 1500, 800, TRUE);

	
	// ###################################################################################
	// ###################### OPEN AND READ CONFIG FILE ##################################
	// ###################################################################################
	// --------------- open file -------------------------------------
	if (argc == 1)
		UTL::fatalError("No input argument has been provided!");
	if (argc > 2)
		UTL::fatalError("Only one arguent mus be provided!");

	std::string configFilePath = argv[1];
	
	std::ifstream configFile(configFilePath);
	if (configFile.is_open() == false)
		UTL::fatalError("Unable to open config File!");
	// --------------- read file --------------------------------------
	// define variables that need to be read
	std::vector<std::string> fuelsNames;  // vector with the names of the mol files
	std::vector<double> Temps;			  // vector of temperatures [K] for the simulations 
	std::vector<double> Press;			  // vector of pressures [atm] for the simulations
	std::string baseMechanismPath;		  // path to base mechanisms' kinetic file
	std::string baseThermoPath;			  // path to base mechanisms' thermo file
	std::string baseNamingPath;			  // path to base mechanisms' glossary
	std::string outFolderPath;			  // path to output folder path
	std::string molFolderPath;			  // path to folder containing the molecules files
	std::string thermoGroupPath;		  // path to the group additivity values csv
	std::string thermoHBIPath;			  // path to the HBI values csv
	std::string rateRulesPath;			  // path to the rate rules csv
	// variables with default value
	int numCores = 2;
	bool reversibeDetailed = false;
	bool generateLumped = false;
	bool useBatch = true;
	double eqRatio = 1.0;
	double tau = 2.0;
	bool lumpWithCoreMech = false;
	bool equilibriumLumping = false;

	if (rdcfg::isKeywordPresent(&configFile, "FUELS"))
		fuelsNames = rdcfg::readStringVec(&configFile, "FUELS");
	else
		UTL::fatalError("In reading config file: No FUELS provided!");

	if (fuelsNames.size() == 0)
		UTL::fatalError("In reading config file: No FUELS provided!");

	if (rdcfg::isKeywordPresent(&configFile, "OUT_FOLDER"))
		outFolderPath = expandEnvVar(rdcfg::readString(&configFile, "OUT_FOLDER"));
	else
		UTL::fatalError("In reading config file: No output folder path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "MOL_FOLDER"))
		molFolderPath = expandEnvVar(rdcfg::readString(&configFile, "MOL_FOLDER"));
	else
		UTL::fatalError("In reading config file: No output folder path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "BASE_MECH"))
		baseMechanismPath = expandEnvVar(rdcfg::readString(&configFile, "BASE_MECH"));
	else
		UTL::fatalError("In reading config file: No base mechanism path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "BASE_THERMO"))
		baseThermoPath = expandEnvVar(rdcfg::readString(&configFile, "BASE_THERMO"));
	else
		UTL::fatalError("In reading config file: No base mechanism's thermodynamic data file path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "BASE_NAMES"))
		baseNamingPath = expandEnvVar(rdcfg::readString(&configFile, "BASE_NAMES"));
	else
		UTL::fatalError("In reading config file: No base mechanism's names file path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "THERMO_GROUP_FILE"))
		thermoGroupPath = expandEnvVar(rdcfg::readString(&configFile, "THERMO_GROUP_FILE"));
	else
		UTL::fatalError("In reading config file: No thermodynamic groups contribution file path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "THERMO_HBI_FILE"))
		thermoHBIPath = expandEnvVar(rdcfg::readString(&configFile, "THERMO_HBI_FILE"));
	else
		UTL::fatalError("In reading config file: No thermodynamic HBI file path defined!");

	if (rdcfg::isKeywordPresent(&configFile, "RATE_RULES_FILE"))
		rateRulesPath = expandEnvVar(rdcfg::readString(&configFile, "RATE_RULES_FILE"));
	else
		UTL::fatalError("In reading config file: No rate rules file path defined!");


	if (rdcfg::isKeywordPresent(&configFile, "REVERSIBLE_DETAILED_MECH"))
	{
		std::string value = rdcfg::readString(&configFile, "REVERSIBLE_DETAILED_MECH");
		if (value == "FALSE")
			reversibeDetailed = false;
		else if (value == "TRUE")
			reversibeDetailed = true;
		else
			UTL::warning("In reading config file: Definition for keyword REVERSIBLE_DETAILED_MECH must be either TRUE or FALSE. Default value FALSE has been set.");
	}

	if (rdcfg::isKeywordPresent(&configFile, "GENERATE_LUMPED"))
	{
		std::string value = rdcfg::readString(&configFile, "GENERATE_LUMPED");
		if (value == "FALSE")
			generateLumped = false;
		else if (value == "TRUE")
			generateLumped = true;
		else
			UTL::warning("In reading config file: Definition for keyword GENERATE_LUMPED must be either TRUE or FALSE. Default value FALSE has been set.");
	}

	if (generateLumped)
	{
		if (rdcfg::isKeywordPresent(&configFile, "NUM_CORES"))
			numCores = rdcfg::readInt(&configFile, "NUM_CORES");

		if (rdcfg::isKeywordPresent(&configFile, "EQ_RATIO"))
			eqRatio = rdcfg::readDouble(&configFile, "EQ_RATIO");

		if (rdcfg::isKeywordPresent(&configFile, "TAU"))
			tau = rdcfg::readDouble(&configFile, "TAU");

		if (rdcfg::isKeywordPresent(&configFile, "LUMPING_T"))
			Temps = rdcfg::readDoubleVec(&configFile, "LUMPING_T");
		else
			UTL::fatalError("In reading config file: If GENERATE_LUMPED is TRUE, LUMPING_T must be provided!");
		if (Temps.size() == 0)
			UTL::fatalError("In reading config file: If GENERATE_LUMPED is TRUE, LUMPING_T must be provided!");


		if (rdcfg::isKeywordPresent(&configFile, "LUMPING_P"))
			Press = rdcfg::readDoubleVec(&configFile, "LUMPING_P");
		else
			UTL::fatalError("In reading config file: If GENERATE_LUMPED is TRUE, LUMPING_P must be provided!");
		if (Press.size() == 0)
			UTL::fatalError("In reading config file: If GENERATE_LUMPED is TRUE, LUMPING_P must be provided!");

		if (rdcfg::isKeywordPresent(&configFile, "SIM_TYPE"))
		{
			std::string value = rdcfg::readString(&configFile, "SIM_TYPE");
			if (value == "BATCH")
			{
				useBatch = true;
				lumpWithCoreMech = true;
			}
			else if (value == "CSTR")
			{
				useBatch = false;
				lumpWithCoreMech = false;
			}
			else if (value == "EQUILIBRIUM")
			{
				equilibriumLumping = true;
			}
			else
				UTL::warning("In reading config file: SIM_TYPE must be either BATCH, CSTR or EQUILIBRIUM. Default value BATCH has been set.");
		}

		if (rdcfg::isKeywordPresent(&configFile, "LUMP_WITH_CORE_MECH"))
		{
			std::string value = rdcfg::readString(&configFile, "LUMP_WITH_CORE_MECH");
			if (value == "FALSE")
			{
				if (useBatch)
				{
					UTL::warning("In reading config file: Keyword LUMP_WITH_CORE_MECH cannot be FALSE when SIM_TYPE is BATCH. TRUE value has been set instead.");
					lumpWithCoreMech = true;
				}
				else
				{
					lumpWithCoreMech = false;
				}
			}
			else if (value == "TRUE")
				lumpWithCoreMech = true;
			else
				UTL::warning("In reading config file: Definition for keyword LUMP_WITH_CORE_MECH must be either TRUE or FALSE. Default value has been set.");

		}
	}


	// print read values to terminal
	UTL::printTitle("VALUES READ FROM CONFIG FILE");
	std::cout << "   FUELS = " << fuelsNames[0] << std::endl;
	for (int i = 1; i < fuelsNames.size(); i++)
		std::cout << "           " << fuelsNames[i] << std::endl;
	std::cout << std::endl;
	std::cout << "   OUT_FOLDER   = " << outFolderPath << std::endl;
	std::cout << std::endl;
	std::cout << "   BASE_MECH   = " << baseMechanismPath << std::endl;
	std::cout << "   BASE_THERMO = " << baseThermoPath << std::endl;
	std::cout << "   BASE_NAMES  = " << baseNamingPath << std::endl;
	std::cout << std::endl;
	std::cout << "   THERMO_GROUP_FILE = " << thermoGroupPath << std::endl;
	std::cout << "   THERMO_HBI_FILE   = " << thermoHBIPath << std::endl;
	std::cout << std::endl;
	std::cout << "   GENRATE_LUMPED = ";
	if (generateLumped)
		std::cout << "TRUE";
	else
		std::cout << "FALSE";
	std::cout << std::endl;
	if (generateLumped)
	{
		std::cout << "   SIM_TYPE  = ";
		if (equilibriumLumping)
		{
			std::cout << "EQUILIBRIUM";
		}
		else
		{
			if (useBatch)
				std::cout << "BATCH";
			else
				std::cout << "CSTR";
		}
		std::cout << std::endl;

		std::cout << "   LUMPING_T = " << Temps[0] << std::endl;
		for (int i = 1; i < Temps.size(); i++)
			std::cout << "               " << Temps[i] << std::endl;
		std::cout << "   LUMPING_P = " << Press[0] << std::endl;
		for (int i = 1; i < Press.size(); i++)
			std::cout << "               " << Press[i] << std::endl;
		std::cout << "   NUM_CORES = " << numCores << std::endl;
	}

	// ############# CHECK IF FILES ARE AVAILABLE ########################################
	if (!UTL::checkFileExistence(baseMechanismPath))
		UTL::fatalError("Base mechanism file is not accessible.");
	if (!UTL::checkFileExistence(baseThermoPath))
		UTL::fatalError("Base mechanism thermodynamic file is not accessible.");
	if (!UTL::checkFileExistence(baseNamingPath))
		UTL::fatalError("Base mechanism naming file is not accessible.");
	if (!UTL::checkFileExistence(thermoGroupPath))
		UTL::fatalError("Thermodynamic group contribution file is not accessible.");
	if (!UTL::checkFileExistence(thermoHBIPath))
		UTL::fatalError("Thermodynamic HBI file is not accessible.");
	if (!UTL::checkFileExistence(rateRulesPath))
		UTL::fatalError("Rate rules file is not accessible.");
	if (!UTL::dirExists(molFolderPath))
		UTL::fatalError("Molecule folder " + molFolderPath + " is not accessible.");

	std::vector<std::string> fuelsPaths(fuelsNames.size());
	for (int i = 0; i < fuelsNames.size(); i++)
	{
		std::string pathToMol = molFolderPath + "\\" + fuelsNames[i] + ".hyd";
		if (!UTL::checkFileExistence(pathToMol))
			UTL::fatalError("Molecule file " + pathToMol + " is not accessible.");
		fuelsPaths[i] = pathToMol;
	}

	// ############# GENERATE REQUIRED FOLDERS ###########################################
	UTL::createDirectory(outFolderPath);

	std::string subMechsFolder = outFolderPath + "\\SUBMECHS";
	UTL::createDirectory(subMechsFolder);

	// ###################################################################################
	// ###################### OPEN AND PARSE FILES #######################################
	// ###################################################################################
	std::cout << std::endl;
	UTL::printTitle("READING FILES");

	// read base mechanism
	CKMechReader baseMech(baseMechanismPath, baseThermoPath, baseNamingPath);
	ChemkinOut chemOut(&baseMech);
	std::cout << "Base mechanism read from : " << baseMechanismPath << std::endl;
	std::cout << "                           " << baseThermoPath << std::endl;
	std::cout << "                           " << baseNamingPath << std::endl << std::endl;

	// read rate rules
	Kinox k(rateRulesPath);
	std::cout << "Rate rules read from : " << rateRulesPath << std::endl << std::endl;

	// read fuels files
	std::cout << "Reading molecules: ";
	std::vector<Molecola> fuels(fuelsNames.size());
	for (int i = 0; i < fuels.size(); i++)
	{
		fuels[i] = Molecola(fuelsPaths[i]);
		std::cout << fuels[i] << "   from " << fuelsPaths[i] << std::endl;
		std::cout << "                   ";
	}
	std::cout << std::endl << std::endl;

	// read thermo files
	ThermoOut thermOut(thermoGroupPath, thermoHBIPath, baseThermoPath, &chemOut);

	// ###################################################################################
	// ###################### INITIALIZE VARIABLES #######################################
	// ###################################################################################
	std::vector<Reaction> totalReactionsList;
	std::vector<Molecola> totalSpeciesList;
	std::vector<LumpedReaction> totalLumpedReactionsList;

	// CREATE SPECIAL MOLECULES
	Molecola N2(1);
	Molecola O2(2);
	Molecola O(3);
	Molecola OH(4);
	Molecola HO2(5);
	Molecola H2O(6);
	Molecola H2O2(7);
	Molecola H(8);
	Molecola H2(9);
	Molecola CO(10);
	Molecola CH3;
	CH3.makeCH3();
	Molecola CH4;
	CH4.makeCH4();
	Molecola C2H5;
	C2H5.makeC2H5();
	Molecola C2H6;
	C2H6.makeC2H6();

	std::vector<Molecola> specialMols = { N2, O2, O, OH, HO2, H2O, H2O2, H, H2,
		CO, CH3, CH4, C2H5, C2H6 };


	// ###################################################################################
	// ####################### GENERATE SUBMECHANISMS ####################################
	// ###################################################################################

	// start time measurement
	std::chrono::steady_clock::time_point startTime = std::chrono::steady_clock::now();
	std::ofstream timeFile;
	timeFile.open(outFolderPath+"\\timeStamps.txt");
	std::chrono::steady_clock::time_point currentTime = std::chrono::steady_clock::now();

	for (int fuelIndex = 0; fuelIndex < fuels.size(); fuelIndex++)
	{
		std::cout << std::endl << std::endl;
		UTL::printTitle(fuelsNames[fuelIndex] + " submechanism");
		std::cout << std::endl;

		// ### GENERATE FOLDERS ###
		std::cout << "Creating output folder ..." << std::endl << std::endl;
		std::string subMechFold = subMechsFolder + "\\" + fuelsNames[fuelIndex];
		UTL::createDirectory(subMechFold);


		Molecola HC = fuels[fuelIndex];

		// ### GENERATE RADICALS ####
		//std::cout << "Generating radicals ..." << std::endl;
		std::vector<Molecola> Rs;
		std::vector<int> cToRMap;
		generateR(HC, &Rs, &cToRMap);
		//std::cout << "   RADICALS:" << std::endl;
		//for (int i = 0; i < Rs.size(); i++)
		//	std::cout << "   " << i + 1 << ")   " << Rs[i] << std::endl;
		//std::cout << std::endl;

		// ### GENERATE ROO ####
		//std::cout << "Generating ROO ..." << std::endl;
		std::vector<Molecola> ROOs;
		std::vector<int> cToROOMap;
		generateROO(HC, &ROOs, &cToROOMap);
		//std::cout << "   ROO:" << std::endl;
		//for (int i = 0; i < ROOs.size(); i++)
		//	std::cout << "   " << i + 1 << ")   " << ROOs[i] << std::endl;
		//std::cout << std::endl;

		// ### GENERATE QOOH ####
		//std::cout << "Generating QOOH ..." << std::endl;
		std::vector<Molecola> QOOHs;
		std::vector<std::vector<int>> cToQOOHMap;
		generateQOOH(HC, &QOOHs, &cToQOOHMap);
		//std::cout << "   QOOH:" << std::endl;
		//for (int i = 0; i < QOOHs.size(); i++)
		//	std::cout << "   " << i + 1 << ")   " << QOOHs[i] << std::endl;
		//std::cout << std::endl;

		// ### GENERATE QOOH ####
		//std::cout << "Generating OOQOOH ..." << std::endl;
		std::vector<Molecola> OOQOOHs;
		std::vector<std::vector<int>> cToOOQOOHMap;
		generateOOQOOH(HC, &OOQOOHs, &cToOOQOOHMap);
		//std::cout << "   OOQOOH:" << std::endl;
		//for (int i = 0; i < OOQOOHs.size(); i++)
		//	std::cout << "   " << i + 1 << ")   " << OOQOOHs[i] << std::endl;
		//std::cout << std::endl;

		// ################ GENERATE REACTIONS ##########################################
		UTL::printEmbeddedString('*', "Generate reactions");

		// ### INITIATION REACTIONS (C-C breakage) ###
		std::vector<Reaction> initReac = initiationReactions(HC, &k);
		//for (auto& rea : initReac)
			//std::cout << rea << std::endl;
		std::cout << "- Initiation reactions added.       ("
			<< initReac.size() << " reactions)" << std::endl;

		// ### H-ABSTRACTION ###
		std::vector<Reaction> hAbsReac;
		std::vector<Reaction> hAbsReacByO2 = hAbstractionReactions(HC, O2, &k);
		std::vector<Reaction> hAbsReacByOH = hAbstractionReactions(HC, OH, &k);
		std::vector<Reaction> hAbsReacByH = hAbstractionReactions(HC, H, &k);
		std::vector<Reaction> hAbsReacByO = hAbstractionReactions(HC, O, &k);
		std::vector<Reaction> hAbsReacByHO2 = hAbstractionReactions(HC, HO2, &k);
		std::vector<Reaction> hAbsReacByCH3 = hAbstractionReactions(HC, CH3, &k);
		std::vector<Reaction> hAbsReacByC2H5 = hAbstractionReactions(HC, C2H5, &k);

		UTL::concatenate(&hAbsReac, &hAbsReacByO2);
		UTL::concatenate(&hAbsReac, &hAbsReacByOH);
		UTL::concatenate(&hAbsReac, &hAbsReacByH);
		UTL::concatenate(&hAbsReac, &hAbsReacByO);
		UTL::concatenate(&hAbsReac, &hAbsReacByHO2);
		UTL::concatenate(&hAbsReac, &hAbsReacByCH3);
		UTL::concatenate(&hAbsReac, &hAbsReacByC2H5);
		//for (auto& rea : hAbsReac)
		//	std::cout << rea << std::endl;
		std::cout << "- H abstraction reactions added.       ("
			<< hAbsReac.size() << " reactions)" << std::endl;

		// ### O2 ADDITION TO R (R+O2->ROO) ###
		std::vector<Reaction> O2addToRReac = O2AdditionToRReactions(Rs, &k);
		//for (auto& rea : O2addToRReac)
		//	std::cout << rea << std::endl;
		std::cout << "- O2 addition to R reactions added.       ("
			<< O2addToRReac.size() << " reactions)" << std::endl;

		// ### O2 ELIMINATION FROM ROO ###
		std::vector<Reaction> O2ElimFromROOReac = O2EliminationFromROOReactions(ROOs, &k);
		//for (auto& rea : O2ElimFromROOReac)
		//	std::cout << rea << std::endl;
		std::cout << "- O2 elimination from ROO reactions added.       ("
			<< O2ElimFromROOReac.size() << " reactions)" << std::endl;

		// ### O2 ADDITION TO QOOH (QOOH+O2->OOQOOH) ###
		std::vector<Reaction> O2addToQOOHReac = O2AdditionToQOOHReactions(QOOHs, &k);
		//for (auto& rea : O2addToQOOHReac)
		//	std::cout << rea << std::endl;
		std::cout << "- O2 addition to QOOH reactions added.       ("
			<< O2addToQOOHReac.size() << " reactions)" << std::endl;

		// ### O2 ELIMINATION FROM OOQOOH ###
		std::vector<Reaction> O2ElimFromOOQOOHReac = O2EliminationFromOOQOOHReactions(OOQOOHs, &k);
		//for (auto& rea : O2ElimFromOOQOOHReac)
		//	std::cout << rea << std::endl;
		std::cout << "- O2 elimination from OOQOOH reactions added.       ("
			<< O2ElimFromOOQOOHReac.size() << " reactions)" << std::endl;

		// ### RADICALS ISOMERIZATION ###
		std::vector<Reaction> RIsomReac = RIsomerizationReaction(Rs, &k);
		if (reversibeDetailed)
		{
			std::vector<Reaction> newVec;
			for (int i = 0; i < RIsomReac.size(); i++)
			{
				bool hasDuplicate = false;
				for (int j = 0; j < i; j++)
				{
					std::vector<Molecola*>reac1=RIsomReac[i].reactantList();
					std::vector<Molecola*>reac2=RIsomReac[j].reactantList();
					std::vector<Molecola*>prod1=RIsomReac[i].productList();
					std::vector<Molecola*>prod2=RIsomReac[j].productList();
					if (*(reac1[0]) == *(prod2[0]) && *(reac2[0]) == *(prod1[0]))
					{
						hasDuplicate = true;
						break;
					}

				}
				if (hasDuplicate == false)
					newVec.push_back(RIsomReac[i]);
			}
			RIsomReac = newVec;
		}
		//for (auto& rea : RIsomReac)
		//	std::cout << rea << std::endl;
		std::cout << "- R isomerization reactions added.       ("
			<< RIsomReac.size() << " reactions)" << std::endl;

		// ### ISOMERIZATION ROO (ROO -> QOOH) ###
		std::vector<Reaction> ROOIsomReac = ROOIsomerizationReactions(ROOs, &k);
		//for (auto& rea : ROOIsomReac)
		//	std::cout << rea << std::endl;
		std::cout << "- ROO isomerization reactions added.       ("
			<< ROOIsomReac.size() << " reactions)" << std::endl;

		// ### ISOMERIZATION QOOH (QOOH -> ROO) ###
		std::vector<Reaction> QOOHIsomReac = QOOHIsomerizationReactions(QOOHs, &k);
		//for (auto& rea : QOOHISomReac)
		//	std::cout << rea << std::endl;
		std::cout << "- QOOH isomerization reactions added.       ("
			<< QOOHIsomReac.size() << " reactions)" << std::endl;

		// ### ISOMERIZATION OOQOOH (OOQOOH -> P(OOH)2) ###
		std::vector<Reaction> OOQOOHIsomReac = OOQOOHIsomerizationReactions(OOQOOHs, &k);
		std::vector<Molecola> POOH2s = getProducts(OOQOOHIsomReac, POOH2_);
		//for (auto& prod : POOH2s)
		//	std::cout << prod << std::endl;
		//for (auto& rea : OOQOOHIsomReac)
		//	std::cout << rea << std::endl;
		std::cout << "- OOQOOH isomerization reactions added.       ("
			<< OOQOOHIsomReac.size() << " reactions)" << std::endl;

		// ### ISOMERIZATION P(OOH)2 (P(OOH)2 -> OOQOOH) ###
		std::vector<Reaction> POOH2IsomReac = POOH2IsomerizationReactions(POOH2s, &k);
		//for (auto& rea : POOH2IsomReac)
		//	std::cout << rea << std::endl;
		std::cout << "- P(OOH)2 isomerization reactions added.       ("
			<< POOH2IsomReac.size() << " reactions)" << std::endl;

		// ### ROO TO OLEFINS (ROO -> OLE + H2O) ###
		std::vector<Reaction> oleFromROOReac = olefinsFromROOReactions(ROOs, &k);
		std::vector<Molecola> OLEs = getProducts(oleFromROOReac, OLE_);
		//for (auto& rea : oleFromROOReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Olefins from ROO reactions added.       ("
			<< oleFromROOReac.size() << " reactions)" << std::endl;
		//for (auto& prod : OLEs)
		//	std::cout << prod << std::endl;

		// ### RADICALS BETA DECOMPOSITION (R -> R' + OLE) ###
		std::vector<Reaction> RBetaDecReac = RBetaDecompositioReactions(Rs, &k);
		//for (auto& rea : RBetaDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Radicals beta decomposition reactions added.       ("
			<< RBetaDecReac.size() << " reactions)" << std::endl;

		// ### OLEFINS FROM R + O2 (R + O2 -> OLE + HO2) ###
		std::vector<Reaction> oleFromRPlusO2reac = olefinsFromRPlusO2Reactions(Rs, &k);
		getAdditionalProducts(&OLEs, oleFromRPlusO2reac, OLE_);
		//for (auto& rea : oleFromRPlusO2reac)
		//	std::cout << rea << std::endl;
		std::cout << "- Olefins from R + O2 reactions added.       ("
			<< oleFromRPlusO2reac.size() << " reactions)" << std::endl;

		// ### DECOMPOSITION BETA-QOOH (BETA-QOOH -> OLE + HO2) ###
		std::vector<Reaction> betaQOOHDecReac = betaQOOHDecompositionReaction(QOOHs, &k);
		getAdditionalProducts(&OLEs, oleFromRPlusO2reac, OLE_);
		//for (auto& rea : betaQOOHDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- beta-QOOH decomposition reactions added.       ("
			<< betaQOOHDecReac.size() << " reactions)" << std::endl;

		// ### DECOMPOSITION GAMMA-QOOH (GAMMA-QOOH -> C'O + OLE'' + OH) ###
		std::vector<Reaction> gammaQOOHDecReac = gammaQOOHDecompositionReaction(QOOHs, &k);
		//for (auto& rea : gammaQOOHDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- gamma-QOOH decomposition reactions added.       ("
			<< gammaQOOHDecReac.size() << " reactions)" << std::endl;

		// ### DECOMPOSITION DELTA-QOOH (DELTA-QOOH -> Q'OOH + OLE'') ###
		std::vector<Reaction> deltaQOOHDecReac = deltaQOOHDecompositionReaction(QOOHs, &k);
		//for (auto& rea : deltaQOOHDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- delta-QOOH decomposition reactions added.       ("
			<< deltaQOOHDecReac.size() << " reactions)" << std::endl;

		// ### Cyclic ethers from QOOH (QOOH -> cEth + OH) ###
		std::vector<Reaction> QOOHToCEthReac = QOOHToCEthReactions(QOOHs, &k);
		std::vector<Molecola> cEths = getProducts(QOOHToCEthReac, cEth_);
		//for (auto& rea : QOOHToCEthReac)
		//	std::cout << rea << std::endl;
		std::cout << "- QOOH to cyclic ethers reactions added.       ("
			<< QOOHToCEthReac.size() << " reactions)" << std::endl;

		// ### Cyclic ethers-OOH from P(OOH)2 (P(OOH)2 -> cEthOOH + OH) ###
		std::vector<Reaction> POOH2ToCEthOOHReac = POOH2ToCEthOOHReactions(POOH2s, &k);
		std::vector<Molecola> cEthOOHs = getProducts(POOH2ToCEthOOHReac, cEthOOH_);
		//for (auto& rea : POOH2ToCEthOOHReac)
		//	std::cout << rea << std::endl;
		std::cout << "- P(OOH)2 to cyclic ethers-OOH reactions added.       ("
			<< POOH2ToCEthOOHReac.size() << " reactions)" << std::endl;

		// ### BETA-P(OOH)2 DECOMPOSITION (BETA-P(OOH)2 -> OLE-OOH + HO2) ###
		std::vector<Reaction> betaPOOH2DecReac = betaPOOH2DecompositionReaction(POOH2s, &k);
		std::vector<Molecola> OLEOOHs = getProducts(betaPOOH2DecReac, oleOOH_);
		//for (auto& rea : betaPOOH2DecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Beta-P(OOH)2 decomposition reactions added.       ("
			<< betaPOOH2DecReac.size() << " reactions)" << std::endl;

		// ### GAMMA-P(OOH)2 DECOMPOSITION (GAMMA-P(OOH)2 -> OLE' + C''O + OH) ###
		std::vector<Reaction> gammaPOOH2DecReac = gammaPOOH2DecompositionReaction(POOH2s, &k);
		//for (auto& rea : gammaPOOH2DecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Gamma-P(OOH)2 decomposition reactions added.       ("
			<< gammaPOOH2DecReac.size() << " reactions)" << std::endl;

		// ### DELTA-P(OOH)2 DECOMPOSITION (DELTAP(OOH)2 -> OLE' + P(OOH)2 + OH) ###
		std::vector<Reaction> deltaPOOH2DecReac = deltaPOOH2DecompositionReaction(POOH2s, &k);
		//for (auto& rea : deltaPOOH2DecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Delta-P(OOH)2 decomposition reactions added.       ("
			<< deltaPOOH2DecReac.size() << " reactions)" << std::endl;

		// ### CYCLICETHER-OOH DECOMPOSITION ###
		std::vector<Reaction> cEthOOHDecReac = cyclicEtherOOHDecompositionReactions(cEthOOHs, &k);
		//for (auto& rea : cEthOOHDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Cyclic ether-OOH decomposition reactions added.       ("
			<< cEthOOHDecReac.size() << " reactions)" << std::endl;

		// ### OOQOOH TO OLE-OOH (OOQOOH -> OLE-OOH + HO2) ###
		std::vector<Reaction> OOQOOHToOLEOOHReac = OOQOOHToOLEOOHReactions(OOQOOHs, &k);
		getAdditionalProducts(&OLEOOHs, OOQOOHToOLEOOHReac, oleOOH_);
		//for (auto& rea : OOQOOHToOLEOOHReac)
		//	std::cout << rea << std::endl;
		std::cout << "- OOQOOH to OLE-OOH reactions added.       ("
			<< OOQOOHToOLEOOHReac.size() << " reactions)" << std::endl;

		// ### OLE-OOH DECOMPOSITION (OLE-OOH -> DEC.PROD. + OH) ###
		std::vector<Reaction> OLEOOHDecReac = OLEOOHDecompositionReactions(OLEOOHs, &k);
		//for (auto& rea : OLEOOHDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- OLE-OOH decomposition reactions added.       ("
			<< OLEOOHDecReac.size() << " reactions)" << std::endl;

		// ### KETOHYDROPEROXYDES FORMATION (OOQOOH -> KHP + OH) ###
		std::vector<Reaction> KHPFormReac = KHPFormationReactions(OOQOOHs, &k);
		std::vector<Molecola> KHPs = getProducts(KHPFormReac, KHP_);
		//for (auto& rea : KHPFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Ketohydroperoxides formation reactions added.       ("
			<< KHPFormReac.size() << " reactions)" << std::endl;

		// ### KETOHYDROPEROXIDES DECOMPOSITION (KHP -> DEC.PROD. + OH) ###
		std::vector<Reaction> KHPDecReac = KHPDecompositionReactions(KHPs, &k);
		//for (auto& rea : KHPDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Ketohydroperoxides decomposition reactions added.       ("
			<< KHPDecReac.size() << " reactions)" << std::endl;

		// ### CYCLIC ETHERS DECOMPOSITION (CETH + OH -> DEC.PRO. + H2O) ###
		std::vector<Reaction> cEthDecReac = cyclicEthersDecompositionReactions(cEths, &k);
		//for (auto& rea : cEthDecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- Cyclic ethers decomposition reactions added.       ("
			<< cEthDecReac.size() << " reactions)" << std::endl;

		// ### ALLYLIC RADICAL FORMATION (OLE + OH -> ALLR + H2O) ###
		std::vector<Reaction> allRadFormRac = allylicRadicalsFormationReactions(OLEs, &k);
		std::vector<Molecola> AllRs = getProducts(allRadFormRac, oleR_);
		//for (auto& prod : AllRs)
		//	std::cout << prod << std::endl;
		//for (auto& rea : allRadFormRac)
		//	std::cout << rea << std::endl;
		std::cout << "- Allylic radicals formation reactions added.       ("
			<< allRadFormRac.size() << " reactions)" << std::endl;

		// ### ALKENYL RO FORMATION (ALLR + HO2 -> ALKRO OH) ###
		std::vector<Reaction> alkROFormReac = alkenylROFormationReactions(AllRs, &k);
		std::vector<Molecola> AlkROs = getProducts(alkROFormReac, alkRO_);
		//for (auto& prod : AlkROs)
		//	std::cout << prod << std::endl;
		//for (auto& rea : alkROFormReac)
		//	std::cout << rea << std::endl;
		std::cout << "- AlkenylRO formation reactions added.       ("
			<< alkROFormReac.size() << " reactions)" << std::endl;

		// ### ALKENYL RO DECOMPOSITION (ALKRO -> DEC.PROD.) ###
		std::vector<Reaction> alkRODecReac = alkenylRODecompositionReactions(AlkROs, &k);
		//for (auto& rea : alkRODecReac)
		//	std::cout << rea << std::endl;
		std::cout << "- AlkenylRO decomposition reactions added.       ("
			<< alkRODecReac.size() << " reactions)" << std::endl;

		// #########################################################################
		// ############### PROCESS SPECIES AND REACTIONS ###########################
		// #########################################################################

		std::cout << std::endl;
		UTL::printEmbeddedString('*', "Processing species and reactions");
		std::vector<Molecola> allSpecies;
		UTL::concatenate(&allSpecies, &specialMols);
		addNewSpecies(&allSpecies, &initReac);
		addNewSpecies(&allSpecies, &hAbsReac);
		addNewSpecies(&allSpecies, &O2addToRReac);
		if(reversibeDetailed == false)
			addNewSpecies(&allSpecies, &O2ElimFromROOReac);
		addNewSpecies(&allSpecies, &O2addToQOOHReac);
		if (reversibeDetailed == false)
			addNewSpecies(&allSpecies, &O2ElimFromOOQOOHReac);
		addNewSpecies(&allSpecies, &RIsomReac);
		addNewSpecies(&allSpecies, &ROOIsomReac);
		if (reversibeDetailed == false)
			addNewSpecies(&allSpecies, &QOOHIsomReac);
		addNewSpecies(&allSpecies, &OOQOOHIsomReac);
		if (reversibeDetailed == false)
			addNewSpecies(&allSpecies, &POOH2IsomReac);
		addNewSpecies(&allSpecies, &oleFromROOReac);
		addNewSpecies(&allSpecies, &RBetaDecReac);
		addNewSpecies(&allSpecies, &oleFromRPlusO2reac);
		addNewSpecies(&allSpecies, &betaQOOHDecReac);
		addNewSpecies(&allSpecies, &gammaQOOHDecReac);
		addNewSpecies(&allSpecies, &deltaQOOHDecReac);
		addNewSpecies(&allSpecies, &QOOHToCEthReac);
		addNewSpecies(&allSpecies, &POOH2ToCEthOOHReac);
		addNewSpecies(&allSpecies, &betaPOOH2DecReac);
		addNewSpecies(&allSpecies, &gammaPOOH2DecReac);
		addNewSpecies(&allSpecies, &deltaPOOH2DecReac);
		addNewSpecies(&allSpecies, &cEthOOHDecReac);
		addNewSpecies(&allSpecies, &OOQOOHToOLEOOHReac);
		addNewSpecies(&allSpecies, &OLEOOHDecReac);
		addNewSpecies(&allSpecies, &KHPFormReac);
		addNewSpecies(&allSpecies, &KHPDecReac);
		addNewSpecies(&allSpecies, &cEthDecReac);
		addNewSpecies(&allSpecies, &allRadFormRac);
		addNewSpecies(&allSpecies, &alkROFormReac);
		addNewSpecies(&allSpecies, &alkRODecReac);

		for (auto& mol : allSpecies)
			mol.fix();

		std::vector<std::string> allNames(allSpecies.size());
		for (int i = 0; i < allSpecies.size(); i++)
		{
			//std::cout << allSpecies[i] << std::endl;
			allNames[i] = chemOut.molToName(allSpecies[i]);
		}

		//for (int i = 0; i < allNames.size(); i++)
		//	for (int j = i + 1; j < allNames.size(); j++)
		//		if(allNames[i] == allNames[j])
		//		{
		//			UTL::error("Two same names found");
		//			std::cout << allNames[i] << "  " << allSpecies[i] << "   " << allSpecies[j] << std::endl;
		//		}


		std::cout << " - List of all species compiled: " << allSpecies.size()
			<< " species." << std::endl;

		int dupReac = 0;
		dupReac += processDuplicateReactions(&initReac);
		dupReac += processDuplicateReactions(&hAbsReac);
		dupReac += processDuplicateReactions(&O2addToRReac);
		dupReac += processDuplicateReactions(&O2ElimFromROOReac);
		dupReac += processDuplicateReactions(&O2addToQOOHReac);
		dupReac += processDuplicateReactions(&O2ElimFromOOQOOHReac);
		dupReac += processDuplicateReactions(&RIsomReac);
		dupReac += processDuplicateReactions(&ROOIsomReac);
		dupReac += processDuplicateReactions(&QOOHIsomReac);
		dupReac += processDuplicateReactions(&OOQOOHIsomReac);
		dupReac += processDuplicateReactions(&POOH2IsomReac);
		dupReac += processDuplicateReactions(&oleFromROOReac);
		dupReac += processDuplicateReactions(&RBetaDecReac);
		dupReac += processDuplicateReactions(&oleFromRPlusO2reac);
		dupReac += processDuplicateReactions(&betaQOOHDecReac);
		dupReac += processDuplicateReactions(&gammaQOOHDecReac);
		dupReac += processDuplicateReactions(&deltaQOOHDecReac);
		dupReac += processDuplicateReactions(&QOOHToCEthReac);
		dupReac += processDuplicateReactions(&POOH2ToCEthOOHReac);
		dupReac += processDuplicateReactions(&betaPOOH2DecReac);
		dupReac += processDuplicateReactions(&gammaPOOH2DecReac);
		dupReac += processDuplicateReactions(&deltaPOOH2DecReac);
		dupReac += processDuplicateReactions(&cEthOOHDecReac);
		dupReac += processDuplicateReactions(&OOQOOHToOLEOOHReac);
		dupReac += processDuplicateReactions(&OLEOOHDecReac);
		dupReac += processDuplicateReactions(&KHPFormReac);
		dupReac += processDuplicateReactions(&KHPDecReac);
		dupReac += processDuplicateReactions(&cEthDecReac);
		dupReac += processDuplicateReactions(&allRadFormRac);
		dupReac += processDuplicateReactions(&alkROFormReac);
		dupReac += processDuplicateReactions(&alkRODecReac);

		std::cout << " - Duplicate reactions merged: " << dupReac <<
			" reactions merged" << std::endl;


		if (reversibeDetailed)
		{
			//for (auto& reac : initReac)
			//	reac.setReversible(true);
			//for (auto& reac : hAbsReac)
			//	reac.setReversible(true);
			for (auto& reac : RIsomReac)
				reac.setReversible(true);
			for (auto& reac : O2addToRReac)
				reac.setReversible(true);
			for (auto& reac : O2addToQOOHReac)
				reac.setReversible(true);
			for (auto& reac : ROOIsomReac)
				reac.setReversible(true);
			for (auto& reac : OOQOOHIsomReac)
				reac.setReversible(true);
			//for (auto& reac : oleFromROOReac)
			//	reac.setReversible(true);
			//for (auto& reac : RBetaDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : oleFromRPlusO2reac)
			//	reac.setReversible(true);
			//for (auto& reac : betaQOOHDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : gammaQOOHDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : deltaQOOHDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : QOOHToCEthReac)
			//	reac.setReversible(true);
			//for (auto& reac : POOH2ToCEthOOHReac)
			//	reac.setReversible(true);
			//for (auto& reac : betaPOOH2DecReac)
			//	reac.setReversible(true);
			//for (auto& reac : gammaPOOH2DecReac)
			//	reac.setReversible(true);
			//for (auto& reac : deltaPOOH2DecReac)
			//	reac.setReversible(true);
			//for (auto& reac : cEthOOHDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : OOQOOHToOLEOOHReac)
			//	reac.setReversible(true);
			//for (auto& reac : OLEOOHDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : KHPFormReac)
			//	reac.setReversible(true);
			//for (auto& reac : KHPDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : cEthDecReac)
			//	reac.setReversible(true);
			//for (auto& reac : allRadFormRac)
			//	reac.setReversible(true);
			//for (auto& reac : alkROFormReac)
			//	reac.setReversible(true);
			//for (auto& reac : alkRODecReac)
			//	reac.setReversible(true);
		}

		std::vector<Reaction> allReactions;
		UTL::concatenate(&allReactions, &initReac);
		UTL::concatenate(&allReactions, &hAbsReac);
		UTL::concatenate(&allReactions, &O2addToRReac);
		if (reversibeDetailed == false)
			UTL::concatenate(&allReactions, &O2ElimFromROOReac);
		UTL::concatenate(&allReactions, &O2addToQOOHReac);
		if (reversibeDetailed == false)
			UTL::concatenate(&allReactions, &O2ElimFromOOQOOHReac);
		UTL::concatenate(&allReactions, &RIsomReac);
		UTL::concatenate(&allReactions, &ROOIsomReac);
		if (reversibeDetailed == false)
			UTL::concatenate(&allReactions, &QOOHIsomReac);
		UTL::concatenate(&allReactions, &OOQOOHIsomReac);
		if (reversibeDetailed == false)
			UTL::concatenate(&allReactions, &POOH2IsomReac);
		UTL::concatenate(&allReactions, &oleFromROOReac);
		UTL::concatenate(&allReactions, &RBetaDecReac);
		UTL::concatenate(&allReactions, &oleFromRPlusO2reac);
		UTL::concatenate(&allReactions, &betaQOOHDecReac);
		UTL::concatenate(&allReactions, &gammaQOOHDecReac);
		UTL::concatenate(&allReactions, &deltaQOOHDecReac);
		UTL::concatenate(&allReactions, &QOOHToCEthReac);
		UTL::concatenate(&allReactions, &POOH2ToCEthOOHReac);
		UTL::concatenate(&allReactions, &betaPOOH2DecReac);
		UTL::concatenate(&allReactions, &gammaPOOH2DecReac);
		UTL::concatenate(&allReactions, &deltaPOOH2DecReac);
		UTL::concatenate(&allReactions, &cEthOOHDecReac);
		UTL::concatenate(&allReactions, &OOQOOHToOLEOOHReac);
		UTL::concatenate(&allReactions, &OLEOOHDecReac);
		UTL::concatenate(&allReactions, &KHPFormReac);
		UTL::concatenate(&allReactions, &KHPDecReac);
		UTL::concatenate(&allReactions, &cEthDecReac);
		UTL::concatenate(&allReactions, &allRadFormRac);
		UTL::concatenate(&allReactions, &alkROFormReac);
		UTL::concatenate(&allReactions, &alkRODecReac);

		std::cout << " - List of reactions created: " << allReactions.size() <<
			" reactions" << std::endl;


		// #########################################################################
		// ####################### PRINT SUBMECHANISM ##############################
		// #########################################################################
		std::cout << std::endl;
		UTL::printEmbeddedString('*', "Printing detailed mechanism");
		std::string detailedMechanismPath = subMechFold + "\\kinetic_detailed.inp";
		std::ofstream detailedMechanism(detailedMechanismPath);
		detailedMechanism << "! " << chemOut.molToName(HC) << " COMBUSTION MECHANISM" << std::endl;
		detailedMechanism << "! " << HC << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "ELEMENTS" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "C" << std::endl;
		detailedMechanism << "H" << std::endl;
		detailedMechanism << "O" << std::endl;
		detailedMechanism << "N" << std::endl;
		detailedMechanism << "AR" << std::endl;
		detailedMechanism << "HE" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "END" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "SPECIES" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "! Generic species" << std::endl;
		//detailedMechanism << "O2     H2    H2O2    H2O   N2   AR    HE   CO   CO2" 
		detailedMechanism << "O2                 H2                 H2O2               H2O"
			<< std::endl;
		detailedMechanism << "N2                 AR                 HE                 CO"
			<< std::endl;
		detailedMechanism << "CO2                HO2                OH                 H"
			<< std::endl;
		detailedMechanism << "CH3                C2H5 " << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "!++++++++++++++++++ " << chemOut.molToName(HC)
			<< " SUBMECHANISM ++++++++++++++++++" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "! Fuel" << std::endl;
		detailedMechanism << chemOut.molToName(HC) << std::endl;
		detailedMechanism << "" << std::endl;
		printSpeciesInFile(&detailedMechanism, Rs, "R alkyl radicals", &chemOut);
		printSpeciesInFile(&detailedMechanism, ROOs, "ROO alkyl peroxy radicals",
			&chemOut);
		printSpeciesInFile(&detailedMechanism, QOOHs, "QOOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, OOQOOHs, "OOQOOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, POOH2s, "P(OOH)2", &chemOut);
		printSpeciesInFile(&detailedMechanism, KHPs, "KHP ketohydroperoxides",
			&chemOut);
		printSpeciesInFile(&detailedMechanism, OLEs, "OLE olefins", &chemOut);
		printSpeciesInFile(&detailedMechanism, cEths, "cEth cyclic ethers",
			&chemOut);
		printSpeciesInFile(&detailedMechanism, cEthOOHs, "CEth-OOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, OLEOOHs, "OLE-OOH ", &chemOut);
		printSpeciesInFile(&detailedMechanism, AllRs, "AllR allylic radicals",
			&chemOut);
		printSpeciesInFile(&detailedMechanism, AlkROs, "AlkRO alkenyl RO",
			&chemOut);
		detailedMechanism << "" << std::endl;

		//std::vector<Molecola> fuel_dec(0);
		std::vector<Molecola> R_dec(0);
		std::vector<Molecola> ROO_dec(0);
		std::vector<Molecola> QOOH_dec(0);
		std::vector<Molecola> OOQOOH_dec(0);
		std::vector<Molecola> OLE_dec(0);
		std::vector<Molecola> CO_dec(0);
		std::vector<Molecola> cEth_dec(0);
		std::vector<Molecola> RO_dec(0);
		std::vector<Molecola> KHP_dec(0);
		std::vector<Molecola> POOH2_dec(0);
		std::vector<Molecola> oleOOH_dec(0);
		std::vector<Molecola> cEthOOH_dec(0);
		std::vector<Molecola> oleR_dec(0);
		std::vector<Molecola> oleCO_dec(0);
		std::vector<Molecola> cEthR_dec(0);
		std::vector<Molecola> cEthCO_dec(0);
		std::vector<Molecola> alkRO_dec(0);
		std::vector<Molecola> special_dec(0);


		for (auto& mol : allSpecies)
		{
			if (mol == CH3 || mol == CH4 || mol == C2H5 || mol == C2H6)
				continue;
			if (mol.parentFuel() == HC)
				continue;
			species typ = mol.kindOfSPecies();
			switch (typ)
			{
				//case fuel_:
				//	fuel_dec.push_back(mol);
				//	break;
			case R_:
				R_dec.push_back(mol);
				break;
			case ROO_:
				ROO_dec.push_back(mol);
				break;
			case QOOH_:
				QOOH_dec.push_back(mol);
				break;
			case OOQOOH_:
				OOQOOH_dec.push_back(mol);
				break;
			case OLE_:
				OLE_dec.push_back(mol);
				break;
			case CO_:
				CO_dec.push_back(mol);
				break;
			case cEth_:
				cEth_dec.push_back(mol);
				break;
			case RO_:
				RO_dec.push_back(mol);
				break;
			case KHP_:
				KHP_dec.push_back(mol);
				break;
			case POOH2_:
				POOH2_dec.push_back(mol);
				break;
			case oleOOH_:
				oleOOH_dec.push_back(mol);
				break;
			case cEthOOH_:
				cEthOOH_dec.push_back(mol);
				break;
			case oleR_:
				oleR_dec.push_back(mol);
				break;
			case oleCO_:
				oleCO_dec.push_back(mol);
				break;
			case cEthR_:
				cEthR_dec.push_back(mol);
				break;
			case cEthCO_:
				cEthCO_dec.push_back(mol);
				break;
			case alkRO_:
				alkRO_dec.push_back(mol);
				break;
			case special_:
				break;
			default:
				UTL::warning("The following decompositino species should not be present");
				std::cout << mol << std::endl;
				break;
			}
		}
		detailedMechanism << "!++++++++++++++++++ DECOMPOSITION PRODUCTS ++++++++++++++++++"
			<< std::endl;
		printSpeciesInFile(&detailedMechanism, R_dec, "R", &chemOut);
		printSpeciesInFile(&detailedMechanism, ROO_dec, "ROO", &chemOut);
		printSpeciesInFile(&detailedMechanism, QOOH_dec, "QOOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, OOQOOH_dec, "OOQOOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, OLE_dec, "OLE", &chemOut);
		printSpeciesInFile(&detailedMechanism, CO_dec, "CO aldehydes/ketones", &chemOut);
		printSpeciesInFile(&detailedMechanism, cEth_dec, "cEth", &chemOut);
		printSpeciesInFile(&detailedMechanism, RO_dec, "RO", &chemOut);
		printSpeciesInFile(&detailedMechanism, KHP_dec, "KHP", &chemOut);
		printSpeciesInFile(&detailedMechanism, POOH2_dec, "POOH2", &chemOut);
		printSpeciesInFile(&detailedMechanism, oleOOH_dec, "OLE-OOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, cEthOOH_dec, "cEth-OOH", &chemOut);
		printSpeciesInFile(&detailedMechanism, oleR_dec, "OLE-R", &chemOut);
		printSpeciesInFile(&detailedMechanism, oleCO_dec, "OLE-CO", &chemOut);
		printSpeciesInFile(&detailedMechanism, cEthR_dec, "cEth-R", &chemOut);
		printSpeciesInFile(&detailedMechanism, cEthCO_dec, "cEth-CO", &chemOut);
		printSpeciesInFile(&detailedMechanism, alkRO_dec, "AlkRO", &chemOut);
		detailedMechanism << "" << std::endl;
		detailedMechanism << "END" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "" << std::endl;
		detailedMechanism << "REACTIONS" << std::endl;
		detailedMechanism << std::endl << "! Initiation reactions:" << std::endl;
		printReactions(&detailedMechanism, initReac, &chemOut);
		detailedMechanism << std::endl << "! H abstraction:" << std::endl;
		printReactions(&detailedMechanism, hAbsReac, &chemOut);
		detailedMechanism << std::endl << "! Oxygen addition to alkyl radicals:"
			<< std::endl;
		printReactions(&detailedMechanism, O2addToRReac, &chemOut);
		if (reversibeDetailed == false)
		{
		detailedMechanism << std::endl << "! Oxygen elimination from ROO:" << std::endl;
		printReactions(&detailedMechanism, O2ElimFromROOReac, &chemOut);
		}
		detailedMechanism << std::endl << "! Oxygen addition to QOOH:" << std::endl;
		printReactions(&detailedMechanism, O2addToQOOHReac, &chemOut);
		if (reversibeDetailed == false)
		{
		detailedMechanism << std::endl << "! Oxygen elimination from OOQOOH:"
			<< std::endl;
		printReactions(&detailedMechanism, O2ElimFromOOQOOHReac, &chemOut);
		}
		detailedMechanism << std::endl << "! Alkyl radicals isomerization:" << std::endl;
		printReactions(&detailedMechanism, RIsomReac, &chemOut);
		detailedMechanism << std::endl << "! ROO isomerization:" << std::endl;
		printReactions(&detailedMechanism, ROOIsomReac, &chemOut);
		if (reversibeDetailed == false)
		{
		detailedMechanism << std::endl << "! QOOH isomerization:" << std::endl;
		printReactions(&detailedMechanism, QOOHIsomReac, &chemOut);
		}
		detailedMechanism << std::endl << "! OOQOOH isomerization:" << std::endl;
		printReactions(&detailedMechanism, OOQOOHIsomReac, &chemOut);
		if (reversibeDetailed == false)
		{
		detailedMechanism << std::endl << "! P(OOH)2 isomerization:" << std::endl;
		printReactions(&detailedMechanism, POOH2IsomReac, &chemOut);
		}
		detailedMechanism << std::endl << "! ROO to olefins:" << std::endl;
		printReactions(&detailedMechanism, oleFromROOReac, &chemOut);
		detailedMechanism << std::endl << "! Alkyl radicals beta decomposition:"
			<< std::endl;
		printReactions(&detailedMechanism, RBetaDecReac, &chemOut);
		detailedMechanism << std::endl << "! Olefins formation from R + O2:"
			<< std::endl;
		printReactions(&detailedMechanism, oleFromRPlusO2reac, &chemOut);
		detailedMechanism << std::endl << "! Beta-QOOH decomposition:" << std::endl;
		printReactions(&detailedMechanism, betaQOOHDecReac, &chemOut);
		detailedMechanism << std::endl << "! Gamma-QOOH decomposition:" << std::endl;
		printReactions(&detailedMechanism, gammaQOOHDecReac, &chemOut);
		detailedMechanism << std::endl << "! Delta-QOOH decomposition:" << std::endl;
		printReactions(&detailedMechanism, deltaQOOHDecReac, &chemOut);
		detailedMechanism << std::endl << "! QOOH to cyclic ethers:" << std::endl;
		printReactions(&detailedMechanism, QOOHToCEthReac, &chemOut);
		detailedMechanism << std::endl << "! P(OOH)2 to cyclic ethers-OOH:"
			<< std::endl;
		printReactions(&detailedMechanism, POOH2ToCEthOOHReac, &chemOut);
		detailedMechanism << std::endl << "! Beta-P(OOH)2 decomposition:" << std::endl;
		printReactions(&detailedMechanism, betaPOOH2DecReac, &chemOut);
		detailedMechanism << std::endl << "! Gamma-P(OOH)2 decomposition:" << std::endl;
		printReactions(&detailedMechanism, gammaPOOH2DecReac, &chemOut);
		detailedMechanism << std::endl << "! Delta-P(OOH)2 decomposition:" << std::endl;
		printReactions(&detailedMechanism, deltaPOOH2DecReac, &chemOut);
		detailedMechanism << std::endl << "! Cyclic ethers-OOH decomposition:"
			<< std::endl;
		printReactions(&detailedMechanism, cEthOOHDecReac, &chemOut);
		detailedMechanism << std::endl << "! OOQOOH to olefins-OOH:" << std::endl;
		printReactions(&detailedMechanism, OOQOOHToOLEOOHReac, &chemOut);
		detailedMechanism << std::endl << "! Olefins-OOH decomposition:" << std::endl;
		printReactions(&detailedMechanism, OLEOOHDecReac, &chemOut);
		detailedMechanism << std::endl << "! Ketohydroperoxides formation:" << std::endl;
		printReactions(&detailedMechanism, KHPFormReac, &chemOut);
		detailedMechanism << std::endl << "! Ketohydroperoxides decomposition:"
			<< std::endl;
		printReactions(&detailedMechanism, KHPDecReac, &chemOut);
		detailedMechanism << std::endl << "! Cyclic ethers decomposition:" << std::endl;
		printReactions(&detailedMechanism, cEthDecReac, &chemOut);
		detailedMechanism << std::endl << "! Allylic radicals formation:" << std::endl;
		printReactions(&detailedMechanism, allRadFormRac, &chemOut);
		detailedMechanism << std::endl << "! Alkenyl RO formation:" << std::endl;
		printReactions(&detailedMechanism, alkROFormReac, &chemOut);
		detailedMechanism << std::endl << "! Alkenyl RO decomposition:" << std::endl;
		printReactions(&detailedMechanism, alkRODecReac, &chemOut);
		detailedMechanism << "" << std::endl;
		detailedMechanism << "END" << std::endl;
		detailedMechanism.close();

		std::cout << " - Kinetic mechanism printed." << std::endl;

		// #########################################################################
		// ####################### PRINT THERMO FILE ###############################
		// #########################################################################

		std::string detailedThermoPath = subMechFold + "\\thermo_detailed.inp";
		std::ofstream detailedThermo(detailedThermoPath);
		detailedThermo << "! THERMODYNAMIC DATA FOR " << chemOut.molToName(HC)
			<< " SUBMECHANISM" << std::endl;
		detailedThermo << "! " << HC << std::endl;
		detailedThermo << std::endl;
		detailedThermo << "THERMO" << std::endl;
		detailedThermo << "300.000  1000.000  5000.000" << std::endl;
		detailedThermo << "AR                G         1  0    0      0G   200.000  6000.00  1000.00      1" << std::endl;
		detailedThermo << "+2.50000000E+00+0.00000000E+00+0.00000000E+00+0.00000000E+00+0.00000000E+00    2" << std::endl;
		detailedThermo << "-7.45375000E+02+4.37967491E+00+2.50000000E+00+0.00000000E+00+0.00000000E+00    3" << std::endl;
		detailedThermo << "+0.00000000E+00+0.00000000E+00-7.45375000E+02+4.37967491E+00+0.00000000E+00    4" << std::endl;
		detailedThermo << "N2                G         2    0    0    0G   200.000  6000.00  1000.00      1" << std::endl;
		detailedThermo << "+2.95257637E+00+1.39690040E-03-4.92631603E-07+7.86010195E-11-4.60755204E-15    2" << std::endl;
		detailedThermo << "-9.23948688E+02+5.87188762E+00+3.53100528E+00-1.23660988E-04-5.02999433E-07    3" << std::endl;
		detailedThermo << "+2.43530612E-09-1.40881235E-12-1.04697628E+03+2.96747038E+00+0.00000000E+00    4" << std::endl;
		detailedThermo << "HE                G        1    0    0    0 G   200.000  6000.00  1000.00      1" << std::endl;
		detailedThermo << "+2.50000000E+00+0.00000000E+00+0.00000000E+00+0.00000000E+00+0.00000000E+00    2" << std::endl;
		detailedThermo << "-7.45375000E+02+9.28723974E-01+2.50000000E+00+0.00000000E+00+0.00000000E+00    3" << std::endl;
		detailedThermo << "+0.00000000E+00+0.00000000E+00-7.45375000E+02+9.28723974E-01+0.00000000E+00    4" << std::endl;
		for (auto& spec : allSpecies)
		{
			if (spec == N2)
				continue;
			detailedThermo << thermOut.NASAOutput(&spec);
		}
		detailedThermo << "END";
		detailedThermo.close();
		std::cout << " - Thermodynamic data file printed." << std::endl;

		UTL::concatenateUnique(&totalReactionsList, &allReactions);
		UTL::concatenateUnique(&totalSpeciesList, &allSpecies);

		currentTime = std::chrono::steady_clock::now();
		timeFile << "Detailed generation for " << chemOut.molToName(HC) << " ended at " <<
			std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count() <<
			" s" << std::endl;

		// #########################################################################
		// ########################### LUMP MECHANISM ##############################
		// #########################################################################
		if (generateLumped)
		{
			UTL::printTitle("GENERATING LUMPED MECHANISM");
			// generate mechanisms
			std::string lumpedMechFilePath = subMechFold + "\\kinetic_lumped.inp";
			std::ofstream lumpedMech(lumpedMechFilePath);
			
			UTL::createDirectory("temp");
			UTL::createDirectory("temp/CKmech");
			std::ofstream kinCK("temp/CKmech/kin.txt");
			kinCK << "ELEMENTS" << std::endl;
			kinCK << "C" << std::endl << "H" << std::endl << "O" << std::endl
				<< "N" << std::endl << "AR" << std::endl << "HE" << std::endl;
			kinCK << "END" << std::endl << std::endl;
			kinCK << "SPECIES" << std::endl;
			kinCK << baseMech.speciesText.str() << std::endl;
			for (auto& name : allNames)
				if (UTL::isPresent(&baseMech.mechSpecies, name) == false)
					kinCK << name << std::endl;
			kinCK << "END" << std::endl << std::endl;
			kinCK << "REACTIONS" << std::endl;
			if (lumpWithCoreMech)
				kinCK << baseMech.reactionsText.str() << std::endl;
			for (auto& reac : allReactions)
			{
				if (isIncluded(reac, &(baseMech.reacs), &chemOut) == false)
					printReaction(&kinCK, reac, &chemOut);
			}
			kinCK << "END" << std::endl << std::endl;
			kinCK.close();

			std::ofstream therCK("temp/CKmech/ther.txt");
			therCK << "THERMO" << std::endl;
			therCK << "300.000  1000.000  5000.000" << std::endl;
			therCK << baseMech.thermoText.str() << std::endl;
			for (auto& spec : allSpecies)
				if (spec.isSpecialMolecule() == 0)
					therCK << thermOut.NASAOutput(&spec);
			kinCK << "END" << std::endl;
			therCK.close();

			// run the simulations
			//Simulation sim("temp/CKmech/kin.txt", "temp/CKmech/ther.txt", Temps, Press,
			//	&allSpecies, &chemOut, numCores);
			Simulation* sim;

			if (equilibriumLumping)
			{
				sim = new Simulation(Temps, &allSpecies, &chemOut, &thermOut, numCores);
			}
			else 
			{
				sim = new Simulation("temp/CKmech/kin.txt", "temp/CKmech/ther.txt", Temps, Press,
					&allSpecies, &chemOut, numCores);
				if (useBatch)
				{
					sim->setBatchParameters(HC, eqRatio, 5000, MAX_SIMULATION_TIME);
				}
				else
				{
					sim->setCSTRParameters(HC, eqRatio, tau, 5000, MAX_SIMULATION_TIME);
				}
			}

			sim->solve();
			std::string resultsDirPath = subMechFold + "\\simResults";
			UTL::createDirectory(resultsDirPath);
			sim->printResults(resultsDirPath);
			
			std::ofstream targetRatesFile(resultsDirPath + "\\targetRates.txt");
			// lump reactions

			if (initReac.size() > 0)
			{
				LumpedReaction initReacLump(initReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(initReacLump);
				lumpedMech << initReacLump.print();
				targetRatesFile << initReacLump.printTargetRates();
			}
			if (hAbsReacByO2.size() > 0)
			{
				LumpedReaction hAbsReacByO2Lump(hAbsReacByO2, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(hAbsReacByO2Lump);
				lumpedMech << hAbsReacByO2Lump.print();
				targetRatesFile << hAbsReacByO2Lump.printTargetRates();
			}
			if (hAbsReacByOH.size() > 0)
			{
				LumpedReaction hAbsReacByOHLump(hAbsReacByOH, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(hAbsReacByOHLump);
				lumpedMech << hAbsReacByOHLump.print();
				targetRatesFile << hAbsReacByOHLump.printTargetRates();
			}
			if (hAbsReacByH.size() > 0)
			{
				LumpedReaction hAbsReacByHLump(hAbsReacByH, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(hAbsReacByHLump);
				lumpedMech << hAbsReacByHLump.print();
				targetRatesFile << hAbsReacByHLump.printTargetRates();
			}
			if (hAbsReacByO.size() > 0)
			{
				LumpedReaction hAbsReacByOLump(hAbsReacByO, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(hAbsReacByOLump);
				lumpedMech << hAbsReacByOLump.print();
				targetRatesFile << hAbsReacByOLump.printTargetRates();
			}
			if (hAbsReacByHO2.size() > 0)
			{
				LumpedReaction hAbsReacByHO2Lump(hAbsReacByHO2, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(hAbsReacByHO2Lump);
				lumpedMech << hAbsReacByHO2Lump.print();
				targetRatesFile << hAbsReacByHO2Lump.printTargetRates();
			}
			if (hAbsReacByCH3.size() > 0)
			{
				LumpedReaction hAbsReacByCH3Lump(hAbsReacByCH3, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(hAbsReacByCH3Lump);
				lumpedMech << hAbsReacByCH3Lump.print();
				targetRatesFile << hAbsReacByCH3Lump.printTargetRates();
			}
			if (hAbsReacByC2H5.size() > 0)
			{
				LumpedReaction hAbsReacByC2H5Lump(hAbsReacByC2H5, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(hAbsReacByC2H5Lump);
				lumpedMech << hAbsReacByC2H5Lump.print();
				targetRatesFile << hAbsReacByC2H5Lump.printTargetRates();
			}
			if (O2addToRReac.size() > 0)
			{
				LumpedReaction O2addToRReacLump(O2addToRReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(O2addToRReacLump);
				lumpedMech << O2addToRReacLump.print();
				targetRatesFile << O2addToRReacLump.printTargetRates();
			}
			if (O2ElimFromROOReac.size() > 0)
			{
				LumpedReaction O2ElimFromROOReacLump(O2ElimFromROOReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(O2ElimFromROOReacLump);
				lumpedMech << O2ElimFromROOReacLump.print();
				targetRatesFile << O2ElimFromROOReacLump.printTargetRates();
			}
			if (O2addToQOOHReac.size() > 0)
			{
				LumpedReaction O2addToQOOHReacLump(O2addToQOOHReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(O2addToQOOHReacLump);
				lumpedMech << O2addToQOOHReacLump.print();
				targetRatesFile << O2addToQOOHReacLump.printTargetRates();
			}
			if (O2ElimFromOOQOOHReac.size() > 0)
			{
				LumpedReaction O2ElimFromOOQOOHReacLump(O2ElimFromOOQOOHReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(O2ElimFromOOQOOHReacLump);
				lumpedMech << O2ElimFromOOQOOHReacLump.print();
				targetRatesFile << O2ElimFromOOQOOHReacLump.printTargetRates();
			}
			if (ROOIsomReac.size() > 0)
			{
				LumpedReaction ROOIsomReacLump(ROOIsomReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(ROOIsomReacLump);
				lumpedMech << ROOIsomReacLump.print();
				targetRatesFile << ROOIsomReacLump.printTargetRates();
			}
			if (QOOHIsomReac.size() > 0)
			{
				LumpedReaction QOOHIsomReacLump(QOOHIsomReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(QOOHIsomReacLump);
				lumpedMech << QOOHIsomReacLump.print();
				targetRatesFile << QOOHIsomReacLump.printTargetRates();
			}
			if (OOQOOHIsomReac.size() > 0)
			{
				LumpedReaction OOQOOHIsomReacLump(OOQOOHIsomReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(OOQOOHIsomReacLump);
				lumpedMech << OOQOOHIsomReacLump.print();
				targetRatesFile << OOQOOHIsomReacLump.printTargetRates();
			}
			if (POOH2IsomReac.size() > 0)
			{
				LumpedReaction POOH2IsomReacLump(POOH2IsomReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(POOH2IsomReacLump);
				lumpedMech << POOH2IsomReacLump.print();
				targetRatesFile << POOH2IsomReacLump.printTargetRates();
			}
			if (oleFromROOReac.size() > 0)
			{
				LumpedReaction oleFromROOReacLump(oleFromROOReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(oleFromROOReacLump);
				lumpedMech << oleFromROOReacLump.print();
				targetRatesFile << oleFromROOReacLump.printTargetRates();
			}
			if (RBetaDecReac.size() > 0)
			{
				LumpedReaction RBetaDecReacLump(RBetaDecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(RBetaDecReacLump);
				lumpedMech << RBetaDecReacLump.print();
				targetRatesFile << RBetaDecReacLump.printTargetRates();
			}
			if (oleFromRPlusO2reac.size() > 0)
			{
				LumpedReaction oleFromRPlusO2reacLump(oleFromRPlusO2reac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(oleFromRPlusO2reacLump);
				lumpedMech << oleFromRPlusO2reacLump.print();
				targetRatesFile << oleFromRPlusO2reacLump.printTargetRates();
			}
			if (betaQOOHDecReac.size() > 0)
			{
				LumpedReaction betaQOOHDecReacLump(betaQOOHDecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(betaQOOHDecReacLump);
				lumpedMech << betaQOOHDecReacLump.print();
				targetRatesFile << betaQOOHDecReacLump.printTargetRates();
			}
			if (gammaQOOHDecReac.size() > 0)
			{
				LumpedReaction gammaQOOHDecReacLump(gammaQOOHDecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(gammaQOOHDecReacLump);
				lumpedMech << gammaQOOHDecReacLump.print();
				targetRatesFile << gammaQOOHDecReacLump.printTargetRates();
			}
			if (deltaQOOHDecReac.size() > 0)
			{
				LumpedReaction deltaQOOHDecReacLump(deltaQOOHDecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(deltaQOOHDecReacLump);
				lumpedMech << deltaQOOHDecReacLump.print();
				targetRatesFile << deltaQOOHDecReacLump.printTargetRates();
			}
			if (QOOHToCEthReac.size() > 0)
			{
				LumpedReaction QOOHToCEthReacLump(QOOHToCEthReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(QOOHToCEthReacLump);
				lumpedMech << QOOHToCEthReacLump.print();
				targetRatesFile << QOOHToCEthReacLump.printTargetRates();
			}
			if (POOH2ToCEthOOHReac.size() > 0)
			{
				LumpedReaction POOH2ToCEthOOHReacLump(POOH2ToCEthOOHReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(POOH2ToCEthOOHReacLump);
				lumpedMech << POOH2ToCEthOOHReacLump.print();
				targetRatesFile << POOH2ToCEthOOHReacLump.printTargetRates();
			}
			if (betaPOOH2DecReac.size() > 0)
			{
				LumpedReaction betaPOOH2DecReacLump(betaPOOH2DecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(betaPOOH2DecReacLump);
				lumpedMech << betaPOOH2DecReacLump.print();
				targetRatesFile << betaPOOH2DecReacLump.printTargetRates();
			}
			if (gammaPOOH2DecReac.size() > 0)
			{
				LumpedReaction gammaPOOH2DecReacLump(gammaPOOH2DecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(gammaPOOH2DecReacLump);
				lumpedMech << gammaPOOH2DecReacLump.print();
				targetRatesFile << gammaPOOH2DecReacLump.printTargetRates();
			}
			if (deltaPOOH2DecReac.size() > 0)
			{
				LumpedReaction deltaPOOH2DecReacLump(deltaPOOH2DecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(deltaPOOH2DecReacLump);
				lumpedMech << deltaPOOH2DecReacLump.print();
				targetRatesFile << deltaPOOH2DecReacLump.printTargetRates();
			}
			if (cEthOOHDecReac.size() > 0)
			{
				LumpedReaction cEthOOHDecReacLump(cEthOOHDecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(cEthOOHDecReacLump);
				lumpedMech << cEthOOHDecReacLump.print();
				targetRatesFile << cEthOOHDecReacLump.printTargetRates();
			}
			if (OOQOOHToOLEOOHReac.size() > 0)
			{
				LumpedReaction OOQOOHToOLEOOHReacLump(OOQOOHToOLEOOHReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(OOQOOHToOLEOOHReacLump);
				lumpedMech << OOQOOHToOLEOOHReacLump.print();
				targetRatesFile << OOQOOHToOLEOOHReacLump.printTargetRates();
			}
			if (OLEOOHDecReac.size() > 0)
			{
				LumpedReaction OLEOOHDecReacLump(OLEOOHDecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(OLEOOHDecReacLump);
				lumpedMech << OLEOOHDecReacLump.print();
				targetRatesFile << OLEOOHDecReacLump.printTargetRates();
			}
			if (KHPFormReac.size() > 0)
			{
				LumpedReaction KHPFormReacLump(KHPFormReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(KHPFormReacLump);
				lumpedMech << KHPFormReacLump.print();
				targetRatesFile << KHPFormReacLump.printTargetRates();
			}
			if (KHPDecReac.size() > 0)
			{
				LumpedReaction KHPDecReacLump(KHPDecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(KHPDecReacLump);
				lumpedMech << KHPDecReacLump.print();
				targetRatesFile << KHPDecReacLump.printTargetRates();
			}
			if (cEthDecReac.size() > 0)
			{
				LumpedReaction cEthDecReacLump(cEthDecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(cEthDecReacLump);
				lumpedMech << cEthDecReacLump.print();
				targetRatesFile << cEthDecReacLump.printTargetRates();
			}
			if (allRadFormRac.size() > 0)
			{
				LumpedReaction allRadFormRacLump(allRadFormRac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(allRadFormRacLump);
				lumpedMech << allRadFormRacLump.print();
				targetRatesFile << allRadFormRacLump.printTargetRates();
			}
			if (alkROFormReac.size() > 0)
			{
				LumpedReaction alkROFormReacLump(alkROFormReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(alkROFormReacLump);
				lumpedMech << alkROFormReacLump.print();
				targetRatesFile << alkROFormReacLump.printTargetRates();
			}
			if (alkRODecReac.size() > 0)
			{
				LumpedReaction alkRODecReacLump(alkRODecReac, sim, &thermOut, &chemOut);
				totalLumpedReactionsList.push_back(alkRODecReacLump);
				lumpedMech << alkRODecReacLump.print();
				targetRatesFile << alkRODecReacLump.printTargetRates();
			}
			targetRatesFile.close();
			
			lumpedMech.close();
			currentTime = std::chrono::steady_clock::now();
			timeFile << "Lumped generation for " << chemOut.molToName(HC) << " ended at " <<
				std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count() <<
				" s" << std::endl;
		}
	}

	// #########################################################################
	// ######################## FILL MISSING REACTIONS #########################
	// #########################################################################
	//std::vector<bool> isSpeciesNonDecomposing(totalSpeciesList.size(), false);
	std::cout << std::endl;
	UTL::printEmbeddedString('#', "PROCESSING COMPLETE MECHANISM");

	int numOfNonConsumingSpecies = 0;
	int initialReacNum = totalReactionsList.size();
	int initialSpeciesNum = totalSpeciesList.size();
	std::vector<Reaction> newlyAddedReactions;
	for (int i = 0; i < totalSpeciesList.size(); i++)
	{
		Molecola spec = totalSpeciesList[i];
		//std::cout << "spec: " << spec << std::endl;
		if (!isThereDecompositionPath(spec, &(baseMech.reacs), &totalReactionsList,
			&chemOut))
		{
			numOfNonConsumingSpecies++;
			std::vector<Reaction> decReac;
			std::vector<Reaction> tempReac;
			std::vector<Molecola> reactantVector = { spec };
			species typ = spec.kindOfSPecies();
			switch (typ)
			{
			case R_:
				decReac = RBetaDecompositioReactions(reactantVector, &k);
				tempReac = olefinsFromRPlusO2Reactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = O2AdditionToRReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				break;
			case ROO_:
				decReac = olefinsFromROOReactions(reactantVector, &k);
				tempReac = ROOIsomerizationReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = O2EliminationFromROOReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				break;
			case QOOH_:
				decReac = betaQOOHDecompositionReaction(reactantVector, &k);
				tempReac = gammaQOOHDecompositionReaction(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = deltaQOOHDecompositionReaction(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = QOOHToCEthReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = O2AdditionToQOOHReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = QOOHIsomerizationReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				break;
			case OOQOOH_:
				decReac = KHPFormationReactions(reactantVector, &k);
				tempReac = OOQOOHToOLEOOHReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = O2EliminationFromOOQOOHReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = OOQOOHIsomerizationReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				break;
			case OLE_:
				decReac = allylicRadicalsFormationReactions(reactantVector, &k);
				break;
			case CO_:
				if (spec.isAldehyde())
					decReac = aldehydesDecompositionReactions(reactantVector, &k);
				else if (spec.isKetone())
				{
					decReac = ketonesDecompositionReactions(reactantVector, &k);
					//std::cout << "Keto found: " << std::endl;
					//for (auto& reac : decReac)
					//	std::cout << "   - " << reac << std::endl;
				}
				break;
			case cEth_:
				decReac = cyclicEthersDecompositionReactions(reactantVector, &k);
				break;
			case RO_:
				UTL::warning("RO found");
				break;
			case KHP_:
				decReac = KHPDecompositionReactions(reactantVector, &k);
				break;
			case POOH2_:
				decReac = betaPOOH2DecompositionReaction(reactantVector, &k);
				tempReac = gammaPOOH2DecompositionReaction(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = deltaPOOH2DecompositionReaction(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = POOH2ToCEthOOHReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				tempReac = POOH2IsomerizationReactions(reactantVector, &k);
				UTL::concatenate(&decReac, &tempReac);
				break;
			case oleOOH_:
				decReac = OLEOOHDecompositionReactions(reactantVector, &k);
				break;
			case cEthOOH_:
				decReac = cyclicEtherOOHDecompositionReactions(reactantVector, &k);
				break;
			case oleR_:
				decReac = alkenylROFormationReactions(reactantVector, &k);
				break;
			case oleCO_:
				if (spec.isOleAldehyde()) 
					decReac = aldehydeOlefinsDecompositionReactions(reactantVector, &k);
				else if (spec.isOleKetone())
				{
					decReac = ketonesOlefinsDecompositionReactions(reactantVector, &k);
				}
				break;
			case cEthR_:
				UTL::warning("cEthR found");
				break;
			case cEthCO_:
				UTL::warning("cEthCO found");
				break;
			case alkRO_:
				decReac = alkenylRODecompositionReactions(reactantVector, &k);
				break;
			case special_:
				UTL::warning("Special found");
				break;
			default:
				UTL::error("unexpected species found");
				break;
			}
			processDuplicateReactions(&decReac);
			//for (auto& reac : decReac)
			//	reac.setReversible(true);
			UTL::concatenate(&totalReactionsList, &decReac);
			addNewSpecies(&totalSpeciesList, &decReac);
			UTL::concatenate(&newlyAddedReactions, &decReac);
		}
	}
	std::cout << " - Fixed " << numOfNonConsumingSpecies
		<< " species without consumption pathways." << std::endl;
	std::cout << "      " << (totalReactionsList.size() - initialReacNum)
		<< " reactions added" << std::endl;
	//std::cout << totalReactionsList.size() << std::endl;
	//std::cout << initialReacNum << std::endl;
	std::cout << "      " << (totalSpeciesList.size() - initialSpeciesNum)
		<< " species added" << std::endl;
	//std::cout << totalSpeciesList.size() << std::endl;
	//std::cout << initialSpeciesNum << std::endl;

	// #########################################################################
	// ################# PRINT COMPLETE DETAILED MECHANISM #####################
	// #########################################################################
	std::ofstream fullDetailedMech(outFolderPath + "\\completeDetailedMech.inp");
	fullDetailedMech << "ELEMENTS" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "C" << std::endl;
	fullDetailedMech << "H" << std::endl;
	fullDetailedMech << "N" << std::endl;
	fullDetailedMech << "O" << std::endl;
	fullDetailedMech << "AR" << std::endl;
	fullDetailedMech << "HE" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "END" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "SPECIES" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "" << std::endl;
	std::vector<std::string> baseMechSpeciesNames = { "O2" , "H2", "H2O2",
		"H2O", "N2", "AR", "HE", "CO", "CO2", "HO2", "OH", "H" };
	UTL::concatenateUnique(&baseMechSpeciesNames, &(baseMech.mechSpecies));
	printSpeciesInFile(&fullDetailedMech, baseMechSpeciesNames,
		"Base mechanism");
	//detect all the parent fuels present in the species
	std::vector<Molecola> allParentFules = fuels;
	for (auto& spec : totalSpeciesList)
	{
		if (spec.isSpecialMolecule() == 0)
			UTL::addUnique(&allParentFules, spec.parentFuel());
	}
	std::vector<std::string> classesNames = {
		"R alkyl radicals",
		"ROO alkyl peroxy radicals",
		"QOOH",
		"OOQOOH",
		"P(OOH)2",
		"KHP ketohydroperoxides",
		"OLE olefins",
		"cEth cyclic ethers",
		"cEth-OOH",
		"OLE-OOH",
		"OLE-R",
		"alkRO alkenyl RO",
		"CO aldehydes/ketones",
		"RO",
		"OLE-CO",
		"cEth-R",
		"cEth-CO" };
	std::vector<species> classes =
	{
		R_,
		ROO_,
		QOOH_,
		OOQOOH_,
		POOH2_,
		KHP_,
		OLE_,
		cEth_,
		cEthOOH_,
		oleOOH_,
		oleR_,
		alkRO_,
		CO_,
		RO_,
		oleCO_,
		cEthR_,
		cEthCO_
	};

	for (auto& parFuel : allParentFules)
	{
		std::string fuelName = chemOut.molToName(parFuel);
		fullDetailedMech << "!";
		for (int i = 0; i < (76 - fuelName.size()) / 2; i++)
			fullDetailedMech << "#";
		fullDetailedMech << " " << fuelName << " ";
		for (int i = 0; i < (76 - fuelName.size()) / 2; i++)
			fullDetailedMech << "#";
		fullDetailedMech << std::endl;
		if (UTL::isPresent(&totalSpeciesList, parFuel)
			&& !UTL::isPresent(&baseMechSpeciesNames, chemOut.molToName(parFuel)))
			printSpeciesInFile(&fullDetailedMech, std::vector<Molecola> {parFuel},
				"fuel", & chemOut);
		for (int i = 0; i < classes.size(); i++)
		{
			std::vector<std::string> tempList;
			for (auto& spec : totalSpeciesList)
			{
				if (spec.isSpecialMolecule() != 0)
					continue;
				if (spec.kindOfSPecies() != classes[i])
					continue;
				if (!(spec.parentFuel() == parFuel))
					continue;
				std::string specName = chemOut.molToName(spec);
				if (!UTL::isPresent(&baseMechSpeciesNames, specName))
					tempList.push_back(specName);
			}
			if (tempList.size() > 0)
			{
				printSpeciesInFile(&fullDetailedMech, tempList,
					classesNames[i]);
			}
		}
	}

	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "END" << std::endl;
	fullDetailedMech << "" << std::endl;

	// print reactions
	fullDetailedMech << "REACTIONS" << std::endl;
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "! BASE MECHANISM START" << std::endl;
	fullDetailedMech << baseMech.reactionsText.str() << std::endl;
	fullDetailedMech << "! BASE MECHANISM END" << std::endl;
	fullDetailedMech << "" << std::endl;

	std::vector<std::string> reactionsLabels = {
		"Initiation reactions",
		"H abstraction",
		"O2 addition to R",
		"O2 elimination from ROO",
		"O2 addition to QOOH",
		"O2 elimination from OOQOOH",
		"Isomerization R",
		"Isomerization ROO",
		"Isomerization QOOH",
		"Isomerization OOQOOH",
		"Isomerization P(OOH)2",
		"ROO to olefins",
		"Radicals beta decomposition",
		"Olefins from R + O2",
		"Decomposition beta-QOOH",
		"Decomposition gamma-QOOH",
		"Decomposition delta-QOOH",
		"Cyclic ethers from QOOH",
		"Ethers-OOH from P(OOH)2",
		"Beta-P(OOH)2 decomposition",
		"Gamma-P(OOH)2 decomposition",
		"Delta P(OOH)2 decomposition",
		"Ether-OOH decomposition",
		"Olefins-OOH from OOQOOH",
		"Ole-OOH decomposition",
		"OOQOOH conversion to ketohydroperoxide",
		"Ketohydroperoxides decomposition",
		"Cyclic ethers decomposition",
		"Allylic radicals formation",
		"Alkenyl RO formation",
		"Alkenyl RO decomposition",
		"Aldehydes decomposition",
		"Olefin aldehydes decomposition",
		"Ketones decomposition",
		"Olefins ketones decomposition"
	};
	for (auto& parFuel : allParentFules)
	{
		std::string fuelName = chemOut.molToName(parFuel);
		fullDetailedMech << "!";
		for (int i = 0; i < (76 - fuelName.size()) / 2; i++)
			fullDetailedMech << "#";
		fullDetailedMech << " " << fuelName << " ";
		for (int i = 0; i < (76 - fuelName.size()) / 2; i++)
			fullDetailedMech << "#";
		fullDetailedMech << std::endl;
		for (auto& label : reactionsLabels)
		{
			std::vector<Reaction> tempList;
			for (auto& reac : totalReactionsList)
			{
				if (reac.reactionLabel() == label)
					if (reac.parentFuel() == parFuel)
						if (isIncluded(reac, &(baseMech.reacs), &chemOut) == false)
							tempList.push_back(reac);
			}
			if (tempList.size() > 0)
			{
				fullDetailedMech << "! " << label << std::endl;
				printReactions(&fullDetailedMech, tempList, &chemOut);
				fullDetailedMech << std::endl;
			}
		}
	}
	//for (auto& reac : totalReactionsList)
	//	if(isIncluded(reac, &(baseMech.reacs), &chemOut) == false)
	//		printReaction(&fullDetailedMech, reac, &chemOut);
	fullDetailedMech << "" << std::endl;
	fullDetailedMech << "END" << std::endl;
	fullDetailedMech.close();

	// print thermo file
	std::ofstream fullDetailedThermo(outFolderPath + "\\completeDetailedThermo.inp");
	fullDetailedThermo << "THERMO" << std::endl;
	fullDetailedThermo << "300.000  1000.000  5000.000" << std::endl;
	fullDetailedThermo << "! BASE MECHANISM START" << std::endl;
	fullDetailedThermo << baseMech.thermoText.str() << std::endl;
	fullDetailedThermo << "! BASE MECHANISM END" << std::endl;
	for (auto& spec : totalSpeciesList)
	{
		if (spec == N2)
			continue;
		fullDetailedThermo << thermOut.NASAOutput(&spec);
	}
	fullDetailedThermo << "END";
	fullDetailedThermo.close();

	// PRINT SPECIES LIST
	std::ofstream speciesListFile(outFolderPath + "\\speciesList.csv");
	speciesListFile << "Mech. name,structure,InChI,class,LumpedSpecies" << std::endl;
	std::vector<std::string> totalSpeciesNamesList(totalSpeciesList.size(), "");
	std::vector<std::string> totalSpeciesInChIsList(totalSpeciesList.size(), "");
	for (int i = 0; i < totalSpeciesList.size(); i++)
	{
		if (totalSpeciesList[i].size() > 4)
		{
			totalSpeciesNamesList[i] = chemOut.molToName(totalSpeciesList[i]);
			totalSpeciesInChIsList[i] = totalSpeciesList[i].inchiName();
			speciesListFile << "\"" << totalSpeciesNamesList[i] << "\","
				<< totalSpeciesList[i] << ",\"" <<
				totalSpeciesInChIsList[i] << "\"," <<
				speciesToText(totalSpeciesList[i].kindOfSPecies()) << "," <<
				chemOut.molToNameLump(totalSpeciesList[i]) << std::endl;
		}
	}

	for (int i = 0; i < totalSpeciesNamesList.size(); i++)
	{
		if (totalSpeciesNamesList[i] == "")
			continue;
		for (int j = i + 1; j < totalSpeciesNamesList.size(); j++)
		{
			if (totalSpeciesNamesList[j] == "")
				continue;
			if (totalSpeciesNamesList[i] == totalSpeciesNamesList[j])
			{
				std::stringstream warningMessage;
				warningMessage << "Two different species found having the same naming: "
					<< std::endl << "    - " << totalSpeciesNamesList[i] << "  " <<
					totalSpeciesList[i] << "   " << totalSpeciesInChIsList[i]
					<< std::endl << "    - " << totalSpeciesNamesList[j] << "  " <<
					totalSpeciesList[j] << "   " << totalSpeciesInChIsList[j] << std::endl;
				warningMessage << "    Probably there is a conflict in the glossary file" << std::endl;
				warningMessage << "    Either remove definition for " << totalSpeciesNamesList[i]
					<< " or add definition for missing InChI." << std::endl << std::endl;
				UTL::warning(warningMessage.str());
			}
		}
	}
	speciesListFile.close();

	chemOut.printLongNameSpeciesMessage();

	// #########################################################################
	// ####################### COMPLETE LUMPED MECHANISM #######################
	// #########################################################################
	if (generateLumped)
	{
		//newlyAddedReactions
		for (auto& parFuel : allParentFules)
		{
			std::vector<Reaction> pertinentReaction;
			for (auto& R : newlyAddedReactions)
			{
				std::vector<Molecola*> reacs = R.reactantList();
				int maxSize = 0;
				Molecola relReac;
				for (auto& reac : reacs)
				{
					if (reac->isSpecialMolecule() == 0 && reac->size() > maxSize)
					{
						maxSize = reac->size();
						relReac = *reac;
					}
				}
				if (relReac.parentFuel() == parFuel)
					pertinentReaction.push_back(R);
			}
			std::vector<std::vector<Reaction>> pertReacByType(reactionsLabels.size());
			for (int i = 0; i < reactionsLabels.size(); i++)
			{
				for (auto& R : pertinentReaction)
					if (R.reactionLabel() == reactionsLabels[i])
						pertReacByType[i].push_back(R);
			}
			for (auto& reacVec : pertReacByType)
			{
				if (reacVec.size() == 0)
					continue;
				LumpedReaction lump(reacVec, &thermOut, Temps, &chemOut);
				totalLumpedReactionsList.push_back(lump);
			}
		}

		// generate species list
		std::vector<std::string> allLumpedNames;
		std::vector<Molecola>    allLumpedMols;
		for (auto& spec : totalSpeciesList)
		{
			if (UTL::addUnique(&allLumpedNames, chemOut.molToNameLump(spec)) == -1)
				allLumpedMols.push_back(spec);
		}

		std::vector<std::string> lumpedSpeciesForMech = { "O2" , "H2", "H2O2",
			"H2O", "N2", "AR", "HE", "CO", "CO2", "HO2", "OH", "H" };
		UTL::concatenateUnique(&lumpedSpeciesForMech, &(baseMech.mechSpecies));
		UTL::concatenateUnique(&lumpedSpeciesForMech, &allLumpedNames);

		std::ofstream fullLumpedMech(outFolderPath + "\\completeLumpedMech.inp");
		fullLumpedMech << "ELEMENTS" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "C" << std::endl;
		fullLumpedMech << "H" << std::endl;
		fullLumpedMech << "N" << std::endl;
		fullLumpedMech << "O" << std::endl;
		fullLumpedMech << "AR" << std::endl;
		fullLumpedMech << "HE" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "END" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "SPECIES" << std::endl;
		fullLumpedMech << "" << std::endl;
		printSpeciesInFile(&fullLumpedMech, lumpedSpeciesForMech, "all");
		fullLumpedMech << "END" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "REACTIONS" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech << "! BASE MECHANISM START" << std::endl;
		fullLumpedMech << baseMech.reactionsText.str() << std::endl;
		fullLumpedMech << "! BASE MECHANISM END" << std::endl;
		fullLumpedMech << "" << std::endl;
		for (auto& R : totalLumpedReactionsList)
			fullLumpedMech << R.print() << std::endl;
		fullLumpedMech << "END" << std::endl;
		fullLumpedMech << "" << std::endl;
		fullLumpedMech.close();

		// print thermo file
		std::ofstream fullLumpedThermo(outFolderPath + "\\completeLumpedThermo.inp");
		fullLumpedThermo << "THERMO" << std::endl;
		fullLumpedThermo << "300.000  1000.000  5000.000" << std::endl;
		fullLumpedThermo << "! BASE MECHANISM START" << std::endl;
		fullLumpedThermo << baseMech.thermoText.str() << std::endl;
		fullLumpedThermo << "! BASE MECHANISM END" << std::endl;
		std::vector<std::string> addedSpecies;
		for (auto& reac : totalLumpedReactionsList)
			if (UTL::addUnique(&addedSpecies, reac.relevantReactantName) == -1)
				fullLumpedThermo << reac.relevantReactantThermo;
		for (int i = 0; i < allLumpedNames.size(); i++)
			if (UTL::isPresent(&addedSpecies, allLumpedNames[i]) == false)
				fullLumpedThermo << thermOut.NASAOutputLumped(allLumpedNames[i],
					{ allLumpedMols[i] }, { 300, 5000 }, { {1},{1} });

		fullLumpedThermo << "END" << std::endl;
		fullLumpedThermo << "" << std::endl;
		fullLumpedThermo.close();
	}

	currentTime = std::chrono::steady_clock::now();
	timeFile << "Mechanism merging ended at " <<
		std::chrono::duration_cast<std::chrono::seconds>(currentTime - startTime).count() <<
		" s" << std::endl;
	timeFile.close();
	return 0;
}






//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~ FUNCTION DEFINITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


void generateR(Molecola HC, std::vector<Molecola>* Rvec, std::vector<int>* cToRMap)
{
	std::vector<Molecola> tempRvec;
	std::vector<int> tempMap(HC.size() + 1);
	for (int i = 1; i <= HC.size(); i++)
	{
		// Generate all the radicals for the 1 groups (carbon)       
		// exept the quaternary carbon                               
		if (HC.tipo(i) != 1) continue;
		if (HC.tipoC(i) == Cq) continue;
		Molecola newrad = HC;
		newrad.tipo(i, 2);
		int isomero = 0;
		for (int j = 0; j < tempRvec.size(); j++) // search if it already exist
		{
			if (newrad == tempRvec[j])
			{
				tempRvec[j].isomeri++;
				tempMap[i] = j; // the radical in the i position is equal to the one in 
				//j position 
				isomero = 1;
			}
		}
		if (isomero == 0)  // if it does not already exist...     
		{
			tempRvec.push_back(newrad);
			tempMap[i] = tempRvec.size() - 1;
		}
	}

	*Rvec = tempRvec;
	*cToRMap = tempMap;
}

void generateROO(Molecola HC, std::vector<Molecola>* ROOvec, std::vector<int>* cToROOMap)
{
	generateR(HC, ROOvec, cToROOMap);
	for (int i = 0; i < ROOvec->size(); i++)
	{
		(*ROOvec)[i].Crad_to_ROO((*ROOvec)[i].trova(2));
	}
}

void generateQOOH(Molecola HC, std::vector<Molecola>* QOOHvec,
	std::vector<std::vector<int>>* cToQOOHMap)
{
	std::vector<Molecola> tempQOOHvec;
	std::vector<std::vector<int>> tempMap(HC.size() + 1);
	for (auto& vec : tempMap)
		vec.resize(HC.size() + 1);

	std::vector<Molecola> ROOvec;
	std::vector<int> cToROOMap;
	generateROO(HC, &ROOvec, &cToROOMap);

	for (int i = 0; i < ROOvec.size(); i++)  // i iterates trough the roo   
	{
		int isomero = 0;
		int quat = 0;
		int dist1[SIZEMAX + 1];
		int dist2[SIZEMAX + 1];
		int dist3[SIZEMAX + 1];
		int dist[SIZEMAX + 1];

		for (int k = 0; k <= (SIZEMAX); k++)
			dist1[k] = dist2[k] = dist3[k] = dist[k] = 0;
		//--------------------------------find the position of the i-th COO*      
		int pos_oo = 1;
		while (ROOvec[i].tipo(pos_oo) != 3) pos_oo++;
		//-----------------find all the postions of all the possible radicals...  
		ROOvec[i].scorri(pos_oo, 1, dist1);
		ROOvec[i].scorri(pos_oo, 2, dist2);
		ROOvec[i].scorri(pos_oo, 3, dist3);
		for (int kk = 1; kk <= SIZEMAX; kk++)
			dist[kk] = dist1[kk] + dist2[kk] + dist3[kk];
		//------------------------------------------- ...and saves them in dist[]   

		for (int pos_rad = 1; pos_rad <= ROOvec[i].size(); pos_rad++)
			if (dist[pos_rad] == 1)					//  iterates trough this COO*  
			{										//  all the positions of the   
				isomero = 0;						//  found radicals              
				quat = 0;
				Molecola temp = ROOvec[i];
				temp.tipo(pos_oo, 4);
				temp.tipo(pos_rad, 2);
				if (temp.tipoC(pos_rad) == Cq) quat = 1;    // if the carbon is quaternary              
				// the next time it has to be skipped      
				for (int j = 0; j < tempQOOHvec.size(); j++)
					if (tempQOOHvec[j] == temp)
					{
						isomero = 1;
						tempQOOHvec[j].isomeri++;
						tempMap[pos_oo][pos_rad] = j;
					};   // end if e for j
				if (!isomero && !quat) {
					tempQOOHvec.push_back(temp);
					tempMap[pos_oo][pos_rad] = tempQOOHvec.size() - 1;
				}
			} // end if and for pos_rad

	}  // end for i 

	*QOOHvec = tempQOOHvec;
	*cToQOOHMap = tempMap;
}

void generateOOQOOH(Molecola HC, std::vector<Molecola>* OOQOOHvec,
	std::vector<std::vector<int>>* cToOOQOOHMap)
{
	generateQOOH(HC, OOQOOHvec, cToOOQOOHMap);
	for (auto& mol : (*OOQOOHvec))
		mol.Crad_to_ROO(mol.trova(2));
}

std::vector<Reaction> initiationReactions(Molecola HC, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola* firstRad;			// array of molecules where the first product is stored
	Molecola* secondRad;		// array of molecules where the second product is stored
	Molecola temp1, temp2;		// molecules where to temporary save the products
	firstRad = new Molecola[HC.size() + 1]; if (!firstRad) exit(200);
	secondRad = new Molecola[HC.size() + 1]; if (!secondRad) exit(201);

	// search for the bonds that can be broken
	std::vector < std::vector<int> > atomsInBond;
	atomsInBond = HC.listOfCCBonds();
	//for (int i = 0; i < atomsInBond.size(); i++)
	//	std::cout << atomsInBond[i][0] << "   " << atomsInBond[i][1] << std::endl;


	int numberOfReactions = 0;
	for (int i = 0; i < atomsInBond.size(); i++)
	{
		bool alreadyPresent = false;	// tells if the an equivalent reaction 
		// is already present
		int isPossible = HC.spezza(atomsInBond[i][0], atomsInBond[i][1], &temp1, &temp2);
		if (isPossible == 0) break;

		for (int j = 1; j <= HC.size(); j++)
		{
			if ((temp1 == firstRad[j] && temp2 == secondRad[j]) || (temp2 == firstRad[j]
				&& temp1 == secondRad[j])) // check if an equivalent reaction is already saved
			{
				firstRad[j].isomeri++;	// the number of equivalent reaction is stored 
				// in the firsRad.isomeri just for convinience
				alreadyPresent = true;
				break;
			}
			int asd = 1;
		}
		if (!alreadyPresent)
		{
			numberOfReactions++;
			firstRad[numberOfReactions] = temp1;
			secondRad[numberOfReactions] = temp2;
		}
	}

	for (int i = 1; i <= numberOfReactions; i++)
	{
		Radicale firstRadicalType = firstRad[i].tipoR(firstRad[i].trova(2));
		Radicale secondRadicalType = secondRad[i].tipoR(secondRad[i].trova(2));

		reactionComment reacomm = k->v_initiation(firstRadicalType, secondRadicalType,
			firstRad[i].isomeri);
		reactions.push_back(Reaction(HC, std::vector<Molecola>{ firstRad[i], secondRad[i] },
			new double[3] { k->A, k->n, k->E }, "Initiation reactions", reacomm));
	}

	return reactions;
}

std::vector<Reaction> hAbstractionReactions(Molecola HC, Molecola absR, Kinox* k)
{
	std::vector<Reaction> reactions;

	if (HC.kindOfSPecies() != fuel_)
	{
		UTL::error("hAbstractionReactions called with HC that is not an hydrocarbon!");
		return std::vector<Reaction> {};
	}
	std::vector<Molecola> Rs;
	std::vector<int> cToRMap;
	generateR(HC, &Rs, &cToRMap);
	std::string nameHrad = k->nameHAbsRad(absR);
	std::string label = "H abstraction";
	Molecola HabsProd = absR;
	if (HabsProd.addH() == 0)
	{
		UTL::error("hAbstractionReactions called with absR that cannot host an hydrogen!");
		return std::vector<Reaction> {};
	}
	for (auto& R : Rs)
	{
		reactionComment reacomm = k->v_h_abstraction(absR, HC.tipoC(R.trova(2)),
			HC.numAbstractableH(R.trova(2)), R.isomeri);
		reactions.push_back(Reaction(std::vector<Molecola>{ HC, absR },
			std::vector<Molecola>{ R, HabsProd }, new double[3] { k->A, k->n, k->E },
			label, reacomm));
	}

	return reactions;
}

std::vector<Reaction> O2AdditionToRReactions(std::vector<Molecola> Rs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& R : Rs)
	{
		if (R.kindOfSPecies() != R_)
		{
			UTL::error("O2AdditionToRReactions called on a non radical species");
			continue;
		}
		Molecola ROO = R;
		ROO.Crad_to_ROO(ROO.trova(2));
		Radicale radicalType = R.tipoR(R.trova(2));
		reactionComment reacomm = k->v_o2_add_r(radicalType);
		reactions.push_back(Reaction(std::vector<Molecola>{ R, O2 }, ROO, new double[3]
			{ k->A, k->n, k->E }, "O2 addition to R", reacomm));
		//chemout.wrireaDetailed( 3, k.A, k.n, k.E, r[i], roo[i]);
	}
	return reactions;
}

std::vector<Reaction> O2EliminationFromROOReactions(std::vector<Molecola> ROOs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& ROO : ROOs)
	{
		if (ROO.kindOfSPecies() != ROO_)
		{
			UTL::error("O2EliminationFromROOReactions called on a non ROO species");
			continue;
		}
		Molecola R = ROO;
		R.COOrad_to_Crad(R.trova(3));
		Radicale radicalType = R.tipoR(R.trova(2));
		reactionComment reacomm = k->v_o2_rem_roo(radicalType);
		reactions.push_back(Reaction(ROO, std::vector<Molecola>{ R, O2 }, new double[3]
			{ k->A, k->n, k->E }, "O2 elimination from ROO", reacomm));
	}
	return reactions;
}

std::vector<Reaction> O2AdditionToQOOHReactions(std::vector<Molecola> QOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& QOOH : QOOHs)
	{
		if (QOOH.kindOfSPecies() != QOOH_)
		{
			UTL::error("O2AdditionToQOOHReactions called on a non QOOH species");
			continue;
		}
		Molecola OOQOOH = QOOH;
		OOQOOH.Crad_to_ROO(OOQOOH.trova(2));
		Radicale radicalType = QOOH.tipoR(QOOH.trova(2));
		reactionComment reacomm = k->v_o2_add_qooh(radicalType);
		reactions.push_back(Reaction(std::vector<Molecola>{ QOOH, O2 }, OOQOOH,
			new double[3] { k->A, k->n, k->E }, "O2 addition to QOOH", reacomm));
	}
	return reactions;
}

std::vector<Reaction> O2EliminationFromOOQOOHReactions(std::vector<Molecola> OOQOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola O2(2);
	for (auto& OOQOOH : OOQOOHs)
	{
		if (OOQOOH.kindOfSPecies() != OOQOOH_)
		{
			UTL::error("O2EliminationFromOOQOOHReactions called on a non OOQOOH species");
			continue;
		}
		Molecola QOOH = OOQOOH;
		QOOH.COOrad_to_Crad(QOOH.trova(3));
		Radicale radicalType = QOOH.tipoR(QOOH.trova(2));
		reactionComment reacomm = k->v_o2_rem_ooqooh(radicalType);
		reactions.push_back(Reaction(OOQOOH, std::vector<Molecola>{ QOOH, O2 },
			new double[3] { k->A, k->n, k->E }, "O2 elimination from OOQOOH", reacomm));
	}
	return reactions;
}

std::vector<Reaction> RIsomerizationReaction(std::vector<Molecola> Rs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& R_reac : Rs)
	{
		if (R_reac.kindOfSPecies() != R_)
		{
			UTL::error("RIsomerizationReaction called on a non R species");
			continue;
		}
		for (int dist = 3; dist < 6; dist++)  // look for 5,6 and 7 atom ring isomerizations
		{
			Anello ring = a5;
			switch (dist)
			{
			case 3:
				ring = a5;
				break;
			case 4:
				ring = a6;
				break;
			case 5:
				ring = a7;
				break;
			default:
				continue;
				break;
			}
			int pos_rad = R_reac.trova(2);
			int trovati[SIZEMAX + 1];
			//-----------------------------isomerization involving 5 atom ring
			R_reac.scorri(pos_rad, dist, trovati);		// flag in the vector "trovati" all the carbons that are at 3 atom distance from the radical
			for (int j = 1; j <= SIZEMAX; j++)
			{
				if (trovati[j] == 1 && R_reac.numAbstractableH(j) != 0)
				{
					Molecola R_prod = R_reac;
					R_prod.Crad_to_C(R_reac.trova(2));
					R_prod.C_to_Crad(j);
					if (R_prod == R_reac)
						continue;
					reactionComment reacomm = k->v_isomerization_r(R_reac.tipoR(pos_rad), R_reac.tipoH(j),
						ring, R_reac.numAbstractableH(j));
					reactions.push_back(Reaction(R_reac, R_prod, new double[3] { k->A, k->n, k->E },
						"Isomerization R", reacomm));
				}
			}

		}
	}
	return reactions;
}

std::vector<Reaction> ROOIsomerizationReactions(std::vector<Molecola> ROOs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& ROO : ROOs)
	{
		if (ROO.kindOfSPecies() != ROO_)
		{
			UTL::error("ROOIsomerizationReactions called on a non ROO species");
			continue;
		}
		int pos_o2 = ROO.trova(3);
		int trovati[SIZEMAX + 1];

		for (int dist = 1; dist <= 3; dist++)
		{
			Anello ringSize;
			switch (dist)
			{
			case 1:
				ringSize = a5;
				break;
			case 2:
				ringSize = a6;
				break;
			case 3:
				ringSize = a7;
				break;
			default:
				break;
			}
			ROO.scorri(pos_o2, dist, trovati);
			for (int j = 1; j <= SIZEMAX; j++)          // j iterates trough the found H
			{
				if (trovati[j] == 1 && ROO.numAbstractableH(j) != 0)
				{
					reactionComment reacomm = k->v_isom_roo(ROO.tipoROO(pos_o2), ROO.tipoH(j),
						ringSize, ROO.numAbstractableH(j));
					int pos_ooh = pos_o2;
					int pos_r = j;
					Molecola QOOH = ROO;
					QOOH.COOrad_to_COOH(pos_o2);
					QOOH.C_to_Crad(j);
					reactions.push_back(Reaction(ROO, QOOH, new double[3] { k->A, k->n, k->E },
						"Isomerization ROO", reacomm));
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> QOOHIsomerizationReactions(std::vector<Molecola> QOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& QOOH : QOOHs)
	{
		if (QOOH.kindOfSPecies() != QOOH_)
		{
			UTL::error("QOOHIsomerizationReactions called on a non QOOH species.");
			continue;
		}
		Anello ring = a5;
		int distance = QOOH.dist(QOOH.trova(2), QOOH.trova(4));
		switch (distance)
		{
		case 1:
			ring = a5;
			break;
		case 2:
			ring = a6;
			break;
		case 3:
			ring = a7;
			break;
		default:
			continue;
			break;
		}
		reactionComment reacomm = k->v_isom_qooh(QOOH.tipoROOH(QOOH.trova(4)),
			QOOH.tipoR(QOOH.trova(2)), ring);
		Molecola ROO = QOOH;
		ROO.Crad_to_C(ROO.trova(2));
		ROO.COOH_to_COOrad(ROO.trova(4));
		reactions.push_back(Reaction(QOOH, ROO, new double[3] { k->A, k->n, k->E },
			"Isomerization QOOH", reacomm));
	}
	return reactions;
}

std::vector<Reaction> OOQOOHIsomerizationReactions(std::vector<Molecola> OOQOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& OOQOOH : OOQOOHs)
	{
		if (OOQOOH.kindOfSPecies() != OOQOOH_)
		{
			UTL::error("OOQOOHIsomerizationReactions called on a non OOQOOH species.");
			continue;
		}
		int j;
		int pos_o2 = OOQOOH.trova(3);
		int trovati[SIZEMAX + 1];
		for (int distance = 1; distance < 4; distance++)
		{
			//std::cout << "  Distance = " << distance << std::endl;
			OOQOOH.scorri(pos_o2, distance, trovati);
			Anello ring = a5;
			switch (distance)
			{
			case 1:
				ring = a5;
				break;
			case 2:
				ring = a6;
				break;
			case 3:
				ring = a7;
				break;
			default:
				continue;
				break;
			}
			for (j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
			{
				if (trovati[j] == 1 && OOQOOH.numAbstractableH(j) != 0
					&& OOQOOH.tipo(j) != 4)
				{
					reactionComment reacomm = k->v_isom_ooqooh(OOQOOH.tipoROO(
						OOQOOH.trova(3)), OOQOOH.tipoH(j), ring,
						OOQOOH.numAbstractableH(j));
					Molecola POOH2 = OOQOOH;
					POOH2.OO_to_OOH(OOQOOH.trova(3));
					POOH2.C_to_Crad(j);

					reactions.push_back(Reaction(OOQOOH, POOH2,
						new double[3] { k->A, k->n, k->E }, "Isomerization OOQOOH",
						reacomm));
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> olefinsFromROOReactions(std::vector<Molecola> ROOs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& ROO : ROOs)
	{
		if (ROO.kindOfSPecies() != ROO_)
		{
			UTL::error("olefinsFromROOReactions called on a non ROO species.");
			continue;
		}
		int pos_o2 = ROO.trova(3);
		int trovati[SIZEMAX + 1];

		ROO.scorri(pos_o2, 1, trovati);	// find the carbon at distance 1 from the carbon 
		// with the OO
		for (int j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
		{
			if (trovati[j] == 1 && ROO.numAbstractableH(j) != 0)
			{
				int pos_r = j;

				if (ROO.tipoROO(pos_o2) == Rp	// Rate rule not available. However it should
					&& ROO.tipoC(pos_r) == Cp)  // happen only for ethane
					continue;

				Molecola OLE = ROO;
				OLE.removeOO(pos_o2);
				OLE.addole(pos_o2, j);
				Molecola HO2(5);

				reactionComment reacomm = k->v_ole_par_roo(ROO.tipoROO(pos_o2),
					ROO.tipoC(pos_r), ROO.numAbstractableH(pos_r));
				reactions.push_back(Reaction(ROO, std::vector<Molecola>{ OLE, HO2},
					new double[3] { k->A, k->n, k->E }, "ROO to olefins", reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> POOH2IsomerizationReactions(std::vector<Molecola> POOH2s, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& POOH2 : POOH2s)
	{
		if (POOH2.kindOfSPecies() != POOH2_)
		{
			UTL::error("POOH2IsomerizationReactions called on a non P(OOH)2 species.");
			continue;
		}
		std::vector<int> posOOH = POOH2.posOOHinPOOH2();
		int posR = POOH2.trova(2);
		for (auto& pos : posOOH)
		{
			int dist = POOH2.dist(posR, pos);
			Anello ring = a5;
			switch (dist)
			{
			case 1:
				ring = a5;
				break;
			case 2:
				ring = a6;
				break;
			case 3:
				ring = a7;
				break;
			default:
				continue;
				break;
			}
			Molecola OOQOOH = POOH2;
			OOQOOH.COOH_to_COOrad(pos);
			OOQOOH.Crad_to_C(posR);
			reactionComment reacomm = k->v_isom_pooh2(POOH2.tipoROOH(pos),
				POOH2.tipoR(posR), ring);
			reactions.push_back(Reaction(POOH2, OOQOOH, new double[3] { k->A, k->n, k->E },
				"Isomerization P(OOH)2", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> RBetaDecompositioReactions(std::vector<Molecola> Rs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& R : Rs)
	{
		if (R.kindOfSPecies() != R_)
		{
			UTL::error("RBetaDecompositioReactions called on a non R species.");
			continue;
		}
		int pos_rad = R.trova(2);
		int posalfa[SIZEMAX + 1];
		R.scorri(pos_rad, 1, posalfa);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (posalfa[j] == 1) // find the C in alpha
			{
				int posbeta[SIZEMAX + 1];
				R.scorri(j, 1, posbeta);  // find the C in beta
				for (int z = 1; z <= SIZEMAX; z++)
				{
					if (posbeta[z] == 1 && z != pos_rad)
					{
						Molecola m1, m2;
						int isPossible = R.spezza(j, z, &m1, &m2);
						if (isPossible == 0)
							continue;

						reactionComment reacomm = k->v_beta_dec_r(R.tipoR(pos_rad),
							m2.tipoR(m2.trova(2)));
						reactions.push_back(Reaction(R,
							std::vector<Molecola>{ m1, m2 },
							new double[3] { k->A, k->n, k->E },
							"Radicals beta decomposition", reacomm));
					}		// end pos beta
				}
			}
		}			// end pos alfa
	}
	return reactions;
}

std::vector<Reaction> olefinsFromRPlusO2Reactions(std::vector<Molecola> Rs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& R : Rs)
	{
		if (R.kindOfSPecies() != R_)
		{
			UTL::error("olefinsFromRPlusO2Reactions called on a non R species.");
			continue;
		}
		int pos_rad = R.trova(2);
		int posalfa[SIZEMAX + 1];
		R.scorri(pos_rad, 1, posalfa);
		Molecola HO2(5);
		Molecola O2(2);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (posalfa[j] == 1 && R.numAbstractableH(j) != 0)
			{
				Molecola OLE = R;
				OLE.Crad_to_C(OLE.trova(2));
				int isPossible = OLE.addole(pos_rad, j);
				if (isPossible == 0)
					continue;
				reactionComment reacomm = k->v_ole_par_r(R.numAbstractableH(j));
				reactions.push_back(Reaction(std::vector<Molecola>{ R, O2 },
					std::vector<Molecola>{ OLE, HO2 }, new double[3] { k->A, k->n, k->E },
					"Olefins from R + O2", reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> betaQOOHDecompositionReaction(std::vector<Molecola> QOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	for (auto& QOOH : QOOHs)
	{
		if (QOOH.kindOfSPecies() != QOOH_)
		{
			UTL::error("betaQOOHDecompositionReaction called on a non QOOH species.");
			continue;
		}
		int posR = QOOH.trova(2);
		int posOOH = QOOH.trova(4);
		int dist = QOOH.dist(posR, posOOH);
		if (dist == 1)
		{
			if (QOOH.tipoR(posR) == Rp
				&& QOOH.tipoROOH(posOOH) == Rp)
				continue;
			Molecola OLE = QOOH.parentFuel();
			OLE.addole(posR, posOOH);

			reactionComment reacomm = k->v_ole_from_beta_qooh(QOOH.tipoROOH(posOOH),
				QOOH.tipoR(posR));
			reactions.push_back(Reaction(QOOH,
				std::vector<Molecola>{ OLE, HO2}, new double[3] { k->A, k->n, k->E },
				"Decomposition beta-QOOH", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> gammaQOOHDecompositionReaction(std::vector<Molecola> QOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& QOOH : QOOHs)
	{
		if (QOOH.kindOfSPecies() != QOOH_)
		{
			UTL::error("gammaQOOHDecompositionReaction called on a non QOOH species.");
			continue;
		}
		int posR = QOOH.trova(2);
		int posOOH = QOOH.trova(4);
		int dist = QOOH.dist(posR, posOOH);
		if (dist == 2)
		{
			int pos_alfa;   // find position where to break bond
			{
				int alfa_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				QOOH.scorri(posOOH, 1, alfa_ooh);
				QOOH.scorri(posR, 1, alfa_rad);
				int j;
				for (j = 1; j <= SIZEMAX; j++)
					if (alfa_ooh[j] == 1 && alfa_rad[j] == 1)
						break;
				pos_alfa = j;
			}
			// break molecule
			Molecola m1, m2;
			int isPossible = QOOH.spezza(posOOH, pos_alfa, &m1, &m2);
			if (isPossible == 0)
				continue;
			int posRprod = m1.trova(2);
			m1.Crad_to_C(posRprod);
			m1.addcheto(posRprod);

			reactionComment reacomm = k->v_ole_from_gamma_qooh(QOOH.tipoROOH(posOOH),
				QOOH.tipoR(posR));
			reactions.push_back(Reaction(QOOH,
				std::vector<Molecola>{m1, m2, OH}, new double[3] { k->A, k->n, k->E },
				"Decomposition gamma-QOOH", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> deltaQOOHDecompositionReaction(std::vector<Molecola> QOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& QOOH : QOOHs)
	{
		if (QOOH.kindOfSPecies() != QOOH_)
		{
			UTL::error("gammaQOOHDecompositionReaction called on a non QOOH species.");
			continue;
		}
		int posR = QOOH.trova(2);
		int posOOH = QOOH.trova(4);
		int dist = QOOH.dist(posR, posOOH);
		if (dist == 3)
		{
			//                 find the position where to break
			int pos_break1;
			{
				int beta_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				QOOH.scorri(posOOH, 2, beta_ooh);
				QOOH.scorri(posR, 1, alfa_rad);
				int j;
				for (j = 1; j <= SIZEMAX; j++)
					if (beta_ooh[j] == 1 && alfa_rad[j] == 1) break;
				pos_break1 = j;
			};
			int pos_break2;
			{
				int alfa_ooh[SIZEMAX + 1];
				int beta_rad[SIZEMAX + 1];
				QOOH.scorri(posOOH, 1, alfa_ooh);
				QOOH.scorri(posR, 2, beta_rad);
				int j;
				for (j = 1; j <= SIZEMAX; j++)
					if (alfa_ooh[j] == 1 && beta_rad[j] == 1) break;
				pos_break2 = j;
			};
			//           generate decomposition products
			Molecola m1, m2;
			int isPossible = QOOH.spezza(pos_break1, pos_break2, &m1, &m2);
			if (isPossible == 0)
				continue;

			reactionComment reacomm = k->v_ole_from_delta_qooh(QOOH.tipoROOH(posOOH),
				QOOH.tipoR(posR));
			reactions.push_back(Reaction(QOOH, std::vector<Molecola>{m1, m2},
				new double[3] { k->A, k->n, k->E },
				"Decomposition delta-QOOH", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> QOOHToCEthReactions(std::vector<Molecola> QOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& QOOH : QOOHs)
	{
		if (QOOH.kindOfSPecies() != QOOH_)
		{
			UTL::error("QOOHToCEthReactions called on a non QOOH species.");
			continue;
		}
		std::string correction = "none";
		int posR = QOOH.trova(2);
		int posOOH = QOOH.trova(4);
		int dist = QOOH.dist(posR, posOOH);
		AnelloO ring;
		if (dist == 1)
			ring = ao3;
		else if (dist == 2)
		{
			ring = ao4;
			int alpha_ooh[SIZEMAX + 1];
			int alpha_rad[SIZEMAX + 1];
			QOOH.scorri(posOOH, 1, alpha_ooh);
			QOOH.scorri(posR, 1, alpha_rad);
			for (int i = 1; i <= SIZEMAX; i++)
			{
				if (alpha_ooh[i] == 1 && alpha_rad[i] == 1)
				{
					if (QOOH.tipoC(i) == Ct)
						correction = "T";
					if (QOOH.tipoC(i) == Cq)
						correction = "Q";
				}
			}
		}
		else if (dist == 3)
		{
			ring = ao5;
			//int alpha_ooh[SIZEMAX + 1];
			//int alpha_rad[SIZEMAX + 1];
			//int beta_ooh[SIZEMAX + 1];
			//int beta_rad[SIZEMAX + 1];
			//QOOH.scorri(posOOH, 1, alpha_ooh);
			//QOOH.scorri(posR,   1, alpha_rad);
			//QOOH.scorri(posOOH, 2, beta_ooh);
			//QOOH.scorri(posR,   2, beta_rad);
			//
			//int pos1 = -1;
			//int pos2 = -1;
			//for (int i = 1; i <= SIZEMAX; i++)
			//{
			//	if (alpha_ooh[i] == 1 && beta_rad[i] == 1)
			//		pos1 = i;
			//	if (beta_ooh[i] == 1 && alpha_rad[i] == 1)
			//		pos2 = i;
			//}
			//if (pos1 != -1 && pos2 != -1)
			//{
			//	if (QOOH.tipoC(pos1) == Cq && QOOH.tipoC(pos2) == Cs)
			//		correction = "QS";
			//	if (QOOH.tipoC(pos1) == Cs && QOOH.tipoC(pos2) == Cq)
			//		correction = "SQ";
			//}
		}
		else
			continue;
		if (QOOH.tipoROOH(posOOH) == Rp && QOOH.tipoR(posR) == Rp && ring == ao3)
			continue;
		Molecola cEth = QOOH.parentFuel();
		cEth.addetero(posR, posOOH);

		reactionComment reacomm = k->v_ether_from_qooh(QOOH.tipoROOH(posOOH),
			QOOH.tipoR(posR), ring, correction);
		//// DEBUG
		//bool toMultiply = false;
		//for (int i = 1; i <= QOOH.size(); i++)
		//{
		//	bool toCheck = false;
		//	if (i == posR || i == posOOH)
		//		toCheck = false;
		//	else if (QOOH.dist(i, posR) + QOOH.dist(i, posOOH) == dist)
		//		toCheck = true;
		//	if (toCheck)
		//		if (QOOH.tipoC(i) == Ct || QOOH.tipoC(i) == Cq)
		//			toMultiply = true;
		//}
		//if(toMultiply)
		//	(*k).A *= 10;
		//// DEBUG
		reactions.push_back(Reaction(QOOH, std::vector<Molecola>{ cEth, OH },
			new double[3] { k->A, k->n, k->E }, "Cyclic ethers from QOOH", reacomm));
	}
	return reactions;
}

std::vector<Reaction> POOH2ToCEthOOHReactions(std::vector<Molecola> POOH2s, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& POOH2 : POOH2s)
	{
		if (POOH2.kindOfSPecies() != POOH2_)
		{
			UTL::error("POOH2ToCEthOOHReactions called on a non P(OOH)2 species.");
			continue;
		}

		std::vector<int> posOOHs = POOH2.posOOHinPOOH2();
		int posR = POOH2.trova(2);
		for (auto& posOOH : posOOHs)
		{
			int dist = POOH2.dist(posR, posOOH);
			AnelloO ring = ao3;
			switch (dist)
			{
			case 1:
				ring = ao3;
				break;
			case 2:
				ring = ao4;
				break;
			case 3:
				ring = ao5;
				break;
			default:
				continue;
				break;
			}
			if (POOH2.tipoROOH(posOOH) == Rp && POOH2.tipoR(posR) == Rp && ring == ao3)
				continue;
			Molecola cEthOOH = POOH2;
			cEthOOH.removeOOH(posOOH);
			cEthOOH.Crad_to_C(posR);
			cEthOOH.addetero(posOOH, posR);
			reactionComment reacomm = k->v_ether_from_pooh2(POOH2.tipoROOH(posOOH),
				POOH2.tipoR(posR), ring);
			reactions.push_back(Reaction(POOH2, std::vector<Molecola>{ cEthOOH, OH },
				new double[3] { k->A, k->n, k->E }, "Ethers-OOH from P(OOH)2", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> betaPOOH2DecompositionReaction(std::vector<Molecola> POOH2s, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	for (auto& POOH2 : POOH2s)
	{
		if (POOH2.kindOfSPecies() != POOH2_)
		{
			UTL::error("betaPOOH2DecompositionReaction called on a non P(OOH)2 species.");
			continue;
		}
		std::vector<int> posOOHs = POOH2.posOOHinPOOH2();
		int posR = POOH2.trova(2);
		for (auto& posOOH : posOOHs)
		{
			if (POOH2.dist(posOOH, posR) == 1)
			{
				if (POOH2.tipoROOH(posOOH) == Rp && POOH2.tipoR(posR) == Rp)
					continue;
				Molecola OLEOOH = POOH2;
				OLEOOH.removeOOH(posOOH);
				OLEOOH.Crad_to_C(posR);
				OLEOOH.addole(posOOH, posR);
				reactionComment reacomm = k->v_pooh2_dec_1(POOH2.tipoROOH(posOOH),
					POOH2.tipoR(posR));
				reactions.push_back(Reaction(POOH2, std::vector<Molecola>{ OLEOOH, HO2 },
					new double[3] { k->A, k->n, k->E }, "Beta-P(OOH)2 decomposition",
					reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> gammaPOOH2DecompositionReaction(std::vector<Molecola> POOH2s, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& POOH2 : POOH2s)
	{
		if (POOH2.kindOfSPecies() != POOH2_)
		{
			UTL::error("gammaPOOH2DecompositionReaction called on a non P(OOH)2 species.");
			continue;
		}
		std::vector<int> posOOHs = POOH2.posOOHinPOOH2();
		int posR = POOH2.trova(2);
		for (auto& posOOH : posOOHs)
		{
			if (POOH2.dist(posOOH, posR) == 2)
			{
				//                  find the position where to break
				int pos_alfa;
				{
					int alfa_ooh[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					POOH2.scorri(posOOH, 1, alfa_ooh);
					POOH2.scorri(posR, 1, alfa_rad);
					int j;
					for (j = 1; j <= SIZEMAX; j++)
						if (alfa_ooh[j] == 1 && alfa_rad[j] == 1) break;
					pos_alfa = j;
				}
				//                             generate the broken molecules
				if (POOH2.tipo(pos_alfa) == 4)
					continue;
				Molecola m1, m2;
				int isPossible = POOH2.spezza(posOOH, pos_alfa, &m1, &m2);
				if (isPossible == 0) continue;

				int pos = m1.trova(2);
				m1.tipo(pos, 1);
				m1.addcheto(pos);
				reactionComment reacomm = k->v_pooh2_dec_2(POOH2.tipoROOH(posOOH),
					POOH2.tipoR(posR));
				reactions.push_back(Reaction(POOH2, std::vector<Molecola>{ m1, m2, OH },
					new double[3] { k->A, k->n, k->E }, "Gamma-P(OOH)2 decomposition",
					reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> deltaPOOH2DecompositionReaction(std::vector<Molecola> POOH2s, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& POOH2 : POOH2s)
	{
		if (POOH2.kindOfSPecies() != POOH2_)
		{
			UTL::error("deltaPOOH2DecompositionReaction called on a non P(OOH)2 species.");
			continue;
		}
		std::vector<int> posOOHs = POOH2.posOOHinPOOH2();
		int posR = POOH2.trova(2);
		for (auto& posOOH : posOOHs)
		{

			int pos_rad = -1;		// find the position where the radical will be located
			{
				int beta_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				POOH2.scorri(posOOH, 1, beta_ooh);
				POOH2.scorri(posR, 2, alfa_rad);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					if (beta_ooh[j] == 1 && alfa_rad[j] == 1)
					{
						pos_rad = j;
						break;
					}
				}
			}
			if (pos_rad == -1)
				continue;

			if (POOH2.tipo(pos_rad) == 4)  // if there is a COOH in this position continue
				continue;

			//                  find the position where to break
			int pos_break1 = -1;
			{
				int beta_ooh[SIZEMAX + 1];
				int alfa_rad[SIZEMAX + 1];
				POOH2.scorri(posOOH, 2, beta_ooh);
				POOH2.scorri(posR, 1, alfa_rad);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					if (beta_ooh[j] == 1 && alfa_rad[j] == 1)
					{
						pos_break1 = j;
						break;
					}
				}
			}
			if (pos_break1 == -1)
				continue;

			int pos_break2 = -1;
			{
				int alfa_ooh[SIZEMAX + 1];
				int beta_rad[SIZEMAX + 1];
				POOH2.scorri(posOOH, 1, alfa_ooh);
				POOH2.scorri(posR, 2, beta_rad);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					if (alfa_ooh[j] == 1 && beta_rad[j] == 1)
					{
						pos_break2 = j;
						break;
					}
				}
			}
			if (pos_break2 == -1)
				continue;

			if (POOH2.tipo(pos_break1) == 4)
				continue;

			//                     break the molecule
			if (pos_break1 == pos_break2)
				continue;
			if (POOH2.areBonded(pos_break1, pos_break2) == 0)
				continue;

			Molecola m1, m2;
			int isPossible = POOH2.spezza(pos_break1, pos_break2, &m1, &m2);
			if (isPossible == 0)
				continue;

			reactionComment reacomm = k->v_pooh2_dec_3(POOH2.tipoROOH(posOOH),
				POOH2.tipoR(posR));
			Reaction reac(POOH2, std::vector<Molecola>{ m1, m2},
				new double[3] { k->A, k->n, k->E }, "Delta P(OOH)2 decomposition", reacomm);
			UTL::addUnique(&reactions, reac);
		}
	}
	return reactions;
}

std::vector<Reaction> cyclicEtherOOHDecompositionReactions(std::vector<Molecola> cEthOOHs,
	Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& cEthOOH : cEthOOHs)
	{
		if (cEthOOH.kindOfSPecies() != cEthOOH_)
		{
			UTL::error("cyclicEtherOOHDecompositionReactions called on a non cEthOOH species.");
			continue;
		}
		int alpha_ooh[SIZEMAX + 1];	// array in which there is 1 if the atom in that 
		// position is in alpha to OOH
		cEthOOH.scorri(cEthOOH.trova(4), 1, alpha_ooh);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (alpha_ooh[j] == 1)
			{
				Molecola m1;
				Molecola m2;
				int isPossible = cEthOOH.spezza(cEthOOH.trova(4), j, &m1, &m2);
				if (isPossible == 0)
					continue;
				//std::cout << m1 << std::endl << m2 << std::endl << std::endl;
				reactionComment reacomm = k->v_etherooh_dec(cEthOOH.tipoROOH(cEthOOH.trova(4)));

				std::vector<Molecola> prods;
				//std::cout << etherOOH[i] << std::endl;
				if (m2.size() == 0)    // if the ring was opened
				{
					m1.Crad_to_keto(cEthOOH.trova(4));
					prods = decomposeLinEthRO(m1);
				}
				else if (m2.kindOfSPecies() == cEthR_)  //do not consider reactions that lead
					//to cyclic ethers CO
				{
					m1.Crad_to_keto(m1.trova(2));
					prods = decomposeCEthR(m2);
					prods.push_back(m1);
				}
				else
				{
					continue;
				}

				prods = fullyDecomposeRO(prods);

				if (prods.size() < 2)
				{
					UTL::error("In cyclicEtherOOHDecompositionReactions less than 2 decomposition products have been found");
					continue;
				}

				prods.push_back(OH);

				reactions.push_back(Reaction(cEthOOH, prods, new double[3] { k->A, k->n, k->E },
					"Ether-OOH decomposition", reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> OOQOOHToOLEOOHReactions(std::vector<Molecola>OOQOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola HO2(5);
	for (auto& OOQOOH : OOQOOHs)
	{
		if (OOQOOH.kindOfSPecies() != OOQOOH_)
		{
			UTL::error("OOQOOHToOLEOOHReactions called on a non OOQOOH species.");
			continue;
		}
		int posOO = OOQOOH.trova(3);
		int trovati[SIZEMAX + 1];

		OOQOOH.scorri(posOO, 1, trovati);	// find the carbon at distance 1 from the carbon with the OO
		for (int j = 1; j <= SIZEMAX; j++)      // j iterates trough the found H
		{
			if (trovati[j] == 1 && OOQOOH.numAbstractableH(j) != 0
				&& OOQOOH.trova(4) != j)
			{
				if (OOQOOH.tipoROO(posOO) == Rp && OOQOOH.tipoC(j) == Cp)
					continue;
				Molecola ole = OOQOOH;
				ole.removeOO(posOO);
				ole.addole(posOO, j);
				reactionComment reacomm = k->v_ole_par_ooqooh(OOQOOH.tipoROO(posOO),
					OOQOOH.tipoC(j), OOQOOH.numAbstractableH(j));
				reactions.push_back(Reaction(OOQOOH, std::vector<Molecola>{ ole, HO2 },
					new double[3] { k->A, k->n, k->E },
					"Olefins-OOH from OOQOOH", reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> OLEOOHDecompositionReactions(std::vector<Molecola> OLEOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& OLEOOH : OLEOOHs)
	{
		if (OLEOOH.kindOfSPecies() != oleOOH_)
		{
			UTL::error("OLEOOHDecompositionReactions called on a non oleOOH species.");
			continue;
		}
		int alfa_ooh[SIZEMAX + 1];					// array in which there is 1 if the atom in that position is in alfta to OOH
		OLEOOH.scorri(OLEOOH.trova(4), 1, alfa_ooh);
		for (int j = 1; j <= SIZEMAX; j++)
		{
			if (alfa_ooh[j] == 1)
			{
				Molecola m1;
				Molecola m2;
				int isPossible = OLEOOH.spezza(OLEOOH.trova(4), j, &m1, &m2);
				if (isPossible == 0)
					continue;
				m1.Crad_to_keto(m1.trova(2));
				reactionComment reacomm = k->v_oleooh_dec(OLEOOH.tipoROOH(OLEOOH.trova(4)));
				reactions.push_back(Reaction(OLEOOH, std::vector<Molecola>{ m1, m2, OH },
					new double[3] { k->A, k->n, k->E }, "Ole-OOH decomposition", reacomm));
			}
		}
	}
	return reactions;
}

std::vector<Reaction> KHPFormationReactions(std::vector<Molecola> OOQOOHs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& OOQOOH : OOQOOHs)
	{
		if (OOQOOH.kindOfSPecies() != OOQOOH_)
		{
			UTL::error("KHPFormationReactions called on a non OOQOOH species.");
			continue;
		}
		int posOO = OOQOOH.trova(3);		// find the position of the oo group
		int posOOH = OOQOOH.trova(4);		// find the position of the ooh group

		int dist = OOQOOH.dist(posOO, posOOH);	// find the distance between the oo 
		//and ooh group
// from the distance find the type of ring formed during the reaction
		Anello ring = a5;
		switch (dist)
		{
		case 1:
			ring = a5;
			break;
		case 2:
			ring = a6;
			break;
		case 3:
			ring = a7;
			break;
		default:		// if the ring is bigger than 7 atoms skip 
			continue;	// since the reaction cannot not happen
			break;
		}

		Molecola KHP = OOQOOH;
		KHP.COOrad_to_COOH(posOO);		// change the oo group in an ooh group
		KHP.tipo(posOOH, 1);					// remove the ooh group and ...
		int isPossible = KHP.addcheto(posOOH);	// ... replace it with the group =o 
		// Check if this step is possible.
		if (isPossible != 0)					// if it is possible save reaction
		{
			reactionComment reacomm = k->v_ooqooh_to_khp(OOQOOH.tipoROO(posOO),
				OOQOOH.tipoROOH(posOOH), ring, OOQOOH.numAbstractableH(posOOH));
			reactions.push_back(Reaction(OOQOOH, std::vector<Molecola>{ KHP, OH },
				new double[3] { k->A, k->n, k->E },
				"OOQOOH conversion to ketohydroperoxide", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> KHPDecompositionReactions(std::vector<Molecola> KHPs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	for (auto& KHP : KHPs)
	{
		if (KHP.kindOfSPecies() != KHP_)
		{
			UTL::error("KHPDecompositionReactions called on a non KHP species.");
			continue;
		}
		int posOOH = KHP.trova(4);		// find the position of the ooh group
		int posO = KHP.trova(7);		// find the position of the =o group
		int dist = KHP.dist(posO, posOOH);	// find the distance between the ooh and =o group
		// from the distance find the type of ring formed during the reaction
		Anello ring = a5;
		switch (dist)
		{
		case 1:
			ring = a5;
			break;
		case 2:
			ring = a6;
			break;
		case 3:
			ring = a7;
			break;
		default:		// if the ring is bigger than 7 atoms skip
			continue;
			break;
		}

		int pos_alpha;
		{
			int alpha_ooh[SIZEMAX + 1];
			int alpha_rad[SIZEMAX + 1];
			KHP.scorri(posO, dist - 1, alpha_ooh);
			KHP.scorri(posOOH, 1, alpha_rad);
			for (int j = 1; j <= SIZEMAX; j++)
			{
				if (alpha_ooh[j] == 1 && alpha_rad[j] == 1)
				{
					pos_alpha = j;
					break;
				}
			}
		}

		Molecola m1, m2;

		int isPossible = KHP.spezza(pos_alpha, posOOH, &m1, &m2);
		if (isPossible == 0)
			continue;
		int pos_rad_m2 = m2.trova(2);	// find the radical on the molecule m2
		m2.tipo(pos_rad_m2, 1);			// remove the radical and..
		m2.addcheto(pos_rad_m2);		// ... replace it with a cheto group

		std::vector<Molecola> products = { m2 };

		if (m1.kindOfSPecies() == RO_)
		{
			std::vector<Molecola> RO_dec_prod = fullyDecomposeRO(m1);
			for (auto& pr : RO_dec_prod)
				products.push_back(pr);
		}
		else
			products.push_back(m1);

		products.push_back(OH);

		reactionComment reacomm = k->v_khp_decomp(KHP.tipoROOH(posOOH), dist);

		reactions.push_back(Reaction(KHP, products, new double[3] { k->A, k->n, k->E },
			"Ketohydroperoxides decomposition", reacomm));
	}
	return reactions;
}

std::vector<Reaction> cyclicEthersDecompositionReactions(std::vector<Molecola> cEths,
	Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola H2O(6);
	Molecola OH(4);
	for (auto& cEth : cEths)
	{
		if (cEth.kindOfSPecies() != cEth_)
		{
			UTL::error("cyclicEthersDecompositionReactions called on a non cyclic ehter species.");
			continue;
		}
		std::vector<int> etherPos = cEth.posEthero();
		int dist = cEth.dist(etherPos[0], etherPos[1]);
		for (int j = 0; j < 2; j++)			// iterates through the two carbon 
			// the oxygen is bonded with
		{
			int posRad = etherPos[j];
			int posCheto = etherPos[1 - j];
			int numH = cEth.numH(posCheto);
			int numH2 = cEth.numH(posRad);

			if (numH == 0 && numH2 == 0) // if there are no abstractable 
				// hydrogens in alpha to the oxygen use the secondary
			{							 // path
				// look for the hydrogens to abstract
				// if there are methyls abstract it from them
				Molecola interm = cEth;
				int alfa_rad[SIZEMAX + 1];
				interm.scorri(posRad, 1, alfa_rad);
				bool normalPathTaken = false;
				int numAbsH = 0;
				int posAbs = 0;
				Carbonio carbType = Cp;
				for (int l = 1; l <= SIZEMAX; l++)
				{
					if (alfa_rad[l] == 1 && !interm.isGroup(l) && interm.tipoC(l) == Cp)
					{
						numAbsH += 3;
						posAbs = l;			// no matter if it is overwritten 
						// we need only one position
					}
				}
				if (posAbs == 0)	// if no methyl was found look for a secondary 
				{					// carbon in alpha
					carbType = Cs;
					for (int l = 1; l <= SIZEMAX; l++)
					{
						if (alfa_rad[l] == 1 && !interm.isGroup(l) && interm.tipoC(l) == Cs)
						{
							numAbsH += 2;
							posAbs = l;			// no matter if it is overwritten 
							// we need only one position
						}
					}
				}

				interm.removeCycEther(posRad, posCheto);
				interm.addole(posAbs, posRad);
				//interm.addcheto(posCheto);

				// find the carbon in alpha to the cheto more close to the double bond
				int alfa_cheto[SIZEMAX + 1];
				interm.scorri(posCheto, 1, alfa_cheto);
				int posAlfaCheto = 0;
				int distAlfaCheto = SIZEMAX;
				for (int l = 1; l <= SIZEMAX; l++)
				{
					if (alfa_cheto[l] == 1)
					{
						if (interm.dist(posAbs, l) < distAlfaCheto)
						{
							distAlfaCheto = interm.dist(posAbs, l);
							posAlfaCheto = l;
						}
					}
				}
				Molecola m1, m2;
				interm.spezza(posCheto, posAlfaCheto, &m1, &m2);
				int posRadM1 = m1.posCrad();
				m1.Crad_to_C(posRadM1);
				m1.addcheto(posRadM1);

				reactionComment reacomm = k->v_h_abstraction(OH, carbType, numAbsH, 1);
				reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
					std::vector<Molecola>{ m1, m2, H2O }, new double[3] { k->A, k->n, k->E },
					"Cyclic ethers decomposition", reacomm));
			}
			else if (numH == 0 && numH2 != 0)
			{
				continue;
			}
			else
			{
				Molecola interm = cEth;
				interm.removeCycEther(posRad, posCheto);
				interm.addcheto(posCheto);
				interm.C_to_Crad(posRad);
				bool reactionHappened = false;
				Molecola m1, m2;
				switch (dist)
				{
				case 1:				// 3 member ring
				{
					// search if there is a carbon in alfa at the cheto and in beta 
					// at the radical
					int alfa_cheto[SIZEMAX + 1];
					int beta_rad[SIZEMAX + 1];
					interm.scorri(posCheto, 1, alfa_cheto);
					interm.scorri(posRad, 2, beta_rad);
					bool normalPathTaken = false;
					for (int l = 1; l <= SIZEMAX; l++)	// normal path: 
					{
						if (alfa_cheto[l] == 1 && beta_rad[l] == 1)
						{
							interm.spezza(posCheto, l, &m1, &m2);
							normalPathTaken = true;
							reactionHappened = true;
							break;
						}
					}
					if (!normalPathTaken)
					{
						for (int l = 1; l <= SIZEMAX; l++)
						{
							if (beta_rad[l] == 1)
							{
								for (int h = 1; h <= interm.size(); h++)
								{
									if (interm.dist(h, l) == 1 && interm.dist(h, posRad) == 1)
									{
										Molecola d1, d2;
										interm.spezza(h, l, &d1, &d2);
										m1 = d1;
										m2 = d2;
										reactionHappened = true;
										break;
									}
									if (reactionHappened)
										break;
								}
							}
						}
					}
				}
				break;
				case 2:				// 4 member ring
				{
					int alfa_cheto[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					interm.scorri(posCheto, 1, alfa_cheto);
					interm.scorri(posRad, 1, alfa_rad);
					for (int l = 1; l <= SIZEMAX; l++)
					{
						if (alfa_cheto[l] == 1 && alfa_rad[l] == 1)
						{
							interm.spezza(posCheto, l, &m1, &m2);
							reactionHappened = true;
							break;
						}
					}
				}
				break;
				case 3:				// 5 member ring
				{
					int alfa_cheto[SIZEMAX + 1];
					int beta_cheto[SIZEMAX + 1];
					int alfa_rad[SIZEMAX + 1];
					int beta_rad[SIZEMAX + 1];
					interm.scorri(posCheto, 1, alfa_cheto);
					interm.scorri(posCheto, 2, beta_cheto);
					interm.scorri(posRad, 1, alfa_rad);
					interm.scorri(posRad, 2, beta_rad);
					for (int l = 1; l <= SIZEMAX; l++)
					{
						if (alfa_cheto[l] == 1 && beta_rad[l] == 1)
						{
							for (int h = 1; h <= SIZEMAX; h++)
							{
								if (beta_cheto[h] == 1 && alfa_rad[h] == 1)
								{
									interm.spezza(h, l, &m1, &m2);
									reactionHappened = true;
									break;
								}
							}
						}
					}
				}
				break;
				default:
					std::cout << "Molecule " << cEth
						<< " has the carbon bonding with the oxygen at a distance of "
						<< dist << ", this is not a valid cyclic ether!" << std::endl;
					break;
				}

				if (reactionHappened)
				{
					std::vector<Molecola> RO_dec_prod;
					bool ism1RO = false;
					bool ism2RO = false;
					if (m1.kindOfSPecies() == RO_)
					{
						RO_dec_prod = decomposeRO(m1);
						ism1RO = true;
					}

					if (m2.kindOfSPecies() == RO_)
					{
						RO_dec_prod = decomposeRO(m2);
						ism2RO = true;
					}

					if (m1.kindOfSPecies() == RO_ && m2.kindOfSPecies() == RO_)
					{
						std::cerr << "ERROR: in decomposition of cyclic ethers,"
							<< " two RO are formed, this should not happen!" << std::endl;
						std::cerr << "Press any key to continue..." << std::endl;
						std::string input_keyboard;
						std::cin >> input_keyboard;
					}

					reactionComment reacomm = k->v_cyc_eth_dec();

					if (ism1RO == true && RO_dec_prod.size() == 2)
					{
						reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
							std::vector<Molecola>{ RO_dec_prod[0], RO_dec_prod[1], m2, H2O },
							new double[3] { k->A, k->n, k->E },
							"Cyclic ethers decomposition", reacomm));
					}
					else if (ism2RO == true && RO_dec_prod.size() == 2)
					{

						reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
							std::vector<Molecola>{ m1, RO_dec_prod[0], RO_dec_prod[1], H2O },
							new double[3] { k->A, k->n, k->E },
							"Cyclic ethers decomposition", reacomm));
					}
					else
					{
						reactions.push_back(Reaction(std::vector<Molecola>{ cEth, OH },
							std::vector<Molecola>{ m1, m2, H2O },
							new double[3] { k->A, k->n, k->E },
							"Cyclic ethers decomposition", reacomm));
					}
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> allylicRadicalsFormationReactions(std::vector<Molecola> OLEs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	Molecola H2O(6);
	for (auto& OLE : OLEs)
	{
		if (OLE.kindOfSPecies() != OLE_)
		{
			UTL::error("allylicRadicalsFormationReactions called on a non olefin species.");
			continue;
		}
		for (int j = 1; j <= OLE.numberOle(); j++)	// iterate through all the double bond
		{											// in the molecule
			std::vector<int> oleCarbons = OLE.posOle(j);
			for (int l = 0; l < 2; l++)		// iterate trhough all (two) carbons 
			{								// involved in the double bond
				int posalfa[SIZEMAX + 1];
				OLE.scorri(oleCarbons[l], 1, posalfa);
				// iterate through all the carbon in alfa at the carbon involved in 
				// the double bond
				for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)
				{
					if (OLE.tipo(h) != 1)
						continue;		// the H can be abstracted only on normal carbons
					if (OLE.isole(h))
						continue;
					int numH = OLE.numH(h);
					if (numH == 0)
						continue; // if there are not hydrogen to abstract skip the reaction

					// save the type of substitution degree of the two carbons 
					// (the carbon we abstract the H from and the carbon invovled in the 
					// double bond in beta position to it) in order to decide if the radical
					// isomerizes or not
					Carbonio alfaC = OLE.tipoC(h);
					Carbonio otherOleC;
					otherOleC = OLE.tipoC(oleCarbons[1 - l]);

					int posRad, pos1Ole, pos2Ole;

					Molecola product = OLE;
					if (carbonioToInt(alfaC) < carbonioToInt(otherOleC))
					{		// radical is going to isomerize
						product.removeOle(oleCarbons[l], oleCarbons[1 - l]);
						product.addole(h, oleCarbons[l]);
						product.C_to_Crad(oleCarbons[1 - l]);
					}
					else	// radical is not going to isomerize
					{
						product.C_to_Crad(h);
					}

					reactionComment reacomm = k->v_allylic_rad_form(alfaC, numH);
					reactions.push_back(Reaction(std::vector<Molecola>{ OLE, OH },
						std::vector<Molecola>{ product, H2O },
						new double[3] { k->A, k->n, k->E },
						"Allylic radicals formation", reacomm));
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> alkenylROFormationReactions(std::vector<Molecola> AllRs, Kinox* k)
{
	std::vector<Reaction> reactions;
	Molecola OH(4);
	Molecola HO2(5);
	for (auto& AllR : AllRs)
	{
		if (AllR.kindOfSPecies() != oleR_)
		{
			UTL::error("alkenylROFormationReactions called on a non oleR species.");
			continue;
		}
		Molecola product = AllR;
		product.Crad_to_COrad(product.trova(2));
		reactionComment reacomm = k->v_alkenyl_ro_form();
		reactions.push_back(Reaction(std::vector<Molecola>{ AllR, HO2 },
			std::vector<Molecola>{ product, OH }, new double[3] { k->A, k->n, k->E },
			"Alkenyl RO formation", reacomm));
	}
	return reactions;
}

/*
std::vector<Reaction> alkenylRODecompositionReactions(std::vector<Molecola> AlkROs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& AlkRO : AlkROs)
	{
		if (AlkRO.kindOfSPecies() != alkRO_)
		{
			UTL::error("alkenylRODecompositionReactions called on a non alkenyl RO species.");
			continue;
		}
		int posRO = AlkRO.trova(9);
		int posalfa[SIZEMAX + 1];
		AlkRO.scorri(posRO, 1, posalfa);		// find the carbons in alpha
		// iterate through all the carbon in alfa at the CO*
		for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)
		{
			if (AlkRO.isole(h))
				continue;	// if it is the carbon involved in the double bond continue
			Carbonio tipoC = AlkRO.tipoC(h);
			Molecola product = AlkRO;
			Molecola m1, m2;
			int isPossible = AlkRO.spezza(posRO, h, &m1, &m2);
			if (isPossible == 0)
			{
				std::cerr << "WARNING: alkenyl RO decomposition failed in breaking the bond,"
					<< "reaction skipped!" << std::endl;
				continue;
			}

			int posRadM1 = m1.trova(2);
			//m1.Crad_to_C(posRadM1);
			//m1.removeAllKeto();
			if (posRadM1 != 0)
				m1.addcheto(posRadM1);
			reactionComment reacomm = k->v_alkenyl_ro_dec(tipoC);
			reactions.push_back(Reaction(AlkRO, std::vector<Molecola>{ m1, m2 },
				new double[3] { k->A, k->n, k->E }, "Alkenyl RO decomposition", reacomm));
		}
	}
	return reactions;
}
*/

std::vector<Reaction> alkenylRODecompositionReactions(std::vector<Molecola> AlkROs, Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& AlkRO : AlkROs)
	{
		if (AlkRO.kindOfSPecies() != alkRO_)
		{
			UTL::error("alkenylRODecompositionReactions called on a non alkenyl RO species.");
			continue;
		}
		int posRO = AlkRO.trova(9);
		int posalfa[SIZEMAX + 1];
		AlkRO.scorri(posRO, 1, posalfa);		// find the carbons in alpha
		// iterate through all the carbon in alfa at the CO*
		std::vector<int> alphaCs;
		for (int h = 1; h <= SIZEMAX; h++)
			if (posalfa[h] == 1)
				alphaCs.push_back(h);
		if (alphaCs.size() == 0)
		{
			//std::cout << "DEBUG: Not found alpha C decomposing AlkRO " << AlkRO << std::endl;
			continue;
		}
		std::vector<int> reactingCs;
		for (auto& ind : alphaCs)
			if (AlkRO.isole(ind) == false)
				reactingCs.push_back(ind);
		if (reactingCs.size() == 0)
			reactingCs.push_back(alphaCs[0]);
		for(auto &h : reactingCs)
		{
			// DEBUG ->
			if (AlkRO.isole(h))
			{
				//std::cout << "DEBUG: The decomposition of this AlkRO would have been skipped " << AlkRO << std::endl;
				int posalfa2[SIZEMAX + 1];
				AlkRO.scorri(h, 1, posalfa2);		// find the carbons in alpha
				std::vector<int> alphaCs2;
				for (int g = 1; g <= SIZEMAX; g++)
					if (posalfa2[g] == 1 && g != posRO)
					{
						alphaCs2.push_back(g);
					}
				for (auto& posA : alphaCs2)
				{
					int posBeta[SIZEMAX + 1];
					AlkRO.scorri(posA, 1, posBeta);		// find the carbons in alpha
					for (int g = 1; g <= SIZEMAX; g++)
					{
						if (posBeta[g] == 1 && g!=h)
						{
							Molecola temp = AlkRO;
							temp.COrad_to_CO(posRO);
							temp.removeOle(posRO, h);
							temp.C_to_Crad(h);
							//std::cout << "asd1 " << temp << std::endl;
							Molecola m1a, m2a;
							int isPossible = temp.spezza(posA, g, &m1a, &m2a);
							//std::cout << isPossible << " " << posA << " " << g << std::endl;
							if (isPossible != 0)
							{
								reactionComment reacomm2("DEGUB");
								//std::cout << "DEBUG: alternative AlkRO decomposition path worked!" << std::endl;
								reactions.push_back(Reaction(AlkRO, std::vector<Molecola>{ m1a, m2a },
									new double[3] { 7.84E+10, 0.588, 33300 }, "Alkenyl RO decomposition", reacomm2));
							}
						}
					}
				}

				//continue;	// if it is the carbon involved in the double bond continue
			}
			// <- DEBUG
			Carbonio tipoC = AlkRO.tipoC(h);
			Molecola product = AlkRO;
			Molecola m1, m2;
			int isPossible = AlkRO.spezza(posRO, h, &m1, &m2);
			if (isPossible == 0)
			{
				//std::cerr << "WARNING: alkenyl RO decomposition failed in breaking the bond with species " << AlkRO
				//	<< ", reaction skipped!" << std::endl;
				continue;
			}

			int posRadM1 = m1.trova(2);
			//m1.Crad_to_C(posRadM1);
			//m1.removeAllKeto();
			if (posRadM1 != 0)
				m1.addcheto(posRadM1);
			reactionComment reacomm = k->v_alkenyl_ro_dec(tipoC);
			reactions.push_back(Reaction(AlkRO, std::vector<Molecola>{ m1, m2 },
				new double[3] { k->A, k->n, k->E }, "Alkenyl RO decomposition", reacomm));
		}
	}
	return reactions;
}

std::vector<Reaction> aldehydesDecompositionReactions(std::vector<Molecola> ALDs, Kinox* k)
{
	Molecola OH(4);
	Molecola H2O(6);
	Molecola CO(10);
	std::vector<Reaction> reactions;
	for (auto& ALD : ALDs)
	{
		if (ALD.isAldehyde() == false)
		{
			UTL::error("aldehydesDecompositionReactions called on a non aldehyde species.");
			continue;
		}
		Molecola m1, m2;
		int posalfa[SIZEMAX + 1];
		ALD.scorri(ALD.posKeto(), 1, posalfa);
		for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)	// find the carbon in alpha
		{
			ALD.spezza(ALD.posKeto(), h, &m1, &m2);
		}

		reactionComment reacomm = k->v_ald_dec();
		reactions.push_back(Reaction(std::vector<Molecola>{ ALD, OH },
			std::vector<Molecola>{ m2, CO, H2O }, new double[3] { k->A, k->n, k->E },
			"Aldehydes decomposition", reacomm));
	}
	return reactions;
}

std::vector<Reaction> aldehydeOlefinsDecompositionReactions(std::vector<Molecola> ALDOLEs, Kinox* k)
{
	Molecola OH(4);
	Molecola H2O(6);
	Molecola CO(10);
	Molecola HCCO(11);
	std::vector<Reaction> reactions;
	for (auto& ALDOLE : ALDOLEs)
	{
		if (ALDOLE.isOleAldehyde() == false)
		{
			UTL::error("aldehydeOlefinsDecompositionReactions called on a non aldehyde olefin species.");
			continue;
		}
		if(ALDOLE.numAbstractableH(ALDOLE.posKeto()) == 1)
		{
			Molecola m1, m2;
			int posalfa[SIZEMAX + 1];
			ALDOLE.scorri(ALDOLE.posKeto(), 1, posalfa);
			bool toSkip = false;
			for (int h = 1; h <= SIZEMAX; h++) if (posalfa[h] == 1)	// find the carbon in alpha
			{
				int isPossible = ALDOLE.spezza(ALDOLE.posKeto(), h, &m1, &m2);
				if (isPossible == 0)
				{
					std::cout << "DEBUG: decomposition ALDOLE did not work for " << ALDOLE << std::endl;
					toSkip = true;
				}
			}
			if (toSkip)
				continue;

			reactionComment reacomm = k->v_ald_dec();
			reactions.push_back(Reaction(std::vector<Molecola>{ ALDOLE, OH },
				std::vector<Molecola>{ m2, CO, H2O }, new double[3] { k->A, k->n, k->E },
				"Olefin aldehydes decomposition", reacomm));
		}
		else if (ALDOLE.isole(ALDOLE.posKeto()))
		{
			int posalpha[SIZEMAX + 1];
			ALDOLE.scorri(ALDOLE.posKeto(), 1, posalpha);
			int posbeta[SIZEMAX + 1];
			ALDOLE.scorri(ALDOLE.posKeto(), 2, posbeta);
			int posgamma[SIZEMAX + 1];
			ALDOLE.scorri(ALDOLE.posKeto(), 3, posgamma);
			Molecola m1, m2;
			for (int gamma = 1; gamma <= SIZEMAX; gamma++) if (posgamma[gamma] == 1 && ALDOLE.numAbstractableH(gamma) > 0)
			{
				Molecola temp = ALDOLE;
				int alfatogamma[SIZEMAX + 1];
				int betatogamma[SIZEMAX + 1];
				temp.scorri(gamma, 1, alfatogamma);
				temp.scorri(gamma, 2, betatogamma);
				for (int alpha = 1; alpha <= SIZEMAX; alpha++) if (posalpha[alpha] == 1 && betatogamma[alpha] == 1)
				{
					for (int beta = 1; beta <= SIZEMAX; beta++) if (posbeta[beta] == 1 && alfatogamma[beta] == 1)
					{
						temp.tipo(gamma, 4);
						int isPossible = temp.spezza(alpha, beta, &m1, &m2);
						m2.removeOle(1);
						m2.removeOle(1);
						int pos1 = m2.trova(4);
						int pos2 = m2.trova(2);
						m2.tipo(pos1, 1);
						m2.tipo(pos2, 1);
						m2.addole(pos1, pos2);
						if (isPossible == 0)
							continue;
						//m1.addole(m1.posKeto(), m1.posCrad());
						std::cout << "DEBUG: ALDOLE dec worked for " << ALDOLE << std::endl;
						std::cout << ALDOLE << "  " << temp << "  " << m1 << "  " << m2 << std::endl;
						reactionComment reacomm = k->v_ald_dec();
						if (m1.numberOfC() == 2)
						{
							reactions.push_back(Reaction(std::vector<Molecola>{ ALDOLE, OH },
								std::vector<Molecola>{ HCCO, m2, H2O }, new double[3] { k->A, k->n, k->E },
								"Olefin aldehydes decomposition", reacomm));
						}
						else if (m1.numberOfC() == 3)
						{
							Molecola C2H3;
							C2H3.makeC2H5();
							C2H3.addole(1, 2);
							reactions.push_back(Reaction(std::vector<Molecola>{ ALDOLE, OH },
								std::vector<Molecola>{ C2H3, CO, m2, H2O }, new double[3] { k->A, k->n, k->E },
								"Olefin aldehydes decomposition", reacomm));
						}
					}
				}
			}


			//for (int h = 1; h <= SIZEMAX; h++) if (posalpha[h] == 1)
			//{
			//	for (int g = 1; g <= SIZEMAX; g++) if (posbeta[g] == 1)
			//	{
			//		Molecola temp = ALDOLE;
			//		for (int t = 1; t <= SIZEMAX; t++) if (posgamma[t] == 1)
			//		{
			//			temp.C_to_Crad(t);
			//			break;
			//		}
			//		int isPossible = temp.spezza(h, g, &m1, &m2);
			//		int posARad[SIZEMAX + 1];
			//		//m2.scorri(m2.posCrad(), 1, posARad);
			//		//for (int t = 1; t <= SIZEMAX; t++) if (posARad[t] == 1 && m2.tipo(t) == 1)
			//		//{
			//		//	int posCrad = m2.posCrad();
			//		//	//std::cout << "m2  " << m2 << "   " << posCrad << "  " << t << std::endl;
			//		//	m2.Crad_to_C(posCrad);
			//		//	m2.addole(posCrad, t);
			//		//}
			//		if (isPossible == 0)
			//			continue;
			//		//m1.addole(m1.posKeto(), m1.posCrad());
			//		std::cout << "DEBUG: ALDOLE dec worked for " << ALDOLE << std::endl;
			//		std::cout << ALDOLE << "  " << m1 << "  " << m2 << std::endl;
			//		reactionComment reacomm = k->v_ald_dec();
			//		reactions.push_back(Reaction(std::vector<Molecola>{ ALDOLE, OH },
			//			std::vector<Molecola>{ HCCO, m2, H2O }, new double[3] { k->A, k->n, k->E },
			//			"Olefin aldehydes decomposition", reacomm));
			//	}
			//}
		}
	}
	return reactions;
}

std::vector<Reaction> ketonesDecompositionReactions(std::vector<Molecola> KETOs,
	Kinox* k)
{
	Molecola OH(4);
	Molecola H2O(6);
	std::vector<Reaction> reactions;
	for (auto& KETO : KETOs)
	{
		if (KETO.isKetone() == false)
		{
			UTL::error("ketonesDecompositionReactions called on a non ketone species.");
			continue;
		}
		for (int i = 1; i <= KETO.size(); i++)
		{
			if (KETO.numAbstractableH(i) > 0)
			{
				int posalpha[SIZEMAX + 1];
				int posbeta[SIZEMAX + 1];
				KETO.scorri(i, 1, posalpha);
				KETO.scorri(i, 2, posbeta);
				for (int j = 1; j <= SIZEMAX; j++)
				{
					for (int h = 1; h <= SIZEMAX; h++)
					{
						if (posalpha[j] == 1 && posbeta[h] == 1)
						{
							if (KETO.dist(j, h) == 1)
							{
								Molecola Ket = KETO;
								Ket.C_to_Crad(i);
								Molecola m1, m2;
								int isPossible = Ket.spezza(j, h, &m1, &m2);
								if (isPossible == 0)
									continue;
								std::vector<Molecola> prod = { m1, m2 };
								prod = fullyDecomposeRO(prod);
								prod.push_back(H2O);
								reactionComment reacomm = k->v_h_abstraction(OH,
									KETO.tipoC(i), KETO.numAbstractableH(i),
									1);
								reactions.push_back(Reaction(std::vector<Molecola>
								{ KETO, OH }, prod,
									new double[3] { k->A, k->n, k->E },
									"Ketones decomposition", reacomm));
							}
						}
					}
				}
			}
		}
	}
	return reactions;
}

std::vector<Reaction> ketonesOlefinsDecompositionReactions(std::vector<Molecola> KETOOLEs,
	Kinox* k)
{
	std::vector<Reaction> reactions;
	for (auto& KETOOLE : KETOOLEs)
	{
		if (KETOOLE.isOleKetone() == false)
		{
			UTL::error("ketonesOlefinsDecompositionReactions called on a non olefin ketone species.");
			continue;
		}
		int posalpha[SIZEMAX + 1];
		int posKeto = KETOOLE.posKeto();
		KETOOLE.scorri(posKeto, 1, posalpha);
		for (int i = 1; i <= SIZEMAX; i++)
		{
			if (posalpha[i] == 1)
			{
				if (KETOOLE.isole(posKeto, i))
					continue;
				Molecola m1, m2;
				KETOOLE.spezza(posKeto, i, &m1, &m2);
				reactionComment reacomm = k->v_initiation(m1.tipoR(m1.trova(2)),
					m2.tipoR(m2.trova(2)), 1);
				std::vector<Molecola> prod = { m1, m2 };
				if (m1.kindOfSPecies() == unidentified_)
					continue;
				if (m2.kindOfSPecies() == unidentified_)
					continue;
				prod = fullyDecomposeRO(prod);
				reactions.push_back(Reaction(std::vector<Molecola>
				{ KETOOLE }, prod,
					new double[3] { k->A, k->n, k->E },
					"Olefins ketones decomposition", reacomm));
			}
		}
	}
	return reactions;
}





std::vector<Molecola> decomposeRO(Molecola m2)
{
	Molecola m3, m4;
	std::vector<Molecola> vec;
	// find first atom of the beta scission
	int found[SIZEMAX + 1];
	int pos2 = m2.posCrad();
	m2.scorri(pos2, 2, found);
	int posBeta1 = 0;
	for (int k = 1; k < SIZEMAX + 1; k++)
	{
		if (found[k] == 1)
		{
			if (m2.dist(m2.posCrad(), k) > 2)
				posBeta1 = k;
		}
	}
	if (posBeta1 != 0) // desired beta scission is possible
	{
		// find second atom of the beta scission
		int found2[SIZEMAX + 1];
		m2.scorri(pos2, 1, found);
		m2.scorri(posBeta1, 1, found2);
		int posBeta2 = 0;
		for (int k = 1; k < SIZEMAX + 1; k++)
			if (found[k] == 1 && found2[k] == 1)
				posBeta2 = k;
		int posRad = m2.posCrad();
		m2.spezza(posBeta1, posBeta2, &m3, &m4);
		vec.push_back(m3);
		vec.push_back(m4);
	}
	else     // desired beta scission is not possible
	{
		m2.scorri(pos2, 1, found);
		int found2[SIZEMAX + 1];
		int posBeta2 = 0;
		for (int k = 0; k < SIZEMAX + 1; k++)
		{
			if (found[k] == 1)
			{
				m2.scorri(k, 1, found2);
				for (int l = 0; l < SIZEMAX + 1; l++)
					if (found2[l] == 1 && l != pos2)
					{
						posBeta1 = k;
						posBeta2 = l;
					}
			}
		}
		if (posBeta1 != 0 && posBeta2 != 0)
		{
			m2.spezza(posBeta1, posBeta2, &m3, &m4);
			vec.push_back(m3);
			vec.push_back(m4);
		}
		else
		{
			vec.push_back(m2);
		}
	}
	//std::cout << m2 << std::endl << m3 << std::endl << m4 << std::endl;

	return vec;
}

std::vector<Molecola> fullyDecomposeRO(Molecola m2) // decompose m2 and its products too (if they are RO) and return the products
{
	std::vector<Molecola> decProds;
	decProds = decomposeRO(m2);
	if (decProds.size() > 1)
	{
		for (int i = 0; i < decProds.size(); i++)
		{
			if (decProds[i].kindOfSPecies() == RO_)
			{
				std::vector<Molecola> secDecProds;
				secDecProds = decomposeRO(decProds[i]);
				decProds.erase(decProds.begin() + i);
				decProds.insert(decProds.end(), secDecProds.begin(), secDecProds.end());
			}
		}
	}
	return decProds;
}

std::vector<Molecola> fullyDecomposeRO(std::vector<Molecola> vec)
{
	std::vector<Molecola> decProds;
	for (int j = 0; j < vec.size(); j++)
	{
		Molecola m2 = vec[j];
		if (m2.kindOfSPecies() == RO_)
		{
			decProds = decomposeRO(m2);
			if (decProds.size() > 1)
			{
				for (int i = 0; i < decProds.size(); i++)
				{
					if (decProds[i].kindOfSPecies() == RO_)
					{
						std::vector<Molecola> secDecProds;
						secDecProds = decomposeRO(decProds[i]);
						decProds.erase(decProds.begin() + i);
						decProds.insert(decProds.end(), secDecProds.begin(), secDecProds.end());
					}
				}
			}
			vec.erase(vec.begin() + j);
			vec.insert(vec.end(), decProds.begin(), decProds.end());
		}
	}

	return vec;
}

std::vector<Molecola> decomposeLinEthRO(Molecola mol)
{
	std::vector<Molecola> vec;
	int dist = mol.dist(mol.posCrad(), mol.trova(8));
	if (dist == 1)
	{
		int found1[SIZEMAX + 1];
		int found2[SIZEMAX + 1];
		mol.scorri(mol.trova(8), 1, found1);
		mol.scorri(mol.posCrad(), 2, found2);
		int pos = 0;
		for (int i = 0; i < SIZEMAX + 1; i++)
			if (found1[i] == 1 && found2[i] == 1)
				pos = i;
		Molecola m1, m2;
		mol.spezza(pos, mol.trova(8), &m1, &m2);
		std::vector<Molecola> prods = decomposeRO(m1);
		vec.push_back(m2);
		for (int i = 0; i < prods.size(); i++)
			vec.push_back(prods[i]);
		//for (int i = 0; i < vec.size(); i++)
		//	std::cout << vec[i] << std::endl;
	}
	else if (dist == 2)
	{
		int found1[SIZEMAX + 1];
		int found2[SIZEMAX + 1];
		mol.scorri(mol.trova(8), 1, found1);
		mol.scorri(mol.posCrad(), 1, found2);
		int pos = 0;
		for (int i = 0; i < SIZEMAX + 1; i++)
			if (found1[i] == 1 && found2[i] == 1)
				pos = i;
		Molecola m1, m2;
		mol.spezza(pos, mol.trova(8), &m1, &m2);
		m2.scorri(m2.posCrad(), 1, found1);
		for (int i = 0; i < SIZEMAX + 1; i++)
			if (found1[i] == 1)
				pos = i;
		m2.removeAtom(m2.posCrad());
		m2.C_to_COrad(pos);
		Molecola m3, m4;
		m2.spezza(m2.trova(9), m2.trova(7), &m3, &m4);
		//std::cout << mol << std::endl << m1 << std::endl << m3 << std::endl << m4 << std::endl;
		vec.push_back(m1);
		vec.push_back(m3);
		vec.push_back(m4);
	}

	return vec;
}

std::vector<Molecola> decomposeCEthR(Molecola m2)
{
	Molecola m3, m4;
	std::vector<Molecola> vec;
	// decompose the radical cyclic ether
	if (m2.posCrad() == m2.posEthero()[0] || m2.posCrad() == m2.posEthero()[1]) // if the radical is on one of the carbons attached to the cyclic O	
	{
		int pos1 = 0;	// oxygen bond that doubles
		int pos2 = 0;	// oxygen bond that breaks
		if (m2.posCrad() == m2.posEthero()[0])
		{
			pos1 = m2.posEthero()[0];
			pos2 = m2.posEthero()[1];
		}
		else
		{
			pos1 = m2.posEthero()[1];
			pos2 = m2.posEthero()[0];
		}

		//TEST

		m2.removeCycEther(1);
		m2.Crad_to_C(pos1);
		if (!m2.ischeto(pos1))
			m2.addcheto(pos1);
		m2.C_to_Crad(pos2);
		std::vector<Molecola> prods = decomposeRO(m2);
		for (int k = 0; k < prods.size(); k++)
			vec.push_back(prods[k]);
	}
	//if (m2.posCrad() == m2.posEthero()[0] || m2.posCrad() == m2.posEthero()[0])
	else
	{
		int found[SIZEMAX + 1];
		int pos1 = 0;		// position of the radical
		int pos2 = 0;		// position of the oxygen bond that breaks
		int pos3 = 0;		// position of the oxygen bond that becomes double bond

		m2.scorri(m2.posCrad(), 1, found);
		for (int k = 0; k < SIZEMAX + 1; k++)
		{
			if (found[k] == 1 && k == m2.posEthero()[0])
			{
				pos1 = m2.posCrad();
				pos2 = m2.posEthero()[0];
				pos3 = m2.posEthero()[1];
			}
			if (found[k] == 1 && k == m2.posEthero()[1])
			{
				pos1 = m2.posCrad();
				pos2 = m2.posEthero()[1];
				pos3 = m2.posEthero()[0];
			}
		}
		if (pos1 != 0 && pos2 != 0 && pos3 != 0)	// radical is in alpha with respect of one carbon bonded to the oxygen
		{
			m2.Crad_to_C(pos1);
			m2.removeCycEther(1);
			m2.addole(pos1, pos2);
			if (!m2.ischeto(pos3))
				m2.addcheto(pos3);
			//Molecola m3, m4;
			m2.scorri(pos3, 1, found);
			int pos4 = 0;
			for (int k = 0; k < SIZEMAX + 1; k++)
			{
				if (found[k] == 1)
				{
					if (pos4 == 0)
						pos4 = k;
					else if (m2.dist(k, pos2) < m2.dist(pos4, pos2))	// select the carbon nearest to the double bond
						pos4 = k;
				}
			}
			m2.spezza(pos3, pos4, &m3, &m4);
			int posRadm3 = m3.posCrad();
			m3.Crad_to_C(posRadm3);
			if (!m3.ischeto(posRadm3))
				m3.addcheto(posRadm3);
			vec.push_back(m3);
			vec.push_back(m4);
		}
		else    // radical is in beta with respect of one of the carbons attached to the oxygen
		{
			m2.scorri(m2.posCrad(), 2, found);
			for (int k = 0; k < SIZEMAX + 1; k++)
			{
				if (found[k] == 1 && k == m2.posEthero()[0])
					pos1 = m2.posEthero()[0];
				if (found[k] == 1 && k == m2.posEthero()[1])
					pos1 = m2.posEthero()[1];
			}
			if (pos1 != 0)
			{
				m2.scorri(m2.posCrad(), 1, found);
				int found2[SIZEMAX + 1];
				m2.scorri(pos1, 1, found2);
				for (int k = 0; k < SIZEMAX + 1; k++)
					if (found[k] == 1 && found2[k] == 1)
						pos2 = k;
				//Molecola m3, m4;
				m2.spezza(pos1, pos2, &m3, &m4);
				//std::cout << m3 << std::endl << m4 << std::endl << std::endl;
				if (m4.size() == 0) // m3 is a linear ether 
				{
					Molecola m5, m6;
					m3.scorri(m3.posCrad(), 2, found);
					m3.scorri(m3.trova(8), 1, found2);
					for (int k = 0; k < SIZEMAX + 1; k++)
						if (found[k] == 1 && found2[k] == 1)
							pos1 = k;
					m3.spezza(pos1, m3.trova(8), &m5, &m6);
					//std::cout << m5 << std::endl << m6 << std::endl;
					vec.push_back(m5);
					vec.push_back(m6);
				}
				else
				{
					std::vector<Molecola> prod = decomposeCEthR(m3);
					vec.push_back(m4);
					vec.push_back(prod[0]);
					if (prod.size() == 2)
						vec.push_back(prod[1]);
				}
			}
		}
	}
	return vec;
}

int getAdditionalProducts(std::vector<Molecola>* molVec,
	std::vector<Reaction> reactions, species kind)
{
	int newElemCount = 0;
	for (auto& reac : reactions)
	{
		std::vector<Molecola*> singProd = reac.productList();
		for (auto& prod : singProd)
		{
			if (prod->kindOfSPecies() == kind)
				if (UTL::addUnique(molVec, *prod) == -1)
					newElemCount++;
		}
	}
	return newElemCount;
}

std::vector<Molecola> getProducts(std::vector<Reaction> reactions, species kind)
{
	std::vector<Molecola> products = {};
	getAdditionalProducts(&products, reactions, kind);
	return products;
}

int carbonioToInt(Carbonio c)
{
	//enum Carbonio { Cp, Cs, Ct, Cq };
	switch (c)
	{
	case Cp:
		return 1;
		break;
	case Cs:
		return 2;
		break;
	case Ct:
		return 3;
		break;
	case Cq:
		return 4;
		break;
	default:
		std::cerr << "ERROR: carbonioToInt called on a non recognized Carbonio: " << c << ". Aborted." << std::endl;
		exit(0);
		break;
	}
}

void addNewSpecies(std::vector<Molecola>* molVec, std::vector<Reaction>* reacVec)
{
	for (auto& reac : *reacVec)
	{
		std::vector<Molecola*> reactants = reac.reactantList();
		for (auto& spec : reactants)
			UTL::addUnique(molVec, *spec);
		std::vector<Molecola*> products = reac.productList();
		for (auto& spec : products)
			UTL::addUnique(molVec, *spec);
	}
}

int processDuplicateReactions(std::vector<Reaction>* reacVec)
{
	std::vector<Reaction> newReac;
	std::vector<bool> toSkip(reacVec->size(), false);
	int mergedReactions = 0;
	for (int i = 0; i < reacVec->size(); i++)
	{
		if (toSkip[i])
			continue;
		Reaction reac = (*reacVec)[i];
		int numDuplicates = 0;
		for (int j = i + 1; j < reacVec->size(); j++)
		{
			if (toSkip[j])
				continue;
			if ((*reacVec)[i] == (*reacVec)[j])
			{
				numDuplicates++;
				toSkip[j] = true;
			}
		}
		if (numDuplicates > 0)
		{
			reac.setMultiplePathMultiplier(numDuplicates + 1);
			reac.setA(reac.A() * double(numDuplicates + 1));
			mergedReactions += numDuplicates;
		}
		newReac.push_back(reac);
	}
	for (int i = 0; i < newReac.size(); i++)
	{
		for (int j = i + 1; j < newReac.size(); j++)
		{
			if (newReac[i].weakEquality(newReac[j]))
			{
				newReac[i].setDuplicate();
				newReac[j].setDuplicate();
			}
		}
	}
	(*reacVec) = newReac;
	return mergedReactions;
}

void printSpeciesInFile(std::ofstream* outfile, std::vector<Molecola> mols,
	std::string label, ChemkinOut* chemOut)
{
	if (mols.size() > 0)
	{
		*outfile << "! " << label << std::endl;
		for (int i = 0; i < mols.size(); i++)
		{
			*outfile << std::left;
			*outfile << std::setw(16) << chemOut->molToName(mols[i])
				<< std::setw(3) << " ";
			if ((i + 1) % 4 == 0 && i != mols.size() - 1)
				*outfile << std::endl;
		}
		*outfile << std::endl;
		*outfile << std::endl;
	}
}

void printSpeciesInFile(std::ofstream* outfile, std::vector<std::string> mols,
	std::string label)
{
	if (mols.size() > 0)
	{
		*outfile << "! " << label << std::endl;
		for (int i = 0; i < mols.size(); i++)
		{
			*outfile << std::left;
			*outfile << std::setw(16) << mols[i]
				<< std::setw(3) << " ";
			if ((i + 1) % 4 == 0 && i != mols.size() - 1)
				*outfile << std::endl;
		}
		*outfile << std::endl;
		*outfile << std::endl;
	}
}

void printReaction(std::ofstream* outfile, Reaction reac, ChemkinOut* chemOut)
{
	//std::cout << reac << std::endl;
	std::vector<Molecola*> reactants = reac.reactantList();
	std::vector<Molecola*> products = reac.productList();

	std::stringstream reactantsSide;
	for (int i = 0; i < reactants.size(); i++)
	{
		reactantsSide << chemOut->molToName(*(reactants[i]));
		if (i != reactants.size() - 1)
			reactantsSide << " + ";
	}

	std::stringstream productsSide;
	for (int i = 0; i < products.size(); i++)
	{
		productsSide << chemOut->molToName(*(products[i]));
		if (i != products.size() - 1)
			productsSide << " + ";
	}

	char cost[50];
	sprintf_s(cost, 50, "   %9.3e %8.3f %8.0f ", reac.A(), reac.n(), reac.E());

	*outfile << std::left << std::setw(23) << reactantsSide.str();
	if(reac.isReversible() == false)
		*outfile << " => ";
	else
		*outfile << " = ";
	*outfile << std::setw(40) << productsSide.str() << " "
		<< std::setw(30) << cost << " "
		<< "!" << reac.printComment() << std::endl;
	if (reac.isDuplicate())
		*outfile << "DUP" << std::endl;

}

void printReactions(std::ofstream* outfile, std::vector<Reaction> reacs,
	ChemkinOut* chemOut)
{
	for (auto& rea : reacs)
		printReaction(outfile, rea, chemOut);
}

bool isThereDecompositionPath(Molecola spec,
	std::vector<baseReaction>* baseMechReacs,
	std::vector<Reaction>* totReactions, ChemkinOut* chemOut)
{
	std::string specName = chemOut->molToName(spec);
	for (auto& baseReac : *baseMechReacs)
	{
		for (auto& reactant : baseReac.reactants)
			if (specName == reactant)
				return true;
		if (baseReac.isReversible)
			for (auto& product : baseReac.products)
				if (specName == product)
					return true;
	}
	for (auto& reaction : *totReactions)
	{
		std::vector<Molecola*> reactants = reaction.reactantList();
		for (auto& reac : reactants)
			if (spec == *reac)
				return true;
		//if (reaction.isReversible())
		//{
		//	std::vector<Molecola*> products = reaction.productList();
		//
		//	for (auto& product : products)
		//		if (spec == *product)
		//			return true;
		//}
	}
	return false;
}

std::string speciesToText(species kind)
{
	switch (kind)
	{
	case fuel_:
		return "HC";
		break;
	case R_:
		return "R";
		break;
	case ROO_:
		return "ROO";
		break;
	case QOOH_:
		return "QOOH";
		break;
	case OOQOOH_:
		return "OOQOOH";
		break;
	case OLE_:
		return "OLE";
		break;
	case CO_:
		return "ALD/KETO";
		break;
	case cEth_:
		return "cEth";
		break;
	case RO_:
		break;
		return "RO";
	case KHP_:
		break;
		return "KHP";
	case ROOH_:
		return "ROO";
		break;
	case POOH2_:
		return "P(OOH)2";
		break;
	case ROOH2_:
		return "R(OOH)2";
		break;
	case oleOOH_:
		return "OLE-OOH";
		break;
	case cEthOOH_:
		return "cEth-OOH";
		break;
	case oleR_:
		return "OLE-R";
		break;
	case oleCO_:
		return "OLE-CO";
		break;
	case cEthR_:
		return "cEth-R";
		break;
	case cEthCO_:
		return "cEth-CO";
		break;
	case lEthRO_:
		return "lEth-RO";
		break;
	case alkRO_:
		return "alkRO";
		break;
	case special_:
		return "special";
		break;
	case unidentified_:
		return "unidentified";
		break;
	}
	return "ERR";
}

bool isIncluded(Reaction reac, std::vector<baseReaction>* reacList,
	ChemkinOut* chemOut)
{
	std::vector<Molecola*> reactants = reac.reactantList();
	std::vector<Molecola*> products = reac.productList();
	std::vector<std::string> reactantsNames;
	std::vector<std::string> productsNames;
	for (auto& spec : reactants)
		reactantsNames.push_back(chemOut->molToName(*spec));
	for (auto& spec : products)
		productsNames.push_back(chemOut->molToName(*spec));

	for (auto& baseReac : *reacList)
	{
		if (baseReac.isReversible == false)
		{
			if (UTL::areEquivalent(&reactantsNames, &(baseReac.reactants))
				&& UTL::areEquivalent(&productsNames, &(baseReac.products)))
				return true;
		}
		else
		{
			if ((UTL::areEquivalent(&reactantsNames, &(baseReac.reactants))
				&& UTL::areEquivalent(&productsNames, &(baseReac.products)))
				|| (UTL::areEquivalent(&productsNames, &(baseReac.reactants))
					&& UTL::areEquivalent(&reactantsNames, &(baseReac.products))))
				return true;
		}
	}
	return false;
}


std::string expandEnvVar(std::string inString)
{
	const char * inCh = inString.c_str();
	char outCh[500];
	ExpandEnvironmentStringsA(inCh, outCh, 500);
	std::string outSt(outCh);
	return outSt;
}