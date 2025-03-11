#pragma once

#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdarg.h>
#include "molec.h"
#include <vector>
#include "RateRule.h"
#include "ReactionComment.h"

//using namespace std;

#ifndef ENUM
#define ENUM
enum Carbonio { Cp, Cs, Ct, Cq };
enum Idrogeno { Hp, Hs, Ht, Hcooh };
enum Radicale { Rpmet, Rpet, Rp, Rs, Rt };
enum Anello { a5, a6, a7 };
enum AnelloO { ao3, ao4, ao5 };
enum Direz { dir, inv };
enum Output { MATRICE, FORMULA };
enum HAbsRad { o2, oh, h, o, ho2, ch3, c2h5 };
#endif

class Kinox
{
	
public:

	RateRule initiation{ 2, true };   // first value is the number of dependencies (or parameters), the second value is wether the rate rule is symmetrical or not
	RateRule hAbstraction{ 2 ,false };
	RateRule isomerizationR{ 3 ,false };
	RateRule betaDecR{ 2 ,false };
	RateRule oleFromR{ 0 ,false };
	RateRule O2AdditionR{ 1 ,false };
	RateRule O2RemovalROO{ 1 ,false };
	RateRule isomROO{ 3 ,false };
	RateRule isomOOQOOH{ 3 ,false };
	RateRule isomQOOH{ 3 ,false };
	RateRule isomPOOH2{ 3 ,false };
	RateRule O2AdditionQOOH{ 1 ,false };
	RateRule O2RemovalOOQOOH{ 1 ,false };
	RateRule OOQOOHToKHP{ 3 ,false };
	RateRule KHPDecomp{ 1 ,false };
	RateRule oleFromROO{ 2 ,false };
	RateRule oleFromBetaQOOH{ 2 ,false };
	RateRule oleFromGammaQOOH{ 2 ,false };
	RateRule oleFromDeltaQOOH{ 2 ,false };
	RateRule POOH2Dec1{ 2 ,false };
	RateRule POOHDec2{ 2 ,false };
	RateRule POOHDec3{ 2 ,false };
	RateRule etherFromQOOH{ 4 ,false };
	RateRule etherFromPOOH2{ 3 ,false };
	RateRule oleFromOOQOOH{ 2 ,false };
	RateRule oleOOHDec{ 1 ,false };
	RateRule etherOOHDec{ 1 ,false };
	RateRule cycEthDec{ 0 ,false };
	RateRule allylicRadForm{ 1 ,false };
	RateRule alkenylROForm{ 0 ,false };
	RateRule alkenylRODec{ 1 ,false };
	RateRule aldDec{ 0 ,false };


public:

	double A;
	double n;
	double E;
	int sizeReactions = 0;
	Kinox(void);    // default constructor
	Kinox(std::string nome);    // constructor
	void leggi(char[80]);

	friend std::ostream& operator<<(std::ostream&, Kinox&);
	friend std::istream& operator>>(std::istream&, Kinox&);

	reactionComment v_initiation(Radicale r1, Radicale r2, int isomers);
	reactionComment v_h_abstraction(Molecola r, Carbonio c, int numH, int isomers);
	reactionComment v_isomerization_r(Radicale r, Idrogeno h, Anello a, int numH);
	reactionComment v_beta_dec_r(Radicale r1, Radicale r2);
	reactionComment v_ole_par_r(int numH);
	reactionComment v_o2_add_r(Radicale r);
	reactionComment v_o2_rem_roo(Radicale r);
	reactionComment v_isom_roo(Radicale r, Idrogeno h, Anello a, int numH);
	reactionComment v_isom_ooqooh(Radicale r, Idrogeno h, Anello a, int numH);
	reactionComment v_isom_qooh(Radicale r1, Radicale r2, Anello a);
	reactionComment v_isom_pooh2(Radicale r1, Radicale r2, Anello a);
	reactionComment v_o2_add_qooh(Radicale r);
	reactionComment v_o2_rem_ooqooh(Radicale r);
	reactionComment v_ooqooh_to_khp(Radicale r1, Radicale r2, Anello a, int numH);
	reactionComment v_khp_decomp(Radicale r, int dist);
	reactionComment v_ole_par_roo(Radicale r1, Carbonio c, int numH);
	reactionComment v_ole_par_ooqooh(Radicale r1, Carbonio c, int numH);
	reactionComment v_ole_from_beta_qooh(Radicale r1, Radicale r2);
	reactionComment v_ole_from_gamma_qooh(Radicale r1, Radicale r2);
	reactionComment v_ole_from_delta_qooh(Radicale r1, Radicale r2);
	reactionComment v_pooh2_dec_1(Radicale r1, Radicale r2);
	reactionComment v_pooh2_dec_2(Radicale r1, Radicale r2);
	reactionComment v_pooh2_dec_3(Radicale r1, Radicale r2);
	reactionComment v_ether_from_qooh(Radicale r1, Radicale r2, AnelloO a, 
		std::string correction);
	reactionComment v_ether_from_pooh2(Radicale r1, Radicale r2, AnelloO a);
	reactionComment v_oleooh_dec(Radicale r);
	reactionComment v_etherooh_dec(Radicale r);
	reactionComment v_cyc_eth_dec();
	reactionComment v_allylic_rad_form(Carbonio c, int numH);
	reactionComment v_alkenyl_ro_form();
	reactionComment v_alkenyl_ro_dec(Carbonio c);
	reactionComment v_ald_dec();

	//void setT(double Temp);
	//void wrireaLump(std::ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	//void wrireaLump1(std::ofstream& stream, Molecola reag, HAbsRad rad, int numC, double A, double n, double E);
	//void wrireaLump3(std::ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	//void wrireaLump4(std::ofstream& stream, int numC, double A, double n, double E);     // R'numC' + O2 -> R'numC'OO
	//void wrireaLump5(std::ofstream& stream, int numC, double A, double n, double E);     // R'numC'OO -> R'numC' + O2
	//void wrireaLump6(std::ofstream& stream, int numC, double A, double n, double E);     // R'numC' + O2 -> OLE'numC' + HO2
	//void wrireaLump7(std::ofstream& stream, int numC, double A, double n, double E);     // R'numC'OO -> Q'numC'OOH
	//void wrireaLump8(std::ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH -> R'numC'OO
	//void wrireaLump9(std::ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH -> ETER'numC' + OH
	//void wrireaLump10(std::ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH -> OLE'numC' + HO2
	//void wrireaLump11(std::ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);
	//void wrireaLump11b(std::ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff1, double A1, double n1, double E1,
	//	std::vector<double> stoicCoeff2, double A2, double n2, double E2);
	//void wrireaLump12(std::ofstream& stream, int numC, double A, double n, double E);     // Q'numC'OOH + O2 -> OOQ'numC'OOH
	//void wrireaLump13(std::ofstream& stream, int numC, double A, double n, double E);     // OOQ'numC'OOH -> Q'numC'OOH + O2
	//void wrireaLump14(std::ofstream& stream, int numC, double A, double n, double E);     // OOQ'numC'OOH -> OQ'numC'OOH + OH
	//void wrireaLump15(std::ofstream& stream, Molecola reag, Molecola prod[], std::vector<double> stoicCoeff, double A, double n, double E);

	std::string nameHAbsRad(Molecola r);				// return a string containing the name of the radical rad
	std::string nameHAbsRadPlusH(Molecola r);			// return a string containing the name of the radical rad when an H is added to it
	void wrirea(std::ofstream& stream, int tiporeaz, double A, double n, double E, Molecola m1, ...);

};




