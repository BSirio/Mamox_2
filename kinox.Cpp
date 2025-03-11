
#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string>
#include <sstream>
#include "kinox.h"
//using namespace std;

//#define FILE_READING_DEBUG
#define R 1.987204


std::istream& eol(std::istream& stream);

Kinox::Kinox(void)
{
	leggi((char*)"kinox.kin");
};

void readParameters(std::ifstream& input, int entries, RateRule& RR_)
{
	std::string cell;
	std::string line;
	std::string identifier;
	std::vector<std::string> dimension(RR_.dimensionality());
	double A;
	double n;
	double E;
	for (int i = 0; i < entries; i++)
	{
		std::getline(input, line);
		std::stringstream linS(line);
		std::getline(linS, cell, ',');
		identifier = cell;
		std::getline(linS, cell, ','); std::getline(linS, cell, ','); // discard two cells

		for (int j = 0; j < RR_.dimensionality(); j++)
		{
			std::getline(linS, cell, ',');
			dimension[j] = cell;
		}
		for (int j = RR_.dimensionality(); j < 4; j++)	//discard remaining cells
			std::getline(linS, cell, ',');

		std::getline(linS, cell, ',');		// get A
		A = std::stod(cell);
		std::getline(linS, cell, ',');		// get n
		n = std::stod(cell);
		std::getline(linS, cell, ',');		// get E
		E = std::stod(cell);

		RR_.addRate(identifier, dimension, A, n, E);
	}
}

Kinox::Kinox(std::string nome)
{
	//leggi(nome);
	std::ifstream input(nome);

	std::string line;
	std::getline(input, line);				// discard first line
	while (std::getline(input, line))		// get all the lines of the file one by one
	{
		std::stringstream lineStream(line);
		std::string index;
		std::getline(lineStream, index, ',');	// get the first cell

		if (index.size() == 0) index = "0";
		std::string cell;

		int value = std::stoi(index);
		switch (value)
		{
		case 1:					// read initiation reaction parameters
		{
			readParameters(input, 15, initiation); 
			break;
		}
		case 2:
		{
			readParameters(input, 21, hAbstraction);
			break;
		}
		case 3:
		{
			readParameters(input, 27, isomerizationR);
			break;
		}
		case 4:
		{
			readParameters(input, 15, betaDecR);
			break;
		}
		case 5:
		{
			readParameters(input, 1, oleFromR);
			break;
		}
		case 6:
		{
			readParameters(input, 3, O2AdditionR);
			break;
		}
		case -6:
		{
			readParameters(input, 3, O2RemovalROO);
			break;
		}
		case 7:
		{
			readParameters(input, 27, isomROO);
			break;
		}
		case -7:
		{
			readParameters(input, 27, isomQOOH);
			break;
		}
		case 8:
		{
			readParameters(input, 3, O2AdditionQOOH);
			break;
		}
		case -8:
		{
			readParameters(input, 3, O2RemovalOOQOOH);
			break;
		}
		case 9:
		{
			readParameters(input, 17, OOQOOHToKHP);
			break;
		}
		case 10:
		{
			readParameters(input, 3, KHPDecomp);
			break;
		}
		case 11:
		{
			readParameters(input, 8, oleFromROO);
			break;
		}
		case 12:
		{
			readParameters(input, 44, etherFromQOOH);
			break;
		}
		case 13:
		{
			readParameters(input, 8, oleFromBetaQOOH);
			std::getline(input, line);
			readParameters(input, 9, oleFromGammaQOOH);
			std::getline(input, line);
			readParameters(input, 9, oleFromDeltaQOOH);
			break;
		}
		case 14:
		{
			readParameters(input, 27, isomOOQOOH);
			break;
		}
		case -14:
		{
			readParameters(input, 27, isomPOOH2);
			break;
		}
		case 15:
		{
			std::vector<std::string> namesRows{ "P", "S", "T" };
			std::vector<std::string> namesColumns{ "P", "S", "T" };
			readParameters(input, 8, oleFromOOQOOH);
			break;
		}
		case 16:
		{
			readParameters(input, 26, etherFromPOOH2);
			break;
		}
		case 17:
		{
			readParameters(input, 8, POOH2Dec1);
			std::getline(input, line);
			readParameters(input, 9, POOHDec2);
			std::getline(input, line);
			readParameters(input, 9, POOHDec3);
			break;
		}
		case 18:
		{
			readParameters(input, 3, oleOOHDec);
			break;
		}
		case 19:
		{
			readParameters(input, 3, etherOOHDec);
			break;
		}
		case 20:
		{
			readParameters(input, 1, cycEthDec);
			break;
		}
		case 21:
		{
			readParameters(input, 3, allylicRadForm);
			break;
		}
		case 22:
		{
			readParameters(input, 1, alkenylROForm);
			break;
		}
		case 23:
		{
			readParameters(input, 4, alkenylRODec);
			break;
		}
		case 24:
		{
			readParameters(input, 1, aldDec);
			break;
		}
		default:
		{
			break;
		}
		}
	}
};

void Kinox::leggi(char nome[80])
{
	std::ifstream myfile(nome);
	if (!myfile)
	{
		std::cerr << "Error in opening the kinetic constants file: " << nome;
		exit(1);
	};
	myfile >> *this;
	myfile.close();
};

std::istream& eol(std::istream& stream)
/////////////////////////////////////////////////////////////
// It is a stream manipulator:
//     stream >> eol      flushes until the end of the line
/////////////////////////////////////////////////////////////
{
	char ch;
	do { stream.get(ch); } while (ch != '\n');
	return stream;
};

std::ostream& operator<<(std::ostream& stream, Kinox& obj)
/////////////////////////////////////////////////////////////
// Write the object into a stream
/////////////////////////////////////////////////////////////
{
	char buff[40];

	return stream;
};

std::istream& operator>>(std::istream& stream, Kinox& obj)
/////////////////////////////////////////////////////////////
// Read the object from a stream
/////////////////////////////////////////////////////////////
{
	char ch;

	return stream;
};

std::string MeEtPST(Radicale r)
{
	switch (r)
	{
	case Rpmet:
		return "me";
		break;
	case Rpet:
		return "et";
		break;
	case Rp:
		return "P";
		break;
	case Rs:
		return "S";
		break;
	case Rt:
		return "T";
		break;
	default:
		std::cout << "Error in MeEtPST: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

std::string MeEtPST_comm(Radicale r)
{
	switch (r)
	{
	case Rpmet:
		return "methyl";
		break;
	case Rpet:
		return "ethyl";
		break;
	case Rp:
		return "primary";
		break;
	case Rs:
		return "secondary";
		break;
	case Rt:
		return "tertiary";
		break;
	default:
		std::cout << "Error in MeEtPST_comm: " << r << " is not a valid value!" << std::endl;
		break;
	}
}


std::string HAbs(HAbsRad r)
{
	switch (r)
	{
	case o2:
		return "O2";
		break;
	case oh:
		return "OH";
		break;
	case h:
		return "H";
		break;
	case o:
		return "O";
		break;
	case ho2:
		return "HO2";
		break;
	case ch3:
		return "CH3";
		break;
	case c2h5:
		return "C2H5";
		break;
	default:
		std::cout << "Error in HAbs: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

//std::string HAbs(Molecola r)
//{
//	Molecola CH3;
//	CH3.makeCH3();
//	Molecola C2H5;
//	C2H5.makeC2H5();
//	
//	if (r.isSpecialMolecule() == 2)
//		return "O2";
//	else if (r.isSpecialMolecule() == 3)
//		return "O";
//	else if (r.isSpecialMolecule() == 4)
//		return "OH";
//	else if (r.isSpecialMolecule() == 5)
//		return "HO2";
//	else if (r.isSpecialMolecule() == 8)
//		return "H";
//	else if (r == CH3)
//		return "CH3";
//	else if (r == C2H5)
//		return "C2H5";
//	else
//		UTL::error("in HAbs: molecule has not been recognized!");
//	return "";
//}

std::string PST(Carbonio c)
{
	switch (c)
	{
	case Cp:
		return "P";
		break;
	case Cs:
		return "S";
		break;
	case Ct:
		return "T";
		break;
	default:
		std::cout << "Error in PST: " << c << " is not a valid value!" << std::endl;
		break;
	}
}


std::string PST_comm(Carbonio c)
{
	switch (c)
	{
	case Cp:
		return "primary";
		break;
	case Cs:
		return "secondary";
		break;
	case Ct:
		return "tertiary";
		break;
	default:
		std::cout << "Error in PST_comm: " << c << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PSTQ(Carbonio c)
{
	switch (c)
	{
	case Cp:
		return "P";
		break;
	case Cs:
		return "S";
		break;
	case Ct:
		return "T";
		break;
	case Cq:
		return "Q";
		break;
	default:
		std::cout << "Error in PSTQ: " << c << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PSTQ_comm(Carbonio c)
{
	switch (c)
	{
	case Cp:
		return "primary";
		break;
	case Cs:
		return "secondary";
		break;
	case Ct:
		return "tertiary";
		break;
	case Cq:
		return "quaternary";
		break;
	default:
		std::cout << "Error in PSTQ_comm: " << c << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PST(Radicale r)
{
	switch (r)
	{
	case Rp:
		return "P";
		break;
	case Rs:
		return "S";
		break;
	case Rt:
		return "T";
		break;
	default:
		std::cout << "Error in PST: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PST_comm(Radicale r)
{
	switch (r)
	{
	case Rp:
		return "primary";
		break;
	case Rs:
		return "secondary";
		break;
	case Rt:
		return "tertiary";
		break;
	default:
		std::cout << "Error in PST_comm: " << r << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PST(Idrogeno h)
{
	switch (h)
	{
	case Hp:
		return "P";
		break;
	case Hs:
		return "S";
		break;
	case Ht:
		return "T";
		break;
	default:
		std::cout << "Error in PST: " << h << " is not a valid value!" << std::endl;
		break;
	}
}

std::string PST_comm(Idrogeno h)
{
	switch (h)
	{
	case Hp:
		return "primary";
		break;
	case Hs:
		return "secondary";
		break;
	case Ht:
		return "tertiary";
		break;
	default:
		std::cout << "Error in PST_comm: " << h << " is not a valid value!" << std::endl;
		break;
	}
}

std::string ring(Anello a)
{
	switch (a)
	{
	case a5:
		return "5";
		break;
	case a6:
		return "6";
		break;
	case a7:
		return "7";
		break;
	default:
		std::cout << "Error in ring: " << a << " is not a valid value!" << std::endl;
		break;
	}
}


std::string ring_comm(Anello a)
{
	switch (a)
	{
	case a5:
		return "5 members";
		break;
	case a6:
		return "6 members";
		break;
	case a7:
		return "7 members";
		break;
	default:
		std::cout << "Error in ring_comm: " << a << " is not a valid value!" << std::endl;
		break;
	}
}

std::string ringO(AnelloO a)
{
	switch (a)
	{
	case ao3:
		return "3";
		break;
	case ao4:
		return "4";
		break;
	case ao5:
		return "5";
		break;
	default:
		std::cout << "Error in ringO: " << a << " is not a valid value!" << std::endl;
		break;
	}
}

std::string ringO_comm(AnelloO a)
{
	switch (a)
	{
	case ao3:
		return "3 members";
		break;
	case ao4:
		return "4 members";
		break;
	case ao5:
		return "5 members";
		break;
	default:
		std::cout << "Error in ringO_comm: " << a << " is not a valid value!" << std::endl;
		break;
	}
}

std::string numCtolabel(int numC)
{
	if (numC == 5)
	{
		return "5";
	}
	if (numC > 5 && numC < 12)
	{
		return "6 to 11";
	}
	if (numC > 11 && numC < 15)
	{
		return "12 to 14";
	}
	if (numC > 14)
	{
		return "15+";
	}
	return "err";
}

reactionComment Kinox::v_initiation(Radicale r1, Radicale r2, int isomers)
{
	std::string RRID;
	std::vector<std::string> param = { MeEtPST(r1) , MeEtPST(r2) };
	std::vector<double> kinPar = initiation.returnRates(param, &RRID);
	A = kinPar[0] * double(isomers);
	n = kinPar[1];
	E = kinPar[2];

	std::string HRComm = "Type radical 1: ";
	HRComm.append(MeEtPST_comm(r1));
	HRComm.append("; Type radical 2: ");
	HRComm.append(MeEtPST_comm(r2));

	reactionComment comment(RRID, HRComm, 1, isomers);
	return comment;
}

reactionComment Kinox::v_h_abstraction(Molecola r, Carbonio c, int numH, int isomers)
{
	std::string RRID;
	std::vector<std::string> param = { nameHAbsRad(r) , PST(c)};
	std::vector<double> kinPar = hAbstraction.returnRates(param, &RRID);
	A = kinPar[0] * double(isomers) * double(numH);
	n = kinPar[1];
	E = kinPar[2];

	std::string comm = "Type of attacked carbon: ";
	comm.append(PST_comm(c));
	comm.append("; correction: ");

	reactionComment comment(RRID, comm, numH, isomers);
	return comment;
}

reactionComment Kinox::v_isomerization_r(Radicale r, Idrogeno h, Anello a, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) , PST(h) , ring(a) };
	std::vector<double> kinPar = isomerizationR.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));
	comm.append("; Type of hydrogen: ");
	comm.append(PST_comm(h));
	comm.append("; Size of ring: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_beta_dec_r(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), MeEtPST(r2) };
	std::vector<double> kinPar = betaDecR.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of decomposing radical: ";
	comm.append(PST_comm(r1));
	comm.append("; Type of expelled radical: ");
	comm.append(MeEtPST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ole_par_r(int numH)
{
	std::string RRID;
	std::vector<std::string> param = { };
	std::vector<double> kinPar = oleFromR.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID, numH);
	return comment;
}

reactionComment Kinox::v_o2_add_r(Radicale r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	std::vector<double> kinPar = O2AdditionR.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_o2_rem_roo(Radicale r)
{
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));
	comm.append("; number of carbons: ");
	
	std::string RRID;
	std::vector<std::string> param = {PST(r)};
	std::vector<double> kinPar = O2RemovalROO.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID, comm);
	return comment;

}


reactionComment Kinox::v_isom_roo(Radicale r, Idrogeno h, Anello a, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r), PST(h), ring(a) };
	std::vector<double> kinPar = isomROO.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OO: ";
	comm.append(PST_comm(r));
	comm.append("; type of H:");
	comm.append(PST_comm(h));
	comm.append("; size of ring: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_isom_ooqooh(Radicale r, Idrogeno h, Anello a, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r), PST(h), ring(a) };
	std::vector<double> kinPar = isomOOQOOH.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of C-OO: ";
	comm.append(PST_comm(r));
	comm.append("; type of H:");
	comm.append(PST_comm(h));
	comm.append("; size of ring: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_isom_qooh(Radicale r1, Radicale r2, Anello a)
{
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; C-H type: ");
	comm.append(PST_comm(r2));
	comm.append("; ring size : ");
	comm.append(ring_comm(a));
	comm.append("; number of carbons: ");
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2), ring(a)};
	std::vector<double> kinPar = isomQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_isom_pooh2(Radicale r1, Radicale r2, Anello a)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1), PST(r2), ring(a) };
	std::vector<double> kinPar = isomPOOH2.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; type C-H: ");
	comm.append(PST_comm(r2));
	comm.append("; ring size: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_o2_add_qooh(Radicale r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	std::vector<double> kinPar = O2AdditionQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_o2_rem_ooqooh(Radicale r)
{
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(r));
	comm.append("; number of carbons: ");
	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	std::vector<double> kinPar = O2RemovalOOQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ooqooh_to_khp(Radicale r1, Radicale r2, Anello a, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) , ring(a) };
	std::vector<double> kinPar = OOQOOHToKHP.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; type C-OOH: ");
	comm.append(PST_comm(r2));
	comm.append("; ring size: ");
	comm.append(ring_comm(a));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_khp_decomp(Radicale r, int dist)
{
	std::string posStr;
	if (dist == 1)
		posStr = "beta";
	if (dist == 2)
		posStr = "gamma";
	if (dist > 2)
		posStr = "other";

	std::string RRID;
	std::vector<std::string> param = { PST(r) , posStr };
	std::vector<double> kinPar = KHPDecomp.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r));
	comm.append("; CO position: ");
	switch (dist) {
	case 1:
		comm.append("beta");
		break;
	case 2:
		comm.append("gamma");
		break;
	default:
		comm.append("other");
		break;
	}

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ole_par_roo(Radicale r1, Carbonio c, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(c) };
	std::vector<double> kinPar = oleFromROO.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; type C-H: ");
	comm.append(PST_comm(c));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_ole_par_ooqooh(Radicale r1, Carbonio c, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(c) };
	std::vector<double> kinPar = oleFromOOQOOH.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OO: ";
	comm.append(PST_comm(r1));
	comm.append("; type C-H: ");
	comm.append(PST_comm(c));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_ole_from_beta_qooh(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = oleFromBetaQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; radical type: ");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ole_from_gamma_qooh(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = oleFromGammaQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];;
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; radical type: ");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ole_from_delta_qooh(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = oleFromDeltaQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; radical type: ");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_pooh2_dec_1(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = POOH2Dec1.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	///return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; radical type: ");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_pooh2_dec_2(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = POOHDec2.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; radical type: ");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_pooh2_dec_3(Radicale r1, Radicale r2)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2) };
	std::vector<double> kinPar = POOHDec3.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; radical type: ");
	comm.append(PST_comm(r2));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ether_from_qooh(Radicale r1, Radicale r2, AnelloO a, 
	std::string correction)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2), ringO(a), correction };
	std::vector<double> kinPar = etherFromQOOH.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; type of radical: ");
	comm.append(PST_comm(r2));
	comm.append("; size of ring: ");
	comm.append(ringO_comm(a));

	reactionComment comment(RRID, comm);
	return comment;
}


reactionComment Kinox::v_ether_from_pooh2(Radicale r1, Radicale r2, AnelloO a)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r1) , PST(r2), ringO(a) };
	std::vector<double> kinPar = etherFromPOOH2.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r1));
	comm.append("; type of radical: ");
	comm.append(PST_comm(r2));
	comm.append("; ring size: ");
	
	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_oleooh_dec(Radicale r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	std::vector<double> kinPar = oleOOHDec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r));
	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_etherooh_dec(Radicale r)
{
	std::string RRID;
	std::vector<std::string> param = { PST(r) };
	std::vector<double> kinPar = etherOOHDec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type C-OOH: ";
	comm.append(PST_comm(r));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_cyc_eth_dec()
{
	std::string RRID;
	std::vector<std::string> param = { };
	std::vector<double> kinPar = cycEthDec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	reactionComment comment(RRID);
	return comment;
}

reactionComment Kinox::v_allylic_rad_form(Carbonio c, int numH)
{
	std::string RRID;
	std::vector<std::string> param = { PST(c) };
	std::vector<double> kinPar = allylicRadForm.returnRates(param, &RRID);
	A = kinPar[0] * double(numH);
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of radical: ";
	comm.append(PST_comm(c));

	reactionComment comment(RRID, comm, numH);
	return comment;
}

reactionComment Kinox::v_alkenyl_ro_form()
{
	std::string RRID;
	std::vector<std::string> param = { };
	std::vector<double> kinPar = alkenylROForm.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];

	reactionComment comment(RRID);
	return comment;
}

reactionComment Kinox::v_alkenyl_ro_dec(Carbonio c)
{
	std::string RRID;
	std::vector<std::string> param = { PSTQ(c) };
	std::vector<double> kinPar = alkenylRODec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	//return A * pow(Temp, n) * exp(-E / R / Temp);
	std::string comm = "Type of carbon: ";
	comm.append(PSTQ_comm(c));

	reactionComment comment(RRID, comm);
	return comment;
}

reactionComment Kinox::v_ald_dec()
{
	std::string RRID;
	std::vector<std::string> param = { };
	std::vector<double> kinPar = aldDec.returnRates(param, &RRID);
	A = kinPar[0];
	n = kinPar[1];
	E = kinPar[2];
	reactionComment comment(RRID);
	return comment;
}





std::string Kinox::nameHAbsRad(Molecola r)
{

	Molecola CH3;
	CH3.makeCH3();
	Molecola C2H5;
	C2H5.makeC2H5();

	if (r.isSpecialMolecule() == 2)
		return "O2";
	else if (r.isSpecialMolecule() == 3)
		return "O";
	else if (r.isSpecialMolecule() == 4)
		return "OH";
	else if (r.isSpecialMolecule() == 5)
		return "HO2";
	else if (r.isSpecialMolecule() == 8)
		return "H";
	else if (r == CH3)
		return "CH3";
	else if (r == C2H5)
		return "C2H5";
	else
		UTL::error("in HAbs: molecule has not been recognized!");
	return "ERR.";
}

std::string Kinox::nameHAbsRadPlusH(Molecola r)
{
	Molecola CH3;
	CH3.makeCH3();
	Molecola C2H5;
	C2H5.makeC2H5();

	if (r.isSpecialMolecule() == 2)
		return "HO2";
	else if (r.isSpecialMolecule() == 3)
		return "OH";
	else if (r.isSpecialMolecule() == 4)
		return "H2O";
	else if (r.isSpecialMolecule() == 5)
		return "H2O2";
	else if (r.isSpecialMolecule() == 8)
		return "H2";
	else if (r == CH3)
		return "CH4";
	else if (r == C2H5)
		return "C2H6";
	else
		UTL::error("in HAbs: molecule has not been recognized!");
	return "ERR.";

}

