#ifndef VD_HEADER

#define VD_HEADER

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/stat.h>
#include <stdlib.h>
#include <vector>
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"

#define _N_AGREGATE 300

//Patch namespace - used to fix known bugs & issues of the software
namespace patch
{
	template < typename T > std::string to_string( const T& n )
	//Bug in the old g++ compilator that crushes when using std::to_string function
	{
		std::ostringstream stm ;
		stm << n ;
		return stm.str() ;
	}
}

///////////////////////////////////////
const int error_code[] = {10,11,12,13,14};
//ERROR CODES:
////    Usage & input files:
////    10 - uncorrect usage of the file (see usage)
////    11 - uncorrect parameter (see usage)
////    12 - no file name after-i option (see usage)
////    13 - file not found
////    14 - uncorrect parameter in GetSlope("string",std::map<...>)

const double Ylimit = 2.;
const double Xlimit = 60.;

#pragma pack(push, 1)
struct vDriftTreeStructure{
	double slave_Y;
	double master_Y;
	double slave_recVDrift;
	int    sectionID;
	unsigned int runNumber;
	unsigned int eventUnixTime;
	unsigned int eventNumber;
};
#pragma pack(pop)

void Usage(char** arg_string)
{
	std::cerr << "\nUsage: " << arg_string[0] << " [options] [argument]" << std::endl;
	std::cerr << "  Options:" << std::endl;
	std::cerr << "    -i                         : use this option to add input file " << std::endl;
	std::cerr << "    -v                         : verbose output " << std::endl;
	std::cerr << "    --root-output              : adds root file to standard output " << std::endl;
	std::cerr << "    --no-factor-multiplication : disable multiplication rule for calibration factor " << std::endl;
	std::cerr << "    --no-TOF-MTPC-factor       : disable TOF-MTPC matching correction " << std::endl;
	std::cerr << "    --check-algorithm          : check the algorithm (multiplication rule & TOF-MTPC correction are off) " << std::endl;
	std::cerr << "  Example: " << std::endl;
	std::cerr << "    ./vDriftAnalyzer -i /home/jdoe/input_file.root\n" << std::endl;
}

void GetCalibTreesFromFile(TFile* inFile, std::map<std::string, TTree*> &trees)
{
	trees["MTPCLvsTOFL"]  = (TTree*) inFile->Get("MTPCLvsTOFL");
	trees["MTPCRvsTOFR"]  = (TTree*) inFile->Get("MTPCRvsTOFR");
	trees["VTPC2vsMTPCL"] = (TTree*) inFile->Get("VTPC2vsMTPCL");
	trees["VTPC2vsMTPCR"] = (TTree*) inFile->Get("VTPC2vsMTPCR");
	trees["VTPC1vsVTPC2"] = (TTree*) inFile->Get("VTPC1vsVTPC2");
}

void ReadBranchesFromTree(TTree* tree, vDriftTreeStructure &data, std::string optSwap = "")
{
	std::string s_name;
	s_name = (std::string) tree->GetName();
	if (optSwap == "swap") std::cout << "ReadBranchesFromTree(): Option 'swap' is ON." << std::endl;

	if (optSwap == ""){
		if (s_name.find("TOF") != std::string::npos){
			tree->SetBranchAddress("tofY", &data.master_Y);
		}
		else{
			tree->SetBranchAddress("master_Y", &data.master_Y);
		}
		tree->SetBranchAddress("slave_Y", &data.slave_Y);
		tree->SetBranchAddress("slave_recVDrift", &data.slave_recVDrift);
	}
	if (optSwap == "swap"){
		if (s_name.find("TOF") != std::string::npos){
			tree->SetBranchAddress("tofY", &data.slave_Y);
			std::cerr << "[WARNING] ReadBranchesFromTree(): Attempt to read drift velocity from TOF!" << std::endl; 
		}
		else{
			tree->SetBranchAddress("master_Y", &data.slave_Y);
			tree->SetBranchAddress("master_recVDrift", &data.slave_recVDrift);
		}
		tree->SetBranchAddress("slave_Y", &data.master_Y);
	}
	tree->SetBranchAddress("sectionId", &data.sectionID);
	tree->SetBranchAddress("runNumber", &data.runNumber);
	tree->SetBranchAddress("eventNumber", &data.eventNumber);
	tree->SetBranchAddress("eventUnixTime", &data.eventUnixTime);
}

inline bool check_file_exists(const std::string& name)
{
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}

void Aggregate(std::string inFileName, std::string outFileName)
{
	std::ifstream infile(inFileName.c_str());
	std::map<long int,double> mTimeVd;
	long int t;
	double   vd;

	while (infile >> t >> vd)
		{
			mTimeVd[t] = vd;
		}
  
}

#endif
