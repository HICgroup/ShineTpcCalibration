#include "vDriftAnalyzer.hpp"

int main(int argc, char** argv)
{
	std::cout << "Starting analyser." << std::endl;

	std::string iFileName;
	TFile* inputFile;
	TFile* f_out;
	bool is_verbose;
	bool is_root_output;
	bool is_no_multiply;
	bool is_no_TOF_factor;
	// long int prev_time;
	std::string s_cut;
	long int prev_time=0;
	is_verbose       = false;
	is_root_output   = false;
	is_no_multiply   = false;
	is_no_TOF_factor = false;
  
	//-------------------------------USAGE-----------------------------------------//
	if (argc<3){
		Usage(argv);
		return error_code[0];
	}

	for(int i=1; i<argc; i++){
		if(std::string(argv[i]) != "-i" &&
		   std::string(argv[i]) != "-v" &&
		   std::string(argv[i]) != "--root-output" &&
		   std::string(argv[i]) != "--no-factor-multiplication" &&
		   std::string(argv[i]) != "--no-TOF-MTPC-factor" &&
		   std::string(argv[i]) != "--check-algorithm"){
			std::cerr << "\n[ERROR]: Unknown parameter " << i << ": " <<  argv[i] << std::endl;
			Usage(argv);
			std::cout << "Analyzer finished with exit code " << error_code[1] << std::endl;
			return error_code[1];
		} else {
			if(std::string(argv[i]) == "-i" && i!=argc-1) {
				iFileName = argv[++i];
				continue;
			}
			if(std::string(argv[i]) == "-i" && i==argc-1) {
				std::cerr << "\n[ERROR]: File name was not specified " << std::endl;
				Usage(argv);
				std::cout << "Analyzer finished with exit code " << error_code[2] << std::endl;
				return error_code[2];
			}
			if(std::string(argv[i]) == "-v") {
				is_verbose = true;
				continue;
			}
			if(std::string(argv[i]) == "--root-output") {
				is_root_output = true;
				continue;
			}
			if(std::string(argv[i]) == "--no-factor-multiplication") {
				is_no_multiply = true;
				continue;
			}
			if(std::string(argv[i]) == "--no-TOF-MTPC-factor") {
				is_no_TOF_factor = true;
				continue;
			}
			if(std::string(argv[i]) == "--check-algorithm") {
				is_no_multiply   = true;
				is_no_TOF_factor = true;
				continue;
			}
		}
	}
	//-------------------------------USAGE-----------(end)-------------------------//

	std::cout << "\nStarting options:" << std::endl;
	std::cout << "\nFile        : " << iFileName.c_str() << std::endl;
	if (is_verbose)  std::cout      << "Verbose mode          : ON" << std::endl;
	if (!is_verbose) std::cout      << "Verbose mode          : OFF" << std::endl;
	if (is_root_output) std::cout   << "ROOT output           : ON" << std::endl;
	if (!is_root_output) std::cout  << "ROOT output           : OFF" << std::endl;
	if (is_no_multiply) std::cout   << "Multiplication factor : OFF - !!! Make sure this is intentional !!!" << std::endl;
	if (is_no_TOF_factor) std::cout << "MTPC-TOF factor       : OFF - !!! Make sure this is intentional !!!" << std::endl;
	std::cout << std::endl;

	if (check_file_exists(iFileName.c_str())){
		inputFile = new TFile(iFileName.c_str(),"ro");
	} else {
		std::cerr << "\n[ERROR]: No file " << iFileName.c_str() << " was found!" << std::endl;
		std::cout << "Analyzer finished with exit code " << error_code[3] << std::endl;
		return error_code[3];
	}
	//-----------------------------------------------------------------------------//

	std::cout << "Getting TTree* from file." << std::endl;

	std::map<std::string, TTree*> dataTrees;

	GetCalibTreesFromFile(inputFile, dataTrees);

	std::map<std::string, vDriftTreeStructure> inData;

	if (is_root_output) f_out = new TFile((iFileName+"_output.root").c_str(),"recreate");

	std::vector<std::string>          v_order;
	std::vector<double>               v_slope;
	std::vector<double>               v_slope_prev;
	std::map<std::string,TGraph*>     m_vDslope;
	std::map<std::string,TGraph*>     m_vDtime;
	std::map<std::string,TGraph*>     m_vDnom;
	std::map<std::string,std::string> filename;
	std::map<std::string, int>        TofFactor;
  
	double v_drift      = 0.;
	double v_drift_average;
	long int time_average=0;
	double v_drift_prev;

	//-----------------------TPC-order-of-calibration------------------------------//
	v_order.push_back("MTPCLvsTOFL");
	v_order.push_back("VTPC2vsMTPCL");
	v_order.push_back("VTPC1vsVTPC2");
	v_order.push_back("MTPCRvsTOFR");
	v_order.push_back("VTPC2vsMTPCR");
	//-----------------------TPC-order-of-calibration-(end)------------------------//
	
	filename["MTPCLvsTOFL"]   = "MTPCL";
	filename["MTPCRvsTOFR"]   = "MTPCR";
	filename["VTPC2vsMTPCL"]  = "VTPC2";
	filename["VTPC2vsMTPCR"]  = "MTPCRfromVTPC2";
	filename["VTPC1vsVTPC2"]  = "VTPC1";

	if (!is_no_TOF_factor) TofFactor["MTPCLvsTOFL"] = 2;
	if (is_no_TOF_factor) TofFactor["MTPCLvsTOFL"] = 1;
	if (!is_no_TOF_factor) TofFactor["MTPCRvsTOFR"] = 2;
	if (is_no_TOF_factor) TofFactor["MTPCRvsTOFR"] = 1;
	TofFactor["VTPC2vsMTPCL"] = 1;
	TofFactor["VTPC2vsMTPCR"] = 1;
	TofFactor["VTPC1vsVTPC2"] = 1;

	if (is_verbose){
		std::cout << "\nCurrent TPC calibration order:" << std::endl;
		for (unsigned int i=0; i<v_order.size();i++){
			std::cout << v_order.at(i) << std::endl;
		}
		std::cout << std::endl;
	}

	for (unsigned int i=0;i<v_order.size();i++){
    
		v_drift_prev = -1.;
    
		if (is_verbose) std::cout << "\n" << v_order.at(i).c_str()<<  ": " << dataTrees[v_order.at(i)]->GetEntriesFast() << " entries" << std::endl;
		if (!is_verbose) std::cout << "\n" << (v_order.at(i)).c_str() << std::endl;
		if (is_root_output) f_out->mkdir(v_order.at(i).c_str());
    
		if (dataTrees[v_order.at(i)]->GetEntriesFast() == 0){
			std::cerr << "\n[WARNING] No data in this tree. Skipping..." << std::endl;
			continue;
		}

		std::cout << "Collecting unbinned data." << std::endl;
    
		if (dataTrees[v_order.at(i)]->GetEntriesFast() == 0){
			std::cerr << "\n[WARNING] Tree has 0 entries. Skipping..." << std::endl;
			continue;
		}

		if (v_order.at(i) != "VTPC2vsMTPCR") ReadBranchesFromTree(dataTrees[v_order.at(i)], inData[v_order.at(i)]);
		if (v_order.at(i) == "VTPC2vsMTPCR") ReadBranchesFromTree(dataTrees[v_order.at(i)], inData[v_order.at(i)], "swap");

		std::map<int,std::vector<double> > vY;
		std::map<int,std::vector<double> > vDY;
		
		for (long int i_tree=0;i_tree<dataTrees[v_order.at(i)]->GetEntriesFast();i_tree++){
			dataTrees[v_order.at(i)]->GetEntry(i_tree);

			if (inData[v_order.at(i)].sectionID == -1) continue;
			
			if ( TMath::Abs(inData[v_order.at(i)].slave_Y-inData[v_order.at(i)].master_Y) <= Ylimit ){
				vY[inData[v_order.at(i)].sectionID].push_back(inData[v_order.at(i)].master_Y);
				vDY[inData[v_order.at(i)].sectionID].push_back( inData[v_order.at(i)].slave_Y-inData[v_order.at(i)].master_Y );
			}

		}

		for (int i_sect=0; i_sect<(int)dataTrees[v_order.at(i)]->GetMaximum("sectionId")+1;i_sect++){ ///<-----ERROR WAS HERE
			TGraph* gr = new TGraph(vY[i_sect].size(),&(vY[i_sect][0]),&(vDY[i_sect][0]));
			gr->SetName(("g_"+v_order.at(i)+"_Id_"+patch::to_string(i_sect)).c_str());
			gr->SetTitle(("graph of dY vs Y from "+v_order.at(i)+" for section "+patch::to_string(i_sect)).c_str());
			TF1* func = new TF1("func","pol1",-Xlimit,Xlimit);
			if (is_verbose) gr->Fit(func,"R");
			if (!is_verbose) gr->Fit(func,"RQ");
			v_slope.push_back(func->GetParameter(1));
			if (is_root_output){
				// std::cout << "  Writing in root file." << std::endl;
				f_out->cd(v_order.at(i).c_str());
				gr->Write();
			}

			delete gr;
			delete func;
			vY[i_sect].clear();
			vDY[i_sect].clear();
		}
		

		std::ofstream                 myfile;
		std::vector<double>           v_drift_avg;
		std::vector<long int>         time_avg;
		std::map<unsigned int,double> vD_Factor;
		std::vector<double>           vtime;
		std::vector<double>           vslope;
		std::vector<double>           vdriftVector;
		std::vector<double>           vdriftNominal;
		std::map<int,std::vector<double> >           slave_Y_fixed;
		std::map<int,std::vector<double> >           v_master_Y;

		std::cout << "Writing txt output file: " << filename[v_order.at(i)] << ".txt" << std::endl;

		myfile.open ((filename[v_order.at(i)]+".txt").c_str());
		//myfile << "# UnixTime[sec] vDrift[cm/usec]\n";

		for (long int i_tree=0;i_tree<dataTrees[v_order.at(i)]->GetEntriesFast();i_tree++){
			dataTrees[v_order.at(i)]->GetEntry(i_tree);

			if (inData[v_order.at(i)].sectionID == -1)                                       continue;
      
			vtime.push_back(inData[v_order.at(i)].eventUnixTime);
			vslope.push_back(v_slope.at(inData[v_order.at(i)].sectionID));
			vdriftNominal.push_back(inData[v_order.at(i)].slave_recVDrift * 1e3);
			
			vD_Factor[inData[v_order.at(i)].eventUnixTime] = 1 / ( 1 + TofFactor[v_order.at(i)] * v_slope.at(inData[v_order.at(i)].sectionID) );
			if (v_order.at(i)!="MTPCLvsTOFL" && v_order.at(i)!="MTPCRvsTOFR" && (!is_no_multiply)){
				for (unsigned int i_factor=0; i_factor<i;i_factor++){
					vD_Factor[inData[v_order.at(i)].eventUnixTime] /= ( 1 + TofFactor[v_order.at(i_factor)] * m_vDslope[v_order.at(i_factor)]->Eval(inData[v_order.at(i)].eventUnixTime) );
				}
			}
      
			v_drift = 1e3*inData[v_order.at(i)].slave_recVDrift * vD_Factor[inData[v_order.at(i)].eventUnixTime];
			vdriftVector.push_back(v_drift);

			if ( TMath::Abs(inData[v_order.at(i)].slave_Y-inData[v_order.at(i)].master_Y) <= Ylimit ){
				slave_Y_fixed[inData[v_order.at(i)].sectionID].push_back( inData[v_order.at(i)].slave_Y * vD_Factor[inData[v_order.at(i)].eventUnixTime] - inData[v_order.at(i)].master_Y );

				// slave_Y_fixed[inData[v_order.at(i)].sectionID].push_back(2*(inData[v_order.at(i)].slave_Y - inData[v_order.at(i)].master_Y ));

				v_master_Y[inData[v_order.at(i)].sectionID].push_back(inData[v_order.at(i)].master_Y);
			}
			
			myfile << inData[v_order.at(i)].eventUnixTime << " " << v_drift << std::endl;
			
			if (is_verbose) std::cout << inData[v_order.at(i)].eventUnixTime << " " << v_drift << " | vD prev " << v_drift_prev << " | sectionId " << inData[v_order.at(i)].sectionID << " slope = " << v_slope.at(inData[v_order.at(i)].sectionID) << " || factor = " << vD_Factor[inData[v_order.at(i)].eventUnixTime] << std::endl;

		}

		std::cout << "Collecting dYvsY with calibrated values for: " << filename[v_order.at(i)] << std::endl;
    
		m_vDslope[v_order.at(i)] = new TGraph(vtime.size(),&vtime[0],&vslope[0]);
		m_vDslope[v_order.at(i)]->SetName(("g_time_slope_"+v_order.at(i)).c_str());
		if (is_root_output){
			// std::cout << "  Writing in root file." << std::endl;
			f_out->cd(v_order.at(i).c_str());
			m_vDslope[v_order.at(i)]->Write();
		}

		for (int i_sect=0; i_sect<(int)dataTrees[v_order.at(i)]->GetMaximum("sectionId")+1;i_sect++){ ///<-----ERROR WAS HERE
			TGraph* gr = new TGraph(slave_Y_fixed[i_sect].size(),&(v_master_Y[i_sect][0]),&(slave_Y_fixed[i_sect][0]));
			gr->SetName(("g_fixed"+v_order.at(i)+"_Id_"+patch::to_string(i_sect)).c_str());
			gr->SetTitle(("graph of fixed dY vs Y from "+v_order.at(i)+" for section "+patch::to_string(i_sect)).c_str());
			TF1* func = new TF1("func","pol1",-Xlimit,Xlimit);
			gr->Fit(func,"RQ");
			if (is_root_output){
				// std::cout << "  Writing in root file." << std::endl;
				f_out->cd(v_order.at(i).c_str());
				gr->Write();
			}

			delete gr;
			delete func;
			slave_Y_fixed[i_sect].clear();
			v_master_Y[i_sect].clear();
		}

		m_vDtime[v_order.at(i)] = new TGraph(vtime.size(),&vtime[0],&vdriftVector[0]);
		m_vDtime[v_order.at(i)]->SetName(("g_vd_calibrated_time_"+v_order.at(i)).c_str());
		if (is_root_output){
			// std::cout << "  Writing in root file." << std::endl;
			f_out->cd(v_order.at(i).c_str());
			m_vDtime[v_order.at(i)]->Write();
		}

		m_vDnom[v_order.at(i)] = new TGraph(vtime.size(),&vtime[0],&vdriftNominal[0]);
		m_vDnom[v_order.at(i)]->SetName(("g_vd_reco_time_"+v_order.at(i)).c_str());
		if (is_root_output){
			// std::cout << "  Writing in root file." << std::endl;
			f_out->cd(v_order.at(i).c_str());
			m_vDnom[v_order.at(i)]->Write();
		}
    
		v_slope_prev.clear();
		v_slope_prev = v_slope;
		v_slope.clear();
		vtime.clear();
		vslope.clear();
		vdriftVector.clear();
		vdriftNominal.clear();
	}

	if (is_root_output) f_out->Close();

	std::cout << "\nAnalyzer finished successfully." << std::endl;
	return 0;
}
