#include <iostream>
#include <cstdio>
#include <TTree.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TF1.h>
#include <TStyle.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLatex.h>

using namespace std;

int validation(char *file="./mergedTrackMatchDump.root"){

	//char *file ="./mergedTrackMatchDump.root"; // place for input file name
	
	//char suffix[10] = "calib";          // Work with calibrated (calib - SHINE method; nominal - LEGACY method) and uncalibrated was done, suffixes were used in output plots names
	//char suffix[10] = "nominal";
	char suffix[10] = "dcs_raw";
	
	char MTPCL_TOFL_name[30], MTPCR_TOFR_name[30], VTPC2_MTPCL_name[30], VTPC2_MTPCR_name[30], VTPC1_VTPC2_name[30], GTPC_VTPC2_name[30], VTPC1_GTPC_name[30], beam_cut_check[50], no_beam_cut[50];
	sprintf(MTPCL_TOFL_name, "./MTPCL_TOFL_%s.png", suffix);
	sprintf(MTPCR_TOFR_name, "./MTPCR_TOFR_%s.png", suffix);
	sprintf(VTPC2_MTPCL_name, "./VTPC2_MTPCL_%s.png", suffix);
	sprintf(VTPC2_MTPCR_name, "./VTPC2_MTPCR_%s.png", suffix);
	sprintf(VTPC1_VTPC2_name, "./VTPC1_VTPC2_%s.png", suffix);
	sprintf(GTPC_VTPC2_name, "./GTPC_VTPC2_%s.png", suffix);
	sprintf(VTPC1_GTPC_name, "./VTPC1_GTPC_%s.png", suffix);
	sprintf(beam_cut_check, "./GTPC_beam_check_%s.png", suffix);
	sprintf(no_beam_cut, "./GTPC_no_beam_%s.png", suffix);
	
	TFile *tree_file = TFile::Open(file);
	if (tree_file == 0){
		cout<<"Tree not opened"<<endl;
     		return 1;
	}
	cout<<"Tree file opened"<<endl;
	tree_file->ls();
	
	TCanvas *window = new TCanvas("window", "window", 5,5,1400,800);
	window->Divide(2,1);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(1);
	int nbin = 200;
	//int beam = 1; //beam size cut [cm]
	double xmin1=-70.0, xmax1=70.0, xmin2=-50.0, xmax2=50.0, ymin1=-70.0, ymax1=70.0, ymin2=-50.0, ymax2=50.0, ymin3=-30.0, ymax3=30.0;
//--------------------- MTPCLvsTOFL ---------------------------------
	TTree *MyTree_MTPCLvsTOFL = new TTree();
	tree_file->GetObject("MTPCLvsTOFL", MyTree_MTPCLvsTOFL);
	if(MyTree_MTPCLvsTOFL == 0){
		cout<<"Tree pointer = NULL"<<endl;
		return 1;
	}
	cout<<"Tree opened"<<endl;
	window->cd(1);
	TH2D *hf1 = new TH2D("hf1","", nbin, xmin1,xmax1,nbin, ymin1,ymax1);
	hf1->SetTitle("MTPCL_y vs TOFL_y");
	hf1->SetYTitle("MTPCL_y");
	hf1->SetXTitle("TOFL_y");	
	MyTree_MTPCLvsTOFL->Draw("slave_Y:tofY>>hf1", "", "colz");

	window->cd(2);
	TH2D *hf2 = new TH2D("hf2", "", nbin, xmin1, xmax1, nbin, ymin3,ymax3);
	hf2->SetTitle("#Delta y vs TOFL_y");
	hf2->SetYTitle("TOFL_y - MTPCL_y");
	hf2->SetXTitle("TOFL_y");
	MyTree_MTPCLvsTOFL->Draw("(tofY-slave_Y):tofY>>hf2", "", "colz");
	TF1 *f1 = new TF1("f1", "[0]+[1]*x", xmin1, xmin2);
	hf2->Fit(f1);
	hf2->Draw("colz");

	window->Print(MTPCL_TOFL_name);

	delete MyTree_MTPCLvsTOFL;
	delete hf1;
	delete hf2;
	delete f1;

//--------------------- MTPCRvsTOFR ---------------------------------
	TTree *MyTree_MTPCRvsTOFR = new TTree();
	tree_file->GetObject("MTPCRvsTOFR", MyTree_MTPCRvsTOFR);
	if(MyTree_MTPCRvsTOFR == 0){
		cout<<"Tree pointer = NULL"<<endl;
		return 1;
	}
	cout<<"Tree opened"<<endl;

	window->cd(1);
	TH2D *hf3 = new TH2D("hf3", "", nbin, xmin1,  xmax1, nbin, ymin1, ymax1);
	hf3->SetTitle("MTPCR_y vs TOFR_y");
	hf3->SetYTitle("MTPCR_y");
	hf3->SetXTitle("TOFR_y");
	MyTree_MTPCRvsTOFR->Draw("slave_Y:tofY>>hf3", "", "colz");	
	
	window->cd(2);
	TH2D *hf4 = new TH2D("hf4", "", nbin, xmin1, xmax1, nbin, ymin3,ymax3);
	hf4->SetTitle("#Delta y vs TOFR_y");
	hf4->SetYTitle("TOFR_y - MTPCR_y");
	hf4->SetXTitle("TOFR_y");
	MyTree_MTPCRvsTOFR->Draw("(tofY-slave_Y):tofY>>hf4", "", "colz");
	TF1 *f2 = new TF1("f2", "[0]+[1]*x", xmin1, xmax1);
	hf4->Fit(f2);
	hf4->Draw("colz");

	window->Print(MTPCR_TOFR_name);

	delete MyTree_MTPCRvsTOFR;
	delete hf3;
	delete hf4;
	delete f2;

//--------------------- VTPC2vsMTPCL ---------------------------------
	TTree *MyTree_VTPC2vsMTPCL = new TTree();
	tree_file->GetObject("VTPC2vsMTPCL", MyTree_VTPC2vsMTPCL);
	if(MyTree_VTPC2vsMTPCL == 0){
		cout<<"Tree pointer = NULL"<<endl;
		return 1;
	}
	cout<<"Tree opened"<<endl;

	window->cd(1);
	TH2D *hf5 = new TH2D("hf5", "", nbin, xmin1, xmax1, nbin, ymin1, ymax1);
	hf5->SetTitle("VTPC2_y vs MTPCL_y");
	hf5->SetYTitle("VTPC2_y");
	hf5->SetXTitle("MTPCL_y");
	MyTree_VTPC2vsMTPCL->Draw("slave_Y:master_Y>>hf5", "", "colz");

	window->cd(2);
	TH2D *hf6 = new TH2D("hf6", "", nbin, xmin1, xmax1, nbin, ymin3, ymax3);
	hf6->SetTitle("#Delta y vs MTPCL_y");
	hf6->SetYTitle("MTPCL_y - VTPC2_y");
	hf6->SetXTitle("MTPCL_y");
	MyTree_VTPC2vsMTPCL->Draw("(master_Y-slave_Y):master_Y>>hf6", "", "colz");
	TF1 *f3 = new TF1("f3", "[0]+[1]*x", xmin1, xmax1);
	hf6->Fit(f3);
	hf6->Draw("colz");
	
	window->Print(VTPC2_MTPCL_name);

	delete MyTree_VTPC2vsMTPCL;
	delete hf5;
	delete hf6;
	delete f3;

//--------------------- VTPC2vsMTPCR ---------------------------------
	TTree *MyTree_VTPC2vsMTPCR = new TTree();
	tree_file->GetObject("VTPC2vsMTPCR", MyTree_VTPC2vsMTPCR);
	if(MyTree_VTPC2vsMTPCR == 0){
		cout<<"Tree pointer = NULL"<<endl;
		return 1;
	}
	cout<<"Tree opened"<<endl;

	window->cd(1);
	TH2D *hf7 = new TH2D("hf7", "", nbin, xmin1, xmax1, nbin, ymin1, ymax1);
	hf7->SetTitle("VTPC2_y vs MTPCR_y");
	hf7->SetYTitle("VTPC2_y");
	hf7->SetXTitle("MTPCR_y");
	MyTree_VTPC2vsMTPCR->Draw("slave_Y:master_Y>>hf7", "", "colz");

	window->cd(2);
	TH2D *hf8 = new TH2D("hf8", "", nbin, xmin1, xmax1, nbin, ymin3, ymax3);
	hf8->SetTitle("#Delta y vs MTPCR_y");
	hf8->SetYTitle("MTPCR_y - VTPC2_y");
	hf8->SetXTitle("MTPCR_y");
	MyTree_VTPC2vsMTPCR->Draw("(master_Y-slave_Y):master_Y>>hf8", "", "colz");
	TF1 *f4 = new TF1("f4", "[0]+[1]*x", xmin1, xmax1);
	hf8->Fit(f4);
	hf8->Draw("colz");
	
	window->Print(VTPC2_MTPCR_name);

	delete MyTree_VTPC2vsMTPCR;
	delete hf7;
	delete hf8;
	delete f4;

//--------------------- VTPC1vsVTPC2 ---------------------------------
	TTree *MyTree_VTPC1vsVTPC2 = new TTree();
	tree_file->GetObject("VTPC1vsVTPC2", MyTree_VTPC1vsVTPC2);
	if(MyTree_VTPC1vsVTPC2 == 0){
		cout<<"Tree pointer = NULL"<<endl;
		return 1;
	}
	cout<<"Tree opened"<<endl;

	window->cd(1);	
	TH2D *hf9 = new TH2D("hf9", "", nbin, xmin2, xmax2, nbin, ymin2, ymax2);
	hf9->SetTitle("VTPC1_y vs VTPC2_y");
	hf9->SetYTitle("VTPC1_y");
	hf9->SetXTitle("VTPC2_y");
	MyTree_VTPC1vsVTPC2->Draw("slave_Y:master_Y>>hf9", "", "colz");

	window->cd(2);
	TH2D *hf10 = new TH2D("hf10", "", nbin, xmin2, xmax2, nbin, ymin3, ymax3);
	hf10->SetTitle("#Delta y vs VTPC2_y");
	hf10->SetYTitle("VTPC2_y - VTPC1_y");
	hf10->SetXTitle("VTPC2_y");
	MyTree_VTPC1vsVTPC2->Draw("(master_Y-slave_Y):master_Y>>hf10", "", "colz");
	TF1 *f5 = new TF1("f5", "[0]+[1]*x", xmin1, xmax1);
	hf10->Fit(f5);
	hf10->Draw("colz");
	
	window->Print(VTPC1_VTPC2_name);

	delete MyTree_VTPC1vsVTPC2;
	delete hf9;
	delete hf10;
	delete f5;

//--------------------- GTPCvsVTPC2 with beam spot cut -------------------
	TTree *MyTree_GTPCvsVTPC2 = new TTree();
	tree_file->GetObject("GTPCvsVTPC2", MyTree_GTPCvsVTPC2);
	if(MyTree_GTPCvsVTPC2 == 0){
		cout<<"Tree pointer = NULL"<<endl;
		return 1;
	}
	cout<<"Tree opened"<<endl;

	window->cd(1);
	TH2D *hf11 = new TH2D("hf11", "", nbin, xmin2, xmax2, nbin, ymin2, ymax2);
	hf11->SetTitle("GTPC_y vs VTPC2_y");
	hf11->SetYTitle("GTPC_y");
	hf11->SetXTitle("VTPC2_y");
	MyTree_GTPCvsVTPC2->Draw("slave_Y:master_Y>>hf11", "!(slave_X>-2 && slave_X<2 && slave_Y>-2 && slave_Y<2) && master_X>-20 && master_X<20 && master_Y>-20 && master_Y<20 && !(master_X>-2 && master_X<2 && master_Y>-2 && master_Y<2)", "colz");

	window->cd(2);
	TH2D *hf12 = new TH2D("hf12", "", nbin, xmin2, xmax2, nbin, ymin3, ymax3);
	hf12->SetTitle("#Delta y vs VTPC2_y");
	hf12->SetYTitle("VTPC2_y - GTPC_y");
	hf12->SetXTitle("VTPC2_y");
	MyTree_GTPCvsVTPC2->Draw("(master_Y-slave_Y):master_Y>>hf12", "!(slave_X>-2 && slave_X<2 && slave_Y>-2 && slave_Y<2) && master_X>-20 && master_X<20 && master_Y>-20 && master_Y<20 && !(master_X>-2 && master_X<2 && master_Y>-2 && master_Y<2)", "colz");
	TF1 *f6 = new TF1("f6", "[0]+[1]*x", xmin2, xmax2);
	hf12->Fit(f6);
	hf12->Draw("colz");
	
	window->Print(GTPC_VTPC2_name);
//----------- beam spot cut -------------------------
	window->cd(1);
	TH2D *hfx = new TH2D("hfx", "", nbin/2, 20, 20, nbin/2, ymin3, ymax3);
	hfx->SetTitle("GTPC_y vs GTPC_x with beam cut");
	hfx->SetYTitle("GTPC_y");
	hfx->SetXTitle("GTPC_x");
	MyTree_GTPCvsVTPC2->Draw("slave_Y:slave_X>>hfx", "!(slave_X>-2 && slave_X<2 && slave_Y>-2 && slave_Y<2) && master_X>-20 && master_X<20 && master_Y>-20 && master_Y<20 && !(master_X>-2 && master_X<2 && master_Y>-2 && master_Y<2)", "colz");
	window->Print(beam_cut_check);
	
//---------- beam spot not cut ----------------------
	window->cd(1);
	TH2D *hfx2 = new TH2D("hfx2", "", nbin/2, 20, 20, nbin/2, ymin3, ymax3);
	hfx2->SetTitle("GTPC_y vs GTPC_x no beam cut");
	hfx2->SetYTitle("GTPC_y");
	hfx2->SetXTitle("GTPC_x");
	MyTree_GTPCvsVTPC2->Draw("slave_Y:slave_X>>hfx2", "", "colz");
	window->cd(2);
	TH2D *hf122 = new TH2D("hf122", "", nbin, xmin2, xmax2, nbin, ymin3, ymax3);
	hf122->SetTitle("#Delta y vs VTPC2_y");
	hf122->SetYTitle("VTPC2_y - GTPC_y");
	hf122->SetXTitle("VTPC2_y");
	MyTree_GTPCvsVTPC2->Draw("(master_Y-slave_Y):master_Y>>hf122", "", "colz");
	hf122->Draw("colz");
	window->Print(no_beam_cut);
	
	delete MyTree_GTPCvsVTPC2;
	delete hf11;
	delete hf12;
	delete hfx;
	delete hfx2;
	delete hf122;
	delete f6;

//--------------------- VTPC1vsGTPC ---------------------------------
	TTree *MyTree_VTPC1vsGTPC = new TTree();
	tree_file->GetObject("VTPC1vsGTPC", MyTree_VTPC1vsGTPC);
	if(MyTree_VTPC1vsGTPC == 0){
		cout<<"Tree pointer = NULL"<<endl;
		return 1;
	}
	cout<<"Tree opened"<<endl;

	window->cd(1);
	TH2D *hf13 = new TH2D("hf13", "", nbin, xmin2, xmax2, nbin, ymin2, ymax2);
	hf13->SetTitle("VTPC1_y vs GTPC_y");
	hf13->SetYTitle("VTPC1_y");
	hf13->SetXTitle("GTPC_y");
	MyTree_VTPC1vsGTPC->Draw("slave_Y:master_Y>>hf13", "", "colz");
	
	window->cd(2);
	TH2D *hf14 = new TH2D("hf14", "", nbin, xmin2, xmax2, nbin, ymin3, ymax3);
	hf14->SetTitle("#Delta y vs GTPC_y");
	hf14->SetYTitle("GTPC_y - VTPC1_y");
	hf14->SetXTitle("GTPC_y");
	MyTree_VTPC1vsGTPC->Draw("(master_Y-slave_Y):master_Y>>hf14", "", "colz");
	TF1 *f7 = new TF1("f7", "[0]+[1]*x", xmin2, xmax2);
	hf14->Fit(f7);
	hf14->Draw("colz");

	window->Print(VTPC1_GTPC_name);

	delete MyTree_VTPC1vsGTPC;
	delete hf13;
	delete hf14;
	delete f7;

	return 0;
}
