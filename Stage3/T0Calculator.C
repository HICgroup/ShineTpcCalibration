/**
  \file
  Analyzer script based on TPC bottom points in order to determine
  globalT0+chamberT0 timing offset, knowing full drift time.
  Based on equation
    (globalT0+chamberT0 + bottomPointRawDriftTime)
      = fullDriftTime
  Here fullDriftTime = driftLength/driftVelocity is assumed to be known, 
  while the bottomPointRawDriftTime is the measured raw drift times of 
  tracks exiting at TPC drift cathode plane, corrected for the known 
  padByPadT0 timing shifts (determined by calibration pulser), and corrected 
  for eventByEventT0 timing shifts which is measured.
  Thus, determines globalT0+chamberT0.

  \author A. Laszlo
  \version $Id:    $
  \date 01 Apr 2014
*/

#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TObject.h>
#include <TClass.h>
#include <TKey.h>
#include <TMath.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;


// Output file name for global calibration parameters.
const string globalT0ValuesFileName = "GlobalT0.xml";
// Output file name for chamber-by-chamber calibration parameters.
const string chamberT0ValuesFileName = "ChamberT0.xml";
// Output file name for fit control plot.
const string t0PlotFileName = "ChamberT0ControlPlots.pdf";
// Reference chamberT0 of VTPC1 in microsecond (fixed by convention).
const double vtpc1ChamberT0InUSec = -0.0946;


/// Structure for bottom point analysis.
struct TPCData {
  // Obligatory header to identify event.
  UInt_t fRunNumber;
  UInt_t fSpillId;
  UInt_t fEventNumber;
  UInt_t fEventUnixTime;
  Int_t fSectionId;
  // Last track cluster location parameters.
  Int_t fLastClusterSectorNumber;
  Int_t fLastClusterPadrowNumber;
  Float_t fLastClusterPadrowXInCM;
  // Last track cluster raw drift time for bottom determination.
  Float_t fLastClusterRawDriftTimeInUSec;
  // Last extrapolated track point raw drift time for bottom determination.
  Float_t fLastPointRawDriftTimeInUSec;
  // Full drift time in usec.
  Float_t fFullDriftTimeInUSec;
};


// Global parameters.
const Float_t histogramMargin = 0.25;
const Float_t refinement = 6.0;
const Int_t nMaxTSlices = 512;
const Float_t minTSlicesInUSec = 0.100;


// Gauss function.
Float_t GaussF(Double_t *x, Double_t *par)
{
  return par[0]*TMath::Exp(-(x[0]-par[1])*(x[0]-par[1])/(2.0*par[2]*par[2]));
}


// Display usage.
void DisplayUsage()
{
  cerr << "\nUsage:\n\tT0Calculator -i rootTTreeFile\n" << endl;
  exit(-1);
}


/// Main function.
int main(int argc, char* argv[])
{
  int exitCode = 0;
  const Float_t driftMeasuringTimeInUSec = nMaxTSlices*minTSlicesInUSec;
  const Float_t tMinInUSec = -histogramMargin*driftMeasuringTimeInUSec;
  const Float_t tMaxInUSec = driftMeasuringTimeInUSec + histogramMargin*driftMeasuringTimeInUSec;
  const Int_t nTBins = lround(nMaxTSlices*refinement*(1.0+histogramMargin+histogramMargin));
  const Float_t tBinWidthInUSec = minTSlicesInUSec/refinement;

  // Parse arguments.
  string inputFileName = "";
  for ( int i=1 ; i<argc ; ++i ) {
    string argValue = argv[i];
    if ( argValue=="-i" ) {
      if ( i+1<argc ) {
        inputFileName = argv[i+1];
        ++i;
      } else {
        cerr << "\nNo file name given with argument -i" << "\n";
        DisplayUsage();
      }
    } else if ( argValue=="-h" ) {
      DisplayUsage();
    } else {
      cerr << "\nInvalid argument " << argv[i] << "\n";
      DisplayUsage();
    }
  }
  if ( inputFileName == "" ) {
    cerr << "\nNo input files specified" << "\n";
    DisplayUsage();
  }

  // Open input file.
  TFile inputFile(inputFileName.c_str());
  vector<TTree*> ttreeArray;
  vector<TPCData> tpcData;
  vector<string> tpcName;
  vector<TH1F*> tpcT0Histogram;
  vector<TF1*> tpcT0Fit;
  vector<double> summedT0Values;
  vector<bool> success;
  int runNumberMin = 0;
  int runNumberMax = 0;

  // Get the TTree objects from input file.
  const TList& listOfKeys = *inputFile.GetListOfKeys();
  TIter nextKey(&listOfKeys);
  TKey* key = 0;
  while ( ( key = (TKey*)nextKey() ) ) {
    const TObject* object = key->ReadObj();
    const string name = key->GetName();
    if ( object->IsA()->InheritsFrom(TTree::Class()) ) {
      TTree* ttree = (TTree*)object;
      ttreeArray.push_back(ttree);
      tpcName.push_back(name);
      cout << "Found TTree: " << name << " in " << inputFileName << endl;
    } else {
      cout << "Found non-TTree object " << name << " in " << inputFileName << endl;
    }
  }
  if ( exitCode ) {
    inputFile.Close();
    return exitCode;
  }
  tpcData.resize(ttreeArray.size());
  for ( unsigned int tpcIndex = 0 ; tpcIndex<tpcData.size() ; ++ tpcIndex ) {
    ttreeArray[tpcIndex]->SetBranchAddress("runNumber", &(tpcData[tpcIndex].fRunNumber));
    ttreeArray[tpcIndex]->SetBranchAddress("spillId", &(tpcData[tpcIndex].fSpillId));
    ttreeArray[tpcIndex]->SetBranchAddress("eventNumber", &(tpcData[tpcIndex].fEventNumber));
    ttreeArray[tpcIndex]->SetBranchAddress("eventUnixTime", &(tpcData[tpcIndex].fEventUnixTime));
    ttreeArray[tpcIndex]->SetBranchAddress("sectionId", &(tpcData[tpcIndex].fSectionId));
    ttreeArray[tpcIndex]->SetBranchAddress("lastClusterSectorNumber", &(tpcData[tpcIndex].fLastClusterSectorNumber));
    ttreeArray[tpcIndex]->SetBranchAddress("lastClusterPadrowNumber", &(tpcData[tpcIndex].fLastClusterPadrowNumber));
    ttreeArray[tpcIndex]->SetBranchAddress("lastClusterPadrowXInCM", &(tpcData[tpcIndex].fLastClusterPadrowXInCM));
    ttreeArray[tpcIndex]->SetBranchAddress("lastClusterRawDriftTimeInUSec", &(tpcData[tpcIndex].fLastClusterRawDriftTimeInUSec));
    ttreeArray[tpcIndex]->SetBranchAddress("lastPointRawDriftTimeInUSec", &(tpcData[tpcIndex].fLastPointRawDriftTimeInUSec));
    ttreeArray[tpcIndex]->SetBranchAddress("fullDriftTimeInUSec", &(tpcData[tpcIndex].fFullDriftTimeInUSec));
  }

  // Process the data to fill histograms.
  tpcT0Histogram.resize(tpcData.size(), 0);
  for ( unsigned int tpcIndex = 0 ; tpcIndex<tpcData.size() ; ++tpcIndex ) {
    tpcT0Histogram[tpcIndex] =
      new TH1F((string("tpcT0_")+tpcName[tpcIndex]).c_str(),
               (string("tpcT0_")+tpcName[tpcIndex]
               + ";TPC globalT0+chamberT0 [us];Entries").c_str(),
               nTBins, driftMeasuringTimeInUSec-tMaxInUSec, driftMeasuringTimeInUSec-tMinInUSec);
    const Int_t nEntries = (Int_t)ttreeArray[tpcIndex]->GetEntries();
    for ( Int_t entryIndex = 0 ; entryIndex<nEntries ; ++entryIndex ) {
      ttreeArray[tpcIndex]->GetEntry(entryIndex);
      if ( runNumberMin == 0 )
        runNumberMin = tpcData[tpcIndex].fRunNumber;
      if ( runNumberMax == 0 )
        runNumberMax = tpcData[tpcIndex].fRunNumber;
      if ( tpcData[tpcIndex].fRunNumber < runNumberMin )
        runNumberMin = tpcData[tpcIndex].fRunNumber;
      if ( runNumberMax < tpcData[tpcIndex].fRunNumber )
        runNumberMax = tpcData[tpcIndex].fRunNumber;
      const Float_t lastPointRawDriftTimeInUSec = tpcData[tpcIndex].fLastPointRawDriftTimeInUSec;
      const Float_t fullDriftTimeInUSec = tpcData[tpcIndex].fFullDriftTimeInUSec;
      tpcT0Histogram[tpcIndex]->Fill(fullDriftTimeInUSec-lastPointRawDriftTimeInUSec);
    }
  }

  // Fit the acquired histograms.
  tpcT0Fit.resize(tpcT0Histogram.size(), 0);
  summedT0Values.resize(tpcT0Histogram.size(), 0.0);
  success.resize(tpcT0Histogram.size(), false);
  for ( unsigned int tpcIndex = 0 ; tpcIndex<tpcT0Histogram.size() ; ++tpcIndex ) {
    TH1F& h = *tpcT0Histogram[tpcIndex];
    const Int_t nBins = h.GetNbinsX();
    Float_t maxValue = h.GetXaxis()->GetBinCenter(1);
    Float_t maxEntries = h.GetBinContent(1);
    for ( Int_t bin = 1 ; bin<=nBins ; ++bin ) {
      const Float_t value = h.GetXaxis()->GetBinCenter(bin);
      const Float_t entries = h.GetBinContent(bin);
      if ( maxEntries < entries ) {
        maxValue = value;
        maxEntries = entries;
      }
    }
    if ( maxEntries > 10 ) {
      const Float_t xMin = maxValue-minTSlicesInUSec*6.0;
      const Float_t xMax = maxValue+minTSlicesInUSec*6.0;
      const Int_t nPars = 3;
      const string funcName = string("Gauss_")+tpcName[tpcIndex];
      tpcT0Fit[tpcIndex] =
        new TF1(funcName.c_str(),
          GaussF,
          xMin,
          xMax,
          nPars);
      tpcT0Fit[tpcIndex]->SetParameter(0, maxEntries);
      tpcT0Fit[tpcIndex]->SetParameter(1, maxValue);
      tpcT0Fit[tpcIndex]->SetParameter(2, minTSlicesInUSec*4.0);
      tpcT0Fit[tpcIndex]->SetParNames("maxEntries", "mean", "sigma");
      h.Fit(funcName.c_str(), "ILL", "", xMin, xMax);
      summedT0Values[tpcIndex] = tpcT0Fit[tpcIndex]->GetParameter(1);
      success[tpcIndex] = true;
      cout << tpcName[tpcIndex] << " "
           << tpcT0Fit[tpcIndex]->GetParameter(1) << " +/- " << tpcT0Fit[tpcIndex]->GetParError(1) << " usec"
           << endl;
    } else {
      summedT0Values[tpcIndex] = 0.0;
      success[tpcIndex] = false;
      cout << tpcName[tpcIndex] << " "
           << 0.0 << " +/- " << 0.0 << " usec"
           << endl;
    }
  }

  // Calculate globalT0 and write out calibration data.
  double globalT0InUSec = 0.0;
  for ( unsigned int tpcIndex = 0 ; tpcIndex<tpcT0Histogram.size() ; ++tpcIndex ) {
    if ( tpcName[tpcIndex] == "VTPC1" && success[tpcIndex] ) {
      globalT0InUSec = summedT0Values[tpcIndex] - vtpc1ChamberT0InUSec;
      break;
    }
  }
  ofstream globalT0ValuesFile;
  globalT0ValuesFile.open(globalT0ValuesFileName.c_str());
  globalT0ValuesFile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
                     << "\n"
                     << "<GlobalT0\n"
                     << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
                     << "  xsi:noNamespaceSchemaLocation=\"[SCHEMAPATH]/TPCGlobalT0_DBFormat.xsd\"\n"
                     << "  t0unit=\"nanosecond\">\n"
                     << "\n"
                     << "  <Entry runStart=\"" << runNumberMin << "\" runStop=\"" << runNumberMax << "\"> " << globalT0InUSec*1000<< " </Entry>\n"
                     << "\n"
                     << "</GlobalT0>\n"
                     << "\n";
  globalT0ValuesFile.close();
  globalT0ValuesFile.clear();

  // Calculate chamberT0 and write out calibration data.
  ofstream chamberT0ValuesFile;
  chamberT0ValuesFile.open(chamberT0ValuesFileName.c_str());
  chamberT0ValuesFile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
                      << "\n"
                      << "<ChamberT0List\n"
                      << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
                      << "  xsi:noNamespaceSchemaLocation=\"[SCHEMAPATH]/TPCChamberT0_DBFormat.xsd\">\n"
                      << "\n";
  chamberT0ValuesFile << "  <ChamberT0 t0unit=\"nanosecond\" runStart=\"" << runNumberMin << "\" runStop=\"" << runNumberMax << "\">\n";
  for ( unsigned int tpcIndex = 0 ; tpcIndex<tpcT0Histogram.size() ; ++tpcIndex ) {
    if ( success[tpcIndex] )
      chamberT0ValuesFile << "    <TPC name=\"" << tpcName[tpcIndex] << "\"> " << (summedT0Values[tpcIndex]-globalT0InUSec)*1000 << " </TPC>\n";
    else
      chamberT0ValuesFile << "    <TPC name=\"" << tpcName[tpcIndex] << "\"> " << 0.0 << " </TPC>\n";
  }
  chamberT0ValuesFile << "  </ChamberT0>\n"
                      << "\n";
  chamberT0ValuesFile << "</ChamberT0List>\n"
                      << "\n";
  chamberT0ValuesFile.close();
  chamberT0ValuesFile.clear();

  // Create control plots.
  gStyle->SetOptStat(kFALSE);
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  c->Clear();
  for ( unsigned int tpcIndex = 0 ; tpcIndex<tpcT0Histogram.size() ; ++tpcIndex ) {
    TH1F& h = *tpcT0Histogram[tpcIndex];
    TF1& f = *tpcT0Fit[tpcIndex];
    h.Draw();
    if ( tpcIndex == 0 )
      c->Print((t0PlotFileName+"(").c_str());
    else if ( tpcIndex+1 == tpcT0Histogram.size() )
      c->Print((t0PlotFileName+")").c_str());
    else
      c->Print(t0PlotFileName.c_str());
  }

  // Close input file.
  inputFile.Close();

  return exitCode;
}
