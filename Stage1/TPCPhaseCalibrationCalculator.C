/**
  \file
  Analyzer script based in order to determine calibration parameters 
  of TPC phase shift. Based on the knowledge that TPC phase TDC has 
  uniform distribution, whose support corresponds to the 
  [0ns, tpcPrimaryClockTimeWidth] time interval.

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
#include <TCanvas.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
using namespace std;


// Output file name for calibration parameters.
const string tpcPhaseCalibValuesFileName = "TPCPhaseCalibration.xml";
// Output file name for fit control plot.
const string tpcPasePlotFileName = "TPCPhaseCalibrationControlPlots.pdf";
const string tpcPaseCanvFileName = "TPCPhaseCalibrationControlPlots.C";

// Global parameters.
const Int_t dynamicalRange = 1023;


// Display usage.
void DisplayUsage()
{
  cerr << "\nUsage:\n\tTPCPhaseCalibrationCalculator -i rootTTreeFile\n" << endl;
  exit(-1);
}


/// Main function.
int main(int argc, char* argv[])
{
  int exitCode = 0;
  // Containers for storing TDC info for TPC phase.
  // Obligatory header to identify event.
  UInt_t runNumber;
  UInt_t spillId;
  UInt_t eventNumber;
  UInt_t eventUnixTime;
  Int_t sectionId;
  // Measured TPC phase in TDC units.
  Int_t tpcPhaseTDC;
  // Primary clock tick width in nanosec.
  Float_t tpcPrimaryClockTimeWidthInNanoSec;

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
  TH1F* tdcHistogram = 0;
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
      cout << "Found TTree: " << name << " in " << inputFileName << endl;
    } else {
      cout << "Found non-TTree object " << name << " in " << inputFileName << endl;
    }
  }
  if ( exitCode ) {
    inputFile.Close();
    return exitCode;
  }
  if ( ttreeArray.size() != 1 ) {
    cerr << "Exactly 1 TTree is expected.\n";
    exitCode = -3;
    inputFile.Close();
    return exitCode;
  }
  ttreeArray[0]->SetBranchAddress("runNumber", &runNumber);
  ttreeArray[0]->SetBranchAddress("spillId", &spillId);
  ttreeArray[0]->SetBranchAddress("eventNumber", &eventNumber);
  ttreeArray[0]->SetBranchAddress("eventUnixTime", &eventUnixTime);
  ttreeArray[0]->SetBranchAddress("sectionId", &sectionId);
  ttreeArray[0]->SetBranchAddress("tpcPhaseTDC", &tpcPhaseTDC);
  ttreeArray[0]->SetBranchAddress("tpcPrimaryClockTimeWidthInNanoSec", &tpcPrimaryClockTimeWidthInNanoSec);

  // Process the data to fill histograms.
  tdcHistogram = new TH1F("tpcPhaseTDCSpectrum",
                          (string("tpcPhaseTDCSpectrum") + ";" + "TPC phase shift [TDC]" + ";" + "Entries").c_str(),
                          dynamicalRange+1, -0.5, dynamicalRange+0.5);
  const Int_t nEntries = (Int_t)ttreeArray[0]->GetEntries();
  for ( Int_t entryIndex = 0 ; entryIndex<nEntries ; ++entryIndex ) {
    ttreeArray[0]->GetEntry(entryIndex);
    if ( entryIndex == 0 ) {
      runNumberMin = runNumber;
    }
    runNumberMax = runNumber;
    tdcHistogram->Fill(tpcPhaseTDC);
  }

  // Process the acquired histogram and calculate calibration data.
  double lowestNonZeroTDC = 0;
  for ( int i=1+1 ; i<=tdcHistogram->GetNbinsX()+1 ; ++i ) {
    const double entries = tdcHistogram->GetBinContent(i);
    if ( entries > 0 ) {
      lowestNonZeroTDC = tdcHistogram->GetBinCenter(i);
    break;
    }
  }
  double highestNonZeroTDC = 0;
  for ( int i=tdcHistogram->GetNbinsX()+1-1-1 ; i>=0 ; --i ) {
    const double entries = tdcHistogram->GetBinContent(i);
    if ( entries > 0 ) {
      highestNonZeroTDC = tdcHistogram->GetBinCenter(i);
      break;
    }
  }

  // Write TDC offset and time span of a TDC count into file.
  ofstream outputFile;
  outputFile.open(tpcPhaseCalibValuesFileName.c_str());
  outputFile << "\n"
             << "<!-- TDCOffset is a dimensionless TDC count, others are in time units. -->\n"
             << "\n"
             << "<TPCPhaseCalibrationList\n"
             << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
             << "  xsi:noNamespaceSchemaLocation=\"[SCHEMAPATH]/TPCPhaseCalibration_DBFormat.xsd\">\n"
             << "\n"
             << "  <TPCPhaseCalibration timeUnit=\"nanosecond\" runStart=\"" << runNumberMin << "\" runStop=\"" << runNumberMax << "\">\n"
             << "    <TDCOffset> " << lowestNonZeroTDC << " </TDCOffset>\n"
             << "    <TDCCountTimeSpan> " << tpcPrimaryClockTimeWidthInNanoSec/(highestNonZeroTDC-lowestNonZeroTDC+1) << " </TDCCountTimeSpan>\n"
             << "    <TimeShift> " << 0.0 << " </TimeShift>\n"
             << "  </TPCPhaseCalibration>\n"
             << "\n"
             << "</TPCPhaseCalibrationList>\n"
             << "\n";
  outputFile.close();
  outputFile.clear();

  // Create control plots.
  gStyle->SetOptStat(kFALSE);
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  c->Clear();
  tdcHistogram->Draw();
  c->Print(tpcPasePlotFileName.c_str());
  c->SaveAs(tpcPaseCanvFileName.c_str());

  // Close input file.
  inputFile.Close();

  return exitCode;
}
