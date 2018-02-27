/**
   \file
  Analyzer script in order to determine calibration parameters 
  MHTDC for MTacc and corresponding S11 matching. It is based on 
  calculation of most probable MTacc - S11 time difference.

   \author A. Laszlo
   \version $Id: MHTDCCalibratorAL.cc 9089 2013-03-26 16:46:48Z munger $
   \date 17 Feb 2012
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
const string mhtdcCalibValuesFileName = "MHTDCCalibration.xml";
// Output file name for fit control plot.
const string mhtdcPlotFileName = "MHTDCCalibrationControlPlots.pdf";
const string mhtdcCanvFileName = "MHTDCCalibrationControlPlots.C";

// Acceptance window of S11 in nanosecond.
const double mhtdcS11TimeWindowInNanoSec = 1.0;


// Display usage.
void DisplayUsage()
{
  cerr << "\nUsage:\n\tMHTDCCalibrationCalculator -i rootTTreeFile\n" << endl;
  exit(-1);
}


/// Main function.
int main(int argc, char* argv[])
{
  int exitCode = 0;
  // Containers for storing TDC info for MHTDC.
  // Obligatory header to identify event.
  UInt_t runNumber;
  UInt_t spillId;
  UInt_t eventNumber;
  UInt_t eventUnixTime;
  Int_t sectionId;
  // Measured MTacc - S11 difference in TDC units.
  Int_t diffMTaccToS11InTDC;
  // TDC unit in nanosec.
  Float_t mhtdcTimeBinWidthInNanoSec;
  // MHTDC dynamic range.
  Int_t mhtdcDynamicRange;

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
  ttreeArray[0]->SetBranchAddress("diffMTaccToS11InTDC", &diffMTaccToS11InTDC);
  ttreeArray[0]->SetBranchAddress("mhtdcTimeBinWidthInNanoSec", &mhtdcTimeBinWidthInNanoSec);
  ttreeArray[0]->SetBranchAddress("mhtdcDynamicRange", &mhtdcDynamicRange);

  // Process the data to fill histograms.
  const Int_t nEntries = (Int_t)ttreeArray[0]->GetEntries();
  for ( Int_t entryIndex = 0 ; entryIndex<nEntries ; ++entryIndex ) {
    ttreeArray[0]->GetEntry(entryIndex);
    if ( entryIndex == 0 ) {
      runNumberMin = runNumber;
      tdcHistogram = new TH1F("mhtdcSpectrumS11",
                              (string("mhtdcSpectrumS11") + ";" + "MTacc - S11 time [nanosecond]" + ";" + "Entries").c_str(),
                              2*mhtdcDynamicRange+1, (-int(mhtdcDynamicRange)-0.5), (int(mhtdcDynamicRange)+0.5));
    }
    runNumberMax = runNumber;
    tdcHistogram->Fill(diffMTaccToS11InTDC);
  }

  // Find the most probable MTacc to S11 correlation.
  double maxEntries = 0.0;
  double maxPosition = 0.0;
  if ( tdcHistogram != 0 ) {
    for ( int i=1+1 ; i<=tdcHistogram->GetNbinsX()+1-1-1 ; ++i ) {
      const double entries = tdcHistogram->GetBinContent(i);
      if ( entries > maxEntries ) {
        const double position = tdcHistogram->GetBinCenter(i);
        maxEntries = entries;
        maxPosition = position;
      }
    }
  }

  // Write calibration data into file.
  ofstream outputFile;
  outputFile.open(mhtdcCalibValuesFileName.c_str());
  outputFile << "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n"
             << "\n"
             << "<MHTDCCalibrationList\n"
             << "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
             << "  xsi:noNamespaceSchemaLocation=\"[SCHEMAPATH]/MHTDCCalibration_DBFormat.xsd\">\n"
             << "\n"
             << "  <MHTDCCalibration timeUnit=\"nanosecond\" runStart=\"" << runNumberMin << "\" runStop=\"" << runNumberMax << "\">\n"
             << "    <MTaccMinusS11Time> " << maxPosition*mhtdcTimeBinWidthInNanoSec << " </MTaccMinusS11Time>\n"
             << "    <S11AcceptanceTimeWindow> " << mhtdcS11TimeWindowInNanoSec << " </S11AcceptanceTimeWindow>\n"
             << "  </MHTDCCalibration>\n"
             << "\n"
             << "</MHTDCCalibrationList>\n"
             << "\n";
  outputFile.close();
  outputFile.clear();

  // Create control plots.
  gStyle->SetOptStat(kFALSE);
  TCanvas* c = new TCanvas("c", "c", 600, 600);
  c->Clear();
  tdcHistogram->Draw();
  c->Print(mhtdcPlotFileName.c_str());
  c->SaveAs(mhtdcCanvFileName.c_str());

  // Close input file.
  inputFile.Close();

  return exitCode;
}
