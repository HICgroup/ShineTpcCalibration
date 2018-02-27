/**
   \file
   implementation of class MHTDCCalibratorAL

   \author A. Laszlo
   \version $Id: MHTDCCalibratorAL.cc 10683 2014-04-22 21:17:27Z laszlo $
   \date 17 Feb 2012
*/


#include "MHTDCCalibratorAL.h"

#include <fwk/CentralConfig.h>
#include <utl/ErrorLogger.h>
#include <utl/UTCDateTime.h>

#include <evt/Event.h>

#include <sstream>
#include <fstream>

#include <TFile.h>


using namespace fwk;
using namespace evt;
using namespace utl;
using namespace std;


namespace MHTDCCalibratorAL {


  /// Init function for MHTDCCalibratorAL
  VModule::EResultFlag
  MHTDCCalibratorAL::Init()
  {
    CentralConfig& cc = CentralConfig::GetInstance();

    Branch topBranch = cc.GetTopBranch("MHTDCCalibratorAL");
    InitVerbosity(topBranch);

    // Get the file name for MHTDC dump output.
    topBranch.GetChild("mhtdcOutputFileName").GetData(fMHTDCOutputFileName);

    // Flag to show that we are not yet initialized.
    fIsInitialized = false;
    fOutputData = 0;

    // Print some message.
    ostringstream info;
    info << " Version: " << GetVersionInfo(VModule::eRevisionNumber) << "\n"
         << "\t * mhtdcOutputFileName = " << fMHTDCOutputFileName << "\n";
    INFO(info);

    return eSuccess;
  }


  /// Process function for MHTDCCalibratorAL
  VModule::EResultFlag
  MHTDCCalibratorAL::Process(Event& event, const AttributeMap& /*attr*/)
  {
    // print some info
    const evt::EventHeader& eventHeader = event.GetEventHeader();
    fRunNumber = eventHeader.GetRunNumber();
    fSpillId = (eventHeader.GetSpillId().IsValid() ? (unsigned int)(eventHeader.GetSpillId()) : 0);
    fEventNumber = eventHeader.GetId();
    fEventUnixTime = UTCDateTime(eventHeader.GetTime()).GetUnixTime().GetSecond();
    fSectionId = -1;
    ostringstream msg;
    msg << "filling MHTDCCalibratorAL tree in run " << fRunNumber
        << " event " << fEventNumber << " ... ";
    INFO(msg);

    // initialize output container in case of first event.
    if ( !fIsInitialized ) {
      // Link fTDCData to a branch of fOutputData tree.
      fOutputData = new TTree("mhtdcData",
                              "MTacc - S11 time differences in TDC units");
      fOutputData->Branch("runNumber", &fRunNumber);
      fOutputData->Branch("spillId", &fSpillId);
      fOutputData->Branch("eventNumber", &fEventNumber);
      fOutputData->Branch("eventUnixTime", &fEventUnixTime);
      fOutputData->Branch("sectionId", &fSectionId);
      fOutputData->Branch("diffMTaccToS11InTDC", &fDiffMTaccToS11InTDC);
      fOutputData->Branch("mhtdcTimeBinWidthInNanoSec", &fMHTDCTimeBinWidthInNanoSec);
      fOutputData->Branch("mhtdcDynamicRange", &fMHTDCDynamicRange);
      fMHTDCTimeBinWidth = det::TimeStructureConst::GetMHTDCTimeBinSpan();
      fMHTDCTimeBinWidthInNanoSec = fMHTDCTimeBinWidth/nanosecond;
      fMHTDCDynamicRange = det::TimeStructureConst::GetMHTDCDynamicRange();
      fIsInitialized = true;
    }

    // get the raw trigger info of the event
    const raw::Trigger& rawTrigger = event.GetRawEvent().GetBeam().GetTrigger();

    // calculate the MHTDC S11 TDC with respect reference (MTacc).
    if ( !rawTrigger.HasTimeStructure(det::TimeStructureConst::eMHTDC, det::TriggerConst::eMTacc) ) {
      ostringstream err;
      err << "MTacc hits not present in MHTDC.\n";
      ERROR(err);
      return eFailure;
    }
    const vector<double>& timeStructureMTacc = rawTrigger.GetTimeStructure(det::TimeStructureConst::eMHTDC, det::TriggerConst::eMTacc);
    if ( timeStructureMTacc.size() != 1 ) {
      ostringstream err;
      err << "Number of MTacc hits (" << timeStructureMTacc.size() << ") in MHTDC is not 1.\n";
      ERROR(err);
      return eFailure;
    }
    const double timeMTacc = timeStructureMTacc[0];
    if ( !rawTrigger.HasTimeStructure(det::TimeStructureConst::eMHTDC, det::TriggerConst::eS1_1) ) {
      ostringstream err;
      err << "S1_1 hits not present in MHTDC.\n";
      ERROR(err);
      return eFailure;
    }
    const vector<double>& timeStructureS11 = rawTrigger.GetTimeStructure(det::TimeStructureConst::eMHTDC, det::TriggerConst::eS1_1);
    if ( timeStructureS11.empty() ) {
      ostringstream err;
      err << "Number of S1_1 hits (" << timeStructureS11.size() << ") in MHTDC is not >= 1.\n";
      ERROR(err);
      return eFailure;
    }
    for ( vector<double>::const_iterator it = timeStructureS11.begin() ; it!= timeStructureS11.end() ; ++it ) {
      const double timeS11 = *it;
      const double timeDiff = timeMTacc - timeS11;

      // save the data to ROOT TTree.
      fDiffMTaccToS11InTDC = lround(timeDiff/fMHTDCTimeBinWidth);
      fOutputData->Fill();
    }

    return eSuccess;
  }


  /// Finish function for MHTDCCalibratorAL
  VModule::EResultFlag
  MHTDCCalibratorAL::Finish()
  {
    // Write output data to file.
    if ( !fMHTDCOutputFileName.empty() ) {
      TFile outputFile(fMHTDCOutputFileName.c_str(), "RECREATE");
      if ( fOutputData != 0 )
        fOutputData->Write();
      outputFile.Close();
      if ( fOutputData !=0 ) {
        delete fOutputData;
        fOutputData = 0;
      }
    }

    INFO("Finished.\n");

    return eSuccess;
  }

}
