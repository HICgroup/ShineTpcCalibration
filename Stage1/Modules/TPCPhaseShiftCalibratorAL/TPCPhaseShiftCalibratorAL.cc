/**
   \file
   implementation of class TPCPhaseShiftCalibratorAL

   \author A. Laszlo
   \version $Id: TPCPhaseShiftCalibratorAL.cc 10692 2014-04-23 14:19:07Z laszlo $
   \date 17 Feb 2012
*/


#include "TPCPhaseShiftCalibratorAL.h"

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


namespace TPCPhaseShiftCalibratorAL {


  /// Init function for TPCPhaseShiftCalibratorAL
  VModule::EResultFlag
  TPCPhaseShiftCalibratorAL::Init()
  {
    CentralConfig& cc = CentralConfig::GetInstance();

    Branch topBranch = cc.GetTopBranch("TPCPhaseShiftCalibratorAL");
    InitVerbosity(topBranch);

    // Get the TPC phase shift calibration file name for output.
    topBranch.GetChild("tpcPhaseOutputFileName").GetData(fTPCPhaseOutputFileName);

    // Flag to show that we are not yet initialized.
    fIsInitialized = false;
    fOutputData = 0;

    // Print some message.
    ostringstream info;
    info << " Version: " << GetVersionInfo(VModule::eRevisionNumber) << "\n"
         << "\t * tpcPhaseOutputFileName = " << fTPCPhaseOutputFileName << "\n";
    INFO(info);

    return eSuccess;
  }


  /// Process function for TPCPhaseShiftCalibratorAL
  VModule::EResultFlag
  TPCPhaseShiftCalibratorAL::Process(Event& event, const AttributeMap& /*attr*/)
  {
    // print some info
    const evt::EventHeader& eventHeader = event.GetEventHeader();
    fRunNumber = eventHeader.GetRunNumber();
    fSpillId = (eventHeader.GetSpillId().IsValid() ? (unsigned int)(eventHeader.GetSpillId()) : 0);
    fEventNumber = eventHeader.GetId();
    fEventUnixTime = UTCDateTime(eventHeader.GetTime()).GetUnixTime().GetSecond();
    fSectionId = -1;
    ostringstream msg;
    msg << "filling TPCPhaseShiftCalibratorAL tree in run " << fRunNumber
        << " event " << fEventNumber << " ... ";
    INFO(msg);

    // initialize output container in case of first event.
    if ( !fIsInitialized ) {
      // Link fTDCData to a branch of fOutputData tree.
      fOutputData = new TTree("tpcPhaseShift",
                              "TPC phase shift measurements with respect to S11 in TDC units");
      fOutputData->Branch("runNumber", &fRunNumber);
      fOutputData->Branch("spillId", &fSpillId);
      fOutputData->Branch("eventNumber", &fEventNumber);
      fOutputData->Branch("eventUnixTime", &fEventUnixTime);
      fOutputData->Branch("sectionId", &fSectionId);
      fOutputData->Branch("tpcPhaseTDC", &fTPCPhaseTDC);
      fOutputData->Branch("tpcPrimaryClockTimeWidthInNanoSec", &fTPCPrimaryClockTimeWidthInNanoSec);
      fTPCPrimaryClockTimeWidthInNanoSec = det::TPCConst::GetTPCPrimaryClockTimeBinSpan()/nanosecond;
      fIsInitialized = true;
    }

    // get the trigger info of the event
    const raw::Trigger& trigger = event.GetRawEvent().GetBeam().GetTrigger();

    // save the TPC phase shift TDC counts
    fTPCPhaseTDC = trigger.GetTDC(det::TriggerConst::eTPCPhase);
    fOutputData->Fill();

    return eSuccess;
  }


  /// Finish function for TPCPhaseShiftCalibratorAL
  VModule::EResultFlag
  TPCPhaseShiftCalibratorAL::Finish()
  {
    // Write output data to file.
    if ( !fTPCPhaseOutputFileName.empty() ) {
      TFile outputFile(fTPCPhaseOutputFileName.c_str(), "RECREATE");
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
