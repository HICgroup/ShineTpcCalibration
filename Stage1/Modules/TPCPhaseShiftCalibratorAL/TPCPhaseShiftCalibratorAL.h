/**
  \file
  Declaration of TPCPhaseShiftCalibratorAL

  \author A. Laszlo
  \version $Id: TPCPhaseShiftCalibratorAL.h 10683 2014-04-22 21:17:27Z laszlo $
  \date 17 Fev 2012
*/

#ifndef _TPCPhaseShiftCalibratorAL_TPCPhaseShiftCalibratorAL_h_
#define _TPCPhaseShiftCalibratorAL_TPCPhaseShiftCalibratorAL_h_

#include <fwk/VModule.h>
#include <string>
#include <TTree.h>

namespace TPCPhaseShiftCalibratorAL {

  /**
    \class TPCPhaseShiftCalibratorAL
    \author A. Laszlo
    \brief Calibrates TPC phase shift TDC scale (phase of S11 to TPC clock)

    \ingroup CalibrationModules
  */

  class TPCPhaseShiftCalibratorAL : public fwk::VModule {

  public:
    fwk::VModule::EResultFlag Init();
    fwk::VModule::EResultFlag Process(evt::Event& event, const utl::AttributeMap& attr);
    fwk::VModule::EResultFlag Finish();

  private:

    // Containers for TDC info for TPC phase in output file.
    // Obligatory header to identify event.
    UInt_t fRunNumber;
    UInt_t fSpillId;
    UInt_t fEventNumber;
    UInt_t fEventUnixTime;
    Int_t fSectionId;
    // Measured TPC phase in TDC units.
    Int_t fTPCPhaseTDC;
    // Primary clock tick width in nanosec.
    Float_t fTPCPrimaryClockTimeWidthInNanoSec;

    // Administration of output file.
    bool fIsInitialized;
    TTree* fOutputData;
    std::string fTPCPhaseOutputFileName;

    REGISTER_MODULE("TPCPhaseShiftCalibratorAL", TPCPhaseShiftCalibratorAL,
                    "$Id: TPCPhaseShiftCalibratorAL.h 10683 2014-04-22 21:17:27Z laszlo $");
  };

}


#endif
