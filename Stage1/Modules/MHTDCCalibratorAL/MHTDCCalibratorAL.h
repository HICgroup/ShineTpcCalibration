/**
  \file
  Declaration of MHTDCCalibratorAL

  \author A. Laszlo
  \version $Id: MHTDCCalibratorAL.h 10683 2014-04-22 21:17:27Z laszlo $
  \date 17 Fev 2012
*/

#ifndef _MHTDCCalibratorAL_MHTDCCalibratorAL_h_
#define _MHTDCCalibratorAL_MHTDCCalibratorAL_h_

#include <fwk/VModule.h>
#include <string>
#include <TTree.h>

namespace MHTDCCalibratorAL {

  /**
    \class MHTDCCalibratorAL
    \author A. Laszlo
    \brief Calculates the time shift and time window of S1_1 MHTDC hit with respect to MTacc.

    \ingroup CalibrationModules
  */

  class MHTDCCalibratorAL : public fwk::VModule {

  public:
    fwk::VModule::EResultFlag Init();
    fwk::VModule::EResultFlag Process(evt::Event& event, const utl::AttributeMap& attr);
    fwk::VModule::EResultFlag Finish();

  private:

    // Containers for storing info for output file.
    // Obligatory header to identify event.
    UInt_t fRunNumber;
    UInt_t fSpillId;
    UInt_t fEventNumber;
    UInt_t fEventUnixTime;
    Int_t fSectionId;
    // Measured MTacc - S11 difference in TDC units.
    Int_t fDiffMTaccToS11InTDC;
    // TDC unit in nanosec.
    Float_t fMHTDCTimeBinWidthInNanoSec;
    // MHTDC dynamic range.
    Int_t fMHTDCDynamicRange;

    // Output file administration.
    bool fIsInitialized;
    TTree* fOutputData;
    double fMHTDCTimeBinWidth;
    std::string fMHTDCOutputFileName;

    REGISTER_MODULE("MHTDCCalibratorAL", MHTDCCalibratorAL,
                    "$Id: MHTDCCalibratorAL.h 10683 2014-04-22 21:17:27Z laszlo $");
  };

}


#endif
