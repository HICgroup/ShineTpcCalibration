/**
  \file
  Declaration of DriftVelocityFromTOFPK

  \author P.Kovesarki
  \version $Id: DriftVelocityFromTOFPK.h 11599 2015-08-31 08:17:16Z laszlo $
  \date 17 Fev 2012
*/

#ifndef _DriftVelocityFromTOFPK_DriftVelocityFromTOFPK_h
#define _DriftVelocityFromTOFPK_DriftVelocityFromTOFPK_h

#include <fwk/VModule.h>
#include <modutils/TrackAutopsy.h>
#include <vector>


namespace DriftVelocityFromTOFPK {

  /**
    \class DriftVelocityFromTOFPK
    \author P. Kovesarki
    \brief Calibrates TPC drift velocity using TPC tracks matched with TOF hits
    \ingroup CalibrationModules
  */


  class DriftVelocityFromTOFPK : public fwk::VModule {

  public:
    fwk::VModule::EResultFlag Init();
    fwk::VModule::EResultFlag Process(evt::Event& event, const utl::AttributeMap& attr);
    fwk::VModule::EResultFlag Finish();

  private:
    /// Enum for Cluster quality.
    enum EClusterQuality {
      eUsedInFit = 0x8
    };

    /// Cut parameters.
    std::vector<unsigned int> fNMinClusters;
    double fLocalXLeftMarginGTPCBeamSpot;
    double fLocalXRightMarginGTPCBeamSpot;
    double fTimeOfFlightMin;
    double fTimeOfFlightMax;
    double fNormalizedChargeMin;
    double fNormalizedChargeMax;

    bool fIsInitialized;
    modutils::TrackAutopsy::TrackMatchData fData; /// Things go here that will be or may be stored into root files
    std::string fTrackMatchOutputFileName;

    UInt_t fRunNumber;
    UInt_t fSpillId;
    UInt_t fEventNumber;
    UInt_t fEventUnixTime;

    REGISTER_MODULE("DriftVelocityFromTOFPK", DriftVelocityFromTOFPK,
                    "$Id: DriftVelocityFromTOFPK.h 11599 2015-08-31 08:17:16Z laszlo $");
  };

}


#endif
