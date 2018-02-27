/**
  \file
  Declaration of TPCBottomPointExtractorAB
  Calculates TPC bottom points in terms of drift time 
  without global and chamber t0 correction (raw drift time). 
  This value is determined for last point of each track presumably exiting 
  at the bottom plane and is dumped to a ROOT TTree.

  \author A. Bogovic, A. Laszlo
  \version $Id: TPCBottomPointExtractorAB.h 10683 2014-04-22 21:17:27Z laszlo $
  \date 17 Fev 2012
*/

#ifndef _TPCBottomPointExtractorAB_TPCBottomPointExtractorAB_h_
#define _TPCBottomPointExtractorAB_TPCBottomPointExtractorAB_h_

#include <fwk/VModule.h>
#include <evt/Index.h>
#include <evt/rec/Cluster.h>
#include <modutils/StraightTrackFitter.h>
#include <modutils/MagneticTrackFitter.h>
#include <string>
#include <vector>
#include <TTree.h>

namespace TPCBottomPointExtractorAB {

  /**
    \class   TPCBottomPointExtractorAB
    \author  A. Bogovic, A.Laszlo
    \brief   Calculates TPC bottom points in terms of drift time 
             without global and chamber t0 correction (raw drift time).

    \ingroup CalibrationModules
  */

  class TPCBottomPointExtractorAB : public fwk::VModule {

  public:
    fwk::VModule::EResultFlag Init();
    fwk::VModule::EResultFlag Process(evt::Event& event, const utl::AttributeMap& attr);
    fwk::VModule::EResultFlag Finish();

  private:

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

    // Output data.
    bool fIsInitialized;
    std::vector<TPCData> fTPCData;
    std::vector<TTree*> fOutputData;
    std::string fTPCBottomPointsOutputFileName;

    // Cut parameters.
    double fLocalXRatioMargin;
    unsigned int fPadrowNumberMargin;
    unsigned int fPadrowNumberMarginGTPC;
    unsigned int fTimeSliceMarginIn512Mode;
    double fLocalXLeftMarginGTPCBeamSpot;
    double fLocalXRightMarginGTPCBeamSpot;
    std::vector<unsigned int> fNMinClustersInChamber;

    // For local refitting.
    std::vector< std::vector< std::pair<evt::Index<evt::rec::Cluster>, double> > > fOrderedClusterIndices;
    std::vector< evt::Index<evt::rec::Cluster> > fClusterIndicesForFit;
    modutils::StraightTrackFitter fStraightTrackFitter;
    modutils::MagneticTrackFitter fMagneticTrackFitter;

    REGISTER_MODULE("TPCBottomPointExtractorAB", TPCBottomPointExtractorAB,
                    "$Id: TPCBottomPointExtractorAB.h 10683 2014-04-22 21:17:27Z laszlo $");
  };

}


#endif
