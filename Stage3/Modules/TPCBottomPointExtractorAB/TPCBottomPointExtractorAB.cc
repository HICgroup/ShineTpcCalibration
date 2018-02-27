/**
   \file
   implementation of class TPCBottomPointExtractorAB

   \author A. Bogovic, A. Laszlo
   \version $Id: TPCBottomPointExtractorAB.cc 11534 2015-07-11 22:41:25Z laszlo $
   \date 17 Feb 2012
*/


#include "TPCBottomPointExtractorAB.h"

#include <fwk/CentralConfig.h>
#include <utl/ErrorLogger.h>
#include <utl/ShineException.h>
#include <utl/UTCDateTime.h>

#include <det/Detector.h>
#include <det/TPCConst.h>
#include <det/TPC.h>
#include <det/TPCChamber.h>
#include <det/TPCSector.h>

#include <evt/Event.h>
#include <evt/RecEvent.h>
#include <evt/rec/Track.h>
#include <evt/rec/Cluster.h>

#include <sstream>
#include <fstream>
#include <cmath>

#include <TFile.h>


using namespace std;
using namespace fwk;
using namespace evt;
using namespace utl;


namespace TPCBottomPointExtractorAB {


  /// Init function for TPCBottomPointExtractorAB
  VModule::EResultFlag
  TPCBottomPointExtractorAB::Init()
  {
    CentralConfig& cc = CentralConfig::GetInstance();

    Branch topBranch = cc.GetTopBranch("TPCBottomPointExtractorAB");
    InitVerbosity(topBranch);

    // Some cut parameters to keep last points away from the sides of chambers.
    topBranch.GetChild("localXRatioMargin").GetData(fLocalXRatioMargin);
    topBranch.GetChild("padRowNumberMargin").GetData(fPadrowNumberMargin);
    topBranch.GetChild("padRowNumberMarginGTPC").GetData(fPadrowNumberMarginGTPC);
    topBranch.GetChild("timeSliceMarginIn512Mode").GetData(fTimeSliceMarginIn512Mode);
    topBranch.GetChild("localXLeftMarginGTPCBeamSpot").GetData(fLocalXLeftMarginGTPCBeamSpot);
    topBranch.GetChild("localXRightMarginGTPCBeamSpot").GetData(fLocalXRightMarginGTPCBeamSpot);

    // Cut on required number of clusters.
    fNMinClustersInChamber.resize((int)det::TPCConst::eNTPCs, 0);
    Branch nMinClustersB = topBranch.GetChild("nMinClusters");
    for ( Branch cB = nMinClustersB.GetFirstChild() ; cB ; ++cB ) {
      const string tpcName = cB.GetName();
      const unsigned int nMinClusters = cB.Get<unsigned int>();
      const det::TPCConst::EId tpcId = det::TPCConst::GetId(tpcName);
      if ( tpcId==det::TPCConst::eUnknown ) {
        ostringstream err;
        err << "Unknown TPC chamber '" << tpcName << "' !";
        ERROR(err);
        throw NonExistentComponentException(err.str());
      }
      const int tpcIndex = int(tpcId)-int(det::TPCConst::eFirstTPC);
      fNMinClustersInChamber.at(tpcIndex) = nMinClusters;
    }

    // Get the TPC bottom point file name for output.
    topBranch.GetChild("tpcBottomPointsOutputFileName").GetData(fTPCBottomPointsOutputFileName);

    // Initialize output container.
    fIsInitialized = false;
    fTPCData.resize((int)det::TPCConst::eNTPCs);
    fOutputData.resize((int)det::TPCConst::eNTPCs, 0);

    // Allocate cluster index containers for ordering.
    fOrderedClusterIndices.resize((int)det::TPCConst::eNTPCs);

    // Print some message.
    ostringstream info;
    info << " Version: " << GetVersionInfo(VModule::eRevisionNumber) << "\n";
    info << "\t * localXRatioMargin = "
         << fLocalXRatioMargin << "\n"
         << "\t * padRowNumberMargin = "
         << fPadrowNumberMargin << "\n"
         << "\t * padRowNumberMarginGTPC = "
         << fPadrowNumberMarginGTPC << "\n"
         << "\t * timeSliceMarginIn512Mode = "
         << fTimeSliceMarginIn512Mode << "\n"
         << "\t * localXLeftMarginGTPCBeamSpot = "
         << (fLocalXLeftMarginGTPCBeamSpot/cm) << " cm" << "\n"
         << "\t * localXRightMarginGTPCBeamSpot = "
         << (fLocalXRightMarginGTPCBeamSpot/cm) << " cm" << "\n";
    info << "\t * nMinClusters = ";
    for ( unsigned int i=0 ; i<fNMinClustersInChamber.size() ; ++i ) {
      info << (i==0 ? string("") : string(",")) << fNMinClustersInChamber[i];
    }
    info << "\n";
    info << "\t * tpcBottomPointsOutputFileName = "
         << fTPCBottomPointsOutputFileName << "\n";
    INFO(info);

    return eSuccess;
  }


  /// Comparison function of cluster index iterators for increasing local Z.
  static
  inline
  bool
  CompareIncreasingLocalZ(const pair<Index<rec::Cluster>, double>& clusterIndexAndZPos1,
                          const pair<Index<rec::Cluster>, double>& clusterIndexAndZPos2)
  {
    return (clusterIndexAndZPos1.second < clusterIndexAndZPos2.second);
  }


  /// Process function for TPCBottomPointExtractorAB
  VModule::EResultFlag
  TPCBottomPointExtractorAB::Process(Event& event, const AttributeMap& /*attr*/)
  {
    // Print some info.
    const evt::EventHeader& eventHeader = event.GetEventHeader();
    const unsigned int runNumber = eventHeader.GetRunNumber();
    const unsigned int spillId = (eventHeader.GetSpillId().IsValid() ? (unsigned int)(eventHeader.GetSpillId()) : 0);
    const unsigned int eventNumber = eventHeader.GetId();
    const unsigned int eventUnixTime = UTCDateTime(eventHeader.GetTime()).GetUnixTime().GetSecond();
    const int sectionId = -1;
    ostringstream msg;
    msg << "filling TPCBottomPointExtractorAB tree in run " << runNumber
        << " event " << eventNumber << " ... ";
    INFO(msg);

    // Get handle of Detector description.
    const det::Detector& detector = det::Detector::GetInstance();
    const det::TPC& detTPC = detector.GetTPC();

    // Initialize output container in case of first event.
    if (!fIsInitialized) {
      // Link fTPCData[tpcIndex] to a branch of fOutputData tree.
      for (unsigned int tpcIndex = 0; tpcIndex<fTPCData.size();
           ++tpcIndex) {
        const det::TPCConst::EId tpcId =
          det::TPCConst::EId(tpcIndex+int(det::TPCConst::eFirstTPC));
        const string& tpcName = det::TPCConst::GetName(tpcId);
        fOutputData[tpcIndex] = new TTree(tpcName.c_str(),
                                          (string("Raw drift times of tracks which supposed to exit at the drift cathode plane in ")+tpcName).c_str());
        fOutputData[tpcIndex]->Branch("runNumber", &(fTPCData[tpcIndex].fRunNumber));
        fOutputData[tpcIndex]->Branch("spillId", &(fTPCData[tpcIndex].fSpillId));
        fOutputData[tpcIndex]->Branch("eventNumber", &(fTPCData[tpcIndex].fEventNumber));
        fOutputData[tpcIndex]->Branch("eventUnixTime", &(fTPCData[tpcIndex].fEventUnixTime));
        fOutputData[tpcIndex]->Branch("sectionId", &(fTPCData[tpcIndex].fSectionId));
        fOutputData[tpcIndex]->Branch("lastClusterSectorNumber", &(fTPCData[tpcIndex].fLastClusterSectorNumber));
        fOutputData[tpcIndex]->Branch("lastClusterPadrowNumber", &(fTPCData[tpcIndex].fLastClusterPadrowNumber));
        fOutputData[tpcIndex]->Branch("lastClusterPadrowXInCM", &(fTPCData[tpcIndex].fLastClusterPadrowXInCM));
        fOutputData[tpcIndex]->Branch("lastClusterRawDriftTimeInUSec", &(fTPCData[tpcIndex].fLastClusterRawDriftTimeInUSec));
        fOutputData[tpcIndex]->Branch("lastPointRawDriftTimeInUSec", &(fTPCData[tpcIndex].fLastPointRawDriftTimeInUSec));
        fOutputData[tpcIndex]->Branch("fullDriftTimeInUSec", &(fTPCData[tpcIndex].fFullDriftTimeInUSec));
      }
      fIsInitialized = true;
    }

    // Get handle of RecEvent.
    const RecEvent& recEvent = event.GetRecEvent();

    using namespace rec;

    const double tbSpan512Mode = det::TPCConst::GetTimeBinSpan(det::TPCConst::e512TBMode);

    modutils::StraightTrackFitResult straightTrackFitResult;
    TrackLocalParameters trackLocalParameters;
    TrackLocalParameters trackLastPointLocalParameters;
    const det::MagneticFieldTracker& magneticFieldTracker = detector.GetMagneticFieldTracker();

    // Loop over tracks, find ones which should exit at the bottom plane and
    // collect the bottom points.
    for (list<Track>::const_iterator trackIt = recEvent.Begin<Track>();
         trackIt != recEvent.End<Track>(); ++ trackIt) {
      const Track& track = *trackIt;
      // Reset number of clusters in each TPC to zero.
      for (int tpcIndex = 0; tpcIndex < int(fOrderedClusterIndices.size());
           ++tpcIndex)
        fOrderedClusterIndices[tpcIndex].clear();

      // Loop over track clusters.
      for (ClusterIndexIterator clusterIndexIt = track.ClustersBegin();
           clusterIndexIt != track.ClustersEnd() ; ++clusterIndexIt) {
        const Index<Cluster>& clusterIndex = *clusterIndexIt;
        const Cluster& cluster = recEvent.Get(clusterIndex);
        const det::TPCConst::EId tpcId = cluster.GetTPCId();
        const int tpcIndex = int(tpcId)-int(det::TPCConst::eFirstTPC);
        const det::TPCChamber& detTPCChamber = detTPC.GetChamber(tpcId);
        const CoordinateSystemPtr& chamberCoordinateSystem =
          detTPCChamber.GetComponentCoordinateSystem();
        const utl::Point& position = cluster.GetPosition();
        const utl::Triple coordinates =
          position.GetCoordinates(chamberCoordinateSystem);
        const double localZ = coordinates.get<2>();
        pair<Index<Cluster>, double> index(clusterIndex, localZ);
        fOrderedClusterIndices[tpcIndex].push_back(index);
      }

      // Order track clusters according to local Z coordinate in each TPC.
      for (int tpcIndex = 0; tpcIndex < int(fOrderedClusterIndices.size());
           ++tpcIndex ) {
        if (fOrderedClusterIndices[tpcIndex].empty())
          continue;
        sort(fOrderedClusterIndices[tpcIndex].begin(),
             fOrderedClusterIndices[tpcIndex].end(), CompareIncreasingLocalZ);
      }

      // Extract last points.
      for ( int tpcIndex = 0 ; tpcIndex < int(fOrderedClusterIndices.size()) ;
            ++tpcIndex ) {
        const vector< pair<evt::Index<Cluster>, double> >& orderedClusterIndexList =
          fOrderedClusterIndices[tpcIndex];
        if ( orderedClusterIndexList.size()<max((unsigned int)4, fNMinClustersInChamber[tpcIndex]) )
          continue;
        const Index<Cluster>& lastClusterIndex =
          fOrderedClusterIndices[tpcIndex].back().first;
        const Cluster& lastCluster = recEvent.Get<Cluster>(lastClusterIndex);
        const det::TPCConst::EId lastClusterTPCId = lastCluster.GetTPCId();
        const unsigned int lastClusterSectorNumber = lastCluster.GetSectorNumber();
        const unsigned int lastClusterPadrowNumber = lastCluster.GetPadrowNumber();
        const Point& lastClusterPosition = lastCluster.GetPosition();
        // Store the local Y coordinate of last point in according tree in terms of drift time.
        const det::TPCChamber& detLastClusterTPCChamber =
          detTPC.GetChamber(lastClusterTPCId);
        const CoordinateSystemPtr& lastClusterTPCCS =
          detLastClusterTPCChamber.GetComponentCoordinateSystem();
        const det::TPCSector& detLastClusterTPCSector = detLastClusterTPCChamber.GetSector(lastClusterSectorNumber);
        const unsigned int lastClusterNPadrows = detLastClusterTPCSector.GetNPadrows();
        const unsigned int padrowNumberMargin = (lastClusterTPCId == det::TPCConst::eGTPC ? fPadrowNumberMarginGTPC : fPadrowNumberMargin);
        // Cut away tracks ending too close to front and back plane.
        if ( lastClusterPadrowNumber<=padrowNumberMargin || lastClusterPadrowNumber+padrowNumberMargin>lastClusterNPadrows )
          continue;
        // Cut away tracks ending too close to amplification (sense-wire) plane.
        const CoordinateSystemPtr& lastClusterTPCSectorCS =
          detLastClusterTPCSector.GetComponentCoordinateSystem();
        const Triple lastClusterSectorCoordinates =
          lastClusterPosition.GetCoordinates(lastClusterTPCSectorCS);
        const double lastClusterSectorY = lastClusterSectorCoordinates.get<1>();
        const double driftVelocity = detLastClusterTPCChamber.GetDriftVelocity();
        const double senseWirePlaneSectorY = detLastClusterTPCSector.GetLocalSenseWirePlaneY();
        const double driftLength = detLastClusterTPCSector.GetDriftVolumeHeight();
        const double fullDriftTime = driftLength/driftVelocity;
        const double globalT0 = detTPC.GetGlobalT0();
        const double lastClusterTPCChamberT0 = detLastClusterTPCChamber.GetChamberT0();
        const double lastClusterRawDriftTime = (senseWirePlaneSectorY-lastClusterSectorY)/driftVelocity - (globalT0+lastClusterTPCChamberT0);
        const int lastClusterRawDriftTimeBinIn512Mode = ::lround((lastClusterRawDriftTime)/tbSpan512Mode);
        if ( lastClusterRawDriftTimeBinIn512Mode<=int(fTimeSliceMarginIn512Mode) )
          continue;
        // Cut away tracks ending too close to side planes.
        const det::TPCPadrow& detLastClusterTPCPadrow = detLastClusterTPCSector.GetPadrow(lastClusterPadrowNumber);
        const CoordinateSystemPtr& lastClusterTPCPadrowCS =
          detLastClusterTPCPadrow.GetComponentCoordinateSystem();
        const int lastClusterPadrowNPads = detLastClusterTPCPadrow.GetNPads();
        const double lastClusterPadrowXFirst = detLastClusterTPCPadrow.GetPadPosition(1).GetCoordinates(lastClusterTPCPadrowCS).get<0>();
        const double lastClusterPadrowXLast = detLastClusterTPCPadrow.GetPadPosition(lastClusterPadrowNPads).GetCoordinates(lastClusterTPCPadrowCS).get<0>();
        const double lastClusterPadrowXMin = min(lastClusterPadrowXFirst, lastClusterPadrowXLast);
        const double lastClusterPadrowXMax = max(lastClusterPadrowXFirst, lastClusterPadrowXLast);
        const double lastClusterPadrowXWidth = ::fabs(lastClusterPadrowXMax-lastClusterPadrowXMin);
        const double lastClusterPadrowXMargin = lastClusterPadrowXWidth*fLocalXRatioMargin;
        const double lastClusterPadrowX = lastClusterPosition.GetCoordinates(lastClusterTPCPadrowCS).get<0>();
        if ( lastClusterPadrowX<lastClusterPadrowXMin+lastClusterPadrowXMargin || lastClusterPadrowXMax-lastClusterPadrowXMargin<lastClusterPadrowX )
          continue;
        // Cut away beam spot y-strip in GTPC to reject on-time and off-time beams.
        if ( lastClusterTPCId == det::TPCConst::eGTPC && fLocalXLeftMarginGTPCBeamSpot<=lastClusterPadrowX && lastClusterPadrowX<=fLocalXRightMarginGTPCBeamSpot )
          continue;
        // Refit local track clusters excluding last cluster to get more accurate last point estimate.
        const Triple lastClusterChamberCoordinates =
          lastClusterPosition.GetCoordinates(lastClusterTPCCS);
        double lastPointChamberX = lastClusterChamberCoordinates.get<0>();
        double lastPointChamberY = lastClusterChamberCoordinates.get<1>();
        const double lastPointChamberZ = lastClusterChamberCoordinates.get<2>();
        fClusterIndicesForFit.clear();
        for ( vector< pair<Index<Cluster>, double> >::const_iterator it=orderedClusterIndexList.begin() ; it!=orderedClusterIndexList.end() ; ++it ) {
          vector< pair<Index<Cluster>, double> >::const_iterator nextIt = it;
          ++nextIt;
          if ( nextIt!=orderedClusterIndexList.end() )
            fClusterIndicesForFit.push_back((*it).first);
        }
        // Fit VTPC1, VTPC2 and GTPC piece with magnetic field tracker.
        if ( lastClusterTPCId==det::TPCConst::eVTPC1 || lastClusterTPCId==det::TPCConst::eVTPC2 || lastClusterTPCId==det::TPCConst::eGTPC ) {
          fMagneticTrackFitter.Init(recEvent, lastClusterTPCCS);
          if ( !fMagneticTrackFitter.KalmanFilter(fClusterIndicesForFit, trackLocalParameters) )
            continue;
          if ( !magneticFieldTracker.TrackToLocalZ(trackLocalParameters, lastPointChamberZ, lastClusterTPCCS, trackLastPointLocalParameters) )
            continue;
          const Triple lastPointChamberCoordinates =
            trackLastPointLocalParameters.GetPosition().GetCoordinates(lastClusterTPCCS);
          lastPointChamberX = lastPointChamberCoordinates.get<0>();
          lastPointChamberY = lastPointChamberCoordinates.get<1>();
        }
        // Fit others with straight track fitter.
        else {
          fStraightTrackFitter.Init(recEvent, lastClusterTPCCS);
          if ( !fStraightTrackFitter.FitStraightLine(fClusterIndicesForFit, straightTrackFitResult) )
            continue;
          const modutils::StraightTrackParameters& straightTrackParameters =
            straightTrackFitResult.fTrackParameters;
          lastPointChamberX = straightTrackParameters.fMX + straightTrackParameters.fNX*lastPointChamberZ;
          lastPointChamberY = straightTrackParameters.fMY + straightTrackParameters.fNY*lastPointChamberZ;
        }
        const Point lastPoint(lastPointChamberX, lastPointChamberY, lastPointChamberZ, lastClusterTPCCS);
        const Triple lastPointSectorCoordinates =
          lastPoint.GetCoordinates(lastClusterTPCSectorCS);
        const double lastPointSectorY = lastPointSectorCoordinates.get<1>();
        const double lastPointRawDriftTime = (senseWirePlaneSectorY-lastPointSectorY)/driftVelocity - (globalT0+lastClusterTPCChamberT0);
        // Store the data for output.
        fTPCData[tpcIndex].fRunNumber = runNumber;
        fTPCData[tpcIndex].fSpillId = spillId;
        fTPCData[tpcIndex].fEventNumber = eventNumber;
        fTPCData[tpcIndex].fEventUnixTime = eventUnixTime;
        fTPCData[tpcIndex].fSectionId = sectionId;
        fTPCData[tpcIndex].fLastClusterSectorNumber = lastClusterSectorNumber;
        fTPCData[tpcIndex].fLastClusterPadrowNumber = lastClusterPadrowNumber;
        fTPCData[tpcIndex].fLastClusterPadrowXInCM = lastClusterPadrowX/cm;
        fTPCData[tpcIndex].fLastClusterRawDriftTimeInUSec = lastClusterRawDriftTime/microsecond;
        fTPCData[tpcIndex].fLastPointRawDriftTimeInUSec = lastPointRawDriftTime/microsecond;
        fTPCData[tpcIndex].fFullDriftTimeInUSec = fullDriftTime/microsecond;
        fOutputData[tpcIndex]->Fill();
      }
    }

    return eSuccess;
  }


  /// Finish function for TPCBottomPointExtractorAB
  VModule::EResultFlag
  TPCBottomPointExtractorAB::Finish()
  {

    // Write output data to file.
    if ( !fTPCBottomPointsOutputFileName.empty() ) {
      TFile outputFile(fTPCBottomPointsOutputFileName.c_str(), "RECREATE");
      for (unsigned int tpcIndex = 0; tpcIndex<fOutputData.size(); ++tpcIndex) {
        if ( fOutputData[tpcIndex] != 0 )
          fOutputData[tpcIndex]->Write();
      }
      outputFile.Close();
      for (unsigned int tpcIndex = 0; tpcIndex<fOutputData.size(); ++tpcIndex) {
        if ( fOutputData[tpcIndex] != 0 ) {
          delete fOutputData[tpcIndex];
          fOutputData[tpcIndex] = 0;
        }
      }
    }

    INFO("Finished.\n");

    return eSuccess;
  }
}
