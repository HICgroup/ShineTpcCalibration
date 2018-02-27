/**
   \file
   implementation of class DriftVelocityFromTOFPK

   \author P. Kovesarki
   \version $Id: DriftVelocityFromTOFPK.cc 11599 2015-08-31 08:17:16Z laszlo $
   \date 6 Aug 2013
*/


#include "DriftVelocityFromTOFPK.h"
#include <fwk/CentralConfig.h>
#include <utl/ErrorLogger.h>

#include <evt/Event.h>
#include <det/Detector.h>
#include <utl/GeometryUtilities.h>
#include <utl/Math.h>
#include <utl/UTCDateTime.h>

#include <sstream>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include "TFile.h"

using namespace fwk;
using namespace evt;
using namespace utl;
using namespace std;


namespace DriftVelocityFromTOFPK {


  /// Init function for DriftVelocityFromTOFPK
  VModule::EResultFlag
  DriftVelocityFromTOFPK::Init()
  {
    CentralConfig& cc = CentralConfig::GetInstance();

    Branch topBranch = cc.GetTopBranch("DriftVelocityFromTOFPK");
    InitVerbosity(topBranch);

    // Get the TPC track length cuts.
    fNMinClusters.resize(modutils::TrackAutopsy::TrackAutopsy::eNTPC, 0);
    topBranch.GetChild("nMinVTPC1Clusters").GetData(fNMinClusters[modutils::TrackAutopsy::TrackAutopsy::eVTPC1]);
    topBranch.GetChild("nMinVTPC2Clusters").GetData(fNMinClusters[modutils::TrackAutopsy::TrackAutopsy::eVTPC2]);
    topBranch.GetChild("nMinMTPCLClusters").GetData(fNMinClusters[modutils::TrackAutopsy::TrackAutopsy::eMTPCL]);
    topBranch.GetChild("nMinMTPCRClusters").GetData(fNMinClusters[modutils::TrackAutopsy::TrackAutopsy::eMTPCR]);
    topBranch.GetChild("nMinGTPCClusters").GetData(fNMinClusters[modutils::TrackAutopsy::TrackAutopsy::eGTPC]);

    // GTPC beam spot strip cut in X to reject off-time particles.
    topBranch.GetChild("localXLeftMarginGTPCBeamSpot").GetData(fLocalXLeftMarginGTPCBeamSpot);
    topBranch.GetChild("localXRightMarginGTPCBeamSpot").GetData(fLocalXRightMarginGTPCBeamSpot);

    // Get the TOF hit cuts.
    topBranch.GetChild("timeOfFlightMin").GetData(fTimeOfFlightMin);
    topBranch.GetChild("timeOfFlightMax").GetData(fTimeOfFlightMax);
    topBranch.GetChild("normalizedChargeMin").GetData(fNormalizedChargeMin);
    topBranch.GetChild("normalizedChargeMax").GetData(fNormalizedChargeMax);

    // Get the file name for track match output.
    topBranch.GetChild("trackMatchOutputFileName").GetData(fTrackMatchOutputFileName);

    // Flag to show that we are not yet initialized.
    fIsInitialized = false;

    // Print some message.
    ostringstream info;
    info << " Version: " << GetVersionInfo(VModule::eRevisionNumber) << "\n";
    info << "\t * nMinVTPC1Clusters = " << fNMinClusters[modutils::TrackAutopsy::TrackAutopsy::eVTPC1] << "\n";
    info << "\t * nMinVTPC2Clusters = " << fNMinClusters[modutils::TrackAutopsy::TrackAutopsy::eVTPC2] << "\n";
    info << "\t * nMinMTPCLClusters = " << fNMinClusters[modutils::TrackAutopsy::TrackAutopsy::eMTPCL] << "\n";
    info << "\t * nMinMTPCRClusters = " << fNMinClusters[modutils::TrackAutopsy::TrackAutopsy::eMTPCR] << "\n";
    info << "\t * nMinGTPCClusters = " << fNMinClusters[modutils::TrackAutopsy::TrackAutopsy::eGTPC] << "\n";
    info << "\t * localXLeftMarginGTPCBeamSpot = " << (fLocalXLeftMarginGTPCBeamSpot/cm) << " cm" << "\n";
    info << "\t * localXRightMarginGTPCBeamSpot = " << (fLocalXRightMarginGTPCBeamSpot/cm) << " cm" << "\n";
    info << "\t * timeOfFlightMin = " << (fTimeOfFlightMin/nanosecond) << " nanosecond" << "\n";
    info << "\t * timeOfFlightMax = " << (fTimeOfFlightMax/nanosecond) << " nanosecond" << "\n";
    info << "\t * normalizedChargeMin = " << fNormalizedChargeMin << "\n";
    info << "\t * normalizedChargeMax = " << fNormalizedChargeMax << "\n";
    info << "\t * trackMatchOutputFileName = " << fTrackMatchOutputFileName << "\n";

    INFO(info);

    return eSuccess;
  }


  /// Process function for DriftVelocityFromTOFPK
  VModule::EResultFlag
  DriftVelocityFromTOFPK::Process(Event& event, const AttributeMap& /*attr*/)
  {
    // print some info
    const evt::EventHeader& eventHeader = event.GetEventHeader();
    const unsigned int runNumber = eventHeader.GetRunNumber();
    const utl::Validated<unsigned int> spillId = eventHeader.GetSpillId();
    const unsigned int eventNumber = eventHeader.GetId();
    const utl::TimeStamp& timeStamp = eventHeader.GetTime();
    fRunNumber = runNumber;
    fSpillId = (spillId.IsValid() ? (unsigned int)spillId : 0);
    fEventNumber = eventNumber;
    fEventUnixTime = UTCDateTime(timeStamp).GetUnixTime().GetSecond();

    det::Detector& detector = det::Detector::GetInstance();
    ostringstream msg;
    msg << "calculating DriftVelocityFromTOFPK values in run " << runNumber
        << " event " << fEventNumber << " ... " << "\n";
    INFO(msg);

    // initialize TTree dump
    if (!fIsInitialized) {
      fData.Init(fEventUnixTime, fEventNumber, fRunNumber, fSpillId);
      fIsInitialized = true;
    }

    // Get handle on global coordinate system in experiment.
    const CoordinateSystemPtr na61CoordSys = detector.GetDetectorCoordinateSystem();
    // Get handle on raw trigger info.
    const raw::Trigger& rawTrigger = event.GetRawEvent().GetBeam().GetTrigger();
    // Get handle on reconstructed event.
    RecEvent& recEvent = event.GetRecEvent();

    // Do here cuts and histograms with rawTrigger, recTrigger, recBeam...
    if (!rawTrigger.IsTrigger(det::TriggerConst::eT2, det::TriggerConst::eAll))
      return eContinueLoop;
    // Check if we have main vertex in recEvent, otherwise take next event.
    if (!recEvent.HasMainVertex())
      return eContinueLoop;
    // Get handle on main vertex.
    const rec::Vertex& mainVertex = recEvent.GetMainVertex();

    modutils::TrackAutopsy::TPCTrackFitter trackSegmentFitter(recEvent, detector.GetMagneticFieldTracker());
    modutils::TrackAutopsy::TrackMatchHandler trackMatchHandler(recEvent, detector.GetMagneticFieldTracker());
    modutils::TrackAutopsy::TrackRegister trackRegister(recEvent, detector.GetMagneticFieldTracker());

    // Get GTPC Sector1 coordinate system for beam spot vetoing.
    const CoordinateSystemPtr gtpcSector1CS = detector.GetTPC().GetChamber(det::TPCConst::eGTPC).GetSector(1).GetComponentCoordinateSystem();

    /// Loop over main vertex tracks.
    const rec::VertexTrackIndexIterator mainVertexTrackBegin = mainVertex.DaughterTracksBegin();
    const rec::VertexTrackIndexIterator mainVertexTrackEnd = mainVertex.DaughterTracksEnd();
    for (rec::VertexTrackIndexIterator mainVertexTrackIt = mainVertexTrackBegin;
         mainVertexTrackIt != mainVertexTrackEnd; ++mainVertexTrackIt) {

      // Get handle of main vertex track.
      const Index<rec::VertexTrack>& mainVertexTrackIndex = *mainVertexTrackIt;
      const rec::VertexTrack& mainVertexTrack = recEvent.Get<rec::VertexTrack>(mainVertexTrackIndex);
      if (mainVertexTrack.GetStatus()!=0)
        continue;
      const Index<rec::Vertex> startVertexIndex = mainVertexTrack.GetStartVertexIndex();
      if ( !recEvent.Has<rec::Vertex>(startVertexIndex) )
        continue;
      const rec::Vertex& startVertex = recEvent.Get(startVertexIndex);
      if ( startVertex.GetIndex()!=mainVertex.GetIndex() )
        continue;

      if (mainVertexTrack.HasTrack()) {
        const Index<rec::Track> trackIndex = mainVertexTrack.GetTrackIndex();
        const rec::Track& track = recEvent.Get<rec::Track>(trackIndex);
        modutils::TrackAutopsy::TrackAutopsy dissectedTrack;
        trackSegmentFitter.SetTarget(dissectedTrack);

        /// Simple loop for sorting the clusters of a track belonging to different detectors
        for (rec::ClusterIndexIterator it = track.ClustersBegin(); it != track.ClustersEnd(); ++it) {
          const Index<rec::Cluster> clusterIndex = *it;
          const rec::Cluster& cluster = recEvent.Get<rec::Cluster>(clusterIndex);
          const int status = cluster.GetStatus();
          if ( !(status&eUsedInFit) )
            continue;
          const det::TPCConst::EId tpcId = cluster.GetTPCId();
          if ( tpcId == det::TPCConst::eGTPC ) {
            const Point& clusterPosition = cluster.GetPosition();
            const Triple clusterPositionCoords = clusterPosition.GetCoordinates(gtpcSector1CS);
            const double clusterLocalX = clusterPositionCoords.get<0>();
            if (fLocalXLeftMarginGTPCBeamSpot<=clusterLocalX && clusterLocalX<=fLocalXRightMarginGTPCBeamSpot)
              continue;
          }
          trackSegmentFitter.StoreCluster(clusterIndex);
        }

        trackSegmentFitter.RefitSegments();
        for(int i=0; i!= modutils::TrackAutopsy::TrackAutopsy::eNTPC; i++)
        {
          if(!dissectedTrack.fTrackSegment[i])
            continue;
          if(dissectedTrack.fClusterIndices[i].size()<fNMinClusters[i])
            continue;
          trackRegister.CreateAndAddTrackToRecEvent(dissectedTrack.fTrackParameters[i], dissectedTrack.fClusterIndices[i]);
        }

        for (evt::rec::TOFMassIndexIterator tofMassIt = mainVertexTrack.TOFMassBegin();
            tofMassIt!= mainVertexTrack.TOFMassEnd(); ++tofMassIt)
        {
          const Index<rec::TOFMass>& tofMassIndex = *tofMassIt;
          const rec::TOFMass& tofMass = recEvent.Get<rec::TOFMass>(tofMassIndex);
          const det::TOFConst::EId tofId = tofMass.GetTOFId();
          if ( tofId!=det::TOFConst::eTOFL && tofId!=det::TOFConst::eTOFR )
            continue;
//          const int scintillatorId = tofMass.GetScintillatorId();
          const double timeOfFlight = tofMass.GetTimeOfFlight();
          const double normalizedCharge = tofMass.GetCharge();
//          const utl::Triple& coords = tofMass.GetPosition().GetCoordinates(na61CoordSys);
//          const double x = coords.get<0>();
//          const double y = coords.get<1>();
//          const double z = coords.get<2>();
//          const utl::Triple& detCoords = detector.GetTOF().GetWall(tofId).GetScintillator(scintillatorId).GetCenterPosition().GetCoordinates(na61CoordSys);
//          const double detX = detCoords.get<0>();
//          const double detY = detCoords.get<1>();
//          const double detZ = detCoords.get<2>();
//          bool isAccepted = false;
          if ( ( fTimeOfFlightMin<=timeOfFlight && timeOfFlight<=fTimeOfFlightMax ) && ( fNormalizedChargeMin<=normalizedCharge && normalizedCharge<=fNormalizedChargeMax ) ) {
            dissectedTrack.fTofMasses.push_back(*tofMassIt);
//            isAccepted = true;
          }
//          cout << "QQRIQ " << det::TOFConst::GetName(tofId) << " " << scintillatorId << " " << (x/cm) << " " << (detX/cm) << " " << (y/cm) << " " << (detY/cm) << " " << (z/cm) << " " << (detZ/cm) << " " << normalizedCharge << " " << (timeOfFlight/nanosecond) << " " << isAccepted << endl;
        }
        trackMatchHandler.DissectedTrackHandler(dissectedTrack, fData);
      }
    }
    return eSuccess;
  }


  /// Finish function for DriftVelocityFromTOFPK
  VModule::EResultFlag
  DriftVelocityFromTOFPK::Finish()
  {
    if (!fTrackMatchOutputFileName.empty()) {
      TFile outputFile(fTrackMatchOutputFileName.c_str(), "RECREATE");
      for(unsigned int i = 0; i != modutils::TrackAutopsy::TrackMatchData::eNTrees; i++)
      {
        fData.fTrackMatchingTrees[i].Write();
      }
      outputFile.Close();
    }

    INFO("Finished.\n");

    return eSuccess;
  }

}
