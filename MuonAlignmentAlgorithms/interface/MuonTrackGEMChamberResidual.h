#ifndef Alignment_MuonAlignmentAlgorithms_MuonTrackGEMChamberResidual_H
#define Alignment_MuonAlignmentAlgorithms_MuonTrackGEMChamberResidual_H
#include "Alignment/MuonAlignmentAlgorithms/interface/MuonChamberResidual.h"

class MuonTrackGEMChamberResidual: public MuonChamberResidual
{
public:
  MuonTrackGEMChamberResidual(edm::ESHandle<GlobalTrackingGeometry> globalGeometry, AlignableNavigator *navigator,
                         DetId chamberId, AlignableDetOrUnitPtr chamberAlignable);
  
  void addResidual(edm::ESHandle<Propagator> prop, const TrajectoryStateOnSurface *tsos, const TrackingRecHit *hit,double, double) override = 0;

  void setSegmentResidual(const reco::MuonChamberMatch *, const reco::MuonSegmentMatch *) override;
};

#endif // Alignment_MuonAlignmentAlgorithms_MuonTrackGEMChamberResidual_H
