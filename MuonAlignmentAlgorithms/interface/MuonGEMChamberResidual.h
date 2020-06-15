#ifndef Alignment_MuonAlignmentAlgorithms_MuonGEMChamberResidual_H
#define Alignment_MuonAlignmentAlgorithms_MuonGEMChamberResidual_H

#include "Alignment/MuonAlignmentAlgorithms/interface/MuonHitsChamberResidual.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"

class MuonGEMChamberResidual: public MuonHitsChamberResidual
{
public:
  MuonGEMChamberResidual(edm::ESHandle<GlobalTrackingGeometry> globalGeometry, AlignableNavigator *navigator,
                         DetId chamberId, AlignableDetOrUnitPtr chamberAlignable);

  void addResidual(edm::ESHandle<Propagator> prop, const TrajectoryStateOnSurface *tsos, const TrackingRecHit *hit,double, double) override;

  void setSegmentResidual(const reco::MuonChamberMatch *, const reco::MuonSegmentMatch *) override {}
};

#endif // Alignment_MuonAlignmentAlgorithms_MuonGEMChamberResidual_H
