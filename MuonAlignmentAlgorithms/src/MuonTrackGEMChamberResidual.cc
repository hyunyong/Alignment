#include "Alignment/MuonAlignmentAlgorithms/interface/MuonTrackGEMChamberResidual.h"


MuonTrackGEMChamberResidual::MuonTrackGEMChamberResidual(edm::ESHandle<GlobalTrackingGeometry> globalGeometry, AlignableNavigator *navigator,
                                                         DetId chamberId, AlignableDetOrUnitPtr chamberAlignable)
  : MuonChamberResidual(globalGeometry, navigator, chamberId, chamberAlignable)
{
  m_type = MuonChamberResidual::kGEM;
  align::GlobalVector zDirection(0., 0., 1.);
  m_sign = m_globalGeometry->idToDet(m_chamberId)->toLocal(zDirection).z() > 0. ? 1. : -1.;
}


void MuonTrackGEMChamberResidual::setSegmentResidual(const reco::MuonChamberMatch *trk, const reco::MuonSegmentMatch *seg)
{
  GEMDetId id(trk->id.rawId());

  GEMSegmentRef segmentGEM = seg->GEMSegmentRef;
  if (segmentGEM.get() != nullptr)
  {
    const GEMSegment* segment = segmentGEM.get();
    m_numHits = segment->nRecHits();
    m_ndof = segment->degreesOfFreedom();
    m_chi2 = segment->chi2();
  }

  align::LocalPoint l_seg(seg->x, seg->y, 0.);
  align::LocalPoint l_trk(trk->x, trk->y, 0.);

  m_residual = trk->x-seg->x;
  m_residual_error = sqrt( pow(trk->xErr, 2) + pow(seg->xErr, 2) );
  m_resslope = trk->dXdZ - seg->dXdZ;
  m_resslope_error = sqrt( pow(trk->dXdZErr, 2) + pow(seg->dXdZErr, 2) );
  
  m_trackx = trk->x;
  m_tracky = trk->y;
  m_trackdxdz = trk->dXdZ;
  m_trackdydz = trk->dYdZ;

  m_segx = seg->x;
  m_segy = seg->y;
  m_segdxdz = seg->dXdZ;
  m_segdydz = seg->dYdZ;

}
