#include "twoTrackResonanceReco.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackSeed.h>
#include <trackbase_historic/TrackAnalysisUtils.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainer.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <ffamodules/CDBInterface.h>

#include <trackreco/ActsPropagator.h>

#include <Math/Vector4D.h>

#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <TParticlePDG.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <vector>

namespace
{
  // Angle between two vectors; atan2 form is accurate for small angles and independent of the vector lengths
  double openingAngle(const Vec3 &a, const Vec3 &b)
  {
    const Vec3 cross{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    return std::atan2(norm(cross), dot(a, b));
  }
}  // namespace

//____________________________________________________________________________..
twoTrackResonanceReco::twoTrackResonanceReco(const std::string &name)
  : SubsysReco(name)
{
}

//____________________________________________________________________________..
twoTrackResonanceReco::~twoTrackResonanceReco() = default;

//____________________________________________________________________________..
int twoTrackResonanceReco::Init(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int twoTrackResonanceReco::InitRun(PHCompositeNode * /*topNode*/)
{
  //m_outfile = new TFile(m_outfile_name.c_str(), "recreate");
  //initTree();

  if (m_use_dEdx_pid)
  {
    init_dEdx_bands();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
bool twoTrackResonanceReco::propagateToPoint(SvtxTrack *track, const Vec3 &target, Vec3 &pos, Vec3 &mom) const
{
  ActsPropagator actsPropagator(m_acts_geometry);

  Acts::Vector3 targetVec(target.x, target.y, target.z);  // cm; ActsPropagator converts internally
  auto surface = actsPropagator.makeVertexSurface(targetVec);

  auto paramsResult = actsPropagator.makeTrackParams(track, m_vertexmap);
  if (!paramsResult.ok())
  {
    return false;
  }

  auto propResult = actsPropagator.propagateTrackFast(paramsResult.value(), surface);
  if (!propResult.ok())
  {
    return false;
  }

  const auto &finalParams = propResult.value().second;
  auto position = finalParams.position(m_acts_geometry->geometry().getGeoContext());
  auto momentum = finalParams.momentum();

  pos = {position.x() / Acts::UnitConstants::cm, position.y() / Acts::UnitConstants::cm, position.z() / Acts::UnitConstants::cm};
  mom = {momentum.x(), momentum.y(), momentum.z()};

  return true;
}

//____________________________________________________________________________..
bool twoTrackResonanceReco::buildSV(SvtxTrack *trackA, SvtxTrack *trackB, Vec3 &vertex, Vec3 &momentumA, Vec3 &momentumB, double &dca) const
{
  Vec3 target{0.5 * (trackA->get_x() + trackB->get_x()),
              0.5 * (trackA->get_y() + trackB->get_y()),
              0.5 * (trackA->get_z() + trackB->get_z())};

  constexpr int kMaxIterations = 10;
  constexpr double kConvergenceTolerance = 1e-4;  // cm

  Vec3 posA, posB;
  for (int iter = 0; iter < kMaxIterations; ++iter)
  {
    if (!propagateToPoint(trackA, target, posA, momentumA) ||
        !propagateToPoint(trackB, target, posB, momentumB))
    {
      return false;
    }

    Vec3 newTarget = 0.5 * (posA + posB);
    double step = norm(newTarget - target);
    target = newTarget;

    if (step < kConvergenceTolerance)
    {
      break;
    }
  }

  vertex = target;
  dca = norm(posA - posB);

  return true;
}

//____________________________________________________________________________..
bool twoTrackResonanceReco::trackToVertexDCA(SvtxTrack *track, const Vec3 &vertex, double &dca) const
{
  Vec3 pos, mom;
  if (!propagateToPoint(track, vertex, pos, mom))
  {
    return false;
  }

  dca = norm(pos - vertex);
  return true;
}

//____________________________________________________________________________..
bool twoTrackResonanceReco::hasSiliconClusters(SvtxTrack *track) const
{
  auto counts = TrackAnalysisUtils::get_cluster_counts(track);
  int mvtxStates = std::get<0>(counts);
  int inttStates = std::get<1>(counts);
  return (mvtxStates > 0 || inttStates > 0);
}

//____________________________________________________________________________..
const SvtxVertex *twoTrackResonanceReco::findMatchingVertex(short int crossing, double referenceZ, int *nMatching) const
{
  const SvtxVertex *best = nullptr;
  double bestMetric = std::numeric_limits<double>::max();
  int n = 0;

  for (auto &iter : *m_vertexmap)
  {
    SvtxVertex *vertex = iter.second;
    if (vertex->get_beam_crossing() != crossing)
    {
      continue;
    }
    ++n;

    double metric;
    if (std::isfinite(referenceZ))
    {
      metric = std::fabs(vertex->get_z() - referenceZ);
    }
    else
    {
      metric = (vertex->get_ndof() > 0) ? vertex->get_chisq() / vertex->get_ndof()
                                         : std::numeric_limits<double>::max();
    }

    if (metric < bestMetric || !best)
    {
      bestMetric = metric;
      best = vertex;
    }
  }

  if (nMatching)
  {
    *nMatching = n;
  }
  return best;
}

//____________________________________________________________________________..
int twoTrackResonanceReco::getNodes(PHCompositeNode *topNode)
{
  m_trackmap = findNode::getClass<SvtxTrackMap>(topNode, m_trackmap_node_name);
  if (!m_trackmap)
  {
    std::cout << PHWHERE << " Could not find " << m_trackmap_node_name << ", aborting event" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_vertexmap = findNode::getClass<SvtxVertexMap>(topNode, m_vertexmap_node_name);
  if (!m_vertexmap)
  {
    std::cout << PHWHERE << " Could not find " << m_vertexmap_node_name << ", aborting event" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_acts_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");
  if (!m_acts_geometry)
  {
    std::cout << PHWHERE << " Could not find ActsGeometry, aborting event" << std::endl;
    return Fun4AllReturnCodes::ABORTEVENT;
  }

  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_cluster_map)
  {
    m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_SEED");
  }
  m_geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");

  if (!m_cluster_map || !m_geom_container)
  {
    m_can_get_dEdx = false;
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void twoTrackResonanceReco::init_dEdx_bands()
{
  std::string dedx_fitparams = CDBInterface::instance()->getUrl("TPC_DEDX_FITPARAM");
  std::cout << PHWHERE << " opening " << dedx_fitparams << std::endl;

  TFile *filefit = TFile::Open(dedx_fitparams.c_str());
  if (!filefit || !filefit->IsOpen())
  {
    std::cerr << PHWHERE << " Error opening the dE/dx calibration file " << dedx_fitparams << std::endl;
    return;
  }

  filefit->GetObject("f_piband", f_pion_plus);
  filefit->GetObject("f_Kband", f_kaon_plus);
  filefit->GetObject("f_pband", f_proton_plus);
  filefit->GetObject("f_piminus_band", f_pion_minus);
  filefit->GetObject("f_Kminus_band", f_kaon_minus);
  filefit->GetObject("f_pbar_band", f_proton_minus);

  m_dEdx_bandMap.insert({-11, f_pion_plus});
  m_dEdx_bandMap.insert({211, f_pion_plus});
  m_dEdx_bandMap.insert({321, f_kaon_plus});
  m_dEdx_bandMap.insert({2212, f_proton_plus});
  m_dEdx_bandMap.insert({11, f_pion_minus});
  m_dEdx_bandMap.insert({-211, f_pion_minus});
  m_dEdx_bandMap.insert({-321, f_kaon_minus});
  m_dEdx_bandMap.insert({-2212, f_proton_minus});

  filefit->Close();
  delete filefit;
}

//____________________________________________________________________________..
double twoTrackResonanceReco::measured_dEdx(SvtxTrack *track) const
{
  TrackSeed *tpcSeed = track->get_tpc_seed();
  if (!tpcSeed)
  {
    return -1;
  }

  float layerThicknesses[4] = {0.0, 0.0, 0.0, 0.0};
  layerThicknesses[0] = m_geom_container->GetLayerCellGeom(7)->get_thickness();
  layerThicknesses[1] = m_geom_container->GetLayerCellGeom(8)->get_thickness();
  layerThicknesses[2] = m_geom_container->GetLayerCellGeom(27)->get_thickness();
  layerThicknesses[3] = m_geom_container->GetLayerCellGeom(50)->get_thickness();

  return TrackAnalysisUtils::calc_dedx(tpcSeed, m_cluster_map, m_acts_geometry, layerThicknesses);
}

//____________________________________________________________________________..
void twoTrackResonanceReco::initTree()
{
  m_outfile->cd();
  m_tree = new TTree("DecayTree", "DecayTree");


  m_tree->Branch("mother_mass", &b_mother_mass, "mother_mass/F");
  m_tree->Branch("mother_pT", &b_mother_pT, "mother_pT/F");
  m_tree->Branch("mother_eta", &b_mother_eta, "mother_eta/F");
  m_tree->Branch("mother_phi", &b_mother_phi, "mother_phi/F");
  m_tree->Branch("mother_DIRA", &b_mother_DIRA, "mother_DIRA/F");
  m_tree->Branch("mother_flight_distance", &b_mother_flight_distance, "mother_flight_distance/F");
  m_tree->Branch("n_matching_vertices", &b_n_matching_vertices, "n_matching_vertices/I");
  m_tree->Branch("mother_PV_DCA", &b_mother_PV_DCA, "mother_PV_DCA/F");
  m_tree->Branch("mother_DIRA_acts", &b_mother_DIRA_acts, "mother_DIRA_acts/F");
  m_tree->Branch("mother_PV_DCA_acts", &b_mother_PV_DCA_acts, "mother_PV_DCA_acts/F");
  m_tree->Branch("mother_opening_angle", &b_mother_opening_angle, "mother_opening_angle/F");
  m_tree->Branch("mother_opening_angle_acts", &b_mother_opening_angle_acts, "mother_opening_angle_acts/F");

  m_tree->Branch("daughter1_mass", &b_daughter1_mass, "daughter1_mass/F");
  m_tree->Branch("daughter1_charge", &b_daughter1_charge, "daughter1_charge/I");
  m_tree->Branch("daughter1_pT", &b_daughter1_pT, "daughter1_pT/F");
  m_tree->Branch("daughter1_eta", &b_daughter1_eta, "daughter1_eta/F");
  m_tree->Branch("daughter1_phi", &b_daughter1_phi, "daughter1_phi/F");
  m_tree->Branch("daughter1_PV_DCA", &b_daughter1_PV_DCA, "daughter1_PV_DCA/F");
  m_tree->Branch("daughter1_dEdx", &b_daughter1_dEdx, "daughter1_dEdx/F");
  m_tree->Branch("daughter1_chi2_per_ndf", &b_daughter1_quality, "daughter1_chi2_per_ndf/F");
  m_tree->Branch("daughter1_p", &b_daughter1_p, "daughter1_p/F");
  m_tree->Branch("daughter1_p_acts", &b_daughter1_p_acts, "daughter1_p_acts/F");
  m_tree->Branch("daughter1_momentum_source", &b_daughter1_momentum_source, "daughter1_momentum_source/I");

  m_tree->Branch("daughter2_mass", &b_daughter2_mass, "daughter2_mass/F");
  m_tree->Branch("daughter2_charge", &b_daughter2_charge, "daughter2_charge/I");
  m_tree->Branch("daughter2_pT", &b_daughter2_pT, "daughter2_pT/F");
  m_tree->Branch("daughter2_eta", &b_daughter2_eta, "daughter2_eta/F");
  m_tree->Branch("daughter2_phi", &b_daughter2_phi, "daughter2_phi/F");
  m_tree->Branch("daughter2_PV_DCA", &b_daughter2_PV_DCA, "daughter2_PV_DCA/F");
  m_tree->Branch("daughter2_dEdx", &b_daughter2_dEdx, "daughter2_dEdx/F");
  m_tree->Branch("daughter2_chi2_per_ndf", &b_daughter2_quality, "daughter2_chi2_per_ndf/F");
  m_tree->Branch("daughter2_p", &b_daughter2_p, "daughter2_p/F");
  m_tree->Branch("daughter2_p_acts", &b_daughter2_p_acts, "daughter2_p_acts/F");
  m_tree->Branch("daughter2_momentum_source", &b_daughter2_momentum_source, "daughter2_momentum_source/I");

  m_tree->Branch("track_to_track_DCA", &b_track_to_track_DCA, "track_to_track_DCA/F");
  m_tree->Branch("both_charge_states_passed", &b_both_charge_states_passed, "both_charge_states_passed/O");

  m_tree->Branch("crossing", &b_crossing, "crossing/S");
  m_tree->Branch("PV_x", &b_PV_x, "PV_x/F");
  m_tree->Branch("PV_y", &b_PV_y, "PV_y/F");
  m_tree->Branch("PV_z", &b_PV_z, "PV_z/F");
  m_tree->Branch("SV_x", &b_SV_x, "SV_x/F");
  m_tree->Branch("SV_y", &b_SV_y, "SV_y/F");
  m_tree->Branch("SV_z", &b_SV_z, "SV_z/F");
}

//____________________________________________________________________________..
void twoTrackResonanceReco::resetBranches()
{
  b_daughter1_mass = 0;
  b_daughter2_mass = 0;
  b_PV_x = b_PV_y = b_PV_z = 0;
  b_SV_x = b_SV_y = b_SV_z = 0;
  b_mother_mass = 0;
  b_mother_pT = 0;
  b_mother_eta = 0;
  b_mother_phi = 0;
  b_mother_DIRA = 0;
  b_mother_flight_distance = 0;
  b_n_matching_vertices = 0;
  b_mother_PV_DCA = 0;
  b_daughter1_charge = 0;
  b_daughter1_pT = 0;
  b_daughter1_eta = 0;
  b_daughter1_phi = 0;
  b_daughter1_PV_DCA = 0;
  b_daughter1_dEdx = -1;
  b_daughter1_quality = -1;
  b_daughter2_charge = 0;
  b_daughter2_pT = 0;
  b_daughter2_eta = 0;
  b_daughter2_phi = 0;
  b_daughter2_PV_DCA = 0;
  b_daughter2_dEdx = -1;
  b_daughter2_quality = -1;
  b_mother_DIRA_acts = 0;
  b_mother_PV_DCA_acts = 0;
  b_mother_opening_angle = 0;
  b_mother_opening_angle_acts = 0;
  b_daughter1_p = b_daughter1_p_acts = 0;
  b_daughter2_p = b_daughter2_p_acts = 0;
  b_daughter1_momentum_source = 0;
  b_daughter2_momentum_source = 0;
  b_track_to_track_DCA = 0;
  b_both_charge_states_passed = false;
}

//____________________________________________________________________________..
int twoTrackResonanceReco::process_event(PHCompositeNode *topNode)
{
  int nodeStatus = getNodes(topNode);
  if (nodeStatus != Fun4AllReturnCodes::EVENT_OK)
  {
    return nodeStatus;
  }

  if (m_trackmap->size() < 2)
  {
    return Fun4AllReturnCodes::EVENT_OK;
  }

  const double daughter1_mass = TDatabasePDG::Instance()->GetParticle(std::abs(m_daughter1_PDGID))->Mass();
  const double daughter2_mass = TDatabasePDG::Instance()->GetParticle(std::abs(m_daughter2_PDGID))->Mass();

  std::vector<SvtxTrack *> goodTracks;
  for (auto &iter : *m_trackmap)
  {
    SvtxTrack *track = iter.second;

    if (!hasSiliconClusters(track))
    {
      continue;
    }

    if (track->get_ndf() == 0 || track->get_chisq()/track->get_ndf() > m_max_track_chi2_per_ndf)
    {
      continue;
    }

    if (track->get_pt() < m_min_daughter_pt)
    {
      continue;
    }

    goodTracks.push_back(track);
  }

  for (std::size_t i = 0; i < goodTracks.size(); ++i)
  {
    for (std::size_t j = i + 1; j < goodTracks.size(); ++j)
    {
      SvtxTrack *trackA = goodTracks[i];
      SvtxTrack *trackB = goodTracks[j];

      if (trackA->get_charge() == trackB->get_charge())
      {
        continue;
      }

      if (trackA->get_crossing() != trackB->get_crossing())
      {
        continue;
      }
      short int crossing = trackA->get_crossing();

      // quick check that this crossing has a vertex at all; the vertex actually used is chosen once the SV is known
      if (!findMatchingVertex(crossing))
      {
        continue;
      }

      Vec3 sv, momentumA, momentumB;
      double daughterDCA;
      if (!buildSV(trackA, trackB, sv, momentumA, momentumB, daughterDCA))
      {
        continue;
      }

      if (daughterDCA > m_track_to_track_DCA_cut)
      {
        continue;
      }

      // The primary vertex is the same-crossing vertex closest in z to the secondary vertex, so a pair is compared with
      // the vertex it plausibly came from, not with whichever vertex of the crossing has the best chi2/ndf.
      int nMatchingVertices = 0;
      const SvtxVertex *primaryVertex = findMatchingVertex(crossing, m_match_vertex_to_sv ? sv.z : std::numeric_limits<double>::quiet_NaN(), &nMatchingVertices);
      if (!primaryVertex)
      {
        continue;
      }

      Vec3 pv{primaryVertex->get_x(), primaryVertex->get_y(), primaryVertex->get_z()};
      Vec3 flight = sv - pv;
      double flightDistance = norm(flight);

      if (flightDistance < m_flight_distance_cut || flightDistance > m_max_flight_distance_cut)
      {
        continue;
      }

      // Momentum used for the cuts, the mass and the stored kinematics: |p| from the TPC seed pT, direction from the
      // silicon (MVTX+INTT) seed evaluated at the secondary vertex. Each piece falls back to the ACTS momentum if it is
      // unavailable (see refineMomentum). The momentum does not depend on the mass hypothesis, so it is computed once
      // per track pair, and the mother direction used for DIRA and the mother IP is built from these same vectors.
      Vec3 refinedMomentumA, refinedMomentumB;
      const int sourceA = refineMomentum(trackA, sv, momentumA, refinedMomentumA);
      const int sourceB = refineMomentum(trackB, sv, momentumB, refinedMomentumB);

      Vec3 motherMomentum = refinedMomentumA + refinedMomentumB;
      if (norm(motherMomentum) <= 0)
      {
        continue;
      }
      Vec3 motherDir = (1.0 / norm(motherMomentum)) * motherMomentum;

      double dira = (flightDistance > 0) ? dot(motherDir, flight) / flightDistance : 0;

      if (dira < m_dira_cut)
      {
        continue;
      }

      Vec3 toVertex = pv - sv;
      double along = dot(toVertex, motherDir);
      double motherIP = norm(toVertex - along * motherDir);

      if (motherIP > m_mother_PV_DCA_cut)
      {
        continue;
      }

      // The same two quantities from the ACTS fit momenta, stored for comparison only (no cut is applied to them)
      double diraActs = 0;
      double motherIPActs = 0;
      {
        const Vec3 motherMomentumActs = momentumA + momentumB;
        const double motherPActs = norm(motherMomentumActs);
        if (motherPActs > 0)
        {
          const Vec3 motherDirActs = (1.0 / motherPActs) * motherMomentumActs;
          diraActs = (flightDistance > 0) ? dot(motherDirActs, flight) / flightDistance : 0;
          const double alongActs = dot(toVertex, motherDirActs);
          motherIPActs = norm(toVertex - alongActs * motherDirActs);
        }
      }

      double pvDcaA, pvDcaB;
      if (!trackToVertexDCA(trackA, pv, pvDcaA) || !trackToVertexDCA(trackB, pv, pvDcaB))
      {
        continue;
      }

      if (std::min(pvDcaA, pvDcaB) < m_daughter_PV_DCA_cut)
      {
        continue;
      }

      double dedxA = m_can_get_dEdx ? measured_dEdx(trackA) : -1;
      double dedxB = m_can_get_dEdx ? measured_dEdx(trackB) : -1;

      bool sameSpecies = (std::abs(m_daughter1_PDGID) == std::abs(m_daughter2_PDGID));
      bool aIsPositive = (trackA->get_charge() > 0);

      std::vector<bool> hypotheses;
      if (sameSpecies)
      {
        hypotheses.push_back(aIsPositive);
      }
      else
      {
        hypotheses.push_back(true);   // daughter1 (pdgID1 species) = trackA
        hypotheses.push_back(false);  // daughter1 (pdgID1 species) = trackB, i.e. the charge conjugate
      }

      // Candidate is declared in the header.
      std::vector<Candidate> candidates;
      for (bool daughter1IsA : hypotheses)
      {
        Candidate c;
        c.daughter1Track = daughter1IsA ? trackA : trackB;
        c.daughter2Track = daughter1IsA ? trackB : trackA;
        const Vec3 &momentum1 = daughter1IsA ? refinedMomentumA : refinedMomentumB;
        const Vec3 &momentum2 = daughter1IsA ? refinedMomentumB : refinedMomentumA;
        const Vec3 &acts_momentum1 = daughter1IsA ? momentumA : momentumB;
        const Vec3 &acts_momentum2 = daughter1IsA ? momentumB : momentumA;
        c.daughter1PvDca = daughter1IsA ? pvDcaA : pvDcaB;
        c.daughter2PvDca = daughter1IsA ? pvDcaB : pvDcaA;
        c.dedx1 = daughter1IsA ? dedxA : dedxB;
        c.dedx2 = daughter1IsA ? dedxB : dedxA;

        c.daughter1Vec = ROOT::Math::PxPyPzMVector(momentum1.x, momentum1.y, momentum1.z, daughter1_mass);
        c.daughter2Vec = ROOT::Math::PxPyPzMVector(momentum2.x, momentum2.y, momentum2.z, daughter2_mass);

        // With each momentum built as |p_TPC| along the silicon-seed direction at the SV, the invariant mass of the pair
        // is exactly m^2 = m1^2 + m2^2 + 2 (E1 E2 - p1 p2 cos(opening angle)), so no separate formula is needed and
        // the mother pT/eta/phi branches stay consistent with the mass.
        c.motherMass = (c.daughter1Vec + c.daughter2Vec).M();

        c.openingAngle = openingAngle(momentum1, momentum2);
        c.openingAngleActs = openingAngle(acts_momentum1, acts_momentum2);
        c.p1 = norm(momentum1);
        c.p2 = norm(momentum2);
        c.p1Acts = norm(acts_momentum1);
        c.p2Acts = norm(acts_momentum2);
        c.source1 = daughter1IsA ? sourceA : sourceB;
        c.source2 = daughter1IsA ? sourceB : sourceA;

        if (c.motherMass >= m_min_mass && c.motherMass <= m_max_mass)
        {
          c.passed = true;

          if (m_use_dEdx_pid)
          {
            int daughter1_signedPDG = static_cast<int>(c.daughter1Track->get_charge()) * std::abs(m_daughter1_PDGID);
            int daughter2_signedPDG = static_cast<int>(c.daughter2Track->get_charge()) * std::abs(m_daughter2_PDGID);

            auto it1 = m_dEdx_bandMap.find(daughter1_signedPDG);
            auto it2 = m_dEdx_bandMap.find(daughter2_signedPDG);
            if (it1 != m_dEdx_bandMap.end() && c.dedx1 > 0)
            {
              double expected = it1->second->Eval(c.daughter1Track->get_p());
              if (std::fabs(c.dedx1 - expected) > m_dEdx_band_width * expected)
              {
                c.passed = false;
              }
            }
            if (it2 != m_dEdx_bandMap.end() && c.dedx2 > 0)
            {
              double expected = it2->second->Eval(c.daughter2Track->get_p());
              if (std::fabs(c.dedx2 - expected) > m_dEdx_band_width * expected)
              {
                c.passed = false;
              }
            }
          }
        }

        candidates.push_back(c);
      }

      bool bothChargeStatesPassed = (!sameSpecies && candidates.size() == 2 && candidates[0].passed && candidates[1].passed);

      for (auto &c : candidates)
      {
        if (!c.passed)
        {
          continue;
        }

        if (candidateCounter < 1)
        {
          m_outfile = new TFile(m_outfile_name.c_str(), "recreate");
          initTree();
        }

        ++candidateCounter;

        resetBranches();

        b_crossing = crossing;

        b_daughter1_mass = daughter1_mass;
        b_daughter2_mass = daughter2_mass;

        b_PV_x = pv.x;
        b_PV_y = pv.y;
        b_PV_z = pv.z;
        b_SV_x = sv.x;
        b_SV_y = sv.y;
        b_SV_z = sv.z;

        auto motherVec = c.daughter1Vec + c.daughter2Vec;
        b_mother_mass = c.motherMass;
        b_mother_pT = motherVec.Pt();
        b_mother_eta = motherVec.Eta();
        b_mother_phi = motherVec.Phi();
        b_mother_DIRA = dira;
        b_mother_flight_distance = flightDistance;
        b_n_matching_vertices = nMatchingVertices;
        b_mother_PV_DCA = motherIP;
        b_mother_DIRA_acts = diraActs;
        b_mother_PV_DCA_acts = motherIPActs;
        b_mother_opening_angle = c.openingAngle;
        b_mother_opening_angle_acts = c.openingAngleActs;

        b_daughter1_charge = c.daughter1Track->get_charge();
        b_daughter1_pT = c.daughter1Track->get_pt();
        b_daughter1_eta = c.daughter1Track->get_eta();
        b_daughter1_phi = c.daughter1Track->get_phi();
        b_daughter1_PV_DCA = c.daughter1PvDca;
        b_daughter1_dEdx = c.dedx1;
        b_daughter1_quality = c.daughter1Track->get_chisq()/c.daughter1Track->get_ndf();
        b_daughter1_p = c.p1;
        b_daughter1_p_acts = c.p1Acts;
        b_daughter1_momentum_source = c.source1;

        b_daughter2_charge = c.daughter2Track->get_charge();
        b_daughter2_pT = c.daughter2Track->get_pt();
        b_daughter2_eta = c.daughter2Track->get_eta();
        b_daughter2_phi = c.daughter2Track->get_phi();
        b_daughter2_PV_DCA = c.daughter2PvDca;
        b_daughter2_dEdx = c.dedx2;
        b_daughter2_quality = c.daughter2Track->get_chisq()/c.daughter2Track->get_ndf();
        b_daughter2_p = c.p2;
        b_daughter2_p_acts = c.p2Acts;
        b_daughter2_momentum_source = c.source2;

        b_track_to_track_DCA = daughterDCA;
        b_both_charge_states_passed = bothChargeStatesPassed;

        m_tree->Fill();
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

bool twoTrackResonanceReco::siliconDirectionAtSV(const SvtxTrack *track, const Vec3 &sv, Vec3 &dir) const
{
  const TrackSeed *siSeed = track->get_silicon_seed();
  if (!siSeed)
  {
    return false;
  }

  // Circle (centre X0,Y0 in the transverse plane), stored azimuth at the point of closest approach to the origin,
  // and polar angle from the dz/dr slope. All come from a fit to the silicon (MVTX + INTT) clusters only.
  const double X0 = siSeed->get_X0();
  const double Y0 = siSeed->get_Y0();
  const double storedPhi = siSeed->get_phi();
  const double theta = siSeed->get_theta();
  if (!std::isfinite(X0) || !std::isfinite(Y0) || !std::isfinite(storedPhi) || !std::isfinite(theta))
  {
    return false;
  }

  // Moving along a circle rotates the tangent by the same angle as the radius vector from the centre, so the
  // tangent at the SV is the radial angle of the SV +/- 90 degrees. The sign (clockwise or counter-clockwise) is the
  // one whose result lies closest to the stored azimuth, which is evaluated close to the SV.
  const double alpha = std::atan2(sv.y - Y0, sv.x - X0);
  const double phiCCW = alpha + M_PI / 2.0;
  const double phiCW = alpha - M_PI / 2.0;
  const double phi = (std::fabs(std::remainder(phiCCW - storedPhi, 2.0 * M_PI)) <
                      std::fabs(std::remainder(phiCW - storedPhi, 2.0 * M_PI)))
                         ? phiCCW
                         : phiCW;

  // The polar angle does not change along a helix
  dir = Vec3{std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)};
  return true;
}

//____________________________________________________________________________..
int twoTrackResonanceReco::refineMomentum(const SvtxTrack *track, const Vec3 &sv, const Vec3 &actsMomentum, Vec3 &out) const
{
  int source = 0;
  const double actsP = norm(actsMomentum);

  // direction: ACTS fit at the SV unless the silicon seed is usable
  Vec3 dir = (actsP > 0) ? (1.0 / actsP) * actsMomentum : Vec3{0, 0, 1};
  if (m_use_silicon_direction)
  {
    Vec3 siliconDir;
    if (siliconDirectionAtSV(track, sv, siliconDir))
    {
      dir = siliconDir;
      source |= 1;
    }
  }

  // magnitude: transverse momentum from the TPC seed, converted to |p| with the polar angle of the direction above
  double p = actsP;
  if (m_use_tpc_momentum)
  {
    const TrackSeed *tpcSeed = track->get_tpc_seed();
    if (tpcSeed)
    {
      const double pt = tpcSeed->get_pt();
      const double sinTheta = std::sqrt(dir.x * dir.x + dir.y * dir.y);
      if (std::isfinite(pt) && pt > 0 && sinTheta > 1e-6)
      {
        p = pt / sinTheta;
        source |= 2;
      }
    }
  }

  out = p * dir;
  return source;
}

//____________________________________________________________________________..
std::array<double,3> twoTrackResonanceReco::unit(const Vec3& v)
{
  const double p = std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
  if (p <= 0)
  {
    return std::array<double,3>{0, 0, 0};
  }
  return std::array<double,3>{v.x/p, v.y/p, v.z/p};
}

//____________________________________________________________________________..
int twoTrackResonanceReco::ResetEvent(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int twoTrackResonanceReco::EndRun(const int /*runnumber*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int twoTrackResonanceReco::End(PHCompositeNode * /*topNode*/)
{
  if (candidateCounter != 0)
  {
    m_outfile->cd();
    m_tree->Write();
    m_outfile->Close();
  }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int twoTrackResonanceReco::Reset(PHCompositeNode * /*topNode*/)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void twoTrackResonanceReco::Print(const std::string &what) const
{
  std::cout << "twoTrackResonanceReco::Print(" << what << ")" << std::endl;
  std::cout << "  daughter PDG IDs: " << m_daughter1_PDGID << ", " << m_daughter2_PDGID << std::endl;
  std::cout << "  mother mass range: [" << m_min_mass << ", " << m_max_mass << "] GeV" << std::endl;
  std::cout << "  daughter DCA cut: " << m_track_to_track_DCA_cut << " cm" << std::endl;
  std::cout << "  flight distance cut: " << m_flight_distance_cut << " cm" << std::endl;
  std::cout << "  mother DIRA cut: " << m_dira_cut << std::endl;
  std::cout << "  mother IP cut: " << m_mother_PV_DCA_cut << " cm" << std::endl;
  std::cout << "  daughter IP cut: " << m_daughter_PV_DCA_cut << " cm" << std::endl;
  std::cout << "  daughter quality cut: " << m_max_track_chi2_per_ndf << std::endl;
  std::cout << "  daughter pT: " << m_min_daughter_pt << " GeV" << std::endl;
  std::cout << "  use dE/dx PID: " << (m_use_dEdx_pid ? "true" : "false") << std::endl;
}
