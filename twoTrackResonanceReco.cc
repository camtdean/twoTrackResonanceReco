#include "twoTrackResonanceReco.h"

#include <fun4all/Fun4AllReturnCodes.h>

#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <phool/phool.h>

#include <trackbase_historic/SvtxTrack.h>
#include <trackbase_historic/SvtxTrackMap.h>
#include <trackbase_historic/TrackAnalysisUtils.h>

#include <globalvertex/SvtxVertex.h>
#include <globalvertex/SvtxVertexMap.h>

#include <trackbase/ActsGeometry.h>
#include <trackbase/TrkrClusterContainer.h>

#include <g4detectors/PHG4TpcGeom.h>
#include <g4detectors/PHG4TpcGeomContainer.h>

#include <ffamodules/CDBInterface.h>

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
constexpr double kFieldConstant = 0.00299792458;  // GeV / (T * cm)

double wrapAngle(double angle)
{
  while (angle > M_PI)
  {
    angle -= 2.0 * M_PI;
  }
  while (angle <= -M_PI)
  {
    angle += 2.0 * M_PI;
  }
  return angle;
}
}

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
  m_outfile = new TFile(m_outfile_name.c_str(), "recreate");
  initTree();

  if (m_use_dEdx_pid)
  {
    init_dEdx_bands();
  }

  return Fun4AllReturnCodes::EVENT_OK;
}

twoTrackResonanceReco::SimpleHelix twoTrackResonanceReco::buildHelix(SvtxTrack *track) const
{
  SimpleHelix helix;

  helix.position = {track->get_x(), track->get_y(), track->get_z()};
  helix.momentum = {track->get_px(), track->get_py(), track->get_pz()};
  helix.charge = track->get_charge();

  double pT = std::hypot(helix.momentum.x, helix.momentum.y);

  helix.rotationSense = (helix.charge * m_Bz >= 0) ? 1.0 : -1.0;
  helix.radius = pT / (kFieldConstant * std::fabs(helix.charge) * std::fabs(m_Bz));

  helix.centerX = helix.position.x + helix.radius * helix.rotationSense * (helix.momentum.y / pT);
  helix.centerY = helix.position.y - helix.radius * helix.rotationSense * (helix.momentum.x / pT);

  return helix;
}

double twoTrackResonanceReco::helixPointDCA(const SimpleHelix &helix, const Vec3 &target, Vec3 &point, Vec3 &direction) const
{
  double toTargetX = target.x - helix.centerX;
  double toTargetY = target.y - helix.centerY;
  double distToCenter = std::hypot(toTargetX, toTargetY);

  if (distToCenter < std::numeric_limits<double>::epsilon())
  {
    distToCenter = std::numeric_limits<double>::epsilon();
  }

  double nx = toTargetX / distToCenter;
  double ny = toTargetY / distToCenter;

  point.x = helix.centerX + helix.radius * nx;
  point.y = helix.centerY + helix.radius * ny;

  double phi0 = std::atan2(helix.position.y - helix.centerY, helix.position.x - helix.centerX);
  double phi = std::atan2(ny, nx);
  double deltaPhi = wrapAngle(phi0 - phi);
  double arcLength = helix.radius * helix.rotationSense * deltaPhi;

  double pT = std::hypot(helix.momentum.x, helix.momentum.y);
  point.z = helix.position.z + (helix.momentum.z / pT) * arcLength;

  direction.x = pT * helix.rotationSense * ny;
  direction.y = -pT * helix.rotationSense * nx;
  direction.z = helix.momentum.z;

  return norm(point - target);
}

double twoTrackResonanceReco::twoHelixDCA(const SimpleHelix &helixA, const SimpleHelix &helixB, Vec3 &vertex, Vec3 &momentumA, Vec3 &momentumB) const
{
  double dx = helixB.centerX - helixA.centerX;
  double dy = helixB.centerY - helixA.centerY;
  double centerDistance = std::hypot(dx, dy);

  if (centerDistance < std::numeric_limits<double>::epsilon())
  {
    return std::numeric_limits<double>::max();
  }

  double nx = dx / centerDistance;
  double ny = dy / centerDistance;

  Vec3 pointA{helixA.centerX + helixA.radius * nx, helixA.centerY + helixA.radius * ny, 0};
  Vec3 pointB{helixB.centerX - helixB.radius * nx, helixB.centerY - helixB.radius * ny, 0};

  double phi0A = std::atan2(helixA.position.y - helixA.centerY, helixA.position.x - helixA.centerX);
  double phiA = std::atan2(ny, nx);
  double arcLengthA = helixA.radius * helixA.rotationSense * wrapAngle(phi0A - phiA);

  double phi0B = std::atan2(helixB.position.y - helixB.centerY, helixB.position.x - helixB.centerX);
  double phiB = std::atan2(-ny, -nx);
  double arcLengthB = helixB.radius * helixB.rotationSense * wrapAngle(phi0B - phiB);

  double pTA = std::hypot(helixA.momentum.x, helixA.momentum.y);
  double pTB = std::hypot(helixB.momentum.x, helixB.momentum.y);

  pointA.z = helixA.position.z + (helixA.momentum.z / pTA) * arcLengthA;
  pointB.z = helixB.position.z + (helixB.momentum.z / pTB) * arcLengthB;

  momentumA = {pTA * helixA.rotationSense * ny, -pTA * helixA.rotationSense * nx, helixA.momentum.z};
  momentumB = {-pTB * helixB.rotationSense * ny, pTB * helixB.rotationSense * nx, helixB.momentum.z};

  vertex = 0.5 * (pointA + pointB);

  return norm(pointA - pointB);
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
const SvtxVertex *twoTrackResonanceReco::findMatchingVertex(short int crossing) const
{
  const SvtxVertex *best = nullptr;
  float bestChi2NDF = std::numeric_limits<float>::max();

  for (auto &iter : *m_vertexmap)
  {
    SvtxVertex *vertex = iter.second;
    if (vertex->get_beam_crossing() != crossing)
    {
      continue;
    }

    float chi2NDF = (vertex->get_ndof() > 0) ? vertex->get_chisq() / vertex->get_ndof()
                                              : std::numeric_limits<float>::max();
    if (chi2NDF < bestChi2NDF)
    {
      bestChi2NDF = chi2NDF;
      best = vertex;
    }
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

  m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER");
  if (!m_cluster_map)
  {
    m_cluster_map = findNode::getClass<TrkrClusterContainer>(topNode, "TRKR_CLUSTER_SEED");
  }
  m_geom_container = findNode::getClass<PHG4TpcGeomContainer>(topNode, "TPCGEOMCONTAINER");
  m_acts_geometry = findNode::getClass<ActsGeometry>(topNode, "ActsGeometry");

  if (!m_cluster_map || !m_geom_container || !m_acts_geometry)
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
  m_tree->Branch("mother_PV_DCA", &b_mother_PV_DCA, "mother_PV_DCA/F");

  m_tree->Branch("daughter1_mass", &b_daughter1_mass, "daughter1_mass/F");
  m_tree->Branch("daughter1_charge", &b_daughter1_charge, "daughter1_charge/I");
  m_tree->Branch("daughter1_pT", &b_daughter1_pT, "daughter1_pT/F");
  m_tree->Branch("daughter1_eta", &b_daughter1_eta, "daughter1_eta/F");
  m_tree->Branch("daughter1_phi_beamline", &b_daughter1_phi_beamline, "daughter1_phi_beamline/F");
  m_tree->Branch("daughter1_PV_DCA", &b_daughter1_PV_DCA, "daughter1_PV_DCA/F");
  m_tree->Branch("daughter1_dEdx", &b_daughter1_dEdx, "daughter1_dEdx/F");

  m_tree->Branch("daughter2_mass", &b_daughter2_mass, "daughter2_mass/F");
  m_tree->Branch("daughter2_charge", &b_daughter2_charge, "daughter2_charge/I");
  m_tree->Branch("daughter2_pT", &b_daughter2_pT, "daughter2_pT/F");
  m_tree->Branch("daughter2_eta", &b_daughter2_eta, "daughter2_eta/F");
  m_tree->Branch("daughter2_phi_beamline", &b_daughter2_phi_beamline, "daughter2_phi_beamline/F");
  m_tree->Branch("daughter2_PV_DCA", &b_daughter2_PV_DCA, "daughter2_PV_DCA/F");
  m_tree->Branch("daughter2_dEdx", &b_daughter2_dEdx, "daughter2_dEdx/F");

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
  b_mother_PV_DCA = 0;
  b_daughter1_charge = 0;
  b_daughter1_pT = 0;
  b_daughter1_eta = 0;
  b_daughter1_phi_beamline = 0;
  b_daughter1_PV_DCA = 0;
  b_daughter1_dEdx = -1;
  b_daughter2_charge = 0;
  b_daughter2_pT = 0;
  b_daughter2_eta = 0;
  b_daughter2_phi_beamline = 0;
  b_daughter2_PV_DCA = 0;
  b_daughter2_dEdx = -1;
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
  std::vector<SimpleHelix> goodHelices;
  for (auto &iter : *m_trackmap)
  {
    SvtxTrack *track = iter.second;

    if (!hasSiliconClusters(track))
    {
      continue;
    }

    goodTracks.push_back(track);
    goodHelices.push_back(buildHelix(track));
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

      const SvtxVertex *primaryVertex = findMatchingVertex(crossing);
      if (!primaryVertex)
      {
        continue;
      }

      Vec3 sv, momentumA, momentumB;
      double daughterDCA = twoHelixDCA(goodHelices[i], goodHelices[j], sv, momentumA, momentumB);

      if (daughterDCA > m_track_to_track_DCA_cut)
      {
        continue;
      }

      Vec3 pv{primaryVertex->get_x(), primaryVertex->get_y(), primaryVertex->get_z()};
      Vec3 flight = sv - pv;
      double flightDistance = norm(flight);

      if (flightDistance < m_flight_distance_cut)
      {
        continue;
      }

      Vec3 motherMomentum = momentumA + momentumB;
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

      Vec3 point, direction;
      double pvDcaA = helixPointDCA(goodHelices[i], pv, point, direction);
      double pvDcaB = helixPointDCA(goodHelices[j], pv, point, direction);

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

      // Candidate is declared in the header, alongside SimpleHelix.
      std::vector<Candidate> candidates;
      for (bool daughter1IsA : hypotheses)
      {
        Candidate c;
        c.daughter1Track = daughter1IsA ? trackA : trackB;
        c.daughter2Track = daughter1IsA ? trackB : trackA;
        const Vec3 &momentum1 = daughter1IsA ? momentumA : momentumB;
        const Vec3 &momentum2 = daughter1IsA ? momentumB : momentumA;
        c.daughter1PvDca = daughter1IsA ? pvDcaA : pvDcaB;
        c.daughter2PvDca = daughter1IsA ? pvDcaB : pvDcaA;
        c.dedx1 = daughter1IsA ? dedxA : dedxB;
        c.dedx2 = daughter1IsA ? dedxB : dedxA;

        c.daughter1Vec = ROOT::Math::PxPyPzMVector(momentum1.x, momentum1.y, momentum1.z, daughter1_mass);
        c.daughter2Vec = ROOT::Math::PxPyPzMVector(momentum2.x, momentum2.y, momentum2.z, daughter2_mass);
        c.motherMass = (c.daughter1Vec + c.daughter2Vec).M();

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
        b_mother_PV_DCA = motherIP;

        b_daughter1_charge = c.daughter1Track->get_charge();
        b_daughter1_pT = c.daughter1Vec.Pt();
        b_daughter1_eta = c.daughter1Track->get_eta();
        b_daughter1_phi_beamline = c.daughter1Track->get_phi();
        b_daughter1_PV_DCA = c.daughter1PvDca;
        b_daughter1_dEdx = c.dedx1;

        b_daughter2_charge = c.daughter2Track->get_charge();
        b_daughter2_pT = c.daughter2Vec.Pt();
        b_daughter2_eta = c.daughter2Track->get_eta();
        b_daughter2_phi_beamline = c.daughter2Track->get_phi();
        b_daughter2_PV_DCA = c.daughter2PvDca;
        b_daughter2_dEdx = c.dedx2;

        b_track_to_track_DCA = daughterDCA;
        b_both_charge_states_passed = bothChargeStatesPassed;

        m_tree->Fill();
      }
    }
  }

  return Fun4AllReturnCodes::EVENT_OK;
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
  m_outfile->cd();
  m_tree->Write();
  m_outfile->Close();
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
  std::cout << "  mother IP cut: " << m_mother_PV_DCA_cut << " cm" << std::endl;
  std::cout << "  mother DIRA cut: " << m_dira_cut << std::endl;
  std::cout << "  daughter IP cut: " << m_daughter_PV_DCA_cut << " cm" << std::endl;
  std::cout << "  use dE/dx PID: " << (m_use_dEdx_pid ? "true" : "false") << std::endl;
  std::cout << "  field strength: " << m_Bz << " T" << std::endl;
}
