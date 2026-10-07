#ifndef TWOTRACKRESONANCERECO_H
#define TWOTRACKRESONANCERECO_H

#include <fun4all/SubsysReco.h>

#include <Math/Vector4D.h>

#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <string>

struct Vec3
{
  double x{0}, y{0}, z{0};
};

inline Vec3 operator+(const Vec3 &a, const Vec3 &b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3 &a, const Vec3 &b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator*(double s, const Vec3 &a) { return {s * a.x, s * a.y, s * a.z}; }
inline double dot(const Vec3 &a, const Vec3 &b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double norm(const Vec3 &a) { return std::sqrt(dot(a, a)); }

class PHCompositeNode;
class SvtxTrack;
class SvtxTrackMap;
class SvtxVertex;
class SvtxVertexMap;
class TrkrClusterContainer;
class PHG4TpcGeomContainer;
class ActsGeometry;
class TFile;
class TTree;
class TF1;

class twoTrackResonanceReco : public SubsysReco
{
 public:
  twoTrackResonanceReco(const std::string &name = "twoTrackResonanceReco");

  ~twoTrackResonanceReco() override;

  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

  void setDaughterPDGIDs(int pdgID1, int pdgID2)
  {
    m_daughter1_PDGID = pdgID1;
    m_daughter2_PDGID = pdgID2;
  }

  void setMotherMassRange(float minMass, float maxMass)
  {
    m_min_mass = minMass;
    m_max_mass = maxMass;
  }

  void setMaxDaughterChi2perNDF(float cut) { m_max_track_chi2_per_ndf = cut; }

  void setMinDaughterPT(float cut) { m_min_daughter_pt = cut; }

  void setDaughterDCACut(float cut) { m_track_to_track_DCA_cut = cut; }

  void setFlightDistanceCut(float cut) { m_flight_distance_cut = cut; }

  // Upper limit on the 3D flight distance (cm). Useful to reject pairs that come from a different vertex than the PV.
  void setMaxFlightDistanceCut(float cut) { m_max_flight_distance_cut = cut; }

  // true (default): the primary vertex is the same-crossing vertex closest in z to the secondary vertex.
  // false: the same-crossing vertex with the best chi2/ndf (previous behaviour).
  void matchVertexToSV(bool use = true) { m_match_vertex_to_sv = use; }

  void setMotherIPCut(float cut) { m_mother_PV_DCA_cut = cut; }

  void setDaughterIPCut(float cut) { m_daughter_PV_DCA_cut = cut; }

  void setDIRACut(float cut) { m_dira_cut = cut; }

  void usePID(bool use = true) { m_use_dEdx_pid = use; }

  void setdEdxBandWidth(float width) { m_dEdx_band_width = width; }

  void setTrackMapName(const std::string &name) { m_trackmap_node_name = name; }
  void setVertexMapName(const std::string &name) { m_vertexmap_node_name = name; }

  void setOutputFileName(const std::string &name) { m_outfile_name = name; }

  // Direction of each daughter at the secondary vertex: true = silicon seed (MVTX+INTT) circle/slope, false = ACTS fit
  void useSiliconDirection(bool use = true) { m_use_silicon_direction = use; }

  // Magnitude of each daughter's momentum: true = TPC seed pT, false = ACTS fit
  void useTpcMomentum(bool use = true) { m_use_tpc_momentum = use; }

  void requireSiliconClusters(bool use = true) { m_require_silicon = use; }

 private:
  struct Candidate
  {
    bool passed{false};
    SvtxTrack *daughter1Track{nullptr};
    SvtxTrack *daughter2Track{nullptr};
    double daughter1PvDca{0};
    double daughter2PvDca{0};
    double dedx1{-1};
    double dedx2{-1};
    ROOT::Math::PxPyPzMVector daughter1Vec;
    ROOT::Math::PxPyPzMVector daughter2Vec;
    double motherMass{0};
    double openingAngle{0};      // angle between the daughters' momenta used for the mass
    double openingAngleActs{0};  // same angle using the ACTS momenta, for comparison
    double p1{0}, p2{0};         // |p| used for the mass
    double p1Acts{0}, p2Acts{0}; // |p| from the ACTS fit, for comparison
    int source1{0}, source2{0};  // bit 0: silicon-seed direction used, bit 1: TPC-seed pT used
  };

  int getNodes(PHCompositeNode *topNode);

  bool hasSiliconClusters(SvtxTrack *track) const;

  // Vertex with the given beam crossing. If referenceZ is finite the one closest in z to referenceZ is returned,
  // otherwise the one with the best chi2/ndf. nMatching, if given, is set to how many vertices have that crossing.
  const SvtxVertex *findMatchingVertex(short int crossing, double referenceZ = std::numeric_limits<double>::quiet_NaN(), int *nMatching = nullptr) const;

  bool propagateToPoint(SvtxTrack *track, const Vec3 &target, Vec3 &pos, Vec3 &mom) const;

  bool buildSV(SvtxTrack *trackA, SvtxTrack *trackB, Vec3 &vertex, Vec3 &momentumA, Vec3 &momentumB, double &dca) const;

  bool trackToVertexDCA(SvtxTrack *track, const Vec3 &vertex, double &dca) const;

  void init_dEdx_bands();

  double measured_dEdx(SvtxTrack *track) const;

  void initTree();

  void resetBranches();

  std::array<double,3> unit(const Vec3& v);

  // Unit vector of the track's direction at the secondary vertex, from the silicon seed (MVTX+INTT hits)
  bool siliconDirectionAtSV(const SvtxTrack *track, const Vec3 &sv, Vec3 &dir) const;

  // Momentum used for the mass: |p| from the TPC seed pT, direction from the silicon seed at the SV.
  // Falls back to the ACTS momentum for whichever piece is unavailable. Returns a bitmask of what was used
  // (bit 0 = silicon direction, bit 1 = TPC pT).
  int refineMomentum(const SvtxTrack *track, const Vec3 &sv, const Vec3 &actsMomentum, Vec3 &out) const;

  SvtxTrackMap *m_trackmap{nullptr};
  SvtxVertexMap *m_vertexmap{nullptr};
  TrkrClusterContainer *m_cluster_map{nullptr};
  PHG4TpcGeomContainer *m_geom_container{nullptr};
  ActsGeometry *m_acts_geometry{nullptr};

  std::string m_trackmap_node_name{"SvtxTrackMap"};
  std::string m_vertexmap_node_name{"SvtxVertexMap"};
  std::string m_outfile_name{"twoTrackResonanceReco.root"};

  int candidateCounter{0};

  int m_daughter1_PDGID{211};
  int m_daughter2_PDGID{-211};

  float m_min_mass{0};
  float m_max_mass{2};

  float m_max_track_chi2_per_ndf{100};
  float m_min_daughter_pt{0.0};

  float m_track_to_track_DCA_cut{999};
  float m_flight_distance_cut{-999};
  float m_max_flight_distance_cut{999};
  bool m_match_vertex_to_sv{true};
  float m_mother_PV_DCA_cut{999};
  float m_dira_cut{-1.1};
  float m_daughter_PV_DCA_cut{-0.1};

  bool m_use_silicon_direction{true};
  bool m_use_tpc_momentum{true};
  bool m_require_silicon{true};

  bool m_use_dEdx_pid{false};
  bool m_can_get_dEdx{true};
  float m_dEdx_band_width{0.2};

  TF1 *f_pion_plus{nullptr};
  TF1 *f_kaon_plus{nullptr};
  TF1 *f_proton_plus{nullptr};
  TF1 *f_pion_minus{nullptr};
  TF1 *f_kaon_minus{nullptr};
  TF1 *f_proton_minus{nullptr};
  std::map<int, TF1 *> m_dEdx_bandMap;

  TFile *m_outfile{nullptr};
  TTree *m_tree{nullptr};

  short b_crossing{0};

  float b_daughter1_mass{0};
  float b_daughter2_mass{0};

  float b_PV_x{0}, b_PV_y{0}, b_PV_z{0};
  float b_SV_x{0}, b_SV_y{0}, b_SV_z{0};

  float b_mother_mass{0};
  float b_mother_pT{0};
  float b_mother_eta{0};
  float b_mother_phi{0};
  float b_mother_DIRA{0};
  float b_mother_flight_distance{0};
  int b_n_matching_vertices{0};
  float b_mother_PV_DCA{0};
  float b_mother_DIRA_acts{0};
  float b_mother_PV_DCA_acts{0};
  float b_mother_opening_angle{0};
  float b_mother_opening_angle_acts{0};

  int b_daughter1_charge{0};
  float b_daughter1_pT{0};
  float b_daughter1_eta{0};
  float b_daughter1_phi{0};
  float b_daughter1_PV_DCA{0};
  float b_daughter1_dEdx{-1};
  float b_daughter1_quality{-1};
  float b_daughter1_p{0};
  float b_daughter1_p_acts{0};
  int b_daughter1_momentum_source{0};

  int b_daughter2_charge{0};
  float b_daughter2_pT{0};
  float b_daughter2_eta{0};
  float b_daughter2_phi{0};
  float b_daughter2_PV_DCA{0};
  float b_daughter2_dEdx{-1};
  float b_daughter2_quality{-1};
  float b_daughter2_p{0};
  float b_daughter2_p_acts{0};
  int b_daughter2_momentum_source{0};

  float b_track_to_track_DCA{0};

  bool b_both_charge_states_passed{false};
};

#endif  // TWOTRACKRESONANCERECO_H
