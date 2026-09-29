#ifndef TWOTRACKRESONANCERECO_H
#define TWOTRACKRESONANCERECO_H

#include <fun4all/SubsysReco.h>

#include <Math/Vector4D.h>

#include <cmath>
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

  void setDaughterDCACut(float cut) { m_track_to_track_DCA_cut = cut; }

  void setFlightDistanceCut(float cut) { m_flight_distance_cut = cut; }

  void setMotherIPCut(float cut) { m_mother_PV_DCA_cut = cut; }

  void usePID(bool use = true) { m_use_dEdx_pid = use; }

  void setdEdxBandWidth(float width) { m_dEdx_band_width = width; }

  void setFieldStrength(float bz) { m_Bz = bz; }

  void setTrackMapName(const std::string &name) { m_trackmap_node_name = name; }
  void setVertexMapName(const std::string &name) { m_vertexmap_node_name = name; }

  void setOutputFileName(const std::string &name) { m_outfile_name = name; }

 private:
  struct SimpleHelix
  {
    Vec3 position;
    Vec3 momentum;
    double charge{0};
    double centerX{0};
    double centerY{0};
    double radius{0};
    double rotationSense{0};
  };

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
  };

  int getNodes(PHCompositeNode *topNode);

  SimpleHelix buildHelix(SvtxTrack *track) const;

  bool hasSiliconClusters(SvtxTrack *track) const;

  const SvtxVertex *findMatchingVertex(short int crossing) const;

  double helixPointDCA(const SimpleHelix &helix, const Vec3 &target, Vec3 &point, Vec3 &direction) const;

  double twoHelixDCA(const SimpleHelix &helixA, const SimpleHelix &helixB,Vec3 &vertex, Vec3 &momentumA, Vec3 &momentumB) const;

  void init_dEdx_bands();

  double measured_dEdx(SvtxTrack *track) const;

  void initTree();

  void resetBranches();

  SvtxTrackMap *m_trackmap{nullptr};
  SvtxVertexMap *m_vertexmap{nullptr};
  TrkrClusterContainer *m_cluster_map{nullptr};
  PHG4TpcGeomContainer *m_geom_container{nullptr};
  ActsGeometry *m_acts_geometry{nullptr};

  std::string m_trackmap_node_name{"SvtxTrackMap"};
  std::string m_vertexmap_node_name{"SvtxVertexMap"};
  std::string m_outfile_name{"twoTrackResonanceReco.root"};

  int m_daughter1_PDGID{211};
  int m_daughter2_PDGID{-211};

  float m_min_mass{0};
  float m_max_mass{999};

  float m_track_to_track_DCA_cut{999};
  float m_flight_distance_cut{-999};
  float m_mother_PV_DCA_cut{999};

  bool m_use_dEdx_pid{false};
  bool m_can_get_dEdx{true};
  float m_dEdx_band_width{0.2};

  float m_Bz{1.4};

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
  float b_mother_PV_DCA{0};

  int b_daughter1_charge{0};
  float b_daughter1_pT{0};
  float b_daughter1_eta{0};
  float b_daughter1_phi_beamline{0};
  float b_daughter1_PV_DCA{0};
  float b_daughter1_dEdx{-1};

  int b_daughter2_charge{0};
  float b_daughter2_pT{0};
  float b_daughter2_eta{0};
  float b_daughter2_phi_beamline{0};
  float b_daughter2_PV_DCA{0};
  float b_daughter2_dEdx{-1};

  float b_track_to_track_DCA{0};

  bool b_both_charge_states_passed{false};
};

#endif  // TWOTRACKRESONANCERECO_H
