#ifndef MACRO_FUN4ALLG4SPHENIX_C
#define MACRO_FUN4ALLG4SPHENIX_C

#include <GlobalVariables.C>

#include "G4Setup_sPHENIX.C" // use default

#include <G4_Global.C>
#include <G4_Input.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_Reco.C>
#include <Trkr_Eval.C>

#include <ffamodules/FlagHandler.h>
#include <ffamodules/HeadReco.h>
#include <ffamodules/SyncReco.h>
#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>

#include <phool/PHRandomSeed.h>
#include <phool/recoConsts.h>

#include <tpctrackreco/TpcCrossingFinder.h>
#include <tpctrackreco/TpcPolyClusterTrkrClusterConverter.h>
#include <tpctrackreco/TpcPolyTrackSeedConverter.h>
#include <tpctrackreco/Tpc_AssembledTrackReco.h>
#include <tpctrackreco/Tpc_ModuleTrackReco.h>
#include <tpctrackreco/Tpc_PolyClusterizer.h>
#include <tpctrackreco/Tpc_PolyTrackReco.h>
#include <tpctrackreco/Tpc_PolyTrackVertexer.h>
#include <trackreco/DSTClusterPruning.h>

#include <Rtypes.h>
#include <TROOT.h>
#include <fstream>

#include <twotrackresonancereco/twoTrackResonanceReco.h>
#include <kfparticle_sphenix/KFParticle_sPHENIX.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libtwoTrackResonanceReco.so)
R__LOAD_LIBRARY(libkfparticle_sphenix.so)

int Fun4All_D0_sim(const int nEvents = 10
                 , const string &outdir = "./"
                 , const int processID = 0
                 , bool doPolytracking = false)
{
  std::stringstream nice_processID;
  nice_processID << std::setw(5) << std::setfill('0') << std::to_string(processID);

  int verbosity = 0;

  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(verbosity);

  PHRandomSeed::Verbosity(1);
  CDBInterface::instance()->Verbosity(1);

  recoConsts *rc = recoConsts::instance();
  //rc->set_IntFlag("RANDOMSEED", 12345678);

  Input::VERBOSITY = 0;

  Input::SIMPLE = true;
  Input::SIMPLE_VERBOSITY = verbosity;

  Input::BEAM_CONFIGURATION = Input::pp_COLLISION;

  Input::DZERO = true;
  Input::DZERO_VERBOSITY = verbosity;

  InputInit();

  if (Input::SIMPLE)
  {
    INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("pi-", 5);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus
                                                                            , PHG4SimpleEventGenerator::Gaus
                                                                            , PHG4SimpleEventGenerator::Gaus);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0.01, 0.01, 5.);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-1, 1);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(1, 10.);
  }

  InputRegister();

  rc->set_IntFlag("RUNNUMBER", 79516); // Use a real run number to get the TPC dead maps

  SyncReco *sync = new SyncReco();
  se->registerSubsystem(sync);

  HeadReco *head = new HeadReco();
  se->registerSubsystem(head);

  FlagHandler *flag = new FlagHandler();
  se->registerSubsystem(flag);

  Enable::QA = false;

  // Heavy-flavor simulation setup
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;

  Enable::PIPE = true;
  Enable::PIPE_ABSORBER = true;

  Enable::MVTX = true;
  Enable::INTT = true;
  Enable::TPC = true;
  Enable::TPC_ABSORBER = false;
  Enable::MICROMEGAS = true;

  Enable::MAGNET = true;
  Enable::MAGNET_ABSORBER = true;

  Enable::PLUGDOOR_ABSORBER = true;

  Enable::BLACKHOLE = true;

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");//CDB::global_tag);
  rc->set_uint64Flag("TIMESTAMP", CDB::timestamp);

  // Initialize the selected subsystems
  G4Init();
  G4Setup();

  //------------------
  // Detector Division
  //------------------
  Mvtx_Cells();
  Intt_Cells();
  TPC_Cells();
  Micromegas_Cells();

  //--------------
  // SVTX tracking
  //--------------
  TrackingInit();
  Mvtx_Clustering();
  Intt_Clustering();
  Micromegas_Clustering();

  if (doPolytracking)
  {
    Tracking_Reco_SiliconSeed_run2pp();  

    auto *converter = new TrackSeedTrackMapConverter("SiliconSeedToSvtxTrackMap");
    converter->setTrackSeedName("SiliconTrackSeedContainer");
    converter->setTrackMapName("SiliconSvtxTrackMap");
    converter->setClusterMapName("TRKR_CLUSTER");
    se->registerSubsystem(converter);  

    auto *finder_svx = new PHSimpleVertexFinder("SiliconVertexFinder");
    finder_svx->Verbosity(0);
    finder_svx->setDcaCut(0.1);
    finder_svx->setTrackPtCut(0.2);
    finder_svx->setBeamLineCut(1);
    finder_svx->setTrackQualityCut(500);
    finder_svx->setNmvtxRequired(3);
    finder_svx->setOutlierPairCut(0.1);
    finder_svx->setTrackMapName("SiliconSvtxTrackMap");
    finder_svx->setVertexMapName("SiliconSvtxVertexMap");
    se->registerSubsystem(finder_svx);  

    se->registerSubsystem(new Tpc_ModuleTrackReco());     // makes TPC_MODULETRACKS
    se->registerSubsystem(new Tpc_AssembledTrackReco());  // makes TPC_ASSEMBLEDTRACKS  

    auto *crossingFinder = new TpcCrossingFinder();
    crossingFinder->Verbosity(0);
    crossingFinder->setIsNonDistortedMC(true);
    crossingFinder->setInputNodeName("TPC_ASSEMBLEDTRACKS");
    crossingFinder->setOutputNodeName("TPC_CROSSING_DECISIONS");
    crossingFinder->setVertexMapNodeName("SiliconSvtxVertexMap");  // optional, configurable
    se->registerSubsystem(crossingFinder);  

    auto *cluster = new Tpc_PolyClusterizer();  // makes TPC_POLYCLUSTERS
    cluster->setIsNonDistortedMC(true);
    cluster->setUseSurveyGeometry(false);
    se->registerSubsystem(cluster);  

    se->registerSubsystem(new Tpc_PolyTrackReco());      // makes TPC_POLYTRACKS
    se->registerSubsystem(new Tpc_PolyTrackVertexer());  // makes TPC_POLYTRACKVERTICES  

    se->registerSubsystem(new TpcPolyTrackSeedConverter());           // converts TPC_POLYTRACKS to TpcTrackSeed
    se->registerSubsystem(new TpcPolyClusterTrkrClusterConverter());  // converts TPC_POLYCLUSTERS to TRKR_CLUSTER  

    Tracking_Reco_TrackMatching_run2pp();
    
    auto *clusterPruner = new DSTClusterPruning("DSTClusterPruning");
    clusterPruner->pruneAllSeeds();
    se->registerSubsystem(clusterPruner);

    Tracking_Reco_TrackFit_run2pp();
    Tracking_Reco_Vertex_run2pp();
  }
  else
  {
    ACTSGEOM::ActsGeomInit();

    auto *tpcclusterizer = new TpcClusterizer;
    tpcclusterizer->Verbosity(verbosity);
    tpcclusterizer->SetDeadChannelMapName("TPC_DEADCHANNELMAP");
    tpcclusterizer->set_do_hit_association(G4TPC::DO_HIT_ASSOCIATION);
    tpcclusterizer->set_min_err_squared(0.000001);
    se->registerSubsystem(tpcclusterizer);  

    auto *tpcclustercleaner = new TpcClusterCleaner;
    tpcclustercleaner->Verbosity(verbosity);
    tpcclustercleaner->set_rphi_error_low_cut(0.001);
    se->registerSubsystem(tpcclustercleaner);

    Tracking_Reco();
  }

  //-----------------
  // Global Vertexing
  //-----------------
  Global_Reco();

  // Heavy-flavor simulation setup
  build_truthreco_tables();

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  std::string output_dir = "./output/";  // Top dir of where the output nTuples will be written
  std::string standard_or_poly = doPolytracking ? "_polyseeding" : "_caseeding";
  std::string header = "output_twoTrackReco_simulation";
  std::string processing_folder = "inReconstruction/";
  std::string trailer = "_" + nice_processID.str() + ".root";

  std::string Dzero_reconstruction_name = "Dzero_reco" + standard_or_poly;  // Used for naming output folder, file and node
  std::string Dzero_output_file_name = header + Dzero_reconstruction_name + trailer;
  std::string Dzero_output_dir = output_dir + Dzero_reconstruction_name + "/";
  std::string Dzero_output_reco_dir = Dzero_output_dir + processing_folder;
  std::string Dzero_output_reco_file = Dzero_output_reco_dir + Dzero_output_file_name;

  std::string makeDirectory = "mkdir -p " + Dzero_output_reco_dir;
  system(makeDirectory.c_str());

  twoTrackResonanceReco* myDzeroReco = new twoTrackResonanceReco("DzeroReco");
  myDzeroReco->setDaughterPDGIDs(321, 211);
  myDzeroReco->setMotherMassRange(1.7, 2.0);
  myDzeroReco->setDaughterDCACut(0.05);
  myDzeroReco->setDIRACut(0.85);
  myDzeroReco->setOutputFileName(Dzero_output_reco_file.c_str());
  se->registerSubsystem(myDzeroReco);


  std::string KFParticle_Dzero_reconstruction_name = "Dzero_reco_KFParticle" + standard_or_poly;  // Used for naming output folder, file and node
  std::string KFParticle_Dzero_output_file_name = header + KFParticle_Dzero_reconstruction_name + trailer;
  std::string KFParticle_Dzero_output_dir = output_dir + KFParticle_Dzero_reconstruction_name + "/";
  std::string KFParticle_Dzero_output_reco_dir = KFParticle_Dzero_output_dir + processing_folder;
  std::string KFParticle_Dzero_output_reco_file = KFParticle_Dzero_output_reco_dir + KFParticle_Dzero_output_file_name;

  makeDirectory = "mkdir -p " + KFParticle_Dzero_output_reco_dir;
  system(makeDirectory.c_str());

  KFParticle_sPHENIX *myDzeroKFParticle = new KFParticle_sPHENIX(KFParticle_Dzero_reconstruction_name);
  myDzeroKFParticle->Verbosity(INT_MAX);
  myDzeroKFParticle->setDecayDescriptor("[D0 -> K^- pi^+]cc");
  myDzeroKFParticle->dontUseGlobalVertex(true);
  myDzeroKFParticle->requireTrackVertexBunchCrossingMatch(true);
  myDzeroKFParticle->constrainToPrimaryVertex();
  myDzeroKFParticle->usePID(false);
  myDzeroKFParticle->allowZeroMassTracks();
  myDzeroKFParticle->magFieldFile("FIELDMAP_TRACKING");
  myDzeroKFParticle->saveOutput(true);

  myDzeroKFParticle->setMinimumTrackPT(0.0);
  myDzeroKFParticle->setMaximumTrackchi2nDOF(100.);
  myDzeroKFParticle->setMinMVTXhits(1);
  myDzeroKFParticle->setMinINTThits(1);
  myDzeroKFParticle->setMinTPChits(0);

  myDzeroKFParticle->setMinimumMass(1.7);
  myDzeroKFParticle->setMaximumMass(2.0);
  myDzeroKFParticle->setMaximumDaughterDCA(0.05);
  myDzeroKFParticle->setMinDIRA(0.85);
  myDzeroKFParticle->setMotherPV_DCA(999);

  myDzeroKFParticle->setOutputName(KFParticle_Dzero_output_reco_file.c_str());
  se->registerSubsystem(myDzeroKFParticle);

  //======================
  // Write the DST
  //======================

  Enable::DSTOUT = false;
  Enable::DSTOUT_COMPRESS = true;
  DstOut::OutputDir = output_dir + "/DSTs/";
  std::string makeDSTDirectory = "mkdir -p " + DstOut::OutputDir;
  system(makeDSTDirectory.c_str());
  DstOut::OutputFile = "DST" + standard_or_poly + trailer;

  if (Enable::DSTOUT)
  {
    std::string FullOutFile = DstOut::OutputDir + "/" + DstOut::OutputFile;

    Fun4AllDstOutputManager *out = new Fun4AllDstOutputManager("DSTOUT", FullOutFile);
    out->StripNode("G4HIT_PIPE");
    out->StripNode("G4HIT_SVTXSUPPORT");
    out->StripNode("G4HIT_PIPE");
    out->StripNode("G4HIT_MVTX");
    out->StripNode("G4HIT_INTT");
    out->StripNode("G4HIT_TPC");
    out->StripNode("G4HIT_MICROMEGAS");
    out->StripNode("TRAINING_HITSET");
    out->StripNode("alignmentTransformationContainer");
    out->StripNode("alignmentTransformationContainerTransient");
    out->StripNode("ActsTrajectories");
    out->StripNode("SvtxAlignmentStateMap");

    if (Enable::DSTOUT_COMPRESS)
    {
        ShowerCompress();
        DstCompress(out);
    }
    se->registerOutputManager(out);
    std::cout << "Saving DST output to: " << FullOutFile << std::endl;
  }
  //-----------------
  // Event processing
  //-----------------

  if (nEvents < 0)
  {
      return 0;
  }

  se->run(nEvents);

  //-----
  // Exit
  //-----

  CDBInterface::instance()->Print(); // print used DB files
  se->End();

  std::ifstream outfile(Dzero_output_reco_file);
  if (outfile.good())
  {
    std::string moveOutput = "mv " + Dzero_output_reco_file + " " + Dzero_output_dir;
    system(moveOutput.c_str());
  }

  std::ifstream outfileKFParticle(KFParticle_Dzero_output_reco_file);
  if (outfileKFParticle.good())
  {
    std::string moveOutput = "mv " + KFParticle_Dzero_output_reco_file + " " + KFParticle_Dzero_output_dir;
    system(moveOutput.c_str());
  }

  std::cout << "All done" << std::endl;
  delete se;

  gSystem->Exit(0);
  return 0;
}
#endif
