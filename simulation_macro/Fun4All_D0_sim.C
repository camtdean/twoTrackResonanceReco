#ifndef MACRO_FUN4ALLG4SPHENIX_C
#define MACRO_FUN4ALLG4SPHENIX_C

#include <GlobalVariables.C>

#include "G4Setup_sPHENIX.C" // use default

#include <G4_Global.C>
#include <G4_Input.C>

#include <Trkr_RecoInit.C>
#include <Trkr_Clustering.C>
#include <Trkr_LaserClustering.C>
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

#include <Rtypes.h>
#include <TROOT.h>
#include <fstream>

#include <twotrackresonancereco/twoTrackResonanceReco.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libtwoTrackResonanceReco.so)

int Fun4All_D0_sim(const int nEvents = 10
                 , const string &outdir = "./"
                 , std::string processID = "0")
{
  Fun4AllServer *se = Fun4AllServer::instance();
  se->Verbosity(0);

  PHRandomSeed::Verbosity(1);
  CDBInterface::instance()->Verbosity(1);

  recoConsts *rc = recoConsts::instance();
  //rc->set_IntFlag("RANDOMSEED", std::stoi(processID));

  Input::VERBOSITY = 0;

  Input::SIMPLE = true;
  //Input::SIMPLE_NUMBER = 2; // if you need 2 of them
  Input::SIMPLE_VERBOSITY = 1;

  Input::BEAM_CONFIGURATION = Input::pp_COLLISION; // Input::AA_COLLISION (default), Input::pA_COLLISION, Input::pp_COLLISION

  // Input::GUN = true;
  // Input::GUN_NUMBER = 3; // if you need 3 of them
  // Input::GUN_VERBOSITY = 1;

  // D0 generator
  Input::DZERO = true;
  Input::DZERO_VERBOSITY = 0;

  InputInit();

  if (Input::SIMPLE)
  {
    INPUTGENERATOR::SimpleEventGenerator[0]->add_particles("pi-", 5);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_function(PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus, PHG4SimpleEventGenerator::Gaus);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_mean(0., 0., 0.);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_vertex_distribution_width(0.01, 0.01, 5.);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_eta_range(-1, 1);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_phi_range(-M_PI, M_PI);
    INPUTGENERATOR::SimpleEventGenerator[0]->set_pt_range(1, 10.);
  }

  // particle gun
  // if you run more than one of these Input::GUN_NUMBER > 1
  // add the settings for other with [1], next with [2]...
  //if (Input::GUN)
  //{
  //    INPUTGENERATOR::Gun[0]->AddParticle("pi-", 0, 1, 0);
  //    INPUTGENERATOR::Gun[0]->set_vtx(0, 0, 0);
  //}

  InputRegister();

  rc->set_IntFlag("RUNNUMBER", 1); //! This need to be set for G4_TrkrSimulation.C TPC()?

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
  Enable::MVTX_CELL = Enable::MVTX && true;
  Enable::MVTX_CLUSTER = Enable::MVTX_CELL && true;

  Enable::INTT = true;
  Enable::INTT_CELL = Enable::INTT && true;
  Enable::INTT_CLUSTER = Enable::INTT_CELL && true;

  Enable::TPC = true;
  Enable::TPC_ABSORBER = false;
  Enable::TPC_CELL = Enable::TPC && true;
  Enable::TPC_CLUSTER = Enable::TPC_CELL && true;

  Enable::MICROMEGAS = true;
  Enable::MICROMEGAS_CELL = Enable::MICROMEGAS && true;
  Enable::MICROMEGAS_CLUSTER = Enable::MICROMEGAS_CELL && true;

  Enable::TRACKING_TRACK = (Enable::MICROMEGAS_CLUSTER && Enable::TPC_CLUSTER && Enable::INTT_CLUSTER && Enable::MVTX_CLUSTER) && true;
  Enable::GLOBAL_RECO = (Enable::MBDFAKE || Enable::MBDRECO || Enable::TRACKING_TRACK) && true;

  Enable::MAGNET = true;
  Enable::MAGNET_ABSORBER = true;

  Enable::PLUGDOOR_ABSORBER = true;

  // new settings using Enable namespace in GlobalVariables.C
  Enable::BLACKHOLE = true;

  Enable::CDB = true;
  rc->set_StringFlag("CDB_GLOBALTAG", CDB::global_tag);
  rc->set_uint64Flag("TIMESTAMP", CDB::timestamp);

  // Initialize the selected subsystems
  G4Init();
  G4Setup();

  //------------------
  // Detector Division
  //------------------

  if (Enable::MVTX_CELL)
      Mvtx_Cells();
  if (Enable::INTT_CELL)
      Intt_Cells();
  if (Enable::TPC_CELL)
      TPC_Cells();
  if (Enable::MICROMEGAS_CELL)
      Micromegas_Cells();

  //--------------
  // SVTX tracking
  //--------------
  if (Enable::TRACKING_TRACK)
  {
      TrackingInit();
  }
  if (Enable::MVTX_CLUSTER)
      Mvtx_Clustering();
  if (Enable::INTT_CLUSTER)
      Intt_Clustering();
  if (Enable::TPC_CLUSTER)
  {
    if (G4TPC::ENABLE_DIRECT_LASER_HITS || G4TPC::ENABLE_CENTRAL_MEMBRANE_HITS)
    {
      TPC_LaserClustering();
    }
    else
    {
      TPC_Clustering();
    }
  }
  if (Enable::MICROMEGAS_CLUSTER)
      Micromegas_Clustering();

  if (Enable::TRACKING_TRACK)
  {
      Tracking_Reco();
  }

  // Heavy-flavor simulation setup
  auto vtxfinder = new PHSimpleVertexFinder;
  vtxfinder->Verbosity(0);
  vtxfinder->setDcaCut(0.5);
  vtxfinder->setTrackPtCut(-99999.);
  vtxfinder->setBeamLineCut(1);
  vtxfinder->setTrackQualityCut(1000000000);
  vtxfinder->setNmvtxRequired(3);
  vtxfinder->setOutlierPairCut(0.1);
  se->registerSubsystem(vtxfinder);

  //-----------------
  // Global Vertexing
  //-----------------

  if (Enable::GLOBAL_RECO && Enable::GLOBAL_FASTSIM)
  {
    std::cout << "You can only enable Enable::GLOBAL_RECO or Enable::GLOBAL_FASTSIM, not both" << std::endl;
    gSystem->Exit(1);
  }
  if (Enable::GLOBAL_RECO)
  {
    Global_Reco();
  }
  else if (Enable::GLOBAL_FASTSIM)
  {
    Global_FastSim();
  }

  // Heavy-flavor simulation setup
  build_truthreco_tables();

  //--------------
  // Set up Input Managers
  //--------------

  InputManagers();

  std::string output_dir = "./";  // Top dir of where the output nTuples will be written
  std::string header = "output_twoTrackReco_Dzero_simulation_";
  std::string processing_folder = "inReconstruction/";
  std::string trailer = "_" + processID + ".root";

  std::string Dzero_reconstruction_name = "Dzero_reco";  // Used for naming output folder, file and node
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


  //======================
  // Write the DST
  //======================

  Enable::DSTOUT = true;
  Enable::DSTOUT_COMPRESS = true;
  DstOut::OutputDir = output_dir;
  DstOut::OutputFile = "DST" + trailer;

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

  std::cout << "All done" << std::endl;
  delete se;
  if (Enable::PRODUCTION)
  {
      Production_MoveOutput();
  }

  gSystem->Exit(0);
  return 0;
}
#endif
