#include <cfloat>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>

#include <GlobalVariables.C>

#include <G4_ActsGeom.C>
#include <Trkr_RecoInit.C>

#include <cdbobjects/CDBTTree.h>

#include <twotrackresonancereco/twoTrackResonanceReco.h>

#include <kfparticle_sphenix/KFParticle_sPHENIX.h>
#include <nnscalemap/nnscalemap.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllInputManager.h>
#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libtwoTrackResonanceReco.so)
R__LOAD_LIBRARY(libkfparticle_sphenix.so)
R__LOAD_LIBRARY(libNNScaleMap.so)

void Fun4All_twoTrackReco(
    const int nEvents = 1500,
    const std::string &inputList = "jobLists/run79516_000.list",
    const int nSkip = 0)
{
  const std::string root_ext = ".root";
  bool isSingleFile = std::filesystem::path(inputList).extension() == ".root";

  std::string inputDST;

  if (isSingleFile)
  {
    inputDST = inputList;
    std::cout << "Single input DST: " << inputDST << std::endl;
  }
  else
  {
    std::ifstream file(inputList.c_str());

    if (!file.is_open())
    {
        std::cerr << "Error: Could not open the file." << std::endl;
    }

    if (std::getline(file, inputDST))
    {
      std::cout << "First line: " << inputDST << std::endl;
    }
    else
    {
      std::cout << "The file is empty." << std::endl;
      exit(1);
    }
  }

  std::pair<int, int> runseg = Fun4AllUtils::GetRunSegment(inputDST);
  int runnumber = runseg.first;
  std::stringstream nice_runnumber;
  nice_runnumber << std::setw(8) << std::setfill('0') << std::to_string(runnumber);

  int segment = runseg.second;
  std::stringstream nice_segment;
  nice_segment << std::setw(5) << std::setfill('0') << std::to_string(segment);

  std::stringstream nice_skip;
  nice_skip << std::setw(5) << std::setfill('0') << std::to_string(nSkip);

  Enable::CDB = true;
  recoConsts *rc = recoConsts::instance();
  rc->set_IntFlag("RUNNUMBER", runnumber);
  rc->set_StringFlag("CDB_GLOBALTAG", "newcdbtag");//2026p003_v001
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  auto *se = Fun4AllServer::instance();
  se->Verbosity(1);

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  TrackingInit();

  Fun4AllInputManager *tracks = new Fun4AllDstInputManager("TrackInputManager");
  if (isSingleFile)
  {
    tracks->AddFile(inputList.c_str());
  }
  else
  {
    tracks->AddListFile(inputList.c_str());
  }
  se->registerInputManager(tracks);

  bool enableLightFlavor = (10000 < segment && segment < 15000) ? true : false; //Dont need a huge amount of V0, and only look at ones that werent trained on
  std::string momReweightFile = "/sphenix/u/cdean/analysis/NNScaleMap/run" + std::to_string(runnumber) + "/kappa_lookup.csv";
  std::string momReweightTrackMap = "SvtxTrackMapMomentumScaleCorrected";

  auto* momScale = new NNScaleMap("NNScaleMap"); 
  momScale->setUseCDB(false);  // NNScaleMap defaults to CDB lookup; use the fixed local file instead
  momScale->setKappaLookupFile(momReweightFile.c_str());
  //momScale->setInputTrackMapName("SvtxTrackMap");
  momScale->setOutputTrackMapName(momReweightTrackMap.c_str());
  momScale->setScaleCovariance(true);
  se->registerSubsystem(momScale);

  //Shared light flavor cuts
  float KS0_mass[2] = {0.4, 0.6};
  float L0_mass[2] = {1.08, 1.15};
  float track_to_track_DCA = 0.1;
  float daughter_PV_DCA = 0.05;
  float daughter_PV_DCA_XY = 0.05;
  float min_flight_distance = 0.2;
  float min_dira = 0.85;

  //Shared D0 cuts
  float D0_mass[2] = {1.7, 2.0};
  float D0_track_to_track_DCA = 0.05;
  float D0_daughter_PV_DCA = 0.002;
  float D0_daughter_PV_DCA_XY = 0.002;
  float D0_min_flight_distance = 0.006;
  float D0_min_dira = 0.95;

  std::string output_dir = "./output/afterMLmomentum/";  // Top dir of where the output nTuples will be written
  std::string header = "output_";
  std::string processing_folder = "inReconstruction/";
  std::string trailer = "_" + nice_runnumber.str() + "_" + nice_segment.str() + "_" + nice_skip.str() + ".root";

 /*
  *  D0 twoTrackReco
  */
  std::string Dzero_reconstruction_name = "Dzero_reco_twoTrackReco";  // Used for naming output folder, file and node
  std::string Dzero_output_file_name = header + Dzero_reconstruction_name + trailer;
  std::string Dzero_output_dir = output_dir + Dzero_reconstruction_name + "/";
  std::string Dzero_output_reco_dir = Dzero_output_dir + processing_folder;
  std::string Dzero_output_reco_file = Dzero_output_reco_dir + Dzero_output_file_name;

  std::string makeDirectory = "mkdir -p " + Dzero_output_reco_dir;
  system(makeDirectory.c_str());

  twoTrackResonanceReco* myDzeroReco = new twoTrackResonanceReco("DzeroReco");
  myDzeroReco->setDaughterPDGIDs(321, 211);
  myDzeroReco->setMotherMassRange(D0_mass[0], D0_mass[1]);
  myDzeroReco->setDaughterDCACut(D0_track_to_track_DCA);
  myDzeroReco->setDaughterIPCut(D0_daughter_PV_DCA);
  myDzeroReco->setFlightDistanceCut(D0_min_flight_distance);
  myDzeroReco->setDIRACut(D0_min_dira);
  myDzeroReco->setOutputFileName(Dzero_output_reco_file.c_str());
  myDzeroReco->setTrackMapName(momReweightTrackMap.c_str());
  myDzeroReco->suseTpcMomentum(false);
  se->registerSubsystem(myDzeroReco);

 /*
  *  D0 KFParticle 
  */
  std::string KFParticle_Dzero_reconstruction_name = "Dzero_reco_KFParticle";  // Used for naming output folder, file and node
  std::string KFParticle_Dzero_output_file_name = header + KFParticle_Dzero_reconstruction_name + trailer;
  std::string KFParticle_Dzero_output_dir = output_dir + KFParticle_Dzero_reconstruction_name + "/";
  std::string KFParticle_Dzero_output_reco_dir = KFParticle_Dzero_output_dir + processing_folder;
  std::string KFParticle_Dzero_output_reco_file = KFParticle_Dzero_output_reco_dir + KFParticle_Dzero_output_file_name;

  makeDirectory = "mkdir -p " + KFParticle_Dzero_output_reco_dir;
  if (enableLightFlavor) system(makeDirectory.c_str());

  KFParticle_sPHENIX *myDzeroKFParticle = new KFParticle_sPHENIX(KFParticle_Dzero_reconstruction_name);
  myDzeroKFParticle->setDecayDescriptor("[D0 -> K^- pi^+]cc");
  myDzeroKFParticle->getDetectorInfo();
  myDzeroKFParticle->get_dEdx_info();
  myDzeroKFParticle->dontUseGlobalVertex();
  myDzeroKFParticle->requireTrackVertexBunchCrossingMatch();
  myDzeroKFParticle->constrainToPrimaryVertex();
  myDzeroKFParticle->usePID(false);
  myDzeroKFParticle->allowZeroMassTracks();
  myDzeroKFParticle->magFieldFile("FIELDMAP_TRACKING");
  myDzeroKFParticle->setTrackMapNodeName(momReweightTrackMap.c_str());
  myDzeroKFParticle->saveOutput();

  myDzeroKFParticle->setMinimumTrackPT(0.0);
  myDzeroKFParticle->setMaximumTrackchi2nDOF(300.);
  myDzeroKFParticle->setMinMVTXhits(1);
  myDzeroKFParticle->setMinINTThits(1);
  myDzeroKFParticle->setMinTPChits(20);
  myDzeroKFParticle->setMinimumTrackPV_DCA(D0_daughter_PV_DCA);
  myDzeroKFParticle->setMinimumTrackPV_DCA_XY(D0_daughter_PV_DCA_XY);

  myDzeroKFParticle->setMinimumMass(D0_mass[0]);
  myDzeroKFParticle->setMaximumMass(D0_mass[1]);
  myDzeroKFParticle->setMaximumDaughterDCA(D0_track_to_track_DCA);
  myDzeroKFParticle->setDecayLengthRange(D0_min_flight_distance, 2);
  myDzeroKFParticle->setMinDIRA(D0_min_dira);
  myDzeroKFParticle->setMotherPV_DCA(999);
  myDzeroKFParticle->setMaximumVertexchi2nDOF(FLT_MAX);

  myDzeroKFParticle->setOutputName(KFParticle_Dzero_output_reco_file.c_str());
  if (enableLightFlavor) se->registerSubsystem(myDzeroKFParticle);

 /*
  *  KS0 twoTrackReco
  */
  std::string Kshort_reconstruction_name = "Kshort_reco_twoTrackReco";  // Used for naming output folder, file and node
  std::string Kshort_output_file_name = header + Kshort_reconstruction_name + trailer;
  std::string Kshort_output_dir = output_dir + Kshort_reconstruction_name + "/";
  std::string Kshort_output_reco_dir = Kshort_output_dir + processing_folder;
  std::string Kshort_output_reco_file = Kshort_output_reco_dir + Kshort_output_file_name;

  makeDirectory = "mkdir -p " + Kshort_output_reco_dir;
  if (enableLightFlavor) system(makeDirectory.c_str());

  twoTrackResonanceReco* myKshortReco = new twoTrackResonanceReco("KshortReco");
  myKshortReco->setMotherMassRange(KS0_mass[0], KS0_mass[1]);
  myKshortReco->setDaughterDCACut(track_to_track_DCA);
  myKshortReco->setDaughterIPCut(daughter_PV_DCA);
  myKshortReco->setFlightDistanceCut(min_flight_distance);
  myKshortReco->setDIRACut(min_dira);
  myKshortReco->setOutputFileName(Kshort_output_reco_file.c_str());
  myKshortReco->setTrackMapName(momReweightTrackMap.c_str());
  myKshortReco->suseTpcMomentum(false);
  if (enableLightFlavor) se->registerSubsystem(myKshortReco);

 /*
  *  KS0 KFParticle 
  */
  std::string KFParticle_Kshort_reconstruction_name = "Kshort_reco_KFParticle";  // Used for naming output folder, file and node
  std::string KFParticle_Kshort_output_file_name = header + KFParticle_Kshort_reconstruction_name + trailer;
  std::string KFParticle_Kshort_output_dir = output_dir + KFParticle_Kshort_reconstruction_name + "/";
  std::string KFParticle_Kshort_output_reco_dir = KFParticle_Kshort_output_dir + processing_folder;
  std::string KFParticle_Kshort_output_reco_file = KFParticle_Kshort_output_reco_dir + KFParticle_Kshort_output_file_name;

  makeDirectory = "mkdir -p " + KFParticle_Kshort_output_reco_dir;
  if (enableLightFlavor) system(makeDirectory.c_str());

  KFParticle_sPHENIX *myKshortKFParticle = new KFParticle_sPHENIX(KFParticle_Kshort_reconstruction_name);
  myKshortKFParticle->setDecayDescriptor("K_S0 -> pi^+ pi^-");
  myKshortKFParticle->getDetectorInfo();
  myKshortKFParticle->get_dEdx_info();
  myKshortKFParticle->dontUseGlobalVertex();
  myKshortKFParticle->requireTrackVertexBunchCrossingMatch();
  myKshortKFParticle->constrainToPrimaryVertex();
  myKshortKFParticle->usePID(false);
  myKshortKFParticle->allowZeroMassTracks();
  myKshortKFParticle->magFieldFile("FIELDMAP_TRACKING");
  myKshortKFParticle->setTrackMapNodeName(momReweightTrackMap.c_str());
  myKshortKFParticle->saveOutput();

  myKshortKFParticle->setMinimumTrackPT(0.0);
  myKshortKFParticle->setMaximumTrackchi2nDOF(300.);
  myKshortKFParticle->setMinMVTXhits(0);
  myKshortKFParticle->setMinINTThits(0);
  myKshortKFParticle->setMinTPChits(0);
  myKshortKFParticle->setMinimumTrackPV_DCA(daughter_PV_DCA);
  myKshortKFParticle->setMinimumTrackPV_DCA_XY(daughter_PV_DCA_XY);

  myKshortKFParticle->setMinimumMass(KS0_mass[0]);
  myKshortKFParticle->setMaximumMass(KS0_mass[1]);
  myKshortKFParticle->setMaximumDaughterDCA(track_to_track_DCA);
  myKshortKFParticle->setDecayLengthRange(min_flight_distance, 40);
  myKshortKFParticle->setMinDIRA(min_dira);
  myKshortKFParticle->setMotherPV_DCA(999);
  myKshortKFParticle->setMaximumVertexchi2nDOF(FLT_MAX);

  myKshortKFParticle->setOutputName(KFParticle_Kshort_output_reco_file.c_str());
  if (enableLightFlavor) se->registerSubsystem(myKshortKFParticle);

 /*
  *  Lambda0 twoTrackReco 
  */
  std::string Lambda0_reconstruction_name = "Lambda0_reco_twoTrackReco";  // Used for naming output folder, file and node
  std::string Lambda0_output_file_name = header + Lambda0_reconstruction_name + trailer;
  std::string Lambda0_output_dir = output_dir + Lambda0_reconstruction_name + "/";
  std::string Lambda0_output_reco_dir = Lambda0_output_dir + processing_folder;
  std::string Lambda0_output_reco_file = Lambda0_output_reco_dir + Lambda0_output_file_name;

  makeDirectory = "mkdir -p " + Lambda0_output_reco_dir;
  //system(makeDirectory.c_str());

  twoTrackResonanceReco* myLambda0Reco = new twoTrackResonanceReco("Lambda0Reco");
  myLambda0Reco->setDaughterPDGIDs(2212, 211);
  myLambda0Reco->setMotherMassRange(L0_mass[0], L0_mass[1]);
  myLambda0Reco->setDaughterDCACut(track_to_track_DCA);
  myLambda0Reco->setDaughterIPCut(daughter_PV_DCA);
  myLambda0Reco->setFlightDistanceCut(min_flight_distance);
  myLambda0Reco->setDIRACut(min_dira);
  myLambda0Reco->setOutputFileName(Lambda0_output_reco_file.c_str());
  myLambda0Reco->requireSiliconClusters(false);
  //se->registerSubsystem(myLambda0Reco);

 /*
  *  Lambda0 KFParticle 
  */
  std::string KFParticle_Lambda0_reconstruction_name = "Lambda0_reco_KFParticle";  // Used for naming output folder, file and node
  std::string KFParticle_Lambda0_output_file_name = header + KFParticle_Lambda0_reconstruction_name + trailer;
  std::string KFParticle_Lambda0_output_dir = output_dir + KFParticle_Lambda0_reconstruction_name + "/";
  std::string KFParticle_Lambda0_output_reco_dir = KFParticle_Lambda0_output_dir + processing_folder;
  std::string KFParticle_Lambda0_output_reco_file = KFParticle_Lambda0_output_reco_dir + KFParticle_Lambda0_output_file_name;

  makeDirectory = "mkdir -p " + KFParticle_Lambda0_output_reco_dir;
  if (enableLightFlavor) system(makeDirectory.c_str());

  KFParticle_sPHENIX *myLambda0KFParticle = new KFParticle_sPHENIX(KFParticle_Lambda0_reconstruction_name);
  myLambda0KFParticle->setDecayDescriptor("[Lambda0 -> proton^+ pi^-]cc");
  myLambda0KFParticle->getDetectorInfo();
  myLambda0KFParticle->get_dEdx_info();
  myLambda0KFParticle->dontUseGlobalVertex();
  myLambda0KFParticle->requireTrackVertexBunchCrossingMatch();
  myLambda0KFParticle->constrainToPrimaryVertex();
  myLambda0KFParticle->usePID(false);
  myLambda0KFParticle->allowZeroMassTracks();
  myLambda0KFParticle->magFieldFile("FIELDMAP_TRACKING");
  myLambda0KFParticle->setTrackMapNodeName(momReweightTrackMap.c_str());
  myLambda0KFParticle->saveOutput();

  myLambda0KFParticle->setMinimumTrackPT(0.0);
  myLambda0KFParticle->setMaximumTrackchi2nDOF(300.);
  myLambda0KFParticle->setMinMVTXhits(0);
  myLambda0KFParticle->setMinINTThits(0);
  myLambda0KFParticle->setMinTPChits(0);
  myLambda0KFParticle->setMinimumTrackPV_DCA(daughter_PV_DCA);
  myLambda0KFParticle->setMinimumTrackPV_DCA_XY(daughter_PV_DCA_XY);

  myLambda0KFParticle->setMinimumMass(mass[0]);
  myLambda0KFParticle->setMaximumMass(mass[1]);
  myLambda0KFParticle->setMaximumDaughterDCA(track_to_track_DCA);
  myLambda0KFParticle->setDecayLengthRange(min_flight_distance, 40);
  myLambda0KFParticle->setMinDIRA(min_dira);
  myLambda0KFParticle->setMotherPV_DCA(999);
  myLambda0KFParticle->setMaximumVertexchi2nDOF(FLT_MAX);

  myLambda0KFParticle->setOutputName(KFParticle_Lambda0_output_reco_file.c_str());
  if (enableLightFlavor) se->registerSubsystem(myLambda0KFParticle);

  se->skip(nSkip);
  se->run(nEvents);
  se->End();
  se->PrintTimer();

  std::ifstream outfileDzero(Dzero_output_reco_file);
  if (outfileDzero.good())
  {
    std::string moveOutput = "mv " + Dzero_output_reco_file + " " + Dzero_output_dir;
    system(moveOutput.c_str());
  }

  std::ifstream outfileDzeroKFParticle(KFParticle_Dzero_output_reco_file);
  if (outfileDzeroKFParticle.good())
  {
    std::string moveOutput = "mv " + KFParticle_Dzero_output_reco_file + " " + KFParticle_Dzero_output_dir;
    system(moveOutput.c_str());
  }

  std::ifstream outfileKshort(Kshort_output_reco_file);
  if (outfileKshort.good())
  {
    std::string moveOutput = "mv " + Kshort_output_reco_file + " " + Kshort_output_dir;
    system(moveOutput.c_str());
  }

  std::ifstream outfileKFParticle(KFParticle_Kshort_output_reco_file);
  if (outfileKFParticle.good())
  {
    std::string moveOutput = "mv " + KFParticle_Kshort_output_reco_file + " " + KFParticle_Kshort_output_dir;
    system(moveOutput.c_str());
  }

  std::ifstream outfileLambda0(Lambda0_output_reco_file);
  if (outfileLambda0.good())
  {
    std::string moveOutput = "mv " + Lambda0_output_reco_file + " " + Lambda0_output_dir;
    system(moveOutput.c_str());
  }

  std::ifstream outfileL0KFParticle(KFParticle_Lambda0_output_reco_file);
  if (outfileL0KFParticle.good())
  {
    std::string moveOutput = "mv " + KFParticle_Lambda0_output_reco_file + " " + KFParticle_Lambda0_output_dir;
    system(moveOutput.c_str());
  }

  delete se;

  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}
