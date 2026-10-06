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

  //Shared cuts
  float mass[2] = {0.4, 0.6};
  float track_to_track_DCA = 0.1;
  float daughter_PV_DCA = 0.05;
  float min_flight_distance = 0.05;
  float min_dira = 0.85;

  std::string output_dir = "./output/openingAngle/";  // Top dir of where the output nTuples will be written
  std::string header = "output_";
  std::string processing_folder = "inReconstruction/";
  std::string trailer = "_" + nice_runnumber.str() + "_" + nice_segment.str() + "_" + nice_skip.str() + ".root";

  std::string Dzero_reconstruction_name = "Dzero_reco_twoTrackReco";  // Used for naming output folder, file and node
  std::string Dzero_output_file_name = header + Dzero_reconstruction_name + trailer;
  std::string Dzero_output_dir = output_dir + Dzero_reconstruction_name + "/";
  std::string Dzero_output_reco_dir = Dzero_output_dir + processing_folder;
  std::string Dzero_output_reco_file = Dzero_output_reco_dir + Dzero_output_file_name;

  std::string makeDirectory = "mkdir -p " + Dzero_output_reco_dir;
  system(makeDirectory.c_str());

  twoTrackResonanceReco* myDzeroReco = new twoTrackResonanceReco("DzeroReco");
  myDzeroReco->setDaughterPDGIDs(321, 211);
  myDzeroReco->setMotherMassRange(1.7, 2.0);
  myDzeroReco->setDaughterDCACut(track_to_track_DCA);
  myDzeroReco->setDaughterIPCut(0.06);
  myDzeroReco->setFlightDistanceCut(0.06);
  myDzeroReco->setDIRACut(min_dira);
  myDzeroReco->setOutputFileName(Dzero_output_reco_file.c_str());
  se->registerSubsystem(myDzeroReco);

  std::string Kshort_reconstruction_name = "Kshort_reco_twoTrackReco";  // Used for naming output folder, file and node
  std::string Kshort_output_file_name = header + Kshort_reconstruction_name + trailer;
  std::string Kshort_output_dir = output_dir + Kshort_reconstruction_name + "/";
  std::string Kshort_output_reco_dir = Kshort_output_dir + processing_folder;
  std::string Kshort_output_reco_file = Kshort_output_reco_dir + Kshort_output_file_name;

  makeDirectory = "mkdir -p " + Kshort_output_reco_dir;
  system(makeDirectory.c_str());

  twoTrackResonanceReco* myKshortReco = new twoTrackResonanceReco("KshortReco");
  myKshortReco->setMotherMassRange(mass[0], mass[1]);
  myKshortReco->setDaughterDCACut(track_to_track_DCA);
  myKshortReco->setDaughterIPCut(daughter_PV_DCA);
  myKshortReco->setFlightDistanceCut(min_flight_distance);
  myKshortReco->setDIRACut(min_dira);
  myKshortReco->setOutputFileName(Kshort_output_reco_file.c_str());
  se->registerSubsystem(myKshortReco);

  std::string KFParticle_Kshort_reconstruction_name = "Kshort_reco_KFParticle";  // Used for naming output folder, file and node
  std::string KFParticle_Kshort_output_file_name = header + KFParticle_Kshort_reconstruction_name + trailer;
  std::string KFParticle_Kshort_output_dir = output_dir + KFParticle_Kshort_reconstruction_name + "/";
  std::string KFParticle_Kshort_output_reco_dir = KFParticle_Kshort_output_dir + processing_folder;
  std::string KFParticle_Kshort_output_reco_file = KFParticle_Kshort_output_reco_dir + KFParticle_Kshort_output_file_name;

  makeDirectory = "mkdir -p " + KFParticle_Kshort_output_reco_dir;
  system(makeDirectory.c_str());

  KFParticle_sPHENIX *myKshortKFParticle = new KFParticle_sPHENIX(KFParticle_Kshort_reconstruction_name);
  myKshortKFParticle->setDecayDescriptor("K_S0 -> pi^+ pi^-");
  myKshortKFParticle->dontUseGlobalVertex(true);
  myKshortKFParticle->requireTrackVertexBunchCrossingMatch(true);
  myKshortKFParticle->constrainToPrimaryVertex();
  myKshortKFParticle->usePID(false);
  myKshortKFParticle->allowZeroMassTracks();
  myKshortKFParticle->magFieldFile("FIELDMAP_TRACKING");
  myKshortKFParticle->saveOutput(true);

  myKshortKFParticle->setMinimumTrackPT(0.0);
  myKshortKFParticle->setMaximumTrackchi2nDOF(100.);
  myKshortKFParticle->setMinMVTXhits(1);
  myKshortKFParticle->setMinINTThits(1);
  myKshortKFParticle->setMinTPChits(0);
  myKshortKFParticle->setMinimumTrackPV_DCA(daughter_PV_DCA);

  myKshortKFParticle->setMinimumMass(mass[0]);
  myKshortKFParticle->setMaximumMass(mass[1]);
  myKshortKFParticle->setMaximumDaughterDCA(track_to_track_DCA);
  myKshortKFParticle->setDecayLengthRange(min_flight_distance, FLT_MAX);
  myKshortKFParticle->setMinDIRA(min_dira);
  myKshortKFParticle->setMotherPV_DCA(999);

  myKshortKFParticle->setOutputName(KFParticle_Kshort_output_reco_file.c_str());
  se->registerSubsystem(myKshortKFParticle);

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

  delete se;

  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}
