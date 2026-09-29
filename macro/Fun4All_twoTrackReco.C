#include <iostream>
#include <fstream>
#include <string>

#include <GlobalVariables.C>

#include <G4_ActsGeom.C>
#include <Trkr_RecoInit.C>

#include <cdbobjects/CDBTTree.h>

#include <twotrackresonancereco/twoTrackResonanceReco.h>

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

void Fun4All_twoTrackReco(
    const int nEvents = 1000,
    const std::string &inputList = "my.list",
    const int nSkip = 0)
{
  std::ifstream file(inputList.c_str());

  if (!file.is_open())
  {
      std::cerr << "Error: Could not open the file." << std::endl;
  }

  std::string inputDST;
  
  if (std::getline(file, inputDST))
  {
    std::cout << "First line: " << inputDST << std::endl;
  }
  else 
  {
    std::cout << "The file is empty." << std::endl;
    exit(1);
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
  rc->set_StringFlag("CDB_GLOBALTAG", "2026p003_v001"); // newcdbtag
  rc->set_uint64Flag("TIMESTAMP", runnumber);
  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");

  auto *se = Fun4AllServer::instance();
  se->Verbosity(1);

  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  TrackingInit();

  Fun4AllInputManager *tracks = new Fun4AllDstInputManager("TrackInputManager");
  tracks->AddListFile(inputList.c_str());
  se->registerInputManager(tracks);

  std::string output_dir = "./";  // Top dir of where the output nTuples will be written
  std::string header = "output_twoTrackReco_";
  std::string processing_folder = "inReconstruction/";
  std::string trailer = "_" + nice_runnumber.str() + "_" + nice_segment.str() + "_" + nice_skip.str() + ".root";

  std::string Kshort_reconstruction_name = "Kshort_reco";  // Used for naming output folder, file and node
  std::string Kshort_output_file_name = header + Kshort_reconstruction_name + trailer;
  std::string Kshort_output_dir = output_dir + Kshort_reconstruction_name + "/";
  std::string Kshort_output_reco_dir = Kshort_output_dir + processing_folder;
  std::string Kshort_output_reco_file = Kshort_output_reco_dir + Kshort_output_file_name;

  std::string makeDirectory = "mkdir -p " + Kshort_output_reco_dir;
  system(makeDirectory.c_str());

  twoTrackResonanceReco* myKshortReco = new twoTrackResonanceReco("KshortReco");
  myKshortReco->setMotherMassRange(0.4, 0.6);
  myKshortReco->setDaughterDCACut(0.1);
  myKshortReco->setFlightDistanceCut(0.8);
  myKshortReco->setOutputFileName(Kshort_output_reco_file.c_str());
  se->registerSubsystem(myKshortReco);

  se->skip(nSkip);
  se->run(nEvents);
  se->End();
  se->PrintTimer();

  std::ifstream outfile(Kshort_output_reco_file);
  if (outfile.good())
  {
    std::string moveOutput = "mv " + Kshort_output_reco_file + " " + Kshort_output_dir;
    system(moveOutput.c_str());
  }

  delete se;

  std::cout << "Finished" << std::endl;
  gSystem->Exit(0);
}