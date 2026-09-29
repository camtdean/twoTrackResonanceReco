#include <iostream>
#include <fstream>
#include <string>

#include <GlobalVariables.C>

//#include <G4_ActsGeom.C>
//#include <G4_Magnet.C>
//#include <QA.C>
//#include <Trkr_Clustering.C>
//#include <Trkr_Reco.C>
//#include <Trkr_RecoInit.C>
//#include <Trkr_TpcReadoutInit.C>

//#include <globalvertex/GlobalVertexReco.h>

#include <cdbobjects/CDBTTree.h>

//#include <tpccalib/PHTpcResiduals.h>

//#include <trackingqa/SiliconSeedsQA.h>
//#include <trackingqa/TpcSeedsQA.h>
//#include <trackingqa/TpcSiliconQA.h>

//#include <trackingdiagnostics/TrackResiduals.h>
//#include <trackingdiagnostics/TrkrNtuplizer.h>

//#include <kfparticle_sphenix/KFParticle_sPHENIX.h>
#include <twotrackresonancereco/twoTrackResonanceReco.h>

#include <ffamodules/CDBInterface.h>

#include <fun4all/Fun4AllDstInputManager.h>
//#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllInputManager.h>
//#include <fun4all/Fun4AllOutputManager.h>
//#include <fun4all/Fun4AllRunNodeInputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <fun4all/Fun4AllUtils.h>

#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)
R__LOAD_LIBRARY(libffamodules.so)
R__LOAD_LIBRARY(libphool.so)
R__LOAD_LIBRARY(libcdbobjects.so)
R__LOAD_LIBRARY(libtwoTrackResonanceReco.so)

void Fun4All_twoTrackReco(
    const int nEvents = 20000,
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

  int rounded_up = 100 * (std::ceil((float) runnumber / 100));
  std::stringstream nice_rounded_up;
  nice_rounded_up << std::setw(8) << std::setfill('0') << std::to_string(rounded_up);

  int rounded_down = 100 * (std::floor((float) runnumber / 100));
  std::stringstream nice_rounded_down;
  nice_rounded_down << std::setw(8) << std::setfill('0') << std::to_string(rounded_down);

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
//  std::string geofile = CDBInterface::instance()->getUrl("Tracking_Geometry");
/*
  TpcReadoutInit(runnumber);
  std::cout << " run: " << runnumber
            << " samples: " << TRACKING::reco_tpc_maxtime_sample
            << " pre: " << TRACKING::reco_tpc_time_presample
            << " vdrift: " << G4TPC::tpc_drift_velocity_reco
            << std::endl;
*/

  // distortion calibration mode
  /*
   * set to true to enable residuals in the TPC with
   * TPC clusters not participating to the ACTS track fit
   */
/*
  G4TRACKING::SC_CALIBMODE = false;
  Enable::MVTX_APPLYMISALIGNMENT = true;
  ACTSGEOM::mvtx_applymisalignment = Enable::MVTX_APPLYMISALIGNMENT;
  TRACKING::streaming_mode = true;
*/
  auto *se = Fun4AllServer::instance();
  se->Verbosity(1);
/*
  Fun4AllRunNodeInputManager *ingeo = new Fun4AllRunNodeInputManager("GeoIn");
  ingeo->AddFile(geofile);
  se->registerInputManager(ingeo);

  TrackingInit();
*/
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

