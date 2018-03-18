void run_eve
(
  TString name = "run3000_s0", 
  TString pathToData = "/Users/ejungwoo/spiritroot/macros/data/",
  //TString parname = "ST.parameters.Commissioning_201604.par",
  TString parname = "ST.parameters.PhysicsRuns_201707.par",
  /*
   * - If dataList is "", deactivate single pad data,
   * - If dataList is set, activate single pad data (independent of reco file).
   *   XXX This may cause serious speed problem if meta data is not set.
   *   (depending on the system) if startEventID is not correct, pad may not match.
  */
  //TString dataList = "list_run3000.txt",
  TString dataList = "",
    Int_t runNo = 3000,
    Int_t startEventID = 0,
   Bool_t useMeta = false,
  TString supplePath = "/data/Q16264/rawdataSupplement"
)
{
  TString spiritroot = TString(gSystem -> Getenv("VMCWORKDIR"))+"/";
  if (pathToData.IsNull())
    pathToData = spiritroot+"macros/data/";

  //TString input     = "/mnt/spirit/analysis/user/leej/SpiRITROOT.trackChargeFix/macros/data/mc0_s0.reco.v1.04.root";

  //TString input     = "/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s200.reco.trackChargeFix.1505.c48e1f8.test.root";
  //TString input     = "/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s200.reco.develop.1520.9312cef.test.root";
  //TString input     = "/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s200.reco.develop.1522.ac6f82e.test.root";
  //TString input     = "/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s0.reco.develop.1525.d4b08aa.test.root";
  //TString input     = "/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2844_s0.reco.develop.1525.d4b08aa.test.root";
  TString input     = "/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s200.reco.develop.1525.d4b08aa.test.root";
  TString output    = "/mnt/spirit/analysis/user/leej/SpiRITROOT.trackChargeFix/macros/data/mc0_s0.eve.v1.04.root";

  //TString input     = "/mnt/spirit/analysis/user/leej/SpiRITROOT.trackChargeFix/macros/data/proton_helix.digi.root";
  //TString output    = "/mnt/spirit/analysis/user/leej/SpiRITROOT.trackChargeFix/macros/data/proton_helix.eve.root";

  //TString input     = "/mnt/spirit/analysis/user/leej/macros/reconstruction/data/run2900_s10.reco.trackChargeFix.1504.9a9e4cc.test.root";//pathToData + name + ".reco.root";
  //TString output    = "/mnt/spirit/analysis/user/leej/macros/reconstruction/data/run2900_s10.reco.trackChargeFix.1504.9a9e4cc.eve.root";//pathToData + name + ".eve.root";

  TString parameter = spiritroot + "parameters/"  + parname;
  TString geomety   = spiritroot + "geometry/geomSpiRIT.man.root";

  FairLogger *logger = FairLogger::GetLogger();
  logger -> SetLogToScreen(true);

  STEveManager *eve = new STEveManager();
  eve -> SetInputFile(input);         // Set input file (string)
  eve -> SetParInputFile(parameter);  // Set parameter file (string)
  eve -> SetOutputFile(output);       // Set output file (string)
  eve -> SetBackgroundColor(kWhite);  // Set background color (Color_t) 
  eve -> SetGeomFile(geomety);        // Set geometry file (string)
  eve -> SetVolumeTransparency(80);   // Set geometry transparency (integer, 0~100)
  eve -> SetViewerPoint(-0.7, 1.1);   // Set camera angle (theta, phi)

  STEveDrawTask *draw = new STEveDrawTask();
  draw -> SetRendering("mc",         false);
  draw -> SetRendering("digi",       false);
  draw -> SetRendering("hit",        false);
  draw -> SetRendering("hitbox",     false);
  draw -> SetRendering("helixhit",   true);
  draw -> SetRendering("helix",      false);
  draw -> SetRendering("cluster",    true);
  draw -> SetRendering("recotrack",  true);

  if (dataList.IsNull() == false) {
    if (useMeta)
      dataList = Form("%s/run_%04d/dataList.txt", supplePath.Data(), runNo);

    STDecoderTask *decoder = new STDecoderTask();
    decoder -> SetUseSeparatedData(true);
    decoder -> SetPersistence(false);
    decoder -> SetUseGainCalibration(false);
    decoder -> SetGGNoiseData("");
    decoder -> SetDataList(dataList);
    decoder -> SetEventID(startEventID);

    if (useMeta) {
      TString metaFile = Form("%s/run_%04d/metadataList.txt", supplePath.Data(), runNo);
      std::ifstream metalistFile(metaFile.Data());
      TString dataFileWithPath;
      for (Int_t iCobo = 0; iCobo < 12; iCobo++) {
        dataFileWithPath.ReadLine(metalistFile);
        dataFileWithPath = Form("%s/run_%04d/%s", supplePath.Data(), runNo, dataFileWithPath.Data());
        decoder -> SetMetaData(dataFileWithPath, iCobo);
      }
    }
    eve -> AddTask(decoder);
  }

  eve -> AddEveTask(draw);
  eve -> Init();
}
