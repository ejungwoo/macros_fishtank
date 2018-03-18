void genfit_reco()
{
  TString input = "/user/leej/macros/charge/data/run2841_s5991.FIX3.root";
  TString tag = "FIX3.genfit";
  TString fSupplePath = "";
  TString fParameterFile = "ST.parameters.JungWoo_201801.par";

  TString pathToData = Form("/user/leej/macros/charge/data/");
  gSystem -> Exec(TString("mkdir -p ") + pathToData);

  TString spiritroot = TString(gSystem -> Getenv("VMCWORKDIR"))+"/";
  TString par = spiritroot+"parameters/"+fParameterFile;
  TString geo = spiritroot+"geometry/geomSpiRIT.man.root";
  TString out = pathToData+"run2841_s5991."+tag+".root"; 
  TString log = pathToData+"run2841_s5991.log"; 

  FairLogger *logger = FairLogger::GetLogger();
  logger -> SetLogToScreen(true);

  FairParAsciiFileIo* parReader = new FairParAsciiFileIo();
  parReader -> open(par);

  FairRunAna* run = new FairRunAna();
  run -> SetGeomFile(geo);
  run -> SetInputFile(input);
  run -> SetOutputFile(out);
  run -> GetRuntimeDb() -> setSecondInput(parReader);

  auto helix = new STHelixTrackingTask();
  helix -> SetPersistence(true);
  helix -> SetClusterPersistence(true);

  auto genfitPID = new STGenfitPIDTask();
  genfitPID -> SetPersistence(true);
  genfitPID -> SetBDCFile("");  
  genfitPID -> SetConstantField();

  run -> AddTask(helix);
  run -> AddTask(genfitPID);

  run -> Init();
  run -> Run(0, 1);

  cout << "Log    : " << log << endl;
  cout << "Input  : " << input << endl;
  cout << "Output : " << out << endl;

  gApplication -> Terminate();
}
