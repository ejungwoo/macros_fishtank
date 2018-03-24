void vertex_summary()
{
  //TString tag = "develop.1522.ac6f82e";
  //TString tag = "develop.1526.c176211";
  TString tag = "develop.1534.2537e2e";

  Int_t runs[] = {2900,2901,2902,2903,2904,2905};
  Int_t numSplits[] = {44, 45, 43, 52, 42, 46};



  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // TPC
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto input_tree = new TChain("cbmsim");
  for (auto iRun = 0; iRun < 5; ++iRun)
    for (auto split = 0; split < numSplits[iRun]; ++split)
      input_tree -> Add(Form("/mnt/spirit/analysis/user/leej/macros/reconstruction/data/run%d_s%d.reco.",runs[iRun],split)+tag+".root");

  TClonesArray *vertexArray = nullptr;
  input_tree -> SetBranchAddress("STVertex", &vertexArray);



  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // OUTPUT
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto fileName = TString("/mnt/spirit/analysis/user/leej/macros/summary/data/summary_vertex.")+tag+".root";
  TFile *output_file = new TFile(fileName,"recreate");
  TTree *output_tree = new TTree("vertex","");

  Double_t vx, vy, vz;
  output_tree -> Branch("vx",      &vx);
  output_tree -> Branch("vy",      &vy);
  output_tree -> Branch("vz",      &vz);



  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // MAKE SUMMARY
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  Long64_t numEvents = input_tree -> GetEntries();
  for (Long64_t entry = 0; entry < numEvents; ++entry)
  {
    if (entry % 5000 == 0)
      cout << "Event " << entry << " / " << numEvents << endl;

    input_tree -> GetEntry(entry);
    if (vertexArray -> GetEntries() != 1)
      continue;

    auto vertexid = -1;
    auto numVertex = vertexArray -> GetEntries();
    STVertex *vertex = nullptr;
    for (auto iv = 0; iv < numVertex; ++iv) {
      vertex = (STVertex *) vertexArray -> At(iv);
      if (vertex -> IsCollisionVertex()) {
        vertexid = iv;
        break;
      }
    }
    if (vertexid == -1)
      continue;

    auto posV = vertex -> GetPos();
    vx = posV.X();
    vy = posV.Y();
    vz = posV.Z();
  }

  output_file -> cd();
  output_tree -> Write();
  cout << "OUTPUT: " << fileName << endl;
  gApplication -> Terminate();
}
