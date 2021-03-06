void make_summary(int iRun = 1)
{
  TString tag = "develop.1534.2537e2e";

  auto GetCobo = [](Double_t z, Double_t x)
  {
    auto layer = (Int_t)(z/12);
    auto row = (Int_t)((x+432)/8);

    if (layer < 28) {
           if (row < 36) return 0;
      else if (row < 54) return 1;
      else if (row < 90) return 3;
      else               return 4;
    }
    else if (layer < 56) {
           if (row < 18) return 1;
      else if (row < 54) return 2;
      else if (row < 72) return 4;
      else               return 5;
    }
    if (layer < 84) {
           if (row < 36) return 6;
      else if (row < 54) return 7;
      else if (row < 90) return 9;
      else               return 10;
    }
    else {
           if (row < 18) return 7;
      else if (row < 54) return 8;
      else if (row < 72) return 10;
      else               return 11;
    }
  };

  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // TPC
  ////////////////////////////////////////////////////////////////////////////////////////////////////////

  auto input_tree = new TChain("cbmsim");
  Int_t runs[] = {2900,2901,2902,2903,2904,2905};
  Int_t numSplits[] = {44, 45, 43, 52, 42, 46};

  for (auto split = 0; split < numSplits[iRun]; ++split)
    input_tree -> Add(Form("/mnt/spirit/analysis/user/leej/macros/reconstruction/data/run%d_s%d.reco.",runs[iRun],split)+tag+".root");

  TClonesArray *clusterArray = nullptr;
  TClonesArray *trackArray = nullptr;
  input_tree -> SetBranchAddress("STHitCluster", &clusterArray);
  input_tree -> SetBranchAddress("STRecoTrack", &trackArray);

  auto fileName = TString(Form("/mnt/spirit/analysis/user/leej/macros/cluster_summary/data/summary_cluster_run%d.",runs[iRun]))+tag+".root";
  TFile *output_file = new TFile(fileName,"recreate");
  TTree *output_tree = new TTree("clusters","");

  Short_t lr, ndf, cobo;
  output_tree -> Branch("lr", &lr);
  output_tree -> Branch("cobo", &cobo);

  Float_t x, z, dy, dxz;
  output_tree -> Branch("x", &x);
  output_tree -> Branch("z", &z);
  output_tree -> Branch("dy", &dy);
  output_tree -> Branch("dxz", &dxz);

  Long64_t numEvents = input_tree -> GetEntries();
  for (Long64_t entry = 0; entry < numEvents; ++entry)
  {
    if (entry % 500 == 0)
      cout << "Event " << entry << " / " << numEvents << endl;

    input_tree -> GetEntry(entry);

    auto ntracks = trackArray ->  GetEntries();
    for (auto trackid = 0; trackid < ntracks; ++trackid)
    {
      auto track = (STRecoTrack *) trackArray -> At(trackid);

      auto dedxArray = track -> GetdEdxPointArray();
      if (dedxArray -> size() < 30)
        continue;
      auto z1 = (dedxArray->at(0)).fPosition.Z();
      auto z2 = (dedxArray->back()).fPosition.Z();
      if (z1 > 60 || z2 < 1284)
        continue;

      for (auto p : *dedxArray)
      {
        auto cluster = (STHitCluster *) clusterArray -> At(p.fClusterID);
        if (cluster -> GetNumHits() < 2)
          continue;

        x = cluster -> GetX();
        z = cluster -> GetZ();
        cobo = GetCobo(z, x);

        auto poca = cluster -> GetPOCA();
        auto pos  = cluster -> GetPosition();
        auto dist = pos - poca;

        if (cluster -> IsLayerCluster()) {
          lr = cluster -> GetLayer();
          dxz = dist.X();
        } else {
          lr = -(cluster -> GetRow());
          dxz = dist.Z();
        }
        dy = dist.Y();

        output_tree -> Fill();
      }
    }
  }

  output_file -> cd();
  output_tree -> Write();
  cout << "OUTPUT: " << fileName << endl;
  gApplication -> Terminate();
}
