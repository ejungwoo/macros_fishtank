void whathappendtoz()
{
  auto file = new TFile("data/summary2.root");
  auto tree = (TTree *) file -> Get("data");

  auto file2 = new TFile("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s200.reco.trackChargeFix.1505.c48e1f8.test.root");
  auto tree2 = (TTree *) file2 -> Get("cbmsim");

  Int_t clusterID;
  Long64_t eventID;
  tree -> SetBranchAddress("eventID",&eventID);
  tree -> SetBranchAddress("clusterID",&clusterID);

  TClonesArray *hitArray = nullptr;
  TClonesArray *clusterArray = nullptr;
  tree2 -> SetBranchAddress("STHit",&hitArray);
  tree2 -> SetBranchAddress("STHitCluster",&clusterArray);

  for (int i = 0; i < 2; ++i)
  {
    tree -> GetEntry(i);
    tree2 -> GetEntry(eventID);

    auto clusternew = new STHitCluster();

    auto cluster = (STHitCluster *) clusterArray -> At(clusterID);
    cout << cluster -> GetRow() << " " << cluster -> GetLayer() << endl;;
    clusternew -> SetLayer(cluster -> GetLayer());
    clusternew -> SetRow(cluster -> GetRow());
    TMatrixD mat = cluster -> GetCovMatrix();
    cout << cluster -> GetX() << " "
      << cluster -> GetY() << " "
      << cluster -> GetZ() << " | "
      << cluster -> GetCharge() << " | "
      << mat(0,0) << " "
      << mat(1,1) << " "
      << mat(2,2) << endl;

    auto ids = cluster -> GetHitIDs();
    for (auto id : *ids) {
      auto hit = (STHit *) hitArray -> At(id);
      cout << hit -> GetX() << " "
        << hit -> GetY() << " "
        << hit -> GetZ() << " | "
        << hit -> GetCharge() << endl;
        clusternew -> AddHit(hit);

    /*
    mat = clusternew -> GetCovMatrix();
    cout << clusternew -> GetX() << " "
      << clusternew -> GetY() << " "
      << clusternew -> GetZ() << " | "
      << clusternew -> GetCharge() << " | "
      << mat(0,0) << " "
      << mat(1,1) << " "
      << mat(2,2) << endl;
      */
    }

    cout << endl;
  }
}
