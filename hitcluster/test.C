void test()
{
  auto tree = new TChain("cbmsim");
  tree -> Add("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2844_s200.reco.develop.1525.d4b08aa.test.root");

  TClonesArray *clusterArray = nullptr;
  TClonesArray *trackArray = nullptr;
  tree -> SetBranchAddress("STHitCluster", &clusterArray);
  tree -> SetBranchAddress("STHelixTrack", &trackArray);

  auto entries = tree -> GetEntries();
  for (auto event = 0; event < 1; ++event)
  {
    tree -> GetEntry(event);
    auto numTracks = trackArray -> GetEntries();
    for (auto iHelix = 0; iHelix < 1; ++iHelix)
    {
      auto helix = (STHelixTrack *) trackArray -> At(iHelix);
      auto ids = helix -> GetClusterIDArray();

      auto maxx = 0.;
      auto maxy = 0.;
      auto maxz = 0.;
      for (auto id : *ids) {
        auto cluster = (STHitCluster *) clusterArray -> At(id);

        auto cov = cluster -> GetCovMatrix();
        auto cx = cov(0,0); if (cx > maxx) maxx = cx;
        auto cy = cov(1,1); if (cy > maxy) maxy = cy;
        auto cz = cov(2,2); if (cz > maxz) maxz = cz;

        auto numHits = cluster -> GetNumHits();
        cout << setw(5) << numHits << setw(15) << cx << setw(15) << cy << setw(15) << cz << endl;
      }
      for (auto id : *ids) {
        auto cluster = (STHitCluster *) clusterArray -> At(id);

        auto cov = cluster -> GetCovMatrix();
        auto cx = (cov(0,0)/maxx+1)/2;
        auto cy = (cov(1,1)/maxy+1)/2;
        auto cz = (cov(2,2)/maxz+1)/2;

        auto numHits = cluster -> GetNumHits();
        cout << setw(5) << numHits << setw(15) << cx << setw(15) << cy << setw(15) << cz << endl;
      }
    }
  }
}
