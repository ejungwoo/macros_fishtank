void track22()
{
  auto file = new TFile("/user/leej/macros/charge/data/run2841_s5991.FIX3.genfit.root");
  auto tree = (TTree *) file -> Get("cbmsim"); 

  TClonesArray *clusterArray = nullptr;
  TClonesArray *helixArray = nullptr;
  TClonesArray *trackArray = nullptr;
  TClonesArray *vertexArray = nullptr;

  tree -> SetBranchAddress("STHitCluster", &clusterArray);
  tree -> SetBranchAddress("STHelixTrack", &helixArray);
  tree -> SetBranchAddress("STRecoTrack", &trackArray);
  tree -> SetBranchAddress("STVertex", &vertexArray);
  tree -> GetEntry(0);

  Int_t numHelices = helixArray -> GetEntries();
  for (auto iHelix = 0; iHelix < numHelices; ++iHelix)
  {
    auto helix = (STHelixTrack *) helixArray -> At(iHelix);
    auto trackID = helix -> GetGenfitID();
    if (trackID != 22)
      continue;
    auto track = (STRecoTrack *) trackArray -> At(trackID);
    helix -> Print();
    track -> Print();

    cout << track -> GetEffCurvature1() << " " << track -> GetEffCurvature2() << " " << track -> GetEffCurvature3() << endl;

    auto dedxArray = track -> GetdEdxPointArray();
    for (auto dedx : *dedxArray) {
      cout << dedx.fClusterID << " " << dedx.fPosition.Z() << " " << dedx.fLength << endl;
    }
  }
}
