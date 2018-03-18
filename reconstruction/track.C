void track()
{
  auto file = new TFile("data/run2900_s0.reco.trackChargeFix.1503.4a82220.root");
  //auto file = new TFile("data/run2900_s0.reco.trackChargeFix.1508.bf6c33b.root");
  auto tree = (TTree *) file -> Get("cbmsim");

  TClonesArray *trackArray = nullptr;
  tree -> SetBranchAddress("STRecoTrack", &trackArray);

  tree -> GetEntry(1);

  auto track = (STRecoTrack *) trackArray -> At(0);
  auto dedxArray = track -> GetdEdxPointArray();
  for (auto dedx : *dedxArray)
    cout << dedx.fClusterID << endl;
    //dedx.Print();
}
