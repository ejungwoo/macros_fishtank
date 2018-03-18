void event()
{
  auto file = new TFile("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s200.reco.develop.1525.d4b08aa.test1.root");
  auto tree =(TTree *) file -> Get("cbmsim");
  TClonesArray *recos = nullptr;
  TClonesArray *vadds = nullptr;
  tree -> SetBranchAddress("STRecoTrack", &recos);
  tree -> SetBranchAddress("VATracks", &vadds);

  auto hist = new TH1D("hist","",100,0,100);
  cout << file -> GetName() << endl;
  for (auto event = 0; event < tree -> GetEntries(); ++event) {
    tree -> GetEntry(event);
    cout << recos -> GetEntries() << " " << vadds -> GetEntries() << " --- " << recos -> GetEntries() - vadds -> GetEntries() << endl;
    hist -> Fill(recos -> GetEntries() - vadds -> GetEntries());
  }
  hist -> Draw();
}
