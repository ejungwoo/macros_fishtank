void make_summary(int iRun = 0)
{
  //TString tag = "develop.1522.ac6f82e";
  TString tag = "develop.1526.c176211";

  Int_t runs[] = {2900,2901,2902,2903,2904,2905};
  Int_t numSplits[] = {44, 45, 43, 52, 42, 46};



  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // TPC
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto input_tree = new TChain("cbmsim");
  for (auto split = 0; split < numSplits[iRun]; ++split)
    input_tree -> Add(Form("/mnt/spirit/analysis/user/leej/macros/reconstruction/data/run%d_s%d.reco.",runs[iRun],split)+tag+".root");

  TClonesArray *vertexArray = nullptr;
  TClonesArray *recoTrackArray = nullptr;
  TClonesArray *vaddTrackArray = nullptr;

  input_tree -> SetBranchAddress("STVertex", &vertexArray);
  input_tree -> SetBranchAddress("STRecoTrack", &recoTrackArray);
  input_tree -> SetBranchAddress("VATracks", &vaddTrackArray);



  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // BEAM
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto input_beam = new TChain("beam");
       input_beam -> Add(Form("/mnt/spirit/analysis/changj/BeamAnalysis/macros/refined/Sn132_all/run%d.refined.root",runs[iRun]));
  bool sigma20;
  input_beam -> SetBranchAddress("sigma20", &sigma20);



  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // OUTPUT
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto fileName = TString(Form("/mnt/spirit/analysis/user/leej/macros/summary/data/summary%d.",runs[iRun]))+tag+".merged.root";
  TFile *output_file = new TFile(fileName,"recreate");
  TTree *output_tree = new TTree("data","");

  bool goodv, goodt, pion;
  Int_t ndf, ntracks;
  Double_t length, vx, vy, vz;

  Int_t pdg;
  Int_t pdg2;

  Double_t dist,  charge,  dedx,  p,  px,  py,  pz,  pocax,  pocay,  pocaz;
  Double_t dist2, charge2, dedx2, p2, px2, py2, pz2, pocax2, pocay2, pocaz2;

  output_tree -> Branch("goodv",   &goodv);
  output_tree -> Branch("goodt",   &goodt);
  output_tree -> Branch("pion",    &pion);
  output_tree -> Branch("vx",      &vx);
  output_tree -> Branch("vy",      &vy);
  output_tree -> Branch("vz",      &vz);
  output_tree -> Branch("ndf",     &ndf);
  output_tree -> Branch("length",  &length);
  output_tree -> Branch("ntracks", &ntracks);
  output_tree -> Branch("charge",  &charge); //reco
  output_tree -> Branch("pdg",     &pdg);
  output_tree -> Branch("dedx",    &dedx);
  output_tree -> Branch("p",       &p);
  output_tree -> Branch("px",      &px);
  output_tree -> Branch("py",      &py);
  output_tree -> Branch("pz",      &pz);
  output_tree -> Branch("dist",    &dist);
  output_tree -> Branch("pocax",   &pocax);
  output_tree -> Branch("pocay",   &pocay);
  output_tree -> Branch("pocaz",   &pocaz);
  output_tree -> Branch("charge2", &charge2); //vadd
  output_tree -> Branch("pdg2",    &pdg2);
  output_tree -> Branch("dedx2",   &dedx2);
  output_tree -> Branch("p2",      &p2);
  output_tree -> Branch("px2",     &px2);
  output_tree -> Branch("py2",     &py2);
  output_tree -> Branch("pz2",     &pz2);
  output_tree -> Branch("dist2",   &dist2);
  output_tree -> Branch("pocax2",  &pocax2);
  output_tree -> Branch("pocay2",  &pocay2);
  output_tree -> Branch("pocaz2",  &pocaz2);


  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  // MAKE SUMMARY
  ////////////////////////////////////////////////////////////////////////////////////////////////////////
  auto f1 = new TF1("f1","[0]*TMath::Landau(x,[1],[2],false)",0,100);

  Long64_t numEvents = input_tree -> GetEntries();
  for (Long64_t entry = 0; entry < numEvents; ++entry)
  {
    if (entry % 5000 == 0)
      cout << "Event " << entry << " / " << numEvents << endl;

    input_tree -> GetEntry(entry);
    input_beam -> GetEntry(entry);

    if (!sigma20)
      continue;

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

    goodv = false;
    if (vz>-13.8&&vz<-10.3 && vx>-15&&vx<15 && vy<-206.06&&vy>-246.06)
      goodv = true;

    ntracks = vaddTrackArray -> GetEntries();
    for (auto trackid2 = 0; trackid2 < ntracks; ++trackid2)
    {
      auto track2 = (STRecoTrack *) vaddTrackArray -> At(trackid2);
      if (track2 -> GetVertexID() != vertexid)
        continue;

      auto trackid = track2 -> GetParentID();
      auto track = (STRecoTrack *) recoTrackArray -> At(trackid);

      auto dedxArray = track -> GetdEdxPointArray();
      ndf = dedxArray -> size();
      if (ndf < 10)
        continue;

      pdg  = STPID::GetPDG(track -> GetPID());
      pdg2 = STPID::GetPDG(track2  -> GetPID());

      auto poca  = track -> GetPOCAVertex() - posV;
      auto poca2 = track2 -> GetPOCAVertex() - posV;

      dist  = poca.Mag();
      dist2 = poca2.Mag();

      pocax  = poca.X();
      pocay  = poca.Y();
      pocaz  = poca.Z();
      pocax2 = poca2.X();
      pocay2 = poca2.Y();
      pocaz2 = poca2.Z();

      length  = (dedxArray->at(ndf-1)).fLength - (dedxArray->at(0)).fLength;

      charge  = track -> GetCharge();
      charge2  = track2  -> GetCharge();

      auto mom  = track -> GetMomentum();
      auto mom2 = track2 -> GetMomentum();

      p  = mom.Mag();
      p2 = mom2.Mag();

      px = mom.X();
      py = mom.Y();
      pz = mom.Z();
      px2 = mom2.X();
      py2 = mom2.Y();
      pz2 = mom2.Z();

      dedx  = track  -> GetdEdxWithCut(0,0.7);
      dedx2 = track2 -> GetdEdxWithCut(0,0.7);

      goodt = false;
      if (ndf>30&&charge!=0&&dist<5)
        goodt = true;

      pion = false;
      if (((p2>310||p2<200)&&dedx2<-0.18*(p2-400)+20) || (p2>200&&p2<310&&dedx2<30))
        pion = true;

      output_tree -> Fill();
    }
  }

  output_file -> cd();
  output_tree -> Write();
  cout << "OUTPUT: " << fileName << endl;
  gApplication -> Terminate();
}

  /*
  TString tag = "trackChargeFix.1499.d50474b";  // before matching helix and reco track index
  TString tag = "trackChargeFix.1500.8a7d00d";  // before format version fix
  TString tag = "trackChargeFix.1501.7e5c17a";  // JungWoo - cut gg
  TString tag = "trackChargeFix.1502.9030344";  // PhysicsRun - use gg
  TString tag = "trackChargeFix.1503.4a82220";  // with parentID
  TString tag = "trackChargeFix.1508.bf6c33b";
  TString tag = "trackChargeFix.1504.9a9e4cc";  // cluster error on=
  */

