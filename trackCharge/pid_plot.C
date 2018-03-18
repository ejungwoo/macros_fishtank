void pid_plot(
  bool drawPIDLine = true,
  STPID::PID selectPID = STPID::kNon, // kPion, kProton, kDeuteron, kTriton, k3He, k4He, KNon
  Int_t nbins = 200,
  Int_t p1 = -500,
  Int_t p2 = 2500,
  Int_t dedx1 = 0,
  Int_t dedx2 = 800
)
{

  gStyle -> SetPadBottomMargin(0.10);
  gStyle -> SetPadLeftMargin(0.11);
  gStyle -> SetPadRightMargin(0.12);
  gStyle -> SetTitleFontSize(0.06);

  auto tree = new TChain("cbmsim");
  //  tree -> Add("./data/run2894_s0.reco_minus1to90.v1.04.root");
  tree -> Add("/user/leej/macros/charge/data/run2841_s0.FIX3.root");


  TClonesArray *trackArray = nullptr;
  TClonesArray *vertexArray = nullptr;

  tree -> SetBranchAddress("STRecoTrack", &trackArray);
  tree -> SetBranchAddress("STVertex", &vertexArray);

  auto hVertexXY = new TH2D("hVertexXY",";x (mm);y (mm)",100,-100,100,100,-300,-150);
  auto hVertexZY = new TH1D("hVertexZY",";z (mm)",200,-100,1344);

  auto histPID = new TH2D("histPID",";p (MeV); dEdx (ADC/mm)",nbins,p1,p2,nbins,dedx1,dedx2);
  histPID -> SetTitle("All events");
  histPID -> SetTitleSize(0.04,"xy");
  histPID -> SetTitleOffset(1.4,"y");
  histPID -> SetTitleOffset(1.1,"x");
  histPID -> GetXaxis() -> CenterTitle();
  histPID -> GetYaxis() -> CenterTitle();
  
  auto histPID2 = new TH2D("histPID2",";p (MeV); dEdx (ADC/mm)",nbins,p1,p2,nbins,dedx1,dedx2);
  histPID2 -> SetTitle("Events with d<10mm");
  histPID2 -> SetTitleSize(0.04,"xy");
  histPID2 -> SetTitleOffset(1.4,"y");
  histPID2 -> SetTitleOffset(1.1,"x");
  histPID2 -> GetXaxis() -> CenterTitle();
  histPID2 -> GetYaxis() -> CenterTitle();

  auto histPID3 = new TH2D("histPID3",";p (MeV); dEdx (ADC/mm)",nbins,p1,p2,nbins,dedx1,dedx2);
  histPID3 -> SetTitle("Events with d>10mm");
  histPID3 -> SetTitleSize(0.04,"xy");
  histPID3 -> SetTitleOffset(1.4,"y");
  histPID3 -> SetTitleOffset(1.1,"x");
  histPID3 -> GetXaxis() -> CenterTitle();
  histPID3 -> GetYaxis() -> CenterTitle();

  auto histDist = new TH1D("histDist","distance (mm)",300,0,100);
  auto histDist2 = new TH1D("histDist2","distance - (-13.3 < Vz <-9.3) (mm)",300,0,100);

  auto hMult = new TH1D("hMult","Multiplicity",100,0,100);
  auto hMult2 = new TH1D("hMult2","Multiplicity - (-13.3 < Vz <-9.3) (mm)",100,0,100);
  
  Int_t feventID = -1;
  Int_t ftrackID = -1;
  Int_t fpid = -1;
  Double_t fp = -999;
  Double_t fdedx = -999;

  auto pidTest = new STPIDTest();

  Int_t nEvents = tree -> GetEntries();
  for (Int_t iEvent = 0; iEvent < nEvents; ++iEvent)
  {
    if (iEvent % 500 == 0)
      cout << "Event " << iEvent << endl;

    tree -> GetEntry(iEvent);

    if (vertexArray -> GetEntries() != 1) continue;
    auto vertex = (STVertex *) vertexArray -> At(0);
    auto posVertex = vertex -> GetPos();

    auto vx = posVertex.X();
    auto vy = posVertex.Y();
    auto vz = posVertex.Z();
    
    hVertexXY -> Fill(posVertex.X(),posVertex.Y());
    hVertexZY -> Fill(posVertex.Z());

    Int_t nTracks = trackArray -> GetEntries();

    hMult->Fill(nTracks);
    if (vz<-9.49569&&vz>-12.80121&&vx>-15&&vx<15&&vy<-206.06&&vy>-246.06)
      hMult2->Fill(nTracks);

    for (auto iTrack = 0; iTrack < nTracks; ++iTrack) {
      auto track = (STRecoTrack *) trackArray -> At(iTrack);
      auto pid = track -> GetPID();

      if (selectPID != STPID::kNon && pid != selectPID) continue;
      if (track -> GetVertexID() < 0) continue;

      auto posTrk = track->GetPosTargetPlane();
      auto distance = sqrt(pow(posTrk.X()-posVertex.X(),2)+pow(posTrk.Y()-posVertex.Y(),2)+pow(posTrk.Z()-posVertex.Z(),2));
      histDist->Fill(distance);

      auto distance2 = sqrt(pow(posTrk.X()-posVertex.X(),2)+pow(posTrk.Y()-posVertex.Y(),2));      
      if (vz<-9.49569&&vz>-12.80121&&vx>-15&&vx<15&&vy<-206.06&&vy>-246.06)
	histDist2->Fill(distance2);	
      
      
      auto dedxArray = track -> GetdEdxPointArray();
      if (dedxArray -> size() < 20) continue;

      auto p = track -> GetCharge() * track -> GetMomentum().Mag();
      auto dedx = track -> GetdEdxWithCut(0,0.7);
      histPID -> Fill(p, dedx);

      if (vz<-9.49569&&vz>-12.80121&&vx>-15&&vx<15&&vy<-206.06&&vy>-246.06){
	if (distance2 <= 10)
	  histPID2 -> Fill(p, dedx);
	else
	  histPID3 -> Fill(p, dedx);
      }
      
    }
  }
  
  auto cvsPID = new TCanvas("cvsPID","",20,20,1200,700);
  histPID -> Draw("colz"); 

  for (auto ipid = 0; ipid < STPID::GetNUMSTPID(); ++ipid) {
    auto pid = static_cast<STPID::PID>(ipid);
    pidTest -> GetdEdxFunction(pid) -> Draw("same");
  }

  auto cvsPID2 = new TCanvas("cvsPID2","",20,20,1200,700);
  histPID2 -> Draw("colz");

  for (auto ipid = 0; ipid < STPID::GetNUMSTPID(); ++ipid) {
    auto pid = static_cast<STPID::PID>(ipid);
    pidTest -> GetdEdxFunction(pid) -> Draw("same");
  }

  auto cvsPID3 = new TCanvas("cvsPID3","",20,20,1200,700);
  histPID3 -> Draw("colz"); 
  
  for (auto ipid = 0; ipid < STPID::GetNUMSTPID(); ++ipid) {
    auto pid = static_cast<STPID::PID>(ipid);
    pidTest -> GetdEdxFunction(pid) -> Draw("same");
  }
  
  auto cvsVXY = new TCanvas("cvsVXY","",20,20,400,300);
  hVertexXY -> SetStats(0);
  hVertexXY -> Draw("colz");
  auto cvsVZY = new TCanvas("cvsVZY","",420,20,400,300);
  cvsVZY -> SetLogy();
  hVertexZY -> SetStats(0);  
  hVertexZY -> Draw("colz");

  auto cvsDist = new TCanvas("cvsDist","",200,200,800,300);
  cvsDist->Divide(2,1);
  cvsDist->cd(1);
  histDist -> Draw();
  cvsDist->cd(2);
  histDist2 -> Draw();

  auto cvsMult = new TCanvas("cvsMult","",20,20,400,300);
  cvsMult->Divide(2,1);
  cvsMult->cd(1);
  hMult->Draw();
  cvsMult->cd(2);
  hMult2->Draw();

}
