#include "style2.h"
using namespace style2;

bool saveFlag = false;
TString fName;

Double_t xx(TVector3 p, bool top = true) { return p.Z(); }
Double_t yy(TVector3 p, bool top = true) { if (top) return p.X(); return p.Y(); }

Int_t colorPositive = kGray;
Int_t colorNegative = kRed;

void drawEvent(TTree *tree, Int_t entry)
{
  TClonesArray *clusterArray = nullptr;
  TClonesArray *helixArray = nullptr;
  TClonesArray *trackArray = nullptr;
  TClonesArray *vertexArray = nullptr;

  tree -> SetBranchAddress("STHitCluster", &clusterArray);
  tree -> SetBranchAddress("STHelixTrack", &helixArray);
  tree -> SetBranchAddress("STRecoTrack", &trackArray);
  tree -> SetBranchAddress("STVertex", &vertexArray);

  tree -> GetEntry(entry);

  vector<TText*> ids_pid;
  auto cvs_pid = c("pid");
  auto hist_pid = new TH2D("pid",fName+" PID;p/Q (MeV/c); dE/dx (ADC/mm)",50,-1000,2000,50,0,800);
  make(hist_pid) -> Draw();
  cvs_pid -> SetGrid();

  auto graph_pid = new TGraph();
       graph_pid -> SetMarkerStyle(24);
       graph_pid -> SetMarkerSize(2);

  for (bool isTopView : {false, true})
  {
    TCanvas *cvs = nullptr;
    TH2D *frame = nullptr;

    if (isTopView) {
      //cvs = c(Form("event%d_top",entry),1000,750);
      cvs = c(Form("event%d_top",entry),900,700);
      frame = new TH2D("frame",fName+" TOP;z;x",100,-50,1344,100,-432,432);
    } else {
      cvs = c(Form("event%d_side",entry),900,700);
      frame = new TH2D("frame2",fName+" SIDE;z;y",100,-50,1344,100,-530,0);
    }

    make(frame) -> Draw();

    vector<TGraph*> vCPoints;
    vector<TGraph*> vLine;
    vector<TText*> ids_eve;

    TGraph *cluster_positive = nullptr;
    TGraph *cluster_negative = nullptr;

    STVertex *vertex = (STVertex *) vertexArray -> At(0);
    if (vertex == nullptr) {
      if (isTopView)
        cout << ">>>>>>>> NO VERTEX" << endl;
      //return;
    }

    TVector3 posVertex;
    if (vertex != nullptr) {
      posVertex = vertex -> GetPos();
      auto vx = posVertex.X();
      auto vy = posVertex.Y();
      auto vz = posVertex.Z();

      if (vz<-9.49569&&vz>-12.80121&&vx>-15&&vx<15&&vy<-206.06&&vy>-246.06) {
        if (isTopView)
          cout << ">>>>>>>> VERTEX NOT IN RANGE" << endl;
        //return;
      }
    }

    Int_t numHelices = helixArray -> GetEntries();
    if (numHelices == 0)
      return;
    for (auto iHelix = 0; iHelix < numHelices; ++iHelix)
    {
      auto helix = (STHelixTrack *) helixArray -> At(iHelix);
      auto trackID = helix -> GetGenfitID();
      //if (helix -> IsPositiveChargeParticle()) cout << trackID << " (+)" << endl;
      //else cout << trackID << " (-)" << endl;
      if (trackID < 0)
        continue;

      auto track = (STRecoTrack *) trackArray -> At(trackID);
      if (track == nullptr) {
        if (isTopView)
          cout << ">>>>>>>> track is nullptr" << endl;
        continue;
      }

      auto charge = track -> GetCharge();
      if (vertex != nullptr) {
        if (track -> GetVertexID() != 0)
          continue;
        auto pocaVertex = track -> GetPOCAVertex();
        auto distVertex = (posVertex - pocaVertex).Mag();
        if (distVertex > 10) {
          continue;
        }
      }

      //auto clusterIDs = track -> GetClusterIDArray();
      auto dedxArray = track -> GetdEdxPointArray();
      if (dedxArray -> size() < 10) {
        //continue;
      }

      if (isTopView) {
        auto x = (track -> GetMomentum()).Mag() / charge;
        auto y = track -> GetdEdxWithCut(0,0.7);
        graph_pid -> SetPoint(graph_pid -> GetN(), x, y);

        auto id_pid = new TText(0,0,Form("%d",trackID));
        id_pid -> SetTextSize(0.025);
        id_pid -> SetTextAlign(22);
        id_pid -> SetX(x);
        id_pid -> SetY(y);
        if (charge < 0) {
          id_pid -> SetTextColor(colorNegative);
          if (isTopView)
            cout << "trackID x y y-Length-In-Period" << endl;
            cout << trackID << " " << x << " " << y << " " << helix -> YLengthInPeriod () << endl;
        }
        ids_pid.push_back(id_pid);
      }

      if (isTopView && (trackID == 0 || trackID == 22)) {
        cout << "-------------------------------------------------------------------------------------------" << endl;
        helix -> DetermineParticleCharge(posVertex);
        cout << "-------------------------------------------------------------------------------------------" << endl;
      }

      auto clusterEnd = (STHitCluster *) clusterArray -> At((dedxArray->at(dedxArray->size()-1)).fClusterID);
      auto positionEnd = clusterEnd -> GetPosition();

      auto id_eve = new TText(0,0,Form("%d",trackID));
      id_eve -> SetTextSize(0.025);
      id_eve -> SetTextAlign(22);
      id_eve -> SetX(xx(positionEnd, isTopView));
      id_eve -> SetY(yy(positionEnd, isTopView));
      if (charge < 0) id_eve -> SetTextColor(colorNegative);
      ids_eve.push_back(id_eve);

      auto clusters = new TGraph();
      //clusters -> SetMarkerStyle(21);
      clusters -> SetMarkerStyle(24);
      clusters -> SetMarkerSize(0.3);
      clusters -> SetMarkerColor(colorPositive);
      clusters -> SetLineColor(colorPositive);
      vCPoints.push_back(clusters);

      vector<STHitCluster*> vCluster;
      for (auto dedx : *dedxArray) {
        auto clusterID = dedx.fClusterID;
        auto cluster = (STHitCluster *) clusterArray -> At(clusterID);
        auto position = cluster -> GetPosition();
        clusters -> SetPoint(clusters -> GetN(), xx(position, isTopView), yy(position, isTopView));
        vCluster.push_back(cluster);
      }

      auto line = new TGraph();
      line -> SetLineColor(colorPositive);
      vLine.push_back(line);

      auto lengthi = 0.;
      auto lengthf = 1.;
      if (charge < 0) {
        lengthi = -1.;
        lengthf = 2.;
      }

      for (Double_t r=lengthi; r<lengthf; r+=0.01) {
        auto position = helix -> InterpolateByRatio(r);
        line -> SetPoint(line -> GetN(), xx(position, isTopView), yy(position, isTopView));
      }

      if (track -> GetCharge() < 0) {
        cluster_negative = clusters;
        clusters -> SetMarkerColor(colorNegative); 
        clusters -> SetLineColor(kBlack); 
        line -> SetLineColor(kBlack); 
        clusters -> SetMarkerSize(0.5); 
      }
      else
        cluster_positive = clusters;
    }

    for (auto line : vLine) line -> Draw("lsame");
    for (auto clusters : vCPoints) clusters -> Draw("psame");
    for (auto t : ids_eve) t -> Draw("same");

    auto legend = new TLegend();
    legend -> AddEntry(cluster_positive,"(+) charged track","lp");
    legend -> AddEntry(cluster_negative,"(-) charged track","lp");
    legend -> AddEntry(ids_eve[0],"No. for track-ID","");
    make(legend) -> Draw();
    if (saveFlag) save(cvs,"png");
    
    if (vertex != nullptr) {
      auto m = new TMarker(xx(posVertex,isTopView),yy(posVertex,isTopView),3);
      m -> SetMarkerSize(2);
      m -> SetMarkerColor(kPink);
      m -> Draw("samep");
    }
  }

  cvs_pid -> cd();
  make(graph_pid) -> Draw("psame");
  for (auto t : ids_pid) t -> Draw("same");
  if (saveFlag) save(cvs_pid);
}

void eve()
{
  gStyle -> SetOptStat(0);

  auto file = new TFile("/user/leej/macros/charge/data/run2841_s5991.FIX3.root");
  auto tree = (TTree *) file -> Get("cbmsim"); 

  int entry = 0;
  saveFlag = false;
  fName = Form("RUN-2841 EVENT-%d",entry); 

  drawEvent(tree, entry);
}
