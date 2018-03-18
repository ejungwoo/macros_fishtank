#include "style2.h"
using namespace style2;

void eve
(
 Int_t runID = 2900,
 Int_t split = 0,
 Int_t entry = 0,
 Int_t trackID = 1
 )
{
  gstat(0);

  //TString name = Form("r%d.s%d.e%d.t%d",runID,split,entry,trackID);
  TString name2 = Form("r%04d_s%02d_e%04d_t%03d",runID,split,entry,trackID);
  TString fileName = Form("/mnt/spirit/analysis/user/leej/reconstruction/data/run%d_s%d.reco.trackChargeFix.1494.0f0b0d1.root",runID,split);

  auto file = new TFile(fileName,"read");
  auto tree = (TTree *) file -> Get("cbmsim");
  TClonesArray *trackArray = new TClonesArray("STRecoTrack",100);
  TClonesArray *vertexArray = new TClonesArray("STVertex",100);
  tree -> SetBranchAddress("STRecoTrack",&trackArray);
  tree -> SetBranchAddress("STVertex",&vertexArray);
  tree -> GetEntry(entry);

  TGraph *graph_selected = new TGraph();
  graph_selected -> SetMarkerStyle(21);
  graph_selected -> SetMarkerSize(0.5);
  graph_selected -> SetMarkerColor(kPink);
  graph_selected -> SetLineColor(kPink);

  TGraph *graph_selected2 = new TGraph();
  graph_selected2 -> SetMarkerStyle(21);
  graph_selected2 -> SetMarkerSize(0.5);
  graph_selected2 -> SetMarkerColor(kPink);
  graph_selected2 -> SetLineColor(kPink);

  TGraph *graph_all = new TGraph();
  graph_all -> SetMarkerStyle(24);
  graph_all -> SetMarkerSize(0.5);
  graph_all -> SetMarkerColor(kGray);
  graph_all -> SetLineColor(kGray);

  TGraph *graph_all2 = new TGraph();
  graph_all2 -> SetMarkerStyle(24);
  graph_all2 -> SetMarkerSize(0.5);
  graph_all2 -> SetMarkerColor(kGray);
  graph_all2 -> SetLineColor(kGray);

  auto vertex = (STVertex *) vertexArray -> At(0);
  TMarker *mVertex = nullptr;
  TMarker *mVertex2 = nullptr;
  if (vertex != nullptr) {
    posVertex = vertex -> GetPos();
    auto vx = posVertex.X();
    auto vy = posVertex.Y();
    auto vz = posVertex.Z();
    mVertex = new TMarker(vz, vx, 3);
    mVertex -> SetMarkerSize(2);
    mVertex -> SetMarkerColor(kBlue);
    mVertex2 = new TMarker(vz, vy, 3);
    mVertex2 -> SetMarkerSize(2);
    mVertex2 -> SetMarkerColor(kBlue);
  }

  auto numTracks = trackArray -> GetEntries();
  for (auto iTrack = 0; iTrack < numTracks; ++iTrack) {

    auto track = (STRecoTrack *) trackArray -> At(iTrack);
    auto dedxArray = track -> GetdEdxPointArray();
    for (auto dedx : *dedxArray) {
      auto pos = dedx.fPosition;
      if (iTrack == trackID)
        graph_selected -> SetPoint(graph_selected->GetN(),pos.Z(),pos.X());
      else
        graph_all -> SetPoint(graph_all->GetN(),pos.Z(),pos.X());
    }
  }

  cc();
  auto frame = new TH2D("framezx",name+";z;x",100,-50,1344,100,-432,432);
  make(frame) -> Draw();
  if (mVertex != nullptr) mVertex -> Draw("psame");
  make(graph_all) -> Draw("psame");
  make(graph_selected) -> Draw("psame");

  cc();
  auto frame2 = new TH2D("framezy",name+" SIDE;z;y",100,-50,1344,100,-530,0);
  make(frame2) -> Draw();
  if (mVertex2 != nullptr) mVertex -> Draw("psame");
  make(graph_all2) -> Draw("psame");
  make(graph_selected2) -> Draw("psame");
}
