#include "style.h"
using namespace style;

void ana()
{
  auto file = new TFile("data/selected.root");
  auto tree = (TTree *) file -> Get("dedx");

  Int_t run, eventid, trackid, mult;
  Double_t mom, dedx;

  tree -> SetBranchAddress("run",&run);
  tree -> SetBranchAddress("eventid",&eventid);
  tree -> SetBranchAddress("trackid",&trackid);
  tree -> SetBranchAddress("mult",&mult);

  tree -> SetBranchAddress("mom",&mom);
  tree -> SetBranchAddress("dedx",&dedx);

  auto hist_dedx = new TH1D("hist","",20,50,110);

  auto ofile = new TFile("dedx.root","recreate");
  auto otree = new TTree("dedx","");

  Int_t n;
  Double_t array[500];

  otree -> Branch("n",&n);
  otree -> Branch("array",&array,"array[n]/D");
  otree -> Branch("mom",&mom);
  otree -> Branch("dedx",&dedx);

  TFile *file_reco = nullptr;
  TTree *reco = nullptr;

  for (auto entry = 0; entry < tree -> GetEntries(); ++entry)
  //for (auto entry = 0; entry < 5; ++entry)
  {
    tree -> GetEntry(entry);

    eventid -= 1;
    Int_t split = eventid/2000;
    eventid = eventid - split*2000;

    hist_dedx -> Fill(dedx);

    TString name = Form("dedxData/run%d_s%d.reco.v1.04.root",run,split);
    if (file_reco == nullptr) {
      file_reco = new TFile(name,"read");
      reco = (TTree *) file_reco -> Get("cbmsim");
    }
    else if (file_reco -> GetName() != name) {
      file_reco -> Close();
      file_reco = new TFile(name,"read");
      reco = (TTree *) file_reco -> Get("cbmsim");
    }

    cout << entry << " / " << tree -> GetEntries() << "   n:" << name << " e" << eventid << " t" << trackid << endl;

    TClonesArray *trackArray = nullptr; 
    reco -> SetBranchAddress("STRecoTrack",&trackArray);
    reco -> GetEntry(eventid);

    auto track = (STRecoTrack *) trackArray -> At(trackid);

    if (track == nullptr)
      continue;

    auto pointArray = track -> GetdEdxPointArray();

    if (pointArray -> size() < 50)
      continue;

    n = pointArray -> size();

    auto i = 0;
    bool a = false;
    for (auto i = 0; i < n; ++i) {
      auto p = pointArray -> at(i);
      array[i] = p.fdE/p.fdx;
    }

    otree -> Fill();
  }

  ofile -> cd();
  otree -> Write();

  cc();
  free(make(hist_dedx)) -> Draw();
}
