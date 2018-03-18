#include "style.h"
using namespace style;

void covariance()
{
  gstat(0);

  auto file = new TFile("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s0.reco.develop.1533.5099ac1.test1.root");
  auto tree = (TTree *) file -> Get("cbmsim");

  TClonesArray *clusterArray = nullptr;
  TClonesArray *helixArray = nullptr;
  TClonesArray *recoArray = nullptr;

  tree -> SetBranchAddress("STHitCluster",&clusterArray);
  tree -> SetBranchAddress("STHelixTrack",&helixArray);
  tree -> SetBranchAddress("STRecoTrack",&recoArray);

  auto histX = new TH1D("histX",";dX",100,-15,5);
  auto histY = new TH1D("histY",";dY",100,-15,5);

  auto histY2 = new TH2D("histY2",";dY;layer",100,-15,15,112,0,112);

  auto numEntries = tree -> GetEntries();
  for (auto event = 0; event < numEntries; ++event) {
    tree -> GetEntry(event);
    auto numReco = recoArray -> GetEntries();
    cout << event << " " << numReco << endl;
    for (auto ireco = 0; ireco < numReco; ++ireco) {
      auto track = (STRecoTrack *) recoArray -> At(ireco);
      auto clusterIDArray = track -> GetClusterIDArray();
      for (auto id : *clusterIDArray) {
        auto cluster = (STHitCluster *) clusterArray -> At(id);
        if (cluster -> IsLayerCluster()) {
          auto poca = cluster -> GetPOCA();
          auto pos  = cluster -> GetPosition();
          auto dist = pos - poca;
          histX -> Fill(dist.X());
          histY -> Fill(dist.Y());
          histY2 -> Fill(dist.Y(),cluster->GetLayer());
        }
      }
    }
  }

  cc(); make(free(histX)) -> Draw();
  cc(); make(free(histY)) -> Draw();

  cc(); make(free(histY2)) -> Draw("colz");
}
