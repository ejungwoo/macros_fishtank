#include "style.h"
using namespace style;

void hit_and_cluster(
    TString branch = "trackChargeFix.1504.9a9e4cc",
    bool saveFigures = false
    )
{
  auto file = new TFile("/mnt/spirit/analysis/user/leej/macros/reconstruction/data/run2900_s10.reco.trackChargeFix.1504.9a9e4cc.test.root");
  auto tree = (TTree *) file -> Get("cbmsim");

  TClonesArray *hits = nullptr;
  TClonesArray *clusters = nullptr;
  tree -> SetBranchAddress("STHit",&hits);
  tree -> SetBranchAddress("STHitCluster",&clusters);

  auto entries = tree -> GetEntries();
  for (auto entry = 0; entry < entries; ++entry)
  {
    tree -> GetEntry(entry);
    cout << "Event: " << entry << endl;
    auto numClusters = clusters -> GetEntries();
    //for (auto iCluster = 0; iCluster < numClusters; ++iCluster) {
    for (auto iCluster = 0; iCluster < 10; ++iCluster) {
      auto cluster = (STHitCluster *) clusters -> At(iCluster);
      auto hitIDs = cluster -> GetHitIDs();
      cout << hitIDs -> size() << endl;
      cout << cluster -> GetCharge() << " " << cluster -> GetX() << " " << cluster -> GetY() << " " << cluster -> GetZ() << endl;
      for (auto id : *hitIDs) {
        auto hit = (STHit *) hits -> At(id);
        cout << "  " << hit -> GetCharge() << " " << hit -> GetX() << " " << hit -> GetY() << " " << hit -> GetZ() << endl;
      }
    }
    return;
  }
}
