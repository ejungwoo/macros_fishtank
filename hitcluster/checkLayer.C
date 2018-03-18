#include "style.h"
using namespace style;

void checkLayer()
{
  auto input_tree = new TChain("cbmsim");
  input_tree -> Add("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s200.reco.trackChargeFix.1505.c48e1f8.test.root");

  TClonesArray *hitArray = nullptr;
  TClonesArray *clusterArray = nullptr;
  input_tree -> SetBranchAddress("STHit", &hitArray);
  input_tree -> SetBranchAddress("STHitCluster", &clusterArray);

  Long64_t countBadLayer = 0;
  Long64_t countBadRow = 0;
  Long64_t countNeither = 0;
  Long64_t countCheckedNeither = 0;

  auto numEvents = input_tree -> GetEntries();
  for (auto eventID = 0; eventID < numEvents; ++eventID)
  {
    input_tree -> GetEntry(eventID);

    auto numCluster = clusterArray -> GetEntries();
    cout << "number of clusters: " << numCluster << endl;
    for (auto clusterID = 0; clusterID < numCluster; ++clusterID)
    {
      auto cluster = (STHitCluster *) clusterArray -> At(clusterID);
      if (cluster -> GetTrackID() < 0)
        continue;
      auto ids = cluster -> GetHitIDs();
      if (ids -> size() < 2)
        continue;

      auto layerCluster = cluster -> IsLayerCluster();
      auto rowCluster   = cluster -> IsRowCluster();

      bool sameLayer = true;
      bool sameRow   = true;

      auto hit = (STHit *) hitArray -> At((*ids)[0]);
      auto layer = hit -> GetLayer();
      auto row   = hit -> GetRow();

      for (auto id : *ids) {
        hit = (STHit *) hitArray -> At(id);
        if (hit -> GetLayer() != layer) { sameLayer = false; }
        if (hit -> GetRow()   != row)   sameRow   = false;
      }

      if ( layerCluster && !sameLayer) ++countBadLayer;
      if ( rowCluster   && !sameRow)    ++countBadRow;
      if (!sameRow      && !sameLayer)  ++countNeither;
      if (!layerCluster && !rowCluster) ++countCheckedNeither;
    }
  }

  cout << "was layerCluster but not on same layer:    " << countBadLayer << endl;
  cout << "was rowCluster   but not on same row:      " << countBadRow << endl;
  cout << "neither same layer nor same row:           " << countCheckedNeither << endl;
  cout << "neither same layer nor same row (checked): " << countNeither << endl;

  gApplication -> Terminate();
}
