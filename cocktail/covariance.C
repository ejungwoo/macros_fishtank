#include "style.h"
using namespace style;

void covariance()
{
  gstat(0);
  zcolor(1);

  auto GetCobo = [](Double_t z, Double_t x)
  {
    auto layer = (Int_t)(z/12);
    auto row = (Int_t)((x+432)/8);

    if (layer < 28) {
           if (row < 36) return 0;
      else if (row < 54) return 1;
      else if (row < 90) return 3;
      else               return 4;
    }
    else if (layer < 56) {
           if (row < 18) return 1;
      else if (row < 54) return 2;
      else if (row < 72) return 4;
      else               return 5;
    }
    if (layer < 84) {
           if (row < 36) return 6;
      else if (row < 54) return 7;
      else if (row < 90) return 9;
      else               return 10;
    }
    else {
           if (row < 18) return 7;
      else if (row < 54) return 8;
      else if (row < 72) return 10;
      else               return 11;
    }
  };

  auto tree = new TChain("cbmsim");
  int runs[] =   {3187, 3191, 3193, 3202, 3203, 3204};
  int splits[] = {   5,    5,    8,    6,    5,    5};

  //auto irun = 0;
  //for (auto irun : {1})
  for (auto irun : {0,1,2})
  //for (auto irun : {3,4,5})
  for (auto split = 0; split < splits[irun]; ++split)
    tree -> Add(Form("data/run%d_s%d.reco.develop.1533.5099ac1.root",runs[irun],split));

  //tree -> Add("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s0.reco.develop.1533.5099ac1.test1.root");
  //tree -> Add("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s1.reco.develop.1533.5099ac1.test1.root");
  //tree -> Add("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s2.reco.develop.1533.5099ac1.test1.root");

  TClonesArray *clusterArray = nullptr;
  TClonesArray *helixArray = nullptr;
  TClonesArray *recoArray = nullptr;

  tree -> SetBranchAddress("STHitCluster",&clusterArray);
  tree -> SetBranchAddress("STHelixTrack",&helixArray);
  tree -> SetBranchAddress("STRecoTrack",&recoArray);

  auto histdX = new TH1D("dx",";dx",500,-8,8);
  auto histdY = new TH1D("dy",";dy",500,-8,8);
  auto histdYvsL = new TH2D("dy_vs_layer",";dY;layer",500,-8,8,112,0,112);
  auto histdYvsX = new TH2D("dy_vs_x",";dY;x",500,-8,8,108,-432,432);

  auto numEntries = tree -> GetEntries();
  for (auto event = 0; event < numEntries; ++event) {
    tree -> GetEntry(event);
    auto numReco = recoArray -> GetEntries();
    for (auto ireco = 0; ireco < numReco; ++ireco) {
      auto track = (STRecoTrack *) recoArray -> At(ireco);

      //auto dedxArray = track -> GetdEdxPointArray();
      //if (dedxArray -> size() < 10) continue;

      auto clusterIDArray = track -> GetClusterIDArray();
      for (auto id : *clusterIDArray) {
        auto cluster = (STHitCluster *) clusterArray -> At(id);

        if (cluster -> IsLayerCluster()) {
          auto poca = cluster -> GetPOCA();
          auto pos  = cluster -> GetPosition();
          auto dist = pos - poca;
          histdX -> Fill(dist.X());

          auto cobo = GetCobo(pos.Z(),pos.X());
          if (cobo == 9)
          {
            histdY -> Fill(dist.Y());
            histdYvsL -> Fill(dist.Y(),cluster->GetLayer());
            histdYvsX -> Fill(dist.Y(),pos.X());
          }

        }
      }
    }
  }

  auto draw = [](TH1 *h) {
    auto c1 = c(h->GetName());
    make(free(h)) -> Draw("colz");
    save(c1);
  };

  draw(histdX);
  draw(histdY);
  draw(histdYvsL);
  draw(histdYvsX);
}
