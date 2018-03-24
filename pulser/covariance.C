#include "style.h"
using namespace style;

void covariance()
{
  gstat(1);
  zcolor(1);

  auto GetCobo = [](Double_t z, Double_t x)
  {
    auto layer = (Int_t)(z/12);
    auto row = (Int_t)((x+432)/8);

    if(row >8) return -1;

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
  tree -> Add("/mnt/spirit/analysis/user/leej/macros/pulser/data/run3236_s0.reco.develop.1534.2537e2e.root");

  TClonesArray *hitArray = nullptr;
  tree -> SetBranchAddress("STHit",&hitArray);

  auto histdY = new TH1D("dy",";dy",500,-65,-55);
  auto histdYvsL = new TH2D("dy_vs_layer",";dY;layer",500,-65,-55,112,0,112);
  auto histdYvsR = new TH2D("dy_vs_x",";dY;x",500,-65,-55,108,0,108);

  auto numEntries = tree -> GetEntries();
  for (auto event = 0; event < numEntries; ++event) {
    tree -> GetEntry(event);
    auto numHits = hitArray -> GetEntries();
    for (auto ihit = 0; ihit < numHits; ++ihit) {
      auto hit = (STHit *) hitArray -> At(ihit);
      auto pos = hit -> GetPosition();
      if (GetCobo(pos.Z(),pos.X()) == 6) {
          histdY -> Fill(pos.Y());
          histdYvsL -> Fill(pos.Y(),hit -> GetLayer());
          histdYvsR -> Fill(pos.Y(),hit -> GetRow());
      }
    }
  }

  auto draw = [](TH1 *h) {
    auto c1 = c(h->GetName());
    make(free(h)) -> Draw("colz");
    //save(c1);
  };

  draw(histdY);
  draw(histdYvsL);
  draw(histdYvsR);
}
