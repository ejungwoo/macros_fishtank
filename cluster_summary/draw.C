#include "style.h"
using namespace style;

void draw()
{
  zcolor(1);

  TString names[] = {TString("/mnt/spirit/analysis/user/leej/macros/cluster_summary/data/summary_cluster_run2900.develop.1534.2537e2e.root"),
    TString("/mnt/spirit/analysis/user/leej/macros/cluster_summary/data/summary_cluster_run2901.develop.1534.2537e2e.root"),
    TString("/mnt/spirit/analysis/user/leej/macros/cluster_summary/data/summary_cluster_cocktail100.develop.1534.2537e2e.root"),
    TString("/mnt/spirit/analysis/user/leej/macros/cluster_summary/data/summary_cluster_cocktail300.develop.1534.2537e2e.root")};

  auto idx = 0;
  for (auto name : names) {
    auto file = new TFile(name);
    auto tree = (TTree *) file -> Get("clusters");
    cc();
    tree -> Draw(Form("z:dy>>histzdy%d(200,-10,10,112,0,1344)",idx),"lr>0","colz");
    ++idx;
  }
}
