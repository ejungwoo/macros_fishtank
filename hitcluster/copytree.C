#include "style.h"
using namespace style;
void copytree()
{
  auto file = new TFile("data/summary.root");
  auto tree = (TTree *) file -> Get("data");
  TCut cut = "layerCluster&&cov.X()<0.02";
  auto ofile = new TFile("data/summary2.root","recreate");
  auto otree = tree -> CopyTree(cut);
  otree -> Write();
  gApplication -> Terminate();
}
