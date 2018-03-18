#include "style.h"
using namespace style;

void draw_plots
(
  //TString branch = "trackChargeFix.1501.7e5c17a",
  TString branch = "trackChargeFix.1508.bf6c33b",
  TString tag = "va"
)
{
  bool save_figures = true;

  TString branchtag = branch+"."+tag;
  TString fileName = TString("data/summary.")+branchtag+".root";

  //auto file = new TFile("data/summary.root");
  auto file = new TFile(fileName,"read");
  auto tree = (TTree *) file -> Get("data");

  TCut cut_good = "goodVertex && distVertex<5 && trackLength>500 && numdEdx>20";
  TCut cut_mpion = cut_good + "pdg==211&&trackCharge<0";
  TCut cut_ppion = cut_good + "trackCharge>0&&p.Mag()<1000&&(dedx<-.1*(p.Mag()-250)+20||dedx<-.5*(p.Mag()-150)+40||dedx<20)&&p.Mag()<600";
  TCut cut_proton = cut_good + "pdg==2212";
  TCut cut_electron = cut_good + "p.Mag()<90&&dedx<30&&trackCharge<0";
  TCut cut_positron = cut_good + "p.Mag()<90&&dedx<30&&trackCharge>0";

  /*
  //auto file_cut = new TFile("data/copy_positive-pion.root","recreate");
  auto file_cut = new TFile("data/copy.root","recreate");
  file_cut -> cd();
  auto tree_cut = tree -> CopyTree(cut);
  tree_cut -> Write();
  */

  gStyle -> SetPalette(kRainBow);
  TCanvas *cvs;

  TCut cuts[] = {cut_good,cut_mpion,cut_ppion,cut_electron,cut_positron};
  Int_t idx = 0;
  for (TCut cut : cuts) {
    auto name = Form("pid%d",idx++) + branchtag;
    cvs = cc(name);
    auto h = make(new TH2D(name,";p (MeV/c); dE/dx (ADC/mm);",500,-800,2000,500,0,1000));
    tree -> Draw(TString("dedx:p.Mag()/trackCharge>>")+name,cut,"colz");
    cvs -> SetLogz();
    if (save_figures) cvs -> SaveAs(TString("figures/goodcut_")+cvs->GetName()+".png");
  }
}
