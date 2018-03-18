#include "style.h"
using namespace style;

void psa()
{
  auto file = new TFile("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s0.reco.develop.1525.d4b08aa.test.root");
  auto tree = (TTree *) file -> Get("cbmsim");

  auto pulse = new STPulse();

  TClonesArray *eventArray = nullptr;
  TClonesArray *hitArray = nullptr;
  tree -> SetBranchAddress("STRawEvent", &eventArray);
  tree -> SetBranchAddress("STHit", &hitArray);

  tree -> GetEntry(0);

  auto event = (STRawEvent *) eventArray -> At(0);

  /*
  Int_t rows[] = {48,48,50};
  Int_t layers[] = {3,10,12};

  for (auto ii = 0; ii < 3; ++ii)
  {
    auto row = rows[ii];
    auto layer = layers[ii];
    */
  for (auto row = 50; row < 53; ++row)
  for (auto layer = 3; layer < 10; ++layer)
  {
    auto ofile = new TFile(Form("psa%d_%d.root",row,layer),"recreate");

    TObjArray pulses;
    auto nHits = hitArray -> GetEntries();
    for (auto iHit = 0; iHit < nHits; ++iHit) {
      auto hit = (STHit *) hitArray -> At(iHit);
      if (hit -> GetRow() == row && hit -> GetLayer() == layer)
        pulses.Add(pulse -> GetPulseFunction(hit));
    }

    auto pad = event -> GetPad(row,layer);
    auto adc = pad -> GetADC();

    auto hist = new TH1D(Form("hist_%d_%d",row,layer),Form("Row=%d Layer=%d;Time Bucket; Charge (ADC)",row,layer),512,0,512);
    hist -> SetLineColor(kGray);
    hist -> SetLineWidth(2);
    for (auto tb = 0; tb < 512; ++tb)
      hist -> SetBinContent(tb+1,adc[tb]);

    auto cvs = c();
    free(make(hist)) -> Draw();
    cout << pulses.GetEntries() << endl;
    for (auto ipulse = 0; ipulse < pulses.GetEntries(); ++ipulse) {
      auto f1 = (TF1 *) pulses[ipulse];
      f1 -> SetLineColor(kPink);
      f1 -> SetLineWidth(1);
      f1 -> SetNpx(5000);
      f1 -> Draw("same");
    }
    auto graph = new TGraph();
    graph -> SetLineColor(kBlack);
    graph -> SetLineWidth(1);

    for (auto tb = 0.; tb < 512; tb+=0.1) {
      auto val = 0;
      for (auto ipulse = 0; ipulse < pulses.GetEntries(); ++ipulse) {
        auto f1 = (TF1 *) pulses[ipulse];
        val += f1 -> Eval(tb);
      }
      graph -> SetPoint(graph->GetN(),tb,val);
    }
    graph -> Draw("lsame");

    ofile -> cd();
    cvs -> Write("cvs");
    hist -> Write("data");
    graph -> Write("fit");
    pulses.Write("hits",TObject::kSingleKey);
  }
}
