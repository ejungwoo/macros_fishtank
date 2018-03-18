#include "../reconstruction/config_draw.h"
#include "style.h"
using namespace style;

void draw_hitcluster_info(
    //TString branch = "trackChargeFix.1504.9a9e4cc",
    TString branch = "trackChargeFix.1505.c48e1f8",
    bool saveFigures = false
    )
{
  zcolor(1);
  //gstat(0);
  auto tail = Shorten(branch);
  TString titles = "";

  vector<TString> whatToPlot;
  whatToPlot.push_back(TString("cov.X()/charge"));
  whatToPlot.push_back(TString("cov.Y()/charge"));
  whatToPlot.push_back(TString("cov.Z()/charge"));

  //vector<TString> whatToPlot;
  //whatToPlot.push_back(TString("cov.X()/100/30"));
  //whatToPlot.push_back(TString("cov.Y()/100/30"));
  //whatToPlot.push_back(TString("cov.Z()/100/30"));

  vector<TString> whatToPlot2;
  whatToPlot2.push_back(TString("cov.X():n"));
  whatToPlot2.push_back(TString("cov.Y():n"));
  whatToPlot2.push_back(TString("cov.Z():n"));

  TCut cut = "";
  //TCut cut = "cov.X()==0";
  //TCut cut = "layerCluster&&n==2";
  //TCut cut = "layerCluster&&lastCluster";//&&n!=1";
  //TCut cut = "cov.X()!=0&&cov.Y()!=0&&cov.Z()!=0";

  cout << "CUT: " << cut << endl;
  auto file = new TFile("data/summary.root");
  cout << file -> GetName() << endl;
  auto tree = (TTree *) file -> Get("data");

  auto c0 = cc("cvs_charge");
  free(make(new TH1D("hist_charge","charge;ADC",200,0,20000)));
  tree -> Draw("charge>>hist_charge");
  //c0 -> SetLogy();

  for (auto what : whatToPlot) {
    auto cname = what+"."+tail; ReName(cname);
    cout << cname << endl;
    auto cvs = cc(cname);
    auto hname = what; ReName(hname);
    auto draw = what + ">>" + hname; ReName(draw, "", true);
    //free(make(new TH1D(hname,hname+titles,200,0,1./10000)));
    free(make(new TH1D(hname,hname+titles,200,0,1)));
    tree -> Draw(draw,cut);
    cvs -> SetLogy();
  }
  return;

  for (auto what : whatToPlot2) {
    auto cname = what+"."+tail; ReName(cname);
    cout << cname << endl;
    auto cvs = cc(cname);
    auto hname = what; ReName(hname);
    auto draw = what + ">>" + hname; ReName(draw, "", true);
    make(new TH2D(hname,hname+titles,20,0,20,200,0,50));
    tree -> Draw(draw,cut,"colz");
    cvs -> SetLogz();
  }
}
