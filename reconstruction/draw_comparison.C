#include "config_draw.h"
#include "style.h"
using namespace style;

void draw_comparison(
    //TString branch = "trackChargeFix.1501.7e5c17a",
    //TString branch = "trackChargeFix.1503.4a82220",
    //TString branch = "trackChargeFix.1508.bf6c33b",
    TString branch = "develop.1522.ac6f82e",
    bool saveFigures = false
    )
{
  //gstat(0);
  auto tail = Shorten(branch);

  vector<TString> whatToPlot;
  whatToPlot.push_back(TString("pocaV_TAG.X()"));
  whatToPlot.push_back(TString("pocaV_TAG.Y()"));
  whatToPlot.push_back(TString("pocaV_TAG.Z()"));
  //whatToPlot.push_back(TString("sqrt(pocaV_TAG.Y()*pocaV_TAG.Y()+pocaV_TAG.Z()*pocaV_TAG.Z())"));
  //whatToPlot.push_back(TString("pocaV_TAG.Mag()"));
  //whatToPlot.push_back(TString("pocaV_TAG.Perp()"));
  //whatToPlot.push_back(TString("pocaV2_TAG.Mag()"));
  //whatToPlot.push_back(TString("pocaV2_TAG.Perp()"));

  vector<TString> whatToPlot2;
  //whatToPlot2.push_back(TString("pocaV_TAG.Y():pocaV_TAG.X()"));
  //whatToPlot2.push_back(TString("pocaV2_TAG.Y():pocaV2_TAG.X()"));

  TString cut_goodV      = "goodVertex";
  TString cut_select     = "goodVertex&&dedx_TAG>50&&dedx_TAG<100&&p_TAG.Mag()>500&&p_TAG.Mag()<500.1&&charge_reco>0";
  TString cut_notDistV   = cut_goodV  + "&&length_TAG>500&&numdEdx>20";
  TString cut_normal     = cut_goodV  + "&&pocaV_reco.Mag()<5&&length_reco>500&&numdEdx>20";
  TString cut_loose      = cut_goodV  + "&&pocaV_TAG.Mag()<10&&length_TAG>300&&numdEdx>10";
  TString cut_mpion      = cut_normal + "&&pdg_TAG==211&&charge_reco<0";
  TString cut_ppion      = cut_normal + "&&charge_reco>0&&p_TAG.Mag()<1000&&(dedx_TAG<-.1*(p_TAG.Mag()-250)+20||dedx_TAG<-.5*(p_TAG.Mag()-150)+40||dedx_TAG<20)&&p_TAG.Mag()<600";
  TString cut_proton     = cut_normal + "&&pdg_TAG==2212";
  TString cut_electron   = cut_normal + "&&p_TAG.Mag()<90&&dedx_TAG<30&&charge_reco<0";
  TString cut_positron   = cut_normal + "&&p_TAG.Mag()<90&&dedx_TAG<30&&charge_reco>0";

  vector<TString> cuts;
  vector<TString> cutNames;
  //cuts.push_back(TString());    cutNames.push_back(TString("noCut"));
  cuts.push_back(cut_select);    cutNames.push_back(TString("selected"));
  cuts.push_back(cut_goodV);    cutNames.push_back(TString("goodV"));
  //cuts.push_back(cut_notDistV); cutNames.push_back(TString("NodVCut"));
  cuts.push_back(cut_normal);   cutNames.push_back(TString("normalCut"));
  cuts.push_back(cut_loose);    cutNames.push_back(TString("looseCut"));
  cuts.push_back(cut_mpion);    cutNames.push_back(TString("mpion"));
  cuts.push_back(cut_ppion);    cutNames.push_back(TString("ppion"));
  cuts.push_back(cut_proton);   cutNames.push_back(TString("proton"));
  cuts.push_back(cut_electron); cutNames.push_back(TString("electron"));
  cuts.push_back(cut_positron); cutNames.push_back(TString("positron"));

  //cuts.push_back(cut_goodV + "&&pocaV_TAG.Z()<.2 && pocaV_TAG.Z()>0 "); cutNames.push_back(TString("range_0<z<0.2"));
  //cuts.push_back(cut_goodV + "&&pocaV_TAG.Z()<.4 && pocaV_TAG.Z()>.2"); cutNames.push_back(TString("range_0.2<z<0.4"));
  //cuts.push_back(cut_goodV + "&&pocaV_TAG.Z()<.6 && pocaV_TAG.Z()>.4"); cutNames.push_back(TString("range_0.4<z<0.6"));
  //cuts.push_back(cut_goodV + "&&pocaV_TAG.Z()<.8 && pocaV_TAG.Z()>.6"); cutNames.push_back(TString("range_0.6<z<0.8"));
  //cuts.push_back(cut_goodV + "&&pocaV_TAG.Z()<1  && pocaV_TAG.Z()>.8"); cutNames.push_back(TString("range_0.8<z<1"));
  //cuts.push_back(cut_goodV + "&&pocaV_TAG.Z()<1.2&& pocaV_TAG.Z()>1"); cutNames.push_back(TString("range_1<z<1.2"));

  //TString tags[] = {TString("vadd")};
  auto numTags = 2;
  TString tags[] = {TString("vadd"), TString("reco")};
  TString options[] = {TString(""), TString("same")};
  Color_t colors[] = {kBlue, kPink};

  auto file = new TFile(TString("data/summary.")+branch+".merged2.root","read");
  auto tree = (TTree *) file -> Get("data");
  file -> Print();
  tree -> Print("toponly");

  gStyle -> SetPalette(kRainBow);
  
  auto draw1 = [tree, branch, cuts, cutNames](TString tag, TString name, TString title, TString formula, Int_t nx, Int_t x1, Int_t x2, bool setSave, Int_t iCut, bool same = false) {
    auto cut = cuts[iCut];
    cut.ReplaceAll("TAG",tag);
    cout << name << " " << tag << " " << " " << cut << " " << formula << endl;
    auto hist = new TH1D(name+"_"+tag+"."+cutNames[iCut],name+";"+title,nx,x1,x2);
    tree -> Project(name+"_"+tag+"."+cutNames[iCut],formula,cut);
    if (tag == "reco")
      hist -> SetLineColor(kPink);
    if (!same) {
      auto cvs = cc(name+"_"+tag+"."+cutNames[iCut]+"."+branch);
      free(make(hist)) -> Draw();
      cvs -> SetLogz();
    } else {
      free(make(hist)) -> Draw("same");
    }
    if (setSave) save((TCanvas *) gPad,"png");
  };

  auto draw2 = [tree, branch, cuts, cutNames](TString tag, TString name, TString title, TString formula, Int_t nx, Int_t x1, Int_t x2, Int_t ny, Int_t y1, Int_t y2, bool setSave, Int_t iCut) {
    auto cut = cuts[iCut];
    cut.ReplaceAll("TAG",tag);
    cout << name << " " << tag << " " << " " << cut << " " << formula << endl;
    auto hist = new TH2D(name+"_"+tag+"."+cutNames[iCut],name+";"+title,nx,x1,x2,ny,y1,y2);
    tree -> Project(name+"_"+tag+"."+cutNames[iCut],formula,cut);
    auto cvs = cc(name+"_"+tag+"."+cutNames[iCut]+"."+branch);
    free(make(hist)) -> Draw("colz");
    cvs -> SetLogz();
    if (setSave)
      save(cvs,"png");
  };

  draw2("vadd","PID","p/Q (MeV/c);dE/dx (ADC/mm)","dedx_vadd:p_vadd.Mag()/charge_reco",400,-800,2000,400,0,1000,true,0);
  draw2("reco","PID","p/Q (MeV/c);dE/dx (ADC/mm)","dedx_reco:p_reco.Mag()/charge_reco",400,-800,2000,400,0,1000,true,0);
  //draw2("vadd","PID","p/Q (MeV/c);dE/dx (ADC/mm)","dedx_vadd:p_vadd.Mag()/charge_reco",400,-800,2000,400,0,1000,true,2);
  //draw2("reco","PID","p/Q (MeV/c);dE/dx (ADC/mm)","dedx_reco:p_reco.Mag()/charge_reco",400,-800,2000,400,0,1000,true,2);
  //draw2("reco","PID","p/Q (MeV/c);dE/dx (ADC/mm)","dedx_reco:p_reco.Mag()/charge_reco",400,-800,2000,400,0,1000,true,4);
  //draw2("reco","PID","p/Q (MeV/c);dE/dx (ADC/mm)","dedx_reco:p_reco.Mag()/charge_reco",400,-800,2000,400,0,1000,true,5);
  return;

  //draw2("vadd","length","reco;vadd","length_vadd:length_reco",200,0,2000,200,0,2000,false,0);
  //return;

  auto iCut = 2;
  {
    draw1("vadd","POCA","","pocaV_vadd.Mag()",400,0,20, false, iCut);
    draw1("reco","POCA","","pocaV_reco.Mag()",400,0,20, false, iCut, true);

    draw1("vadd","Length","","length_vadd",200,0,2000, false, iCut);
    draw1("reco","Length","","length_reco",200,0,2000, false, iCut, true);

    draw1("vadd","num_dEdx","","numdEdx",200,0,200, false, iCut);

  }
  return;

  /*
  for (auto what : whatToPlot) {
    for (auto iCut = 0; iCut < cuts.size(); ++iCut) {
      auto cname = what+cutNames[iCut]+"."+tail; ReName(cname, "all");
      auto cvs = c(cname);
      auto lg = new TLegend();
      for (auto iTag = 0; iTag < numTags; ++iTag) {
        auto tag = tags[iTag];
        auto hname = what + "_" + cutNames[iCut]; ReName(hname, tag);
        auto draw = what + ">>" + hname; ReName(draw, tag, true);
        auto cut = cuts[iCut]; ReName(cut, tag, true);
        auto entries = tree -> GetEntries(cut);
        cout << "H1: " << draw << " WITH_CUT " << cut << ":  " << entries << endl;

        auto h1 = make(new TH1D(hname,hname+";distance from vertex (mm);",200,-10,10));
        h1 -> SetLineColor(colors[iTag]);
        h1 -> SetMinimum(0);
        lg -> AddEntry(h1,tag+Form(": %lld",entries),"l");
        tree -> Draw(draw,cut,options[iTag]);
      }
      make(lg) -> Draw();
      cvs -> SetGrid();
      if (saveFigures) save(cvs,"png");
    }
  }

  for (auto what : whatToPlot2) {
    for (auto iCut = 0; iCut < cuts.size(); ++iCut) {
      for (auto iTag = 0; iTag < numTags; ++iTag) {
        auto tag = tags[iTag];
        auto cname = what+cutNames[iCut]+"_"+tag+"."+tail; ReName(cname, "all");
        auto cvs = cc(cname);
        auto hname = what + "_" + cutNames[iCut]; ReName(hname, tag);
        auto draw = what + ">>" + hname; ReName(draw, tag, true);
        auto cut = cuts[iCut]; ReName(cut, tag, true);
        auto entries = tree -> GetEntries(cut);
        cout << "H2: " << draw << " WITH_CUT " << cut << ":  " << entries << endl;

        auto h2 = make(new TH2D(hname,hname+";dx (mm); dy (mm)",200,0,10,200,0,10));
        tree -> Draw(draw,cut,"colz");
        if (saveFigures) save(cvs,"png");
      }
    }
  }
  */

  /*
  auto nameCut = TString("copyCut");
  auto fileCut = new TFile(TString("data/")+nameCut+".root","recreate");
  auto treeCut = tree -> CopyTree(cut);
       treeCut -> Write();
   */
}
