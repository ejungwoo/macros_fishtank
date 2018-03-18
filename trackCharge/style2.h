namespace style2
{
  TCanvas *c    (TString name="",double w=0,double h=0);
  TCanvas *cc   (TString name="",double w=0,double h=0);

   TGraph *make (TGraph  *graph);
      TH1 *make (TH1     *h);
      TF1 *make (TF1     *f);
  TLegend *make (TLegend *legend);
  TCanvas *make (TCanvas *cvs);

     void save  (TCanvas *cvs,TString format="pdf");

  /********************************************************/

  void init();

     int fICvs=0;
     int fWCvs=680;
     int fHCvs=550;

  double fWUnit=0.006;
  double fHUnit=0.035;

  double fWStat=0.25;
  double fHStat=0.18;
  double fXStat=0.95;
  double fYStat=0.88;

  double fTMargin=0.12;
  double fBMargin=0.10;
  double fLMargin=0.12;
  double fRMargin=0.055;
  double fRMarginH1=0.03;
  double fRMarginH2=0.155;

  double fMainTitleSize=0.05;
  double fAxisLabelSize=0.03;
  double fAxisTitleSize=0.04;
  double fXTitleOffset=1.15;
  double fYTitleOffset=1.35;
};

void style2::init() {
  fXStat=1.-fRMargin;
  fYStat=1.-fTMargin;
  gStyle->SetStatStyle(0);
  gStyle->SetStatW(fWStat);
  gStyle->SetStatH(fHStat);
  gStyle->SetStatX(fXStat);
  gStyle->SetStatY(fYStat);
  gStyle->SetTitleFontSize(fMainTitleSize);
}

TCanvas *style2::c(TString name,double w,double h) {
  fRMargin=fRMarginH1;
  init();
  if(w==0) w=fWCvs;
  if(h==0) h=fHCvs;
  if(name.IsNull()) name=Form("canvas-%d",fICvs);
  auto cvs=new TCanvas(name,name,(fICvs+1)*20,(fICvs+1)*20,w,h);
  make(cvs);
  return cvs;
}

TCanvas *style2::cc(TString name,double w,double h) {
  fRMargin=fRMarginH2;
  init();
  if(w==0) w=fWCvs;
  if(h==0) h=fHCvs;
  if(name.IsNull()) name=Form("canvas-%d",fICvs);
  auto cvs=new TCanvas(name,name,(fICvs+1)*20,(fICvs+1)*20,w,h);
  make(cvs);
  return cvs;
}

TCanvas *style2::make(TCanvas *cvs) {
  cvs->SetMargin(fLMargin,fRMargin,fBMargin,fTMargin);
  ++fICvs;
  return cvs;
}

TGraph *style2::make(TGraph *graph) {
  init();
  if(graph->GetMarkerStyle()==1)
    graph->SetMarkerStyle(20);
  make(graph->GetHistogram());
  return graph;
}

TH1 *style2::make(TH1 *h) {
  init();
  h->GetXaxis()->CenterTitle();
  h->GetXaxis()->SetTitleOffset(fXTitleOffset);
  h->GetXaxis()->SetTitleSize(fAxisTitleSize);
  h->GetXaxis()->SetLabelSize(fAxisLabelSize);
  h->GetYaxis()->CenterTitle();
  h->GetYaxis()->SetTitleOffset(fYTitleOffset);
  h->GetYaxis()->SetTitleSize(fAxisTitleSize);
  h->GetYaxis()->SetLabelSize(fAxisLabelSize);
  h->GetZaxis()->SetLabelSize(fAxisLabelSize);
  return h;
}

TLegend *style2::make(TLegend *legend) {
  init();
  auto y2=1.-fTMargin;
  if(gStyle->GetOptStat() != 0) y2=gStyle->GetStatY()-gStyle->GetStatH();
  auto lmax=0;
  TIter next(legend->GetListOfPrimitives());
  while (auto make=(TLegendEntry *) next()) {
    auto l=TString(make->GetLabel()).Length();
    if(lmax<l) lmax=l;
  }
  
  lmax+=3;
  auto y1=y2-legend->GetNRows()*fHUnit;
  auto x2=1.-fRMargin;
  auto x1=x2-lmax*fWUnit;
  if(x1>x2-fWStat) x1=x2-fWStat;
  legend->SetX1(x1);
  legend->SetX2(x2);
  legend->SetY1(y1);
  legend->SetY2(y2);
  legend->SetFillStyle(1001);
  return legend;
}

TF1 *style2::make(TF1 *f) {
  f -> SetNpx(1000);
  return f;
}

void style2::save(TCanvas *cvs,TString format) {
  TString path=TString(gSystem->Getenv("PWD"))+"/figures/";
  gSystem->Exec(TString("mkdir -p ")+path);
  TString head=cvs->GetName();
  TString name=head+"."+format;
  int version=0;
  while (1) {
    char* c=gSystem->Which(path.Data(),name.Data());
    if(TString(c).IsNull()) {
      cvs->SaveAs(path+name);
      break;
    }
    name=head+"."+TString::Itoa(++version,36)+"."+format;
  } 
}
