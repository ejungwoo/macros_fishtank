#include "style.h"
using namespace style;

void draw()
{
  zcolor(1);

  TString tag = "develop.1522.ac6f82e";
  //TString tag = "develop.1526.c176211";

  auto tree = new TChain("data");
  for (auto run : {2900, 2901})
    tree -> Add(Form("/mnt/spirit/analysis/user/leej/macros/summary/data/summary%04d.%s.merged.root",run,tag.Data()));

  auto draw2 = [tree,tag](TCut cut, TString name, TString title, TString formula, Int_t nx, Int_t x1, Int_t x2, Int_t ny, Int_t y1, Int_t y2, bool saveFlag)
  {
    cout << name << " << cut:" << TString(cut) << ", formula:" << formula << endl;

    auto hist = new TH2D(name,title,nx,x1,x2,ny,y1,y2);
    tree -> Project(name,formula,cut);

    auto c2 = cc(name+"."+tag);
    free(make(hist)) -> Draw("colz");
    c2 -> SetLogz();

    if (saveFlag) save(c2,"png");
  };

  int nx=400, x1=-800, x2=3000;
  int ny=400, y1=   0, y2=1000;

  //TCut cutvrtx  = "vz>-13.8&&vz<-10.3&&vx>-15&&vx<15&&vy<-206.06&&vy>-246.06";
  //TCut cuttrck  = "ndf>30&&charge!=0&&dist<5";

  TCut cutv= "goodv";
  TCut cutt= "goodt";
  TCut cutall = cutv + cutt;
  TCut cutpi  = TCut("pion") + cutall;
  TCut cutpip = TCut("charge2>0") + cutpi;
  TCut cutpin = TCut("charge2<0") + cutpi;

  draw2(cutall, "PID",                             "PID;p/Q (MeV/c);dE/dx (ADC/mm)","dedx:p/charge",  nx,x1,x2,ny,y1,y2,1);
  draw2(cutall, "PID+Vertex",            "PID (+Vertex);p/Q (MeV/c);dE/dx (ADC/mm)","dedx2:p2/charge",nx,x1,x2,ny,y1,y2,1);
  draw2(cutpip, "PID+Vertex_pi+","PID (+Vertex) #pi^{+};p/Q (MeV/c);dE/dx (ADC/mm)","dedx2:p2/charge",nx,x1,x2,ny,y1,y2,1);
  draw2(cutpin, "PID+Vertex_pi-","PID (+Vertex) #pi^{-};p/Q (MeV/c);dE/dx (ADC/mm)","dedx2:p2/charge",nx,x1,x2,ny,y1,y2,1);
}
