void draw2()
{
  //gStyle -> SetOptStat(0);
  gStyle -> SetPalette(55);
  gStyle -> SetNumberContours(100);
  gStyle -> SetPadRightMargin(0.12);
  gStyle -> SetPadTopMargin(0.05);
  gStyle -> SetPadBottomMargin(0.15);
  gStyle -> SetPadLeftMargin(0.13);
  gStyle -> SetTitleOffset(1.15, "x");
  gStyle -> SetTitleOffset(1.05, "y");
  gStyle -> SetTitleSize(0.06, "x");
  gStyle -> SetTitleSize(0.06, "y");
  gStyle -> SetTitleSize(0.06, "z");
  gStyle -> SetLabelSize(0.06, "x");
  gStyle -> SetLabelSize(0.06, "y");
  gStyle -> SetLabelSize(0.06, "z");
  gStyle -> SetLegendTextSize(0.06);
  gStyle -> SetOptLogz();

  auto tree = new TChain("data");
  tree -> Add("/mnt/spirit/analysis/user/leej/macros/summary/data/summary2900.develop.1522.ac6f82e.fit.root");
  tree -> Add("/mnt/spirit/analysis/user/leej/macros/summary/data/summary2901.develop.1522.ac6f82e.fit.root");
  //tree -> Add("/mnt/spirit/analysis/user/leej/macros/summary/data/summary2900.develop.1526.c176211.merged.root");
  //tree -> Add("/mnt/spirit/analysis/user/leej/macros/summary/data/summary2901.develop.1526.c176211.merged.root");
  //tree -> Add("/mnt/spirit/analysis/user/leej/macros/summary/data/summary2902.develop.1522.ac6f82e.fit.root");
  //tree -> Add("/mnt/spirit/analysis/user/leej/macros/summary/data/summary2903.develop.1522.ac6f82e.fit.root");
  //tree -> Add("/mnt/spirit/analysis/user/leej/macros/summary/data/summary2904.develop.1522.ac6f82e.fit.root");
  //tree -> Add("/mnt/spirit/analysis/user/leej/macros/summary/data/summary2905.develop.1522.ac6f82e.fit.root");

  auto draw2 = [tree](TCut cut, TString name, TString title, TString formula, Int_t nx, Int_t x1, Int_t x2, Int_t ny, Int_t y1, Int_t y2, bool saveFlag)
  {
    cout << name << " << cut:" << TString(cut) << ", formula:" << formula << endl;

    auto hist = new TH2D(name,title,nx,x1,x2,ny,y1,y2);
    tree -> Project(name,formula,cut);

    auto c2 = new TCanvas("cvs", "", 840, 500);
    hist -> Draw("colz");
    c2 -> SetLogz();

    //if (saveFlag) save(c2,"png");
  };

  int nx=400, x1=-800, x2=3000;
  int ny=400, y1=   0, y2=1000;

  //TCut cutvrtx  = "vz>-13.8&&vz<-10.3&&vx>-15&&vx<15&&vy<-206.06&&vy>-246.06";
  //TCut cuttrck  = "ndf>30&&charge!=0&&dist<5";

  TCut cutv= "goodv";
  TCut cutt= "ndf>30&&charge!=0&&dist<5";
  TCut cutall = cutv + cutt;
  TCut cutpi  = TCut("pion") + cutall;
  TCut cutpip = TCut("charge2>0") + cutpi;
  TCut cutpin = TCut("charge2<0") + cutpi;

  //draw2(cutall, "PID",                             "PID;p/Q (MeV/c);dE/dx (ADC/mm)","dedx:p/charge",  nx,x1,x2,ny,y1,y2,0);
  draw2(cutall, "PID+Vertex",            "PID (+Vertex);p/Q (MeV/c);dE/dx (ADC/mm)","dedx2:p2/charge",nx,x1,x2,ny,y1,y2,0);
  //draw2(cutpip, "PID+Vertex_pi+","PID (+Vertex) #pi^{+};p/Q (MeV/c);dE/dx (ADC/mm)","dedx2:p2/charge",nx,x1,x2,ny,y1,y2,0);
  //draw2(cutpin, "PID+Vertex_pi-","PID (+Vertex) #pi^{-};p/Q (MeV/c);dE/dx (ADC/mm)","dedx2:p2/charge",nx,x1,x2,ny,y1,y2,0);
}
