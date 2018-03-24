#include <cmath>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>
#include <unistd.h>

#include        "TROOT.h"
#include        "TAttText.h"
#include        "TAxis.h"
#include        "TCanvas.h"
#include        "TChain.h"
#include        "TCut.h"
#include        "TF1.h"
#include        "TFile.h"
#include        "TGraph.h"
#include        "TGraphAsymmErrors.h"
#include        "TGraphErrors.h"
#include        "TH1.h"
#include        "TH2.h"
#include        "THistPainter.h"
#include        "TCutG.h"
#include        "TKey.h"
#include        "TLatex.h"
#include        "TLegend.h"
#include        "TMath.h"
#include        "TMatrixD.h"
#include        "TMinuit.h"
#include        "TMultiGraph.h"
#include        "TNtuple.h"
#include        "TPave.h"
#include        "TPaveText.h"
#include        "TPoint.h"
#include        "TRandom.h"
#include        "TRandom3.h"
#include        "TRint.h"
#include        "TStyle.h"
#include        "TString.h"
#include        "TTree.h"
#include        "TH1F.h"
#include        "TH2F.h"
#include        "TProfile.h"

using namespace std;

Int_t RunNumber=999;
Int_t run2run[3000];

Double_t nnp_traw_to_tcal(Int_t,Int_t,Int_t);
Double_t start_traw_to_tcal(Int_t,Int_t);

Double_t T0 = 750.;
Double_t VT0 = 1102.;

const Double_t Distance = 856.105 + 280; // Mizuki's photogrametry to the blue frame + to the 1st plane surface
const Double_t vDistance = 856.105 + 280 - 32.5; // wrong 812.3;  // EOS experiment distance 32.5 cm
const Double_t th0 = 29.579 * 3.14159265 / 180.;

const Double_t c=29.979245; //cm/ns

const int num_planes=4;
const int num_bars=400;

const Int_t EhitRange=400;  // was 600 

Double_t pedestal1[num_bars],pedestal2[num_bars];

Double_t vped1[8]={0.82,0.57,0.42,0.87,0.76,0.68,1.41,1.35};
Double_t vped2[8]={1.00,1.11,0.82,1.18,1.23,0.78,0.78,0.58};
//Double_t vtoff[8]={142.7,139.4,137.7,140.7,162.8,164.7,163.8,142.3};

Double_t vtoff[8]={142.84+0.0404,
		   139.49+0.1797,
		   137.56+0.1781,
		   140.71-0.0091,
		   162.86+0.0071,
		   164.88-0.1448,
		   163.66+0.1192,
		   142.24+0.0975};

Double_t vtoffcorr[8]={2.2406,
		       2.3823,
		       2.2606,
		       1.9740,
		       2.9501,
		       2.8075,
		       2.4053,
		       0};

//Double_t vtdiff[8]={42.02,39.22,45.62,48.72,-9.08,-2.38,10.32,-46.08};

Double_t vtdiff[8]={42.335-0.0011,
		    39.541+0.0125,
		    45.840-0.0088,
		    48.954+0.0248,
		    -9.097+0.0270,
		    -2.366-0.0075,
		    10.329+0.0261,
		    -46.305+0.0316};

//Double_t vtdiff[8]={0,0,0,0,0,0,0,0};

//Double_t vtdiffcorr[8]={};
Double_t vtdiffcorr[8]={1.408+0.1121,
			1.109-0.0250,
			1.365-0.0043,
			1.868+0.0465,
			-0.720-0.0437,
			0.383+0.014,
			-2.091-0.0413,
			0};

struct TCluster {
  Int_t dim;
  Double_t xf;
  Double_t yf;
  Double_t zf;
  Double_t tf;
  Double_t xl;
  Double_t yl;
  Double_t zl;
  Double_t tl;
  Double_t ax;
  Double_t ay;
  Double_t etot;
};
TCluster clusters[num_bars+1];

const int ENTRY_TABLE_LEN = 100;
struct TCal {
  Int_t raw;
  Double_t ns;
};
struct TCalTable {
  TCal table[ENTRY_TABLE_LEN];
  size_t table_len;
};
TCalTable tcal_start[6 + 1];
TCalTable tcal_nnp[num_bars + 1][6 + 1];
Double_t tdiff[num_bars+1],tsync[num_bars+1],vscint[num_bars+1];
Double_t ediff[num_bars+1],esync[num_bars+1],att[num_bars+1];
Double_t tsynccorr[num_bars+1]={ };

Double_t xwalk[10],awalk[10],bwalk[10];

Double_t wlk(Double_t x);
Double_t QDC(Double_t channel);

TCutG *helium, *protons, *deuterons, *tritons, *z1, *vdoublescut, *fastz1, *discont1, *discon1a, *dyvsYcut, *dxvsXcut;
TCutG *protons_tpc, *deuterons_tpc, *tritons_tpc, *YvsXtpccut, *XvsZtpccut, *dyvsdxtpccut, *dyvsdztpccut;
void opencuts() {
  TFile *f = new TFile("helium.root");
  helium = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("protons.root");
  protons = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("deuterons.root");
  deuterons = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("tritons.root");
  tritons = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("Z1.root");
  z1 = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("vdoubles.root");
  vdoublescut = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("fastZ1.root");
  fastz1 = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("discont1.root");
  discont1 = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("discon1a.root");
  discon1a = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("dyvsYcut.root");
  dyvsYcut = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("dxvsXcut1.root");
  dxvsXcut = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("protons_tpc.root");
  protons_tpc = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("deuterons_tpc.root");
  deuterons_tpc = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("tritons_tpc.root");
  tritons_tpc = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("YvsXtpccut.root");
  YvsXtpccut = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("XvsZtpccut.root");
  XvsZtpccut = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("dyvsdxtpccut.root");
  dyvsdxtpccut = (TCutG*)f->Get("CUTG");
  f->Close();
  f = new TFile("dyvsdztpccut.root");
  dyvsdztpccut = (TCutG*)f->Get("CUTG");
  f->Close();
}

TH1F* hMult = new TH1F("hMult","Multiplicity",100,-0.5,99.5);
TH1F* hMultPart = new TH1F("hMultPart","Multiplicity partial hits",100,-0.5,99.5);
TH1F* hMultPmt = new TH1F("hMultPmt","Multiplicity pmt",100,-0.5,99.5);
TH1F* hMultV = new TH1F("hMultV","Multiplicity vertical plane",100,-0.5,99.5);
TH1F* hMultH = new TH1F("hMultH","Multiplicity horizontal plane",100,-0.5,99.5);

TH1F* hMultv = new TH1F("hMultv","Veto time Multiplicity",10,-0.5,9.5);
TH1F* hMultvq = new TH1F("hMultvq","Veto charge Multiplicity",10,-0.5,9.5);
TH2F* hVdoubles = new TH2F("hVdoubles","Veto mult in bar vs bar",8,0.5,8.5,10,-0.5,9.5);
TH2F* hVdoublesc = new TH2F("hVdoublesc","Veto mult in bar vs bar c",8,0.5,8.5,10,-0.5,9.5);

TH2F* hEvsBar = new TH2F("hEvsBar","Energy vs Bars",num_bars,0.5,num_bars+0.5,2000,0,200);

TH2F* hEvsBarcosm = new TH2F("hEvsBarcosm","Energy vs Bars cosm",num_bars,0.5,num_bars+0.5,2000,0,200);
TProfile* hEvsBarcosmprof = new TProfile("hEvsBarcosmprof","Energy vs Bars cosm",num_bars,0.5,num_bars+0.5,0,2000);
TH1F* hBars = new TH1F("hBars","Bars",num_bars,0.5,num_bars+0.5);

TH2F* hT1vsBar = new TH2F("hT1vsBar","Time1 vs Bars",num_bars,0.5,num_bars+0.5,1000,-2000,1000);
TH2F* hT2vsBar = new TH2F("hT2vsBar","Time2 vs Bars",num_bars,0.5,num_bars+0.5,1200,-2000,1000);
TH2F* hE1vsBar = new TH2F("hE1vsBar","Energy1 vs Bars",num_bars,0.5,num_bars+0.5,2000,0,2000);
TH2F* hE2vsBar = new TH2F("hE2vsBar","Energy2 vs Bars",num_bars,0.5,num_bars+0.5,2000,0,2000);
TH2F* hT3vsBar = new TH2F("hT3vsBar","Time3 vs Bars",num_bars,0.5,num_bars+0.5,80,-10,1600);
TH2F* hT5vsBar = new TH2F("hT5vsBar","Time5 vs Bars",num_bars,0.5,num_bars+0.5,1000,-10,40);
TH2F* hT4vsBar = new TH2F("hT4vsBar","Time4 vs Bars",num_bars,0.5,num_bars+0.5,80,-10,1600);
TH2F* hT6vsBar = new TH2F("hT6vsBar","Time6 vs Bars",num_bars,0.5,num_bars+0.5,1000,-10,40);

TH2F* hTofvsEhitall = new TH2F("hTofvsEhitall","Tof vs Ehit",1000,0.,500.,3000,0.,300.);
TH2F* hTofvsEhitallVeto = new TH2F("hTofvsEhitallVeto","Tof vs Ehit Veto",1000,0.,500.,3000,0.,300.);
TH2F* hTofvsEhitallVetoq = new TH2F("hTofvsEhitallVetoq","Tof vs Ehit Vetoq",1000,0.,500.,3000,0.,300.);
TH2F* hTofvsErawall = new TH2F("hTofvsErawall","Tof vs Eraw",300,0.,1600.,100,30.,44.);
TH2F* hTofvsErawallcorr = new TH2F("hTofvsErawallcorr","Tof vs Eraw",300,0.,1600.,400,20.,60.);
TH1F* hTof[num_planes*2+1];
TH1F* hTofc[num_planes*2+1];
TH1F* hTofch[num_planes*2+1];
TH2F* hTofvsEhit[num_planes*2+1];
TH2F* hTofvsEhitv[num_planes*2+1];
TH2F* hTofvsEhitnv[num_planes*2+1];
TH2F* hTofvsEhitvq[num_planes*2+1];
TH2F* hTofvsBar = new TH2F("hTofvsBar","Tof vs Bars",num_bars,0.5,num_bars+0.5,3000,0,300);
TH2F* hTofvsBarVeto = new TH2F("hTofvsBarVeto","Tof vs Bars veto",num_bars,0.5,num_bars+0.5,3000,0,300);
TH2F* hTofvsBarVetoq = new TH2F("hTofvsBarVetoq","Tof vs Bars vetoq",num_bars,0.5,num_bars+0.5,3000,0,300);
TH2F* hTofcvsBar = new TH2F("hTofcvsBar","Tof vs Bars",num_bars,0.5,num_bars+0.5,3000,0,300);
TH2F* hTofvsZall = new TH2F("hTofvsZall","Tof vs Z beam",200,0.,40.,3000,0,300);
TH2F* hTofvsX = new TH2F("hTofvsX","Tof vs X beam",200,-150.,150.,3000,0,300);
TH2F* hTofvsY = new TH2F("hTofvsY","Tof vs Y beam",200,-150.,150.,3000,0,300);
TH2F* hTofvsPath = new TH2F("hTofvsPath","Tof vs Path beam",2000,820.,900.,2000,0,200);
TH2F* hVelvsBar = new TH2F("hVelvsBar","Velocity vs Bars",num_bars,0.5,num_bars+0.5,2000,0,50);
TH2F* hVelvsBarVeto = new TH2F("hVelvsBarVeto","Velocity vs Bars veto",num_bars,0.5,num_bars+0.5,2000,0,50);
TH2F* hVelvsBarVetoq = new TH2F("hVelvsBarVetoq","Velocity vs Bars vetoq",num_bars,0.5,num_bars+0.5,2000,0,50);
TH2F* hTofvsRun = new TH2F("hTofvsRun","Tof vs run ",2000,2000.5,4000.5,3000,0.,300.);

TH2F* hVTofvsEhitall = new TH2F("hVTofvsEhitall","Veto Tof vs Ehitall",500,0.,50.,1000,0.,250.);
TH2F* hVTofvsEhitStop = new TH2F("hVTofvsEhitStop","Veto Tof vs Ehit stopped",500,0.,50.,1000,0.,250.);
TH2F* hVTofvsEhitMatch = new TH2F("hVTofvsEhitMatch","Veto Tof vs Ehit match",500,0.,50.,1000,0.,250.);
TH2F* hVTofvsEhitMatch2 = new TH2F("hVTofvsEhitMatch2","Veto Tof vs Ehit match2",500,0.,50,1000,0.,250.);
TH2F* hVTofvsEhitNoMatch = new TH2F("hVTofvsEhitNoMatch","Veto Tof vs Ehit no match",500,0.,50.,1000,0.,250.);
TH2F* hVTofvsEhitMulti = new TH2F("hVTofvsEhitMulti","Veto Tof vs Ehit mult",500,0.,50.,1000,0.,250.);
TH2F* hVEhitvsEhit = new TH2F("hVEhitvsEhit","Veto E vs E",500,0.,1000,500,0.,50.);
TH2F* hVEhitvsEhitch = new TH2F("hVEhitvsEhitch","Veto E vs E charged",500,0.,1000,500,0.,50.);
TH2F* hVTofvsBar = new TH2F("hVTofvsBar","VTof vs Vbar",8,0.5,8.5,1000,0.,250.);
TH2F* hVTofvsBargamma = new TH2F("hVTofvsBargamma","VTof vs Vbar gamma",8,0.5,8.5,1000,0.,250.);
TH2F* hVTofvsX = new TH2F("hVTofvsX","VTof vs X",500,-170,170,2000,0.,150.);
TH2F* hVTofvsXneigh = new TH2F("hVTofvsXneigh","VTof vs X neighbours",500,-170,170,2000,0.,150.);
TH2F* hVTofvsXgamma = new TH2F("hVTofvsXgamma","VTof vs X gamma",500,-170,170,2000,0.,150.);

TH2F* hVT1vsVT2 = new TH2F("hVT1vsVT2","VT1 vs VT2",500,-250,250,500,-250,250);

TH2F* hVXvsX[9];
TH2F* hVXvsXax[9];
TH2F* hVXvsXcut[9];
TH2F* hVXvsXch[9];
TH2F* hVTdiffvsY[9];
TH2F* hVTdiffvsYneigh[9];
TH2F* hVTdiffvsYcut[9];
TH2F* hVTdiffvsYch[9];

TH2F* hVTofvsEhit[9];

TH2F* hDYvsY = new TH2F("hDYvsY","Veto DY vs Y",1000,-170.,170.,500,-50.,50.);
TH2F* hDXvsX = new TH2F("hDXvsX","Veto DX vs X",1000,-170.,170.,50,-50.,50.);
TH2F* hDYvsZ = new TH2F("hDYvsZ","Veto DY vs Z",200,0.,40.,500,-50.,50.);
TH2F* hDXvsZ = new TH2F("hDXvsZ","Veto Dx vs Z",200,0.,40.,50,-50.,50.);

TH2F* hdxvsX = new TH2F("hdxvsX","Tpc dx vs X",1000,-170.,170.,1000,-1000.,1000.);
TH2F* hdyvsY = new TH2F("hdyvsY","Tpc dy vs Y",1000,-170.,170.,1000,-500.,500.);
TH2F* hdxvsXn = new TH2F("hdxvsXn","Tpc dx vs X n",1000,-170.,170.,1000,-1000.,1000.);
TH2F* hdyvsYn = new TH2F("hdyvsYn","Tpc dy vs Y n",1000,-170.,170.,1000,-500.,500.);
TH2F* hdxvsXc = new TH2F("hdxvsXc","Tpc dx vs X c",1000,-170.,170.,1000,-1000.,1000.);
TH2F* hdyvsYc = new TH2F("hdyvsYc","Tpc dy vs Y c",1000,-170.,170.,1000,-500.,500.);

TH2F* hdyvsdx = new TH2F("hdyvsdx","Tpc dy vs dx",1000,-1.2,1.2,1000,-1.2,1.2);
TH2F* hdyvsdz = new TH2F("hdyvsdz","Tpc dy vs dz",1000,-1.2,1.2,1000,-1.2,1.2);
TH2F* hdxvsdz = new TH2F("hdxvsdz","Tpc dx vs dz",1000,-1.2,1.2,1000,-1.2,1.2);
TH2F* hdyvsdxc = new TH2F("hdyvsdxc","Tpc dy vs dx c",1000,-1.2,1.2,1000,-1.2,1.2);
TH2F* hdxvsdzc = new TH2F("hdxvsdzc","Tpc dx vs dz",1000,-1.2,1.2,1000,-1.2,1.2);
TH2F* hdyvsdzc = new TH2F("hdyvsdzc","Tpc dy vs dz",1000,-1.2,1.2,1000,-1.2,1.2);
TH2F* hdyvsdxn = new TH2F("hdyvsdxn","Tpc dy vs dx n",1000,-1.2,1.2,1000,-1.2,1.2);
TH2F* hdyvsdzn = new TH2F("hdyvsdzn","Tpc dy vs dz n",1000,-1.2,1.2,1000,-1.2,1.2);
TH2F* hdxvsdzn = new TH2F("hdxvsdzn","Tpc dx vs dz n",1000,-1.2,1.2,1000,-1.2,1.2);

TH2F* hdxvsMomc = new TH2F("hdxvsmomc","Tpc dx vs mom c",1000,0.,6000.,1000,-1.2,1.2);
TH2F* hdyvsMomc = new TH2F("hdyvsmomc","Tpc dy vs mom c",1000,0.,6000.,1000,-1.2,1.2);
TH2F* hdzvsMomc = new TH2F("hdzvsmomc","Tpc dz vs mom c",1000,0.,6000.,1000,-1.2,1.2);

TH2F* hYvsXTpc=new TH2F("hYvsXTpc","Y vs X tpc",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hYvsZTpc=new TH2F("hYvsZTpc","Y vs Z tpc",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hXvsZTpc=new TH2F("hXvsZTpc","X vs Z tpc",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hYvsXTpcc=new TH2F("hYvsXTpcc","Y vs X tpc c",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hYvsZTpcc=new TH2F("hYvsZTpcc","Y vs Z tpc c",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hXvsZTpcc=new TH2F("hXvsZTpcc","X vs Z tpc c",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hYvsXTpcn=new TH2F("hYvsXTpcn","Y vs X tpc n",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hYvsZTpcn=new TH2F("hYvsZTpcn","Y vs Z tpc n",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hXvsZTpcn=new TH2F("hXvsZTpcn","X vs Z tpc n",1000,-3000.,3000.,1000,-3000.,3000.);

TH2F* hXvsZTpccz1=new TH2F("hXvsZTpccz1","X vs Z tpc z1",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hXvsZTpccp=new TH2F("hXvsZTpccp","X vs Z tpc p",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hXvsZTpccd=new TH2F("hXvsZTpccd","X vs Z tpc d",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hXvsZTpcct=new TH2F("hXvsZTpcct","X vs Z tpc t",1000,-3000.,3000.,1000,-3000.,3000.);
TH2F* hXvsZTpcche=new TH2F("hXvsZTpcche","X vs Z tpc he",1000,-3000.,3000.,1000,-3000.,3000.);

TH1F* hDirTpc=new TH1F("hDirTpc","direction tpc",720,-120.,120.);
TH1F* hDirTpcc=new TH1F("hDirTpcc","direction tpc c",720,-120.,120.);
TH1F* hDirTpcn=new TH1F("hDirTpcn","direction tpc n",720,-120.,120.);

TH2F* hdxvsdist = new TH2F("hdxvsdist","Tpc dx vs dist",1000,-10.,100.,1000,-3000.,3000.);
TH2F* hdyvsdist = new TH2F("hdyvsdist","Tpc dy vs dist",1000,-10.,100.,1000,-3000.,3000.);

TH2F* hTofvsMom = new TH2F("hTofvsMom","Tof vs Mom",1000,-3000.,3000.,1000,0.,200.);
TH2F* hTofvsMomc = new TH2F("hTofvsMomc","Tof vs Mom c",1000,-3000.,3000.,1000,0.,200.);

TH2F* hVTofvsMomc = new TH2F("hVTofvsMomc","VTof vs Mom c",1000,-3000.,3000.,1000,0.,200.);
TH2F* hDedxvsVehit = new TH2F("hDedxvsVehit","dedx vs n V etot",1000,0.,50.,1000,0.,1000.);

TH2F* hDedxvsMom = new TH2F("hDedxvsMom","dedx vs Mom",1000,-3000.,3000.,1000,0.,1000.);
TH2F* hDedxvsMomc = new TH2F("hDedxvsMomc","dedx vs Mom c",1000,-3000.,3000.,1000,0.,1000.);
TH2F* hDedxvsMomn = new TH2F("hDedxvsMomn","dedx vs Mom n",1000,-3000.,3000.,1000,0.,1000.);
TH2F* hDedxvsEhitClust = new TH2F("hDedxvsEhitClust","dedx vs n etot",1000,0.,1000.,1000,0.,1000.);
TH2F* hDedxvsEhitClustc = new TH2F("hDedxvsEhitClustc","dedx vs n etot c",1000,0.,1000.,1000,0.,1000.);

TH2F* hDedxvsMomcz1 = new TH2F("hDedxvsMomcz1","dedx vs Mom c z1",1000,-3000.,3000.,1000,0.,1000.);
TH2F* hDedxvsMomcp = new TH2F("hDedxvsMomcp","dedx vs Mom c p",1000,-3000.,3000.,1000,0.,1000.);
TH2F* hDedxvsMomcd = new TH2F("hDedxvsMomcd","dedx vs Mom c d",1000,-3000.,3000.,1000,0.,1000.);
TH2F* hDedxvsMomct = new TH2F("hDedxvsMomct","dedx vs Mom c t",1000,-3000.,3000.,1000,0.,1000.);
TH2F* hDedxvsMomche = new TH2F("hDedxvsMomche","dedx vs Mom c he",1000,-3000.,3000.,1000,0.,1000.);

TH2F* hXYtpc = new TH2F("hXYtpc","Y vs X tpc tracks",680,-170,170,680,-170,170); 
TH2F* hXYtpcc = new TH2F("hXYtpcc","Y vs X tpc tracks c",680,-170,170,680,-170,170); 

TH2F* hNmultvsTrmult=new TH2F("hNmultvsTrmult","NMultiplicity vs tpc mult",150,-0.5,149.5,40,-0.5,39.5);
TH2F* hCmultvsTrmult=new TH2F("hCmultvsTrmult","CMultiplicity vs tpc mult",150,-0.5,149.5,40,-0.5,39.5);
TH2F* hVmultvsTrmult=new TH2F("hVmultvsTrmult","VMultiplicity vs tpc mult",150,-0.5,149.5,40,-0.5,39.5);

TH1F* hTrmult=new TH1F("hTrmult","tpc mult",150,-0.5,149.5);

TH2F* hThetavsPhiTpc=new TH2F("hThetavsPhiTpc","theta vs phi tpc",360,0.,360,180,0.,180.);
TH2F* hThetavsPhiTpcn=new TH2F("hThetavsPhiTpcn","theta vs phi tpc n",360,0.,360,180,0.,180.);

TH2F* hVzvsPhiTpc=new TH2F("hVzvsPhiTpc","Vz vs phi tpc",360,0.,360,400,-600.,1200.);
TH2F* hVzvsPhiTpcc=new TH2F("hVzvsPhiTpcc","Vz vs phi tpc c",360,0.,360,400,-600.,1200.);

TH2F* hDYvsDX = new TH2F("hDYvsDX","DY vs DX",400,-100.,100.,400,-100.,100.);
TH2F* hDYvsDXcut = new TH2F("hDYvsDXcut","DY vs DX neutron",400,-100.,100.,400,-100.,100.);
TH2F* hDYvsDXacut = new TH2F("hDYvsDXacut","DY vs DX charged",400,-100.,100.,400,-100.,100.);
TH2F* hDYvsDXcorr = new TH2F("hDYvsDXcorr","DY vs DX",400,-100.,100.,400,-100.,100.);
TH2F* hDYvsVdoubles = new TH2F("hDYvsVdoubles","DY vs Vdoubles",10,-0.5,9.5,400,-100.,100.);

TH2F* hTofdiffvsBar = new TH2F("hTofdiffvsBar","Tn - Tv vs bar neutrons",400,0.5,400.5,1000,-30.,30.);
TH2F* hTofdiffvsVbar = new TH2F("hTofdiffvsVbar","Tn - Tv vs Vbar neutrons",8,0.5,8.5,1000,-30.,30.);
TH2F* hTofdiffvsBarch = new TH2F("hTofdiffvsBarch","Tn - Tv vs bar charged",400,0.5,400.5,1000,-30.,30.);
TH2F* hTofdiffvsVbarch = new TH2F("hTofdiffvsVbarch","Tn - Tv vs Vbar charged",8,0.5,8.5,1000,-30.,30.);
TH2F* hVetoEraw = new TH2F("hVetoEraw", "Veto Eraw",8,0.5,8.5,500,0,6000);
TH2F* hVetoEcal = new TH2F("hVetoEcal", "Veto Ecal",8,0.5,8.5,500,0,50);
TH2F* hVetoEcal1 = new TH2F("hVetoEcal1", "Veto Ecal1",8,0.5,8.5,500,0,50);
TH2F* hVetoEcal2 = new TH2F("hVetoEcal2", "Veto Ecal2",8,0.5,8.5,500,0,50);
TH2F* hVetoTcal1 = new TH2F("hVetoTcal1", "Veto Tcal1",8,0.5,8.5,500,0,250);
TH2F* hVetoTcal2 = new TH2F("hVetoTcal2", "Veto Tcal2",8,0.5,8.5,500,0,250);

TH2F* hVTdiffEcal1[9];
TH2F* hVTdiffEcal2[9];
TH2F* hTofdiffcorrvsEvntBar[9];

TH2F* hEvsBarVeto = new TH2F("hEvsBarVeto","Energy vs Bars Veto",num_bars,0.5,num_bars+0.5,2000,0,200);
TH2F* hEvsBarVetoq = new TH2F("hEvsBarVetoq","Energy vs Bars Vetoq",num_bars,0.5,num_bars+0.5,2000,0,200);
TH2F* hVetovsRun = new TH2F("hVetovsRun","Veto Edep vs run ",2000,2000.5,4000.5,3000,0.,50.);
TH2F* hVTofvsRun = new TH2F("hVTofvsRun","Veto Tof vs run ",2000,2000.5,4000.5,3000,0.,250.);

TH2F* hTdiffvsBar = new TH2F("hTdiffvsBar","Tdiff vs Bars",num_bars,0.5,num_bars+0.5,8000,-400,400);
TH2F* hEdiffvsBar = new TH2F("hEdiffvsBar","Ediff vs Bars",num_bars,0.5,num_bars+0.5,1000,-10.,10.);

TH2F* htt1_tcal1 = new TH2F("htt1_tcal1","tt1 - tcal1",num_bars,0.5,num_bars+0.5,1000,-10,10);
TH2F* htt2_tcal2 = new TH2F("htt2_tcal2","tt2 - tcal2",num_bars,0.5,num_bars+0.5,1000,-10,10);

TH2F* hE1vsPos[num_bars+1];
TH2F* hE2vsPos[num_bars+1];

TH1F* hClusterMult=new TH1F("hClusterMult","Cluster multiplicity",40,-0.5,39.5);
TH1F* hClusterMultAfter=new TH1F("hClusterMultAfter","Cluster multiplicity after",40,-0.5,39.5);
TH1F* hClusterDim=new TH1F("hClusterDim","Cluster dimension",50,0.,50.);
TH2F* hTofvsEhitClust = new TH2F("hTofvsEhitClust","Tof vs Ehit clusters",500,0.,1000.,3000,0.,300.);
TH2F* hTofvsEhitClustBack = new TH2F("hTofvsEhitClustBack","Tof vs Ehit clusters backscattering",500,0.,1000.,3000,0.,300.);
TH2F* hFirstZ = new TH2F("hFirstZ","Tof vs Z clusters",500,0.,40.,3000,0.,300.);
TH2F* hTofvsEhitClustVeto = new TH2F("hTofvsEhitClustVeto","Tof vs Ehit neutrons",500,0.,1000.,3000,0.,300.);
TH2F* hTofvsEhitClustVetoStrong = new TH2F("hTofvsEhitClustVetoStrong","Tof vs Ehit neutrons strong",500,0.,1000.,3000,0.,300.);
TH2F* hTofvsEhitClustNoVeto = new TH2F("hTofvsEhitClustNoVeto","Tof vs Ehit charged",500,0.,1000.,3000,0.,300.);
TH2F* hTofvsEhitClustNoVeto2 = new TH2F("hTofvsEhitClustNoVeto2","Tof vs Ehit charged 2",500,0.,1000.,3000,0.,300.);
TH2F* hTofvsBarClustVeto = new TH2F("hTofvsBarClustVeto","Tof vs Bar clusters",num_bars,0.5,num_bars+0.5,3000,0.,300.);
TH2F* hTofvsBarClustNoVeto = new TH2F("hTofvsBarClustNoVeto","Tof vs Bar charged clusters",num_bars,0.5,num_bars+0.5,3000,0.,300.);
TH2F* hEvsBarClustVeto = new TH2F("hEvsBarClustVeto","E vs Bar neutrons",num_bars,0.5,num_bars+0.5,500,0.,1000.);
TH2F* hEvsBarClustNoVeto = new TH2F("hEvsBarClustNoVeto","E vs Bar charged",num_bars,0.5,num_bars+0.5,500,0.,1000.);
TH2F* hTofvsEhitClustv[num_planes*2+1];
TH2F* hTofvsEhitClustnv[num_planes*2+1];
TH2F* hTofvsEhitClustp[num_planes*2+1];
TH2F* hTofvsEhitClustvstrong[num_planes*2+1];

TH2F* hBaracut = new TH2F("hBaracut","y vs bar acut",400,0.5,400.5,680,-170,170); 

TH2F* hTdiffClust = new TH2F("hTdiffClust","tdiff in cluster",500,-20,20,500,-20,20);

TH1F* hClusterMultNeutrons=new TH1F("hClusterMultNeutrons","Cluster multiplicity",40,-0.5,39.5);
TH1F* hClusterMultCharged=new TH1F("hClusterMultCharged","Cluster multiplicity",40,-0.5,39.5);
TH2F* hClNovsClNo=new TH2F("hClNovsClNo","Cluster mult veto vs neuland",40,-0.5,39.5,40,-0.5,39.5);

TH2F* hClvsTotEne=new TH2F("hClvsTotEne","Cluster multiplicity vs Total energy",500,0,500,40,-0.5,39.5);

TH2F* hEvent[20000];
TH2F* hEvntXY[20000];
TH2F* hEvntXZ[20000];
TH2F* hEvntYZ[20000];
TH1F* hClnovsEvnt = new TH1F("hClnovsEvnt","Clno vs Evnt",20000,.5,20000.5);

TH2F* hDax = new TH2F("hDax","ax difference vs clust dim",50,0.5,50.5,1000,-200,200);
TH2F* hDay = new TH2F("hDay","ay difference vs clust dim",50,0.5,50.5,1000,-200,200);
TH2F* hDayvsDax = new TH2F("hDayvsDax","ay difference vs ax difference",1000,-4,4,1000,-4,4);
TH2F* haxvsax = new TH2F("haxvsax","ax vs ax",1000,-200,200,1000,-200,200);
TH2F* hayvsax = new TH2F("hayvsax","ay vs ax",1000,-200,200,1000,-200,200);
TH2F* hDirdiff = new TH2F("hDirdiff","direction difference",500,-0.2,0.2,500,-0.2,0.2);
TH2F* hDirdiff2 = new TH2F("hDirdiff2","direction difference after",500,-0.2,0.2,500,-0.2,0.2);

TH2F* haycorr = new TH2F("haycorr","ay difference",1000,-200.,200.,1000,-200,200);
TH2F* haycorr1 = new TH2F("haycorr1","ay difference",1000,-200.,200.,1000,-200,200);
TH2F* haycorr2 = new TH2F("haycorr2","ay difference",1000,-200.,200.,1000,-200,200);
TH2F* hayvsay = new TH2F("hayvsay","ay difference",200,-200.,200.,200,-200,200);

TH2F* hclxz = new TH2F("hclxz","xz",400,-200.,200.,400,-200,200);
TH2F* hclyz = new TH2F("hclyz","yz",1000,-10.,20.,1000,-20,50);

TH2F* hClustCorrel_xz = new TH2F("hClustCorrel_xz","clustcorrel xz",500,-50.,50.,500,-200,200);
TH2F* hClustCorrel_yz = new TH2F("hClustCorrel_yz","clustcorrel yz",500,-50.,50.,500,-200,200);
TH2F* hClustCorrel_tz = new TH2F("hClustCorrel_tz","clustcorrel tz",500,-50.,50.,500,-20,20);
TH2F* hClustCorrel_xx = new TH2F("hClustCorrel_xx","clustcorrel xx",500,-150.,150.,500,-200,200);
TH2F* hClustCorrel_yy = new TH2F("hClustCorrel_yy","clustcorrel yy",500,-150.,150.,500,-200,200);
TH2F* hClustCorrel_tt = new TH2F("hClustCorrel_tt","clustcorrel tt",500,0.,150.,500,-20,20);

TH1F* hTofdiff=new TH1F("hTofdiff","Tn - Tv",4000,-200,200);
TH1F* hTofdiffch=new TH1F("hTofdiffn","Tn - Tv ch",4000,-200,200);
TH2F* hTofdiffvsEcal=new TH2F("hTofdiffvsEcal","Tn - Tv",500,0,50,1000,-200,200);
TH2F* hTofdiffvsEhit=new TH2F("hTofdiffvsEhit","Tn - Tv",1000,0,500,1000,-200,200);
TH2F* hVTofvsTof=new TH2F("hVTofvsTof","VTof vs Tof",1000,0,300,500,0,250);
TH2F* hTofdiffvsEcalch=new TH2F("hTofdiffvsEcalch","Tn - Tv vs Ecal ch",500,0,50,1000,-200,200);
TH2F* hTofdiffvsEhitch=new TH2F("hTofdiffvsEhitch","Tn - Tv vs Ehit ch",1000,0,500,1000,-200,200);

TH2F* hVTofvsTofch=new TH2F("hVTofvsTofch","VTof vs Tof ch",1000,0,300,500,0,250);

TH2F* hVVelvsVel = new TH2F("hVVelvsVel", "VVel vs Vel",500,0,50,500,0,50);
TH2F* hVVelvsVelch = new TH2F("hVVelvsVelch", "VVel vs Vel ch",500,0,50,500,0,50);

TH2F* hPlane[2*num_planes+1];
TH2F* hPlaneBack[2*num_planes+1];
TH2F* hPlanecut[2*num_planes+1];
TH2F* hPlaneacut[2*num_planes+1];

TH2F* hXY[2*num_planes+1];
TH2F* hXYcut[2*num_planes+1];
TH2F* hXYacut[2*num_planes+1];

TH2F* hXZclust = new TH2F("hXZclust","XY hits cluster",8,0.,40.,170,-170,170); 
TH2F* hYZclust = new TH2F("hYZclust","XY hits cluster",8,0.,40.,170,-170,170); 

TH1F* hTstart = new TH1F("hTstart","Time start",2000,-42000,42000);

void loop(TChain *fChain, TChain *fChainx)
{
  char name1[60]; 
  char name2[60];

  for(int i=0;i<20000;i++){
      sprintf(name1,"hEvent%d",i);
      sprintf(name2,"X-Y X-Z Y-Z %d",i);     
      hEvent[i] = new TH2F(name1,name2, 73,0.5,73.5,73,0.5,73.5); 

      sprintf(name1,"hEvntXY%d",i);
      sprintf(name2,"Y-X %d",i);     
      hEvntXY[i] = new TH2F(name1,name2, 216,-135.,135.,216,-135.,135.); 
      sprintf(name1,"hEvntXZ%d",i);
      sprintf(name2,"Z-X %d",i);     
      hEvntXZ[i] = new TH2F(name1,name2, 216,-135.,135.,40,-5.,45.); 
      sprintf(name1,"hEvntYZ%d",i);
      sprintf(name2,"Y-Z %d",i);     
      hEvntYZ[i] = new TH2F(name1,name2, 40,-5.,45.,216,-135.,135.); 
  }

  for(int i=1;i<=num_planes*2;i++){
      sprintf(name1,"hTof%d",i);
      sprintf(name2,"Tof%d",i);     
      hTof[i] = new TH1F(name1,name2, 8000,-200.,200.); 
      sprintf(name1,"hTofc%d",i);
      sprintf(name2,"Tofc%d",i);     
      hTofc[i] = new TH1F(name1,name2, 8000,-200.,200.); 
      sprintf(name1,"hTofch%d",i);
      sprintf(name2,"Tofch%d",i);     
      hTofch[i] = new TH1F(name1,name2, 8000,-200.,200.); 
  }
  for(int i=1;i<=num_planes*2;i++){
      sprintf(name1,"hTofvsEhit%d",i);
      sprintf(name2,"Tof vs Ehit %d",i);     
      hTofvsEhit[i] = new TH2F(name1,name2, 200,0.,EhitRange,3000,0.,300.); 
  }
  for(int i=1;i<=num_planes*2;i++){
      sprintf(name1,"hTofvsEhitv%d",i);
      sprintf(name2,"Tof vs Ehit %d",i);     
      hTofvsEhitv[i] = new TH2F(name1,name2, 200,0.,EhitRange,3000,0.,300.); 
      sprintf(name1,"hTofvsEhitnv%d",i);
      sprintf(name2,"Tof vs Ehit charged %d",i);     
      hTofvsEhitnv[i] = new TH2F(name1,name2, 200,0.,EhitRange,3000,0.,300.); 
  }
  for(int i=1;i<=num_planes*2;i++){
      sprintf(name1,"hTofvsEhitvq%d",i);
      sprintf(name2,"Tof vs Ehit %d",i);     
      hTofvsEhitvq[i] = new TH2F(name1,name2, 200,0.,EhitRange,3000,0.,300.); 
  }
  for(int i=1;i<=num_planes*2;i++){
      sprintf(name1,"hTofvsEhitClustp%d",i);
      sprintf(name2,"Tof vs Ehit clusters %d",i);     
      hTofvsEhitClustp[i] = new TH2F(name1,name2, 200,0.,EhitRange,3000,0.,300.); 
      sprintf(name1,"hTofvsEhitClustv%d",i);
      sprintf(name2,"Tof vs Ehit clusters %d",i);     
      hTofvsEhitClustv[i] = new TH2F(name1,name2, 200,0.,EhitRange,3000,0.,300.); 
      sprintf(name1,"hTofvsEhitClustnv%d",i);
      sprintf(name2,"Tof vs Ehit charged clusters %d",i);     
      hTofvsEhitClustnv[i] = new TH2F(name1,name2, 200,0.,EhitRange,3000,0.,300.); 
      sprintf(name1,"hTofvsEhitClustvstrong%d",i);
      sprintf(name2,"Tof vs Ehit clusters veto strong %d",i);     
      hTofvsEhitClustvstrong[i] = new TH2F(name1,name2, 200,0.,EhitRange,3000,0.,300.); 
  }
  for(int i=1;i<=2*num_planes;i++){
      sprintf(name1,"XY%i",i);
      sprintf(name2,"Xhit vs Yhit  %i",i);     
      hXY[i] = new TH2F(name1,name2,680,-170,170,680,-170,170); 
      sprintf(name1,"XYcut%i",i);
      sprintf(name2,"Xhit vs Yhit cut  %i",i);     
      hXYcut[i] = new TH2F(name1,name2,680,-170,170,680,-170,170); 
      sprintf(name1,"XYacut%i",i);
      sprintf(name2,"Xhit vs Yhit acut  %i",i);     
      hXYacut[i] = new TH2F(name1,name2,680,-170,170,680,-170,170); 
      sprintf(name1,"plane%i",i);
      sprintf(name2,"Xhit vs Yhit clust %i",i);     
      hPlane[i] = new TH2F(name1,name2,680,-170,170,680,-170,170); 
      sprintf(name1,"planeback%i",i);
      sprintf(name2,"Xhit vs Yhit clust %i",i);     
      hPlaneBack[i] = new TH2F(name1,name2,170,-170,170,170,-170,170); 
  }
  for(int i=1;i<=2*num_planes;i++){
      sprintf(name1,"plane%icut",i);
      sprintf(name2,"Xhit vs Yhit clust %i",i);     
      hPlanecut[i] = new TH2F(name1,name2,680,-170,170,680,-170,170); 
  }
  for(int i=1;i<=2*num_planes;i++){
      sprintf(name1,"plane%iacut",i);
      sprintf(name2,"Xhit vs Yhit clust %i",i);     
      hPlaneacut[i] = new TH2F(name1,name2,680,-170,170,680,-170,170); 
  }
  for(int i=1;i<=8;i++){
    sprintf(name1,"VTdiffEcal1%i",i);
    sprintf(name2,"VEcal1 vs position %i",i);     
    hVTdiffEcal1[i] = new TH2F(name1,name2,200,-30,30,500,0,50);
    sprintf(name1,"VTdiffEcal2%i",i);
    sprintf(name2,"VEcal2 vs position %i",i);     
    hVTdiffEcal2[i] = new TH2F(name1,name2,200,-30,30,500,0,50);
    sprintf(name1,"VTdiffvsY%i",i);
    sprintf(name2,"VTdiff vs Y %i",i);     
    hVTdiffvsY[i] = new TH2F(name1,name2,500,-170.,170.,500,-100,100); // check overlapping !!!!!!
    sprintf(name1,"VTdiffvsYneigh%i",i);
    sprintf(name2,"VTdiff vs Y neighbours %i",i);     
    hVTdiffvsYneigh[i] = new TH2F(name1,name2,500,-170.,170.,500,-100,100); // check overlapping !!!!!!
    sprintf(name1,"VTdiffvsYcut%i",i);
    sprintf(name2,"VTdiff vs Y neutrons %i",i);     
    hVTdiffvsYcut[i] = new TH2F(name1,name2,500,-170.,170.,500,-100,100);
    sprintf(name1,"VTdiffvsYch%i",i);
    sprintf(name2,"VTdiff vs Y charged %i",i);     
    hVTdiffvsYch[i] = new TH2F(name1,name2,500,-170.,170.,500,-100,100);
    sprintf(name1,"hVTofvsEhit%i",i);
    sprintf(name2,"VTof vs Ehit %i",i);     
    hVTofvsEhit[i] = new TH2F(name1,name2,500,0.,50.,500,0.,250.);
    sprintf(name1,"hVXvsX%i",i);
    sprintf(name2,"VX vs X %i",i);     
    hVXvsX[i] = new TH2F(name1,name2,500,-170.,170.,10,-157.5,157.5);
    sprintf(name1,"hVXvsXax%i",i);
    sprintf(name2,"VX vs X ax %i",i);     
    hVXvsXax[i] = new TH2F(name1,name2,500,-170.,170.,10,-157.5,157.5);
    sprintf(name1,"hVXvsXcut%i",i);
    sprintf(name2,"VX vs X cut %i",i);     
    hVXvsXcut[i] = new TH2F(name1,name2,500,-170.,170.,10,-157.5,157.5);
    sprintf(name1,"hVXvsXch%i",i);
    sprintf(name2,"VX vs X ch %i",i);     
    hVXvsXch[i] = new TH2F(name1,name2,500,-170.,170.,10,-157.5,157.5);
  }

  Int_t iHit,iBar,jBar;

 // Declaration of leaf types
 Int_t           coin_ch;
 Int_t           f3_traw0;
 Int_t           f3_traw1;
 Int_t           f3_qraw0;
 Int_t           f3_qraw1;
 Double_t        f3_time0;
 Double_t        f3_time1;
 Double_t        f3_time;
 Double_t        f3_timediff;
 Int_t           f5_traw0;
 Int_t           f5_traw1;
 Int_t           f5_qraw0;
 Int_t           f5_qraw1;
 Double_t        f5_time0;
 Double_t        f5_time1;
 Double_t        f5_time;
 Double_t        f5_timediff;
 Int_t           f7_traw0;
 Int_t           f7_traw1;
 Int_t           f7_qraw0;
 Int_t           f7_qraw1;
 Double_t        f7_time0;
 Double_t        f7_time1;
 Double_t        f7_time;
 Double_t        f7_timediff;
 Double_t        f13_1time;
 Double_t        f13_2time;
 Double_t        sbt_t0;
 Double_t        sbt1_t0;
 Double_t        sbt2_t0;
 Double_t        sbt2_qaveraw;
 Double_t        sbt1_qaveraw;
 Double_t        t_0slw;
 Double_t        sbt1_tlraw;
 Double_t        sbt1_trraw;
 Double_t        sbt2_tlraw;
 Double_t        sbt2_trraw;
 Int_t           nl_tms1;
 Int_t           nl_tms2;
 Int_t           nl_tms3;
 Int_t           nl_nbar;
 Int_t           nl_bar[400];   //[nl_nbar]
 Double_t        nl_qraw0[400];   //[nl_nbar]
 Double_t        nl_qraw1[400];   //[nl_nbar]
 Double_t        nl_qcal0[400];   //[nl_nbar]
 Double_t        nl_qcal1[400];   //[nl_nbar]
 Int_t           nl_t1[400];   //[nl_nbar]
 Int_t           nl_t2[400];   //[nl_nbar]
 Int_t           nl_t3[400];   //[nl_nbar]
 Int_t           nl_t4[400];   //[nl_nbar]
 Int_t           nl_t5[400];   //[nl_nbar]
 Int_t           nl_t6[400];   //[nl_nbar]
 Double_t        nl_traw0[400];   //[nl_nbar]
 Double_t        nl_traw1[400];   //[nl_nbar]
 Double_t        nl_tcal0[400];   //[nl_nbar]
 Double_t        nl_tcal1[400];   //[nl_nbar]
 Double_t        nl_x[400];   //[nl_nbar]
 Double_t        nl_y[400];   //[nl_nbar]
 Double_t        nl_z[400];   //[nl_nbar]
 Bool_t          nl_fired0[400];   //[nl_nbar]
 Bool_t          nl_fired1[400];   //[nl_nbar]
 Bool_t          nl_fired[400];   //[nl_nbar]
 Int_t           nl_mul;
 Int_t           nl_mulv;
 Int_t           nl_mulv1;
 Int_t           nl_mulv3;
 Int_t           nl_mulv5;
 Int_t           nl_mulv7;
 Int_t           nl_mulh;
 Int_t           nl_mulh0;
 Int_t           nl_mulh2;
 Int_t           nl_mulh4;
 Int_t           nl_mulh6;
 Int_t           nlv_nbar;
 Int_t           nlv_bar[9];   //[nlv_nbar]
 Double_t        nlv_qraw0[9];   //[nlv_nbar]
 Double_t        nlv_qraw1[9];   //[nlv_nbar]
 Double_t        nlv_qcal0[9];   //[nlv_nbar]
 Double_t        nlv_qcal1[9];   //[nlv_nbar]
 Double_t        nlv_traw0[9];   //[nlv_nbar]
 Double_t        nlv_traw1[9];   //[nlv_nbar]
 Double_t        nlv_tcl0[9];   //[nlv_nbar]
 Double_t        nlv_tcl1[9];   //[nlv_nbar]
 Double_t        nlv_tcal0[9];   //[nlv_nbar]
 Double_t        nlv_tcal1[9];   //[nlv_nbar]
 Double_t        nlv_x[9];   //[nlv_nbar]
 Double_t        nlv_y[9];   //[nlv_nbar]
 Double_t        nlv_z[9];   //[nlv_nbar]
 Bool_t          nlv_fired[9];   //[nlv_nbar]
 Int_t           nlv_mul;
 
 // List of branch
 TBranch        *b_coin_ch;   //!
 TBranch        *b_f3_traw0;   //!
 TBranch        *b_f3_traw1;   //!
 TBranch        *b_f3_qraw0;   //!
 TBranch        *b_f3_qraw1;   //!
 TBranch        *b_f3_time0;   //!
 TBranch        *b_f3_time1;   //!
 TBranch        *b_f3_time;   //!
 TBranch        *b_f3_timediff;   //!
 TBranch        *b_f5_traw0;   //!
 TBranch        *b_f5_traw1;   //!
 TBranch        *b_f5_qraw0;   //!
 TBranch        *b_f5_qraw1;   //!
 TBranch        *b_f5_time0;   //!
 TBranch        *b_f5_time1;   //!
 TBranch        *b_f5_time;   //!
 TBranch        *b_f5_timediff;   //!
 TBranch        *b_f7_traw0;   //!
 TBranch        *b_f7_traw1;   //!
 TBranch        *b_f7_qraw0;   //!
 TBranch        *b_f7_qraw1;   //!
 TBranch        *b_f7_time0;   //!
 TBranch        *b_f7_time1;   //!
 TBranch        *b_f7_time;   //!
 TBranch        *b_f7_timediff;   //!
 TBranch        *b_f13_1time;   //!
 TBranch        *b_f13_2time;   //!
 TBranch        *b_sbt_t0;   //!
 TBranch        *b_sbt1_t0;   //!
 TBranch        *b_sbt2_t0;   //!
 TBranch        *b_sbt2_qaveraw;   //!
 TBranch        *b_sbt1_qaveraw;   //!
 TBranch        *b_t_0slw;   //!
 TBranch        *b_sbt1_tlraw;   //!
 TBranch        *b_sbt1_trraw;   //!
 TBranch        *b_sbt2_tlraw;   //!
 TBranch        *b_sbt2_trraw;   //!
 TBranch        *b_nl_tms1;   //!
 TBranch        *b_nl_tms2;   //!
 TBranch        *b_nl_tms3;   //!
 TBranch        *b_nl_nbar;   //!
 TBranch        *b_nl_bar;   //!
 TBranch        *b_nl_qraw0;   //!
 TBranch        *b_nl_qraw1;   //!
 TBranch        *b_nl_qcal0;   //!
 TBranch        *b_nl_qcal1;   //!
 TBranch        *b_nl_t1;   //!
 TBranch        *b_nl_t2;   //!
 TBranch        *b_nl_t3;   //!
 TBranch        *b_nl_t4;   //!
 TBranch        *b_nl_t5;   //!
 TBranch        *b_nl_t6;   //!
 TBranch        *b_nl_traw0;   //!
 TBranch        *b_nl_traw1;   //!
 TBranch        *b_nl_tcal0;   //!
 TBranch        *b_nl_tcal1;   //!
 TBranch        *b_nl_x;   //!
 TBranch        *b_nl_y;   //!
 TBranch        *b_nl_z;   //!
 TBranch        *b_nl_fired0;   //!
 TBranch        *b_nl_fired1;   //!
 TBranch        *b_nl_fired;   //!
 TBranch        *b_nl_mul;   //!
 TBranch        *b_nl_mulv;   //!
 TBranch        *b_nl_mulv1;   //!
 TBranch        *b_nl_mulv3;   //!
 TBranch        *b_nl_mulv5;   //!
 TBranch        *b_nl_mulv7;   //!
 TBranch        *b_nl_mulh;   //!
 TBranch        *b_nl_mulh0;   //!
 TBranch        *b_nl_mulh2;   //!
 TBranch        *b_nl_mulh4;   //!
 TBranch        *b_nl_mulh6;   //!
 TBranch        *b_nlv_nbar;   //!
 TBranch        *b_nlv_bar;   //!
 TBranch        *b_nlv_qraw0;   //!
 TBranch        *b_nlv_qraw1;   //!
 TBranch        *b_nlv_qcal0;   //!
 TBranch        *b_nlv_qcal1;   //!
 TBranch        *b_nlv_traw0;   //!
 TBranch        *b_nlv_traw1;   //!
 TBranch        *b_nlv_tcl0;   //!
 TBranch        *b_nlv_tcl1;   //!
 TBranch        *b_nlv_tcal0;   //!
 TBranch        *b_nlv_tcal1;   //!
 TBranch        *b_nlv_x;   //!
 TBranch        *b_nlv_y;   //!
 TBranch        *b_nlv_z;   //!
 TBranch        *b_nlv_fired;   //!
 TBranch        *b_nlv_mul;   //!
 
 fChain->SetBranchAddress("coin_ch", &coin_ch, &b_coin_ch);
 fChain->SetBranchAddress("f3_traw0", &f3_traw0, &b_f3_traw0);
 fChain->SetBranchAddress("f3_traw1", &f3_traw1, &b_f3_traw1);
 fChain->SetBranchAddress("f3_qraw0", &f3_qraw0, &b_f3_qraw0);
 fChain->SetBranchAddress("f3_qraw1", &f3_qraw1, &b_f3_qraw1);
 fChain->SetBranchAddress("f3_time0", &f3_time0, &b_f3_time0);
 fChain->SetBranchAddress("f3_time1", &f3_time1, &b_f3_time1);
 fChain->SetBranchAddress("f3_time", &f3_time, &b_f3_time);
 fChain->SetBranchAddress("f3_timediff", &f3_timediff, &b_f3_timediff);
 fChain->SetBranchAddress("f5_traw0", &f5_traw0, &b_f5_traw0);
 fChain->SetBranchAddress("f5_traw1", &f5_traw1, &b_f5_traw1);
 fChain->SetBranchAddress("f5_qraw0", &f5_qraw0, &b_f5_qraw0);
 fChain->SetBranchAddress("f5_qraw1", &f5_qraw1, &b_f5_qraw1);
 fChain->SetBranchAddress("f5_time0", &f5_time0, &b_f5_time0);
 fChain->SetBranchAddress("f5_time1", &f5_time1, &b_f5_time1);
 fChain->SetBranchAddress("f5_time", &f5_time, &b_f5_time);
 fChain->SetBranchAddress("f5_timediff", &f5_timediff, &b_f5_timediff);
 fChain->SetBranchAddress("f7_traw0", &f7_traw0, &b_f7_traw0);
 fChain->SetBranchAddress("f7_traw1", &f7_traw1, &b_f7_traw1);
 fChain->SetBranchAddress("f7_qraw0", &f7_qraw0, &b_f7_qraw0);
 fChain->SetBranchAddress("f7_qraw1", &f7_qraw1, &b_f7_qraw1);
 fChain->SetBranchAddress("f7_time0", &f7_time0, &b_f7_time0);
 fChain->SetBranchAddress("f7_time1", &f7_time1, &b_f7_time1);
 fChain->SetBranchAddress("f7_time", &f7_time, &b_f7_time);
 fChain->SetBranchAddress("f7_timediff", &f7_timediff, &b_f7_timediff);
 fChain->SetBranchAddress("f13_1time", &f13_1time, &b_f13_1time);
 fChain->SetBranchAddress("f13_2time", &f13_2time, &b_f13_2time);
 fChain->SetBranchAddress("sbt_t0", &sbt_t0, &b_sbt_t0);
 fChain->SetBranchAddress("sbt1_t0", &sbt1_t0, &b_sbt1_t0);
 fChain->SetBranchAddress("sbt2_t0", &sbt2_t0, &b_sbt2_t0);
 fChain->SetBranchAddress("sbt2_qaveraw", &sbt2_qaveraw, &b_sbt2_qaveraw);
 fChain->SetBranchAddress("sbt1_qaveraw", &sbt1_qaveraw, &b_sbt1_qaveraw);
 fChain->SetBranchAddress("t_0slw", &t_0slw, &b_t_0slw);
 fChain->SetBranchAddress("sbt1_tlraw", &sbt1_tlraw, &b_sbt1_tlraw);
 fChain->SetBranchAddress("sbt1_trraw", &sbt1_trraw, &b_sbt1_trraw);
 fChain->SetBranchAddress("sbt2_tlraw", &sbt2_tlraw, &b_sbt2_tlraw);
 fChain->SetBranchAddress("sbt2_trraw", &sbt2_trraw, &b_sbt2_trraw);
 fChain->SetBranchAddress("nl_tms1", &nl_tms1, &b_nl_tms1);
 fChain->SetBranchAddress("nl_tms2", &nl_tms2, &b_nl_tms2);
 fChain->SetBranchAddress("nl_tms3", &nl_tms3, &b_nl_tms3);
 fChain->SetBranchAddress("nl_nbar", &nl_nbar, &b_nl_nbar);
 fChain->SetBranchAddress("nl_bar", nl_bar, &b_nl_bar);
 fChain->SetBranchAddress("nl_qraw0", nl_qraw0, &b_nl_qraw0);
 fChain->SetBranchAddress("nl_qraw1", nl_qraw1, &b_nl_qraw1);
 fChain->SetBranchAddress("nl_qcal0", nl_qcal0, &b_nl_qcal0);
 fChain->SetBranchAddress("nl_qcal1", nl_qcal1, &b_nl_qcal1);
 fChain->SetBranchAddress("nl_traw0", nl_traw0, &b_nl_traw0);
 fChain->SetBranchAddress("nl_traw1", nl_traw1, &b_nl_traw1);
 fChain->SetBranchAddress("nl_tcal0", nl_tcal0, &b_nl_tcal0);
 fChain->SetBranchAddress("nl_tcal1", nl_tcal1, &b_nl_tcal1);
 fChain->SetBranchAddress("nl_x", nl_x, &b_nl_x);
 fChain->SetBranchAddress("nl_y", nl_y, &b_nl_y);
 fChain->SetBranchAddress("nl_z", nl_z, &b_nl_z);
 fChain->SetBranchAddress("nl_fired0", nl_fired0, &b_nl_fired0);
 fChain->SetBranchAddress("nl_fired1", nl_fired1, &b_nl_fired1);
 fChain->SetBranchAddress("nl_fired", nl_fired, &b_nl_fired);
 fChain->SetBranchAddress("nl_mul", &nl_mul, &b_nl_mul);
 fChain->SetBranchAddress("nl_mulv", &nl_mulv, &b_nl_mulv);
 fChain->SetBranchAddress("nl_mulv1", &nl_mulv1, &b_nl_mulv1);
 fChain->SetBranchAddress("nl_mulv3", &nl_mulv3, &b_nl_mulv3);
 fChain->SetBranchAddress("nl_mulv5", &nl_mulv5, &b_nl_mulv5);
 fChain->SetBranchAddress("nl_mulv7", &nl_mulv7, &b_nl_mulv7);
 fChain->SetBranchAddress("nl_mulh", &nl_mulh, &b_nl_mulh);
 fChain->SetBranchAddress("nl_mulh0", &nl_mulh0, &b_nl_mulh0);
 fChain->SetBranchAddress("nl_mulh2", &nl_mulh2, &b_nl_mulh2);
 fChain->SetBranchAddress("nl_mulh4", &nl_mulh4, &b_nl_mulh4);
 fChain->SetBranchAddress("nl_mulh6", &nl_mulh6, &b_nl_mulh6);
 fChain->SetBranchAddress("nlv_nbar", &nlv_nbar, &b_nlv_nbar);
 fChain->SetBranchAddress("nlv_bar", nlv_bar, &b_nlv_bar);
 fChain->SetBranchAddress("nlv_qraw0", nlv_qraw0, &b_nlv_qraw0);
 fChain->SetBranchAddress("nlv_qraw1", nlv_qraw1, &b_nlv_qraw1);
 fChain->SetBranchAddress("nlv_qcal0", nlv_qcal0, &b_nlv_qcal0);
 fChain->SetBranchAddress("nlv_qcal1", nlv_qcal1, &b_nlv_qcal1);
 fChain->SetBranchAddress("nl_t1", nl_t1, &b_nl_t1);
 fChain->SetBranchAddress("nl_t2", nl_t2, &b_nl_t2);
 fChain->SetBranchAddress("nl_t3", nl_t3, &b_nl_t3);
 fChain->SetBranchAddress("nl_t4", nl_t4, &b_nl_t4);
 fChain->SetBranchAddress("nl_t5", nl_t5, &b_nl_t5);
 fChain->SetBranchAddress("nl_t6", nl_t6, &b_nl_t6);
 fChain->SetBranchAddress("nlv_traw0", nlv_traw0, &b_nlv_traw0);
 fChain->SetBranchAddress("nlv_traw1", nlv_traw1, &b_nlv_traw1);
 fChain->SetBranchAddress("nlv_tcl0", nlv_tcl0, &b_nlv_tcl0);
 fChain->SetBranchAddress("nlv_tcl1", nlv_tcl1, &b_nlv_tcl1);
 fChain->SetBranchAddress("nlv_tcal0", nlv_tcal0, &b_nlv_tcal0);
 fChain->SetBranchAddress("nlv_tcal1", nlv_tcal1, &b_nlv_tcal1);
 fChain->SetBranchAddress("nlv_x", nlv_x, &b_nlv_x);
 fChain->SetBranchAddress("nlv_y", nlv_y, &b_nlv_y);
 fChain->SetBranchAddress("nlv_z", nlv_z, &b_nlv_z);
 fChain->SetBranchAddress("nlv_fired", nlv_fired, &b_nlv_fired);
 fChain->SetBranchAddress("nlv_mul", &nlv_mul, &b_nlv_mul);

  //#####################################################################################

 Int_t           run;
 Int_t           eventid;
 Double_t        dedx;
 Double_t        mom;
 Int_t           charge;
 Int_t           pid;
 Int_t           stpid;
 Double_t        ndf;
 Double_t        dx;
 Double_t        dy;
 Double_t        dz;
 Int_t           vid;
 Int_t           parentvid;
 Double_t        vx;
 Double_t        vy;
 Double_t        vz;
 Double_t        pocavx;
 Double_t        pocavy;
 Double_t        pocavz;
 Double_t        tpcdist;
 Bool_t          sigma10;
 Bool_t          sigma15;
 Bool_t          sigma20;
 Bool_t          sigma20z;
 Int_t           tpcmult;
 Double_t        pt;
 Double_t        pzCM;
 Double_t        KECM;
 Double_t        rapL;
 Double_t        rapCM;
 Double_t        rapCMNorm;
 Double_t        projx;
 Double_t        projy;
 Double_t        projz;
 Double_t        phiL;
 Double_t        thetaL;
 Double_t        phiCM;
 Double_t        thetaCM;
 Bool_t          proton;
 Bool_t          pip;
 Bool_t          pim;
 Bool_t          pipbg;
 Double_t        ea;
 
 TBranch        *b_run;   //!
 TBranch        *b_eventid;   //!
 TBranch        *b_dedx;   //!
 TBranch        *b_mom;   //!
 TBranch        *b_charge;   //!
 TBranch        *b_pid;   //!
 TBranch        *b_stpid;   //!
 TBranch        *b_ndf;   //!
 TBranch        *b_dx;   //!
 TBranch        *b_dy;   //!
 TBranch        *b_dz;   //!
 TBranch        *b_vid;   //!
 TBranch        *b_parentvid;   //!
 TBranch        *b_vx;   //!
 TBranch        *b_vy;   //!
 TBranch        *b_vz;   //!
 TBranch        *b_pocavx;   //!
 TBranch        *b_pocavy;   //!
 TBranch        *b_pocavz;   //!
 TBranch        *b_tpcdist;   //!
 TBranch        *b_sigma10;   //!
 TBranch        *b_sigma15;   //!
 TBranch        *b_sigma20;   //!
 TBranch        *b_sigma20z;   //!
 TBranch        *b_tpcmult;   //!
 TBranch        *b_pt;   //!
 TBranch        *b_pzCM;   //!
 TBranch        *b_KECM;   //!
 TBranch        *b_rapL;   //!
 TBranch        *b_rapCM;   //!
 TBranch        *b_rapCMNorm;   //!
 TBranch        *b_projx;   //!
 TBranch        *b_projy;   //!
 TBranch        *b_projz;   //!
 TBranch        *b_phiL;   //!
 TBranch        *b_thetaL;   //!
 TBranch        *b_phiCM;   //!
 TBranch        *b_thetaCM;   //!
 TBranch        *b_proton;   //!
 TBranch        *b_pip;   //!
 TBranch        *b_pim;   //!
 TBranch        *b_pipbg;   //!
 TBranch        *b_ea;   //!
 
 fChainx->SetBranchAddress("run", &run, &b_run);
 fChainx->SetBranchAddress("eventid", &eventid, &b_eventid);
 fChainx->SetBranchAddress("dedx", &dedx, &b_dedx);
 fChainx->SetBranchAddress("mom", &mom, &b_mom);
 fChainx->SetBranchAddress("charge", &charge, &b_charge);
 fChainx->SetBranchAddress("pid", &pid, &b_pid);
 fChainx->SetBranchAddress("stpid", &stpid, &b_stpid);
 fChainx->SetBranchAddress("ndf", &ndf, &b_ndf);
 fChainx->SetBranchAddress("dx", &dx, &b_dx);
 fChainx->SetBranchAddress("dy", &dy, &b_dy);
 fChainx->SetBranchAddress("dz", &dz, &b_dz);
 fChainx->SetBranchAddress("vid", &vid, &b_vid);
 fChainx->SetBranchAddress("parentvid", &parentvid, &b_parentvid);
 fChainx->SetBranchAddress("vx", &vx, &b_vx);
 fChainx->SetBranchAddress("vy", &vy, &b_vy);
 fChainx->SetBranchAddress("vz", &vz, &b_vz);
 fChainx->SetBranchAddress("pocavx", &pocavx, &b_pocavx);
 fChainx->SetBranchAddress("pocavy", &pocavy, &b_pocavy);
 fChainx->SetBranchAddress("pocavz", &pocavz, &b_pocavz);
 fChainx->SetBranchAddress("dist", &tpcdist, &b_tpcdist);
 fChainx->SetBranchAddress("sigma10", &sigma10, &b_sigma10);
 fChainx->SetBranchAddress("sigma15", &sigma15, &b_sigma15);
 fChainx->SetBranchAddress("sigma20", &sigma20, &b_sigma20);
 fChainx->SetBranchAddress("sigma20z", &sigma20z, &b_sigma20z);
 fChainx->SetBranchAddress("mult", &tpcmult, &b_tpcmult);
 fChainx->SetBranchAddress("pt", &pt, &b_pt);
 fChainx->SetBranchAddress("pzCM", &pzCM, &b_pzCM);
 fChainx->SetBranchAddress("KECM", &KECM, &b_KECM);
 fChainx->SetBranchAddress("rapL", &rapL, &b_rapL);
 fChainx->SetBranchAddress("rapCM", &rapCM, &b_rapCM);
 fChainx->SetBranchAddress("rapCMNorm", &rapCMNorm, &b_rapCMNorm);
 fChainx->SetBranchAddress("projx", &projx, &b_projx);
 fChainx->SetBranchAddress("projy", &projy, &b_projy);
 fChainx->SetBranchAddress("projz", &projz, &b_projz);
 fChainx->SetBranchAddress("phiL", &phiL, &b_phiL);
 fChainx->SetBranchAddress("thetaL", &thetaL, &b_thetaL);
 fChainx->SetBranchAddress("phiCM", &phiCM, &b_phiCM);
 fChainx->SetBranchAddress("thetaCM", &thetaCM, &b_thetaCM);
 fChainx->SetBranchAddress("proton", &proton, &b_proton);
 fChainx->SetBranchAddress("pip", &pip, &b_pip);
 fChainx->SetBranchAddress("pim", &pim, &b_pim);
 fChainx->SetBranchAddress("pipbg", &pipbg, &b_pipbg);
 fChainx->SetBranchAddress("ea", &ea, &b_ea);
   
 int  plane;
 
 Double_t t1[num_bars+1], t2[num_bars+1], t3[num_bars+1], t4[num_bars+1], t5[num_bars+1], t6[num_bars+1], e1[num_bars+1], e2[num_bars+1];
 int fired[num_bars+1];
 Double_t pos[num_bars+1];
 Double_t thit[num_bars+1];
 Double_t ehit[num_bars+1];
 Double_t xhit[num_bars+1],yhit[num_bars+1],zhit[num_bars+1];
 Double_t tof[num_bars+1];
 Double_t tofc[num_bars+1];
 Double_t tofnowalk[num_bars+1];

 Int_t clbar[num_bars+1];
 Int_t cl[num_bars+1][num_bars+1];

 Double_t tt1[num_bars+1], tt2[num_bars+1];
 Double_t tt1nowalk[num_bars+1], tt2nowalk[num_bars+1];
 Double_t ee1[num_bars+1], ee2[num_bars+1];
 int mult, multpart, multpmt, multV, multH;

 int veto, vetoq;
 int vmult, vqmult;

 Double_t tcal1[num_bars+1], tcal2[num_bars+1];

 int plmultV[num_planes];
 int plmultH[num_planes];
 
 Double_t dist = -1.e20;
 Double_t tstart;
 Double_t s, v[num_bars+1], beta;

 TRandom3 rnd;
 
 Double_t a,b,path,sx, sy, sz, sx2, sy2, sz2, sxy, syz, sxz;
 int n;
 Double_t ah,bh,pathh,sxh, syh, sx2h, sy2h, sxyh;
 int nh;
 
 Long64_t nentries = fChain->GetEntries();

 printf("get entries %lld, %lld\n", nentries, fChainx->GetEntries());
 
 //nentries = 100;
 
 //TH2F* hBarvsEvnt = new TH2F("hBarvsEvnt","Bar vs event",1000,0,nentries,num_bars,0.5,num_bars+0.5);
 //TH2F* hBarvsEvnt2 = new TH2F("hBarvsEvnt2","Bar vs event 2",1000,0,nentries,num_bars,0.5,num_bars+0.5);
 TH2F* hTofvsEvnt = new TH2F("hTofvsEvnt","Tof vs event ",1000,0,nentries,3000,0.,300.);
 TH2F* hEvntvsRun = new TH2F("hEvntvsRun","Evnt vs Run ",2000,2000.5,4000.5,1000,0,nentries);
 TH2F* hVetovsEvnt = new TH2F("hVetovsEvnt","Veto Edep vs event ",1000,0,nentries,3000,0.,50.);
 TH2F* hVTofvsEvnt = new TH2F("hVTofvsEvnt","Veto Tof vs event ",1000,0,nentries,3000,0.,250.);
 TH2F* hVTofvsEvntgamma = new TH2F("hVTofvsEvntgamma","Veto Tof vs event gamma",1000,0,nentries,3000,0.,250.);
 TH2F* hTofdiffvsEvnt = new TH2F("hTofdiffvsEvnt","Tn - Tv vs event ",1000,0,nentries,1000,-100.,100.);
 TH2F* hTofdiffcorrvsEvnt = new TH2F("hTofdiffcorrvsEvnt","Tncorr - Tv vs event ",1000,0,nentries,1000,-100.,100.);
 TH2F* hTstartvsEvnt = new TH2F("hTstartvsEvnt","Tstart vs event ",1000,0,nentries,2000,1000.,1300.);

 TH2F* hDXvsEvntcut = new TH2F("hDXvsEvntcut","DX vs event ",1000,0,nentries,400,-100.,100.);
 TH2F* hDYvsEvntcut = new TH2F("hDYvsEvntcut","DY vs event ",1000,0,nentries,400,-100.,100.);

  for(int i=1;i<=8;i++){
    sprintf(name1,"TofdiffcorrvsEvntBar%i",i);
    sprintf(name2,"Tncorr - Tv vs Evnt bar%i",i);     
    hTofdiffcorrvsEvntBar[i] = new TH2F(name1,name2,1000,0,nentries,1000,-100.,100.);
  }

 TH2F* hClmultvsEvnt = new TH2F("hClmultvsEvnt","Clmult vs event ",1000,0,nentries,50,-0.5,49.5);
 TH2F* hNmultvsEvnt = new TH2F("hNmultvsEvnt","Nmult vs event ",1000,0,nentries,50,-0.5,49.5);
 TH2F* hChmultvsEvnt = new TH2F("hChmultvsEvnt","Chmult vs event ",1000,0,nentries,50,-0.5,49.5);

 // for(int i=0;i<num_bars+1;i++){
 //   sprintf(name1,"T3vsEvntc%i",i);
 //   sprintf(name2,"T3 of tac channel %i",i);     
 //   hT3vsEvnt[i] = new TH2F(name1,name2,100,0,nentries,500,0.,2000.);
 //   sprintf(name1,"T4vsEvntc%i",i);
 //   sprintf(name2,"T4 of tac channel %i",i);     
 //   hT4vsEvnt[i] = new TH2F(name1,name2,1000,0,nentries,100,0.,2000.);
 // }
 // for(int i=0;i<num_bars+1;i++){
 //   sprintf(name1,"hPed1vsEvntc%i",i);
 //   sprintf(name2,"Ped1 of tac channel %i",i);     
 //   hPed1vsEvnt[i] = new TH2F(name1,name2,100,0,nentries,200,-10.,200.);
 // }
 // for(int i=0;i<num_bars+1;i++){
 //   sprintf(name1,"hPed2vsEvntc%i",i);
 //   sprintf(name2,"Ped2 of tac channel %i",i);     
 //   hPed2vsEvnt[i] = new TH2F(name1,name2,100,0,nentries,200,-10.,200.);
 // }

 Int_t runold=0;
 Int_t evt_cnt=0;
 Int_t evtcnt=0;

 //char katname[100];
 //char kyoname[100];
 //char ampname[100];

 double totdetene=0.;

 Int_t tpc_cnt = 0;
 Int_t tpc_track_cnt = 0;
 //UShort_t tpc_tracks[1000000]={0}; 
 
 Int_t current_run;
 Int_t current_event;
 Int_t event_offset = 0;

 fChainx->GetEvent(tpc_cnt);
 current_run = run;
 current_event = eventid + event_offset;

 
 for (Long64_t jentry=0; jentry<nentries;jentry++) {

   fChain->GetEvent(jentry);

   //printf("----------------\n");	  
   if ((float(jentry)/20000.)==int(jentry/20000)) {
     cout << "event: " << jentry << " of " << nentries << endl;
   }

   //if (nl_mul==0&&nlv_mul==0) continue;

   RunNumber = run2run[fChain->GetTreeNumber()];

   //printf("run %d    tpcrun %d\n", RunNumber, run); 
   
   evt_cnt++;
   if(RunNumber!=runold) {
     evt_cnt=1;
     runold=RunNumber;
   }

   if (RunNumber>=3077) {
     for (int i=79;i<=94;i++) {
       tsynccorr[i]=-0.18;
     }
   } else {
     for (int i=79;i<=94;i++) {
       tsynccorr[i]=0.0;
     }
   }

   VT0 = 1159.3;

   if (RunNumber>=0&&RunNumber<10) T0 = 756.5-7.;
   if (RunNumber>=10&&RunNumber<1736) T0 = 756.5;
   if (RunNumber>=1736&&RunNumber<1745) T0 = 756.5+40.-185.;
   if (RunNumber>=1745&&RunNumber<1747) T0 = 756.5;
   if (RunNumber>=1747&&RunNumber<1781) T0 = 756.5+40.-185.;
   if (RunNumber>=1781&&RunNumber<1785) T0 = 756.5;
   if (RunNumber>=1785&&RunNumber<1787) T0 = 756.5-5.;
   if (RunNumber>=1787&&RunNumber<1801) T0 = 756.5;
   if (RunNumber>=1801&&RunNumber<1809) T0 = 756.5+40.-185.;
   if (RunNumber>=1809&&RunNumber<1811) T0 = 756.5-5.;
   if (RunNumber>=1811&&RunNumber<1818) T0 = 756.5+40.-185.;
   if (RunNumber>=1818&&RunNumber<1847) T0 = 756.5-5.;
   if (RunNumber>=1847&&RunNumber<1882) T0 = 756.5+40.-270.;
   if (RunNumber>=1882&&RunNumber<2180) T0 = 756.5+40.-340.;
   if (RunNumber>=2180&&RunNumber<2210) T0 = 360.2;
   if (RunNumber>=2210&&RunNumber<2839) T0 = 560.2;

   // EOS2
   if (RunNumber>=2839&&RunNumber<2908) T0 = 584.5;
   if (RunNumber>=2908&&RunNumber<3010) T0 = 584.5;
   if (RunNumber>=3010&&RunNumber<3013) {
     T0 = 584.7;
     VT0 = 1159.5;
   }
   if (RunNumber>=3013&&RunNumber<3014) {
     T0 = 584.3;
     VT0 = 1159.1;
   }
   if (RunNumber>=3014&&RunNumber<3115) {
     T0 = 584.8;
     VT0 = 1159.6;
   }

   if (RunNumber>=3037&&RunNumber<3040) {
     T0 = 584.5;
     VT0 = 1159.0;
   }

   if (RunNumber>=3040&&RunNumber<3055) {  // no runs
     T0 = 584.5;
     VT0 = 1159.0;
   }

   if (RunNumber>=3056&&RunNumber<3105) {
     T0 = 584.6;
     VT0 = 1159.25;
   }

   if (RunNumber>=3105&&RunNumber<3113) {
     T0 = 584.6;
     VT0 = 1159.2;
   }

   if (RunNumber>=3113&&RunNumber<3115) {  // almost no events
     T0 = 584.5;
     VT0 = 1159.0;
   }

   if (RunNumber>=3115&&RunNumber<3126) {  // almost no events
     T0 = 607.4;
     VT0 = 1181.9;
   }

   if (RunNumber>=3127&&RunNumber<3137) {
     T0 = 607.4;
     VT0 = 1181.9;
   }

   if (RunNumber>=3137) {
     T0 = 582.7;          // 582.7
     VT0 = 1157.3;
   }

   T0 = T0 + 9.16 + 0.2 - 0.15;
 
   mult=0;
   multpart=0;
   multpmt=0;
   multV=0;
   multH=0;

   vmult=0;
   vqmult=0;
   
   for (int p=0;p<num_planes;p++) {
     plmultV[p]=0;
     plmultH[p]=0;
   }
   
   for (int i=0;i<=num_bars;i++) {
     t1[i]=0.;
     t2[i]=0.;
     t3[i]=0.;
     t4[i]=0.;
     t5[i]=0.;
     t6[i]=0.;
     e1[i]=0.;
     e2[i]=0.;
     fired[i]=0;
     fired[i]=0;
     
     tt1[i]=0.;
     tt2[i]=0.;
     tt1nowalk[i]=0.;
     tt2nowalk[i]=0.;
     ee1[i]=0.;
     ee2[i]=0.;
     ehit[i]=0.;
     thit[i]=-1000.;
     pos[i]=-1000.;
     tof[i]=-1000.;
     tofc[i]=-1000.;
     tofnowalk[i]=-1000.;

     tcal1[i]=0.;
     tcal2[i]=0.;

     clbar[i]=0; 
     for(Int_t j=1;j<=num_bars;j++){
       cl[i][j]=0;
     }
   }
   
   Double_t tms1 = start_traw_to_tcal(2, nl_tms1);
   Double_t tms2 = nl_tms2 * 24.9982;
   Double_t tms3 = start_traw_to_tcal(6, nl_tms3);
   tstart = tms1 + tms2 - tms3;

   hTstartvsEvnt->Fill(jentry,tstart);

   Int_t vdoubles[9];
   Int_t vneighbour[9];

   for (iHit=0;iHit<nlv_nbar;iHit++) {
     if (nlv_bar[iHit]<8) {

       if (nlv_qraw0[iHit]>5000) nlv_tcal0[iHit]=(nlv_traw0[iHit]+3.2)/13.;
       if (nlv_qraw1[iHit]>5000) nlv_tcal1[iHit]=(nlv_traw1[iHit]+13.2)/13.;

       if (RunNumber<1906) {
	 nlv_tcal0[iHit] = 4.*nlv_tcal0[iHit] - 120.;
	 nlv_tcal1[iHit] = 4.*nlv_tcal1[iHit] - 120.;
       }
       nlv_qcal0[iHit] = nlv_qcal0[iHit] - vped1[nlv_bar[iHit]];
       nlv_qcal1[iHit] = nlv_qcal1[iHit] - vped2[nlv_bar[iHit]];

       nlv_tcal0[iHit] = nlv_tcal0[iHit] - vtoff[nlv_bar[iHit]] -  
	 vtdiff[nlv_bar[iHit]]/2.-tstart + VT0;
       nlv_tcal1[iHit] = nlv_tcal1[iHit] - vtoff[nlv_bar[iHit]] +
	 vtdiff[nlv_bar[iHit]]/2.-tstart + VT0;

       hVetovsRun->Fill(RunNumber,sqrt(nlv_qcal0[iHit]*nlv_qcal1[iHit]));
       hVetovsEvnt->Fill(jentry,sqrt(nlv_qcal0[iHit]*nlv_qcal1[iHit]));

       vdoubles[nlv_bar[iHit]+1]=0;
       vneighbour[nlv_bar[iHit]+1]=0;

     }
   }

   for (iHit=0;iHit<nl_nbar;iHit++) {
     iBar=nl_bar[iHit]+1;
     t1[iBar] = nnp_traw_to_tcal(iBar, 1, nl_t1[iHit]);
     t2[iBar] = nnp_traw_to_tcal(iBar, 2, nl_t2[iHit]);
     t3[iBar] = nl_t3[iHit] * 24.9982;
     t4[iBar] = nl_t4[iHit] * 24.9982;
     t5[iBar] = nnp_traw_to_tcal(iBar, 5, nl_t5[iHit]);
     t6[iBar] = nnp_traw_to_tcal(iBar, 6, nl_t6[iHit]);

     t1[iBar] = t1[iBar] + t3[iBar] - t5[iBar];
     t2[iBar] = t2[iBar] + t4[iBar] - t6[iBar];
 
     tcal1[iBar] = nl_tcal0[iHit];
     tcal2[iBar] = nl_tcal1[iHit];

     e1[iBar] = QDC(nl_qraw0[iHit]-pedestal1[iBar]);
     e2[iBar] = QDC(nl_qraw1[iHit]-pedestal2[iBar]);
 
     //ee1[iBar] = nl_qcal0[iHit];
     //ee2[iBar] = nl_qcal1[iHit];
     if (nl_fired0[iHit]==1) fired[iBar]++;
     if (nl_fired1[iHit]==1) fired[iBar]++;
   }

   veto = 0;
   vetoq= 0;

   for(int j=0;j<nlv_nbar;j++) {
     if (nlv_bar[j]==8) continue;
     
     if((nlv_qcal0[j]>0.1||nlv_qcal0[j]!=nlv_qcal0[j])||
	(nlv_qcal1[j]>0.1||nlv_qcal1[j]!=nlv_qcal1[j])) {
       vetoq=1;
       vqmult++;
     }

     if(nlv_tcal0[j]>0&&nlv_tcal1[j]>0) {
       veto=1;
       vmult++;
       for(int i=0;i<nlv_nbar;i++) {
	 if (nlv_bar[i]!=8&&fabs(fabs(nlv_bar[i]-nlv_bar[j])-1)<0.5&&nlv_tcal0[i]>0&&nlv_tcal1[i]>0) {
	   vneighbour[nlv_bar[i]+1]=1;
	   vneighbour[nlv_bar[j]+1]=1;
	 }
       }
     }
   }

   for(int j=0;j<nlv_nbar;j++) {
     if (nlv_bar[j]==8) continue;
     if (nlv_bar[j]<7&&vneighbour[nlv_bar[j]+1]==1&&vneighbour[nlv_bar[j]+2]==1) {
       nlv_tcal0[j] = nlv_tcal0[j] - vtdiffcorr[nlv_bar[j]]/2.-vtoffcorr[nlv_bar[j]];
       nlv_tcal1[j] = nlv_tcal1[j] + vtdiffcorr[nlv_bar[j]]/2.-vtoffcorr[nlv_bar[j]];
     }
     if (nlv_bar[j]>0&&vneighbour[nlv_bar[j]+1]==1&&vneighbour[nlv_bar[j]]==1) {
       nlv_tcal0[j] = nlv_tcal0[j] - vtdiffcorr[nlv_bar[j]-1]/2.-vtoffcorr[nlv_bar[j]-1];
       nlv_tcal1[j] = nlv_tcal1[j] + vtdiffcorr[nlv_bar[j]-1]/2.-vtoffcorr[nlv_bar[j]-1];
     }

     hVTofvsEhitall->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
     hVTofvsEhit[nlv_bar[j]+1]->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
     hVTofvsRun->Fill(RunNumber,(nlv_tcal0[j]+nlv_tcal1[j])/2.);
     hVTofvsEvnt->Fill(jentry,(nlv_tcal0[j]+nlv_tcal1[j])/2.);
     hVTofvsBar->Fill(nlv_bar[j],(nlv_tcal0[j]+nlv_tcal1[j])/2.);

   }


   Float_t xcl, ycl, zcl, tcl;
   Int_t clno=0;
   Int_t first_cl=1;
   Int_t olddim, oldcl=0, newclno=0;
   Int_t not_in_cluster;
   
   for (iBar=1;iBar<=num_bars;iBar++) {
     multpmt+=fired[iBar];
     if (fired[iBar]>0) {
       multpart++;
       //hBarvsEvnt->Fill(jentry,iBar);
     }
     if (fired[iBar]==2) {
       mult++;
       multV = multV + ((iBar-1)/50)%2;
       multH = multH + ((iBar-1)/50+1)%2;
       
       plane = (iBar-1)/100;
       
       plmultV[plane] = plmultV[plane] + ((iBar-1)/50)%2;
       plmultH[plane] = plmultH[plane] + ((iBar-1)/50+1)%2;
       
       //hBarvsEvnt2->Fill(jentry,iBar);
       
       if (e1[iBar]<0) printf("kaj je ovo?\n");
       
       hBars->Fill(iBar);
       hT1vsBar->Fill(iBar,t1[iBar]);
       hT2vsBar->Fill(iBar,t2[iBar]);
       hT3vsBar->Fill(iBar,t3[iBar]);
       hT5vsBar->Fill(iBar,t5[iBar]);
       hT4vsBar->Fill(iBar,t4[iBar]);
       hT6vsBar->Fill(iBar,t6[iBar]);
       hE1vsBar->Fill(iBar,e1[iBar]);
       hE2vsBar->Fill(iBar,e2[iBar]);
       
       //hT3vsEvnt[iBar]->Fill(jentry,t3[iBar]);
       //hT4vsEvnt[iBar]->Fill(jentry,t4[iBar]);
       //hPed1vsEvnt[iBar]->Fill(jentry,e1[iBar]);
       //hPed2vsEvnt[iBar]->Fill(jentry,e2[iBar]);

       // sync level

       ee1[iBar]=e1[iBar]*ediff[iBar]*esync[iBar];
       ee2[iBar]=e2[iBar]/ediff[iBar]*esync[iBar];
       
       tt1[iBar]=t1[iBar]-tdiff[iBar]/2.-tsync[iBar]-tsynccorr[iBar]+wlk(e1[iBar]);
       tt2[iBar]=t2[iBar]+tdiff[iBar]/2.-tsync[iBar]-tsynccorr[iBar]+wlk(e2[iBar]);
       tt1nowalk[iBar]=t1[iBar]-tdiff[iBar]/2.-tsync[iBar];
       tt2nowalk[iBar]=t2[iBar]+tdiff[iBar]/2.-tsync[iBar];

       htt1_tcal1->Fill(iBar,tt1[iBar]-tstart-tcal1[iBar]);
       htt2_tcal2->Fill(iBar,tt2[iBar]-tstart-tcal2[iBar]);

       // hit level
       pos[iBar]=(tt1[iBar]-tt2[iBar])*vscint[iBar];
       thit[iBar]=(tt1[iBar]+tt2[iBar])/2.;
 
       // PMT saturation
       ee1[iBar] = ee1[iBar]/(1-0.0110*ee1[iBar]); // 0.0096
       ee2[iBar] = ee2[iBar]/(1-0.0110*ee2[iBar]);
       ehit[iBar]=sqrt(ee1[iBar]*ee2[iBar]);

       if (ehit[iBar]<0) ehit[iBar]=10000.;

       if (fabs(pos[iBar])>500.) continue;

       hEvsBar->Fill(iBar,ehit[iBar]);

       if (((iBar-1)/50)%2==1) {
	 xhit[iBar] = (iBar-76-plane*100)*5.+2.5 + rnd.Uniform(5.) - 2.5;
	 yhit[iBar] = pos[iBar];
	 zhit[iBar] = 5. + plane*10.+2.5 + rnd.Uniform(5.) - 2.5;
       } else {
	 xhit[iBar] = pos[iBar];
	 yhit[iBar] = (iBar-26-plane*100)*5.+2.5 + rnd.Uniform(5.) - 2.5;
	 zhit[iBar] = plane*10.+2.5 + rnd.Uniform(5.) - 2.5;
       }

       //thit[iBar]=thit[iBar]+0.2*((int)(zhit[iBar]/2.5-1.)); // added delay to avoid time resolution problem ???

       hXY[(iBar-1)/50+1]->Fill(xhit[iBar],yhit[iBar]);

       s = sqrt((Distance + zhit[iBar])*(Distance + zhit[iBar]) +
		(xhit[iBar])*(xhit[iBar]) + yhit[iBar]*yhit[iBar]);

       tof[iBar] = thit[iBar] - tstart + T0;
       tofc[iBar] = tof[iBar] / s * Distance;
       tofnowalk[iBar] =  ((tt1nowalk[iBar]+tt2nowalk[iBar])/2. - tstart + T0) / s * Distance;
      
       hTstart->Fill(tstart);
      
       v[iBar] = s/tof[iBar];
       beta = v[iBar]/c;
       
       hTofvsBar->Fill(iBar,tof[iBar]);
       hTofcvsBar->Fill(iBar,tofc[iBar]);
       hVelvsBar->Fill(iBar,v[iBar]);

       //hTofvsX->Fill(xhit[iBar],tof[iBar]);
       //hTofvsY->Fill(yhit[iBar],tof[iBar]);

       hTofvsEhit[(iBar-1)/50+1]->Fill(ehit[iBar],tofc[iBar]);
       hTofvsEhitall->Fill(ehit[iBar],tofc[iBar]);

       hTdiffvsBar->Fill(iBar,tt1[iBar]-tt2[iBar]);
      
       hTofvsEvnt->Fill(jentry,tofc[iBar]);
       hTofvsRun->Fill(RunNumber,tof[iBar]);
       hEvntvsRun->Fill(RunNumber,jentry);

       hTofvsPath->Fill(s,tof[iBar]);
       
       if(ehit[iBar]>0.){
	 hTof[(iBar-1)/50+1]->Fill(tof[iBar]);
	 hTofc[(iBar-1)/50+1]->Fill(tofc[iBar]);
	 if(ehit[iBar]>10.0) hTofch[(iBar-1)/50+1]->Fill(tofc[iBar]);
	 hEdiffvsBar->Fill(iBar,log(ee1[iBar]/ee2[iBar]));
	 hTofvsZall->Fill(zhit[iBar],tof[iBar]);
       }

        // walk correction calibration
       hTofvsErawall->Fill(e1[iBar],tofnowalk[iBar]);
       hTofvsErawall->Fill(e2[iBar],tofnowalk[iBar]);
       hTofvsErawallcorr->Fill(e1[iBar],tofc[iBar]);
       hTofvsErawallcorr->Fill(e2[iBar],tofc[iBar]);

       //printf("%d %lf\n",iBar,ehit[iBar],tof[iBar]);

       if (vetoq>=0) { // ******************* CLUSTERS ********************
	 
	 if(fabs(yhit[iBar])<170.&&fabs(xhit[iBar])<170.&&v[iBar]>0.) { //&&v[iBar]<27) {  // y 95
	   
	   not_in_cluster=0;
	   
	   for(Int_t i=1;i<=clno;i++) { // existing clusters 
	     //printf("cluster %d  dim  %d \n",i, cl[i][0]);
	     for(Int_t j=1;j<=cl[i][0];j++) { // hits in clusters
	       if (cl[i][j]==iBar) continue;  // don't take the same bar if in cluster
	       
	       if (clbar[iBar]==i) continue; // if bar is in the cluster???
	       
	       //printf("    hit %d clbar %d\n",cl[i][j], clbar[cl[i][j]]);
	       
	       xcl = xhit[cl[i][j]];
	       ycl = yhit[cl[i][j]];
	       zcl = zhit[cl[i][j]];
	       tcl = thit[cl[i][j]];
	       
	       if (fabs(xhit[iBar]-xcl)<7.5&&
		   fabs(yhit[iBar]-ycl)<7.5&&
		   fabs(zhit[iBar]-zcl)<7.5) {
		   //&&fabs(thit[iBar]-tcl-0.4)<0.6) {

		 hTdiffClust->Fill(zhit[iBar]-zcl,thit[iBar]-tcl);
		 
		 not_in_cluster=0;
		 
		 if (clbar[iBar]==0) { // not yet in cluster
		   //printf("not yet\n");
		   cl[i][0]++;
		   //printf("cl %d dim  %d\n",i,cl[i][0]);
		   cl[i][cl[i][0]]=iBar;
		   clbar[iBar]=i;
		   newclno = clno;
		 } else if (clbar[iBar]!=i) { // already in a different cluster clbar[iBar]
		   //printf("already\n");
		   
		   oldcl = clbar[iBar];
		   olddim=cl[oldcl][0];
		   
		   //printf("oldcl   %d   dim   %d\n", oldcl, olddim);
		   
		   for(Int_t k=1;k<=cl[i][0];k++) {
		     
		     //printf("cluster %d  hit %d  bar  %d\n", i, k, cl[i][k]);
		     
		     cl[oldcl][olddim+k] = cl[i][k];
		     clbar[cl[i][k]] = oldcl;
		     
		     cl[i][k] = 0;
		     
		   }
		   cl[oldcl][0] = olddim+cl[i][0];
		   cl[i][0]=0;
		   newclno = clno - 1;
		   
		   if (i<clno) {
		     //printf("nije zadnji cluster\n");
		     for (Int_t l=i+1;l<=clno;l++) {
		       cl[l-1][0]=cl[l][0];
		       for (Int_t m=1;m<=cl[l][0];m++) {
			 cl[l-1][m]=cl[l][m];
			 cl[l][m]=0;
			 clbar[cl[l-i][m]]=l-1;
		       }
		       cl[l][0]=0;
		     }
		   }
		 }
	       } else {
		 //printf("not neighbour -> next hit\n");
		 if (clbar[iBar]==0) not_in_cluster = 1;
	       }
	     }
	   }
	   
	   if (not_in_cluster==1) {
	     newclno = clno+1;
	     //printf("new cluster  %d\n", newclno);
	     cl[newclno][0]=1;
	     cl[newclno][1]=iBar;
	     clbar[iBar]=newclno;
	   }
	   
	   clno=newclno;
	   if (first_cl==1) {
	     //printf("first iBar %d\n", iBar);
	     clno=1;
	     cl[1][0]=1;
	     cl[1][1]=iBar;
	     clbar[iBar]=1;
	     first_cl=0;
	   }
	 }
       }   // ******************* CLUSTERS ********************
       
       //     for(Int_t i=1;i<=clno;i++) {
       //       printf("cluster %d\n",i);
       //       for(Int_t j=1;j<=cl[i][0];j++) {
       // 	printf("       hit %d\n",cl[i][j]);
       // 	printf("       clbar %d\n",clbar[cl[i][j]]);
       //       }
       //     }
       
     }
   }
   
   hClusterMult->Fill(clno);
   hClmultvsEvnt->Fill(jentry,clno);
   int nclno=0, cclno=0;

   double axfit = -1000., ayfit = -1000.;
   
   if (clno>0) {
     for (Int_t i=1;i<=clno;i++) {
       clusters[i].dim=cl[i][0];

       hClusterDim->Fill(cl[i][0]);
       clusters[i].etot=0.;

       clusters[i].tf = 10000.;
       clusters[i].tl = 0.;
 
       sx=0.;
       sy=0.;
       sz=0.;
       sx2=0.;
       sy2=0.;
       sz2=0.;
       sxz=0.;
       syz=0.;
 
       for (Int_t j=1;j<=clusters[i].dim;j++) {
	 //if (ehit[cl[i][j]]>5.)
	 clusters[i].etot = clusters[i].etot + ehit[cl[i][j]];

	 sx = sx + xhit[cl[i][j]];
	 sy = sy + yhit[cl[i][j]]; 
	 sz = sz + zhit[cl[i][j]];
	 sx2 = sx2 + xhit[cl[i][j]]*xhit[cl[i][j]];
	 sy2 = sy2 + yhit[cl[i][j]]*yhit[cl[i][j]];
	 sz2 = sz2 + zhit[cl[i][j]]*zhit[cl[i][j]];
	 sxz = sxz + xhit[cl[i][j]]*zhit[cl[i][j]];
	 syz = syz + yhit[cl[i][j]]*zhit[cl[i][j]];	 
	 
	 // the first hit in the cluster
	 if (tof[cl[i][j]]<clusters[i].tf) {
	   clusters[i].tf=tof[cl[i][j]];
	   clusters[i].xf=xhit[cl[i][j]];
	   clusters[i].yf=yhit[cl[i][j]];
	   clusters[i].zf=zhit[cl[i][j]];
	 }
	 // the last hit in the cluster
	 if (tof[cl[i][j]]>clusters[i].tl) {
	   clusters[i].tl=tof[cl[i][j]];
	   clusters[i].xl=xhit[cl[i][j]];
	   clusters[i].yl=yhit[cl[i][j]];
	   clusters[i].zl=zhit[cl[i][j]];
	 }	    
       }

       // "direction" of the cluster
       if(clusters[i].dim>1&&clusters[i].zl!=clusters[i].zf) {

	 // from linear fit
	 axfit = (clusters[i].dim*sxz-sx*sz)/(clusters[i].dim*sz2-sz*sz);  //  X-Z plane
	 ayfit = (clusters[i].dim*syz-sy*sz)/(clusters[i].dim*sz2-sz*sz); //  Y-Z plane

	 // from the first and the last hit in the cluster only
	 clusters[i].ax = (clusters[i].xl-clusters[i].xf)/(clusters[i].zl-clusters[i].zf);  //  X-Z plane
   	 clusters[i].ay = (clusters[i].yl-clusters[i].yf)/(clusters[i].zl-clusters[i].zf); //  Y-Z plane

       } else { // only one hit in the cluster or zl == zf
	 clusters[i].ax = clusters[i].xf/(Distance+clusters[i].zf);
	 clusters[i].ay = clusters[i].yf/(Distance+clusters[i].zf);
       }

       // difference of two methods
       if (clusters[i].dim>2&&clusters[i].zl!=clusters[i].zf) {
	 hDayvsDax->Fill(clusters[i].ax - axfit, clusters[i].ay - ayfit);
       }
     } 
     
     TCluster clhelp;
     
     // sort clusters in time
     for (Int_t j=clno;j>=2;j--) {
       for (Int_t i=1;i<j;i++) {
     	 if(clusters[i].tf>clusters[i+1].tf) {
     	   clhelp = clusters[i];
     	   clusters[i] = clusters[i+1];
     	   clusters[i+1] = clhelp;
     	 }
       }
     }
     
     // discontinous clusters
     for (Int_t i=1;i<clno;i++) {
       for (Int_t j=i+1;j<=clno;j++) {
	 
	 hDirdiff->Fill(clusters[i].xf/(Distance+clusters[i].zf)-clusters[j].xf/(Distance+clusters[j].zf),
			clusters[i].yf/(Distance+clusters[i].zf)-clusters[j].yf/(Distance+clusters[j].zf));
	 
	 //if(fabs(clusters[i].xl/(Distance+clusters[i].zl)-clusters[j].xf/(Distance+clusters[j].zf))<0.02&&  //0.02
	 //   fabs(clusters[i].yl/(Distance+clusters[i].zl)-clusters[j].yf/(Distance+clusters[j].zf))<0.02&&  //0.02
	 if(sqrt(pow((clusters[i].xl/(Distance+clusters[i].zl)-clusters[j].xf/(Distance+clusters[j].zf)),2.)+
	 	 pow((clusters[i].yl/(Distance+clusters[i].zl)-clusters[j].yf/(Distance+clusters[j].zf)),2.))<0.03&&  //0.02
	    fabs(clusters[i].zl-clusters[j].zf)<22.5) {  // 12.5 or 17.5
	   
	   clusters[i].etot=clusters[i].etot+clusters[j].etot;
	   clusters[i].dim=clusters[i].dim+clusters[j].dim;
	   clusters[i].tl=clusters[j].tl;
	   clusters[i].xl=clusters[j].xl;
	   clusters[i].yl=clusters[j].yl;
	   clusters[i].zl=clusters[j].zl;
	   
	   for (Int_t k=j;k<=clno;k++) clusters[k]=clusters[k+1];
	   j=j-1;
	   clno=clno-1;
	 }
       }
       if(clusters[i].dim>1) {
	 clusters[i].ax = (clusters[i].xl-clusters[i].xf)/(clusters[i].zl-clusters[i].zf);  //  X-Z plane
   	 clusters[i].ay = (clusters[i].yl-clusters[i].yf)/(clusters[i].zl-clusters[i].zf); //  Y-Z plane
       } else {
	 clusters[i].ax = clusters[i].xf/(Distance+clusters[i].zf);
	 clusters[i].ay = clusters[i].yf/(Distance+clusters[i].zf);
       }
     }

     for (Int_t i=1;i<clno;i++) {
       for (Int_t j=i+1;j<=clno;j++) {
	 hDirdiff2->Fill(clusters[i].xf/(Distance+clusters[i].zf)-clusters[j].xf/(Distance+clusters[j].zf),
			clusters[i].yf/(Distance+clusters[i].zf)-clusters[j].yf/(Distance+clusters[j].zf));

	 hClustCorrel_xz->Fill(clusters[j].zf-clusters[i].zf,clusters[j].xf-clusters[i].xf);
	 hClustCorrel_yz->Fill(clusters[j].zf-clusters[i].zf,clusters[j].yf-clusters[i].yf);
	 hClustCorrel_tz->Fill(clusters[j].zf-clusters[i].zf,clusters[j].tf-clusters[i].tf);

	 hClustCorrel_xx->Fill(clusters[i].xf,clusters[j].xf-clusters[i].xf);
	 hClustCorrel_yy->Fill(clusters[i].yf,clusters[j].yf-clusters[i].yf);
	 hClustCorrel_tt->Fill(clusters[i].tf,clusters[j].tf-clusters[i].tf);
	 

       }
     }
     
     Double_t st0, sp0, ct0, cp0;
     Double_t stp, spp, ctp, cpp;
     Double_t stn, spn, ctn, cpn;
     Double_t sintp, costp, sintn, costn;
     Double_t axn, ayn;

     // elastic scattering, still not done and applied
     for (Int_t i=1;i<clno;i++) {
       
       st0 = sqrt(clusters[i].xf*clusters[i].xf + clusters[i].yf*clusters[i].yf)/
	 sqrt(clusters[i].xf*clusters[i].xf + clusters[i].yf*clusters[i].yf +
	      (Distance+clusters[i].zf)*(Distance+clusters[i].zf));
       ct0 = sqrt(1-st0*st0);
       
       sp0 = sin(atan2(clusters[i].xf,clusters[i].yf));
       cp0 = cos(atan2(clusters[i].xf,clusters[i].yf));
       
       stp = sqrt((clusters[i].ax*clusters[i].ax+clusters[i].ay*clusters[i].ay)/
		  (clusters[i].ax*clusters[i].ax+clusters[i].ay*clusters[i].ay+1.));
       ctp = sqrt(1-stp*stp);
       
       spp = sin(atan2(clusters[i].ax,clusters[i].ay));
       cpp = sin(atan2(clusters[i].ax,clusters[i].ay));
       
       for (Int_t j=i+1;j<=clno;j++) {
	 
	 axn = (clusters[j].xf-clusters[i].xf)/(clusters[j].zf-clusters[i].zf);
	 ayn = (clusters[j].yf-clusters[i].yf)/(clusters[j].zf-clusters[i].zf);
	 
	 stn = sqrt((axn*axn+ayn*ayn)/
		    (axn*axn+ayn*ayn+1.));
	 ctn = sqrt(1-stn*stn);
	 
	 spn = sin(atan2(axn,ayn));
	 cpn = sin(atan2(axn,ayn));
	 
	 costp = st0*cp0*stp*cpp+st0*sp0*stp*spp+ct0*ctp;
	 costn = st0*cp0*stn*cpn+st0*sp0*stn*spn+ct0*ctn;
	 
	 sintp = sqrt(1-costp*costp);
	 sintn = sqrt(1-costn*costn);
       }
     } 

     Int_t doubleflag=0;
     
     // VETO multiple hits
     for(int j=0;j<nlv_nbar;j++) {
       if (nlv_bar[j]==8) continue;

       if (vdoublescut->IsInside(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.)) doubleflag=1;
      
       for (Int_t i=1;i<=clno;i++) {
	 if (fabs(clusters[i].xf*vDistance/(Distance+clusters[i].zf)-(3.5-nlv_bar[j])*31.5)<16.&&
	     //if (fabs(clusters[i].ax*(vDistance-Distance-clusters[i].zf)+clusters[i].xf-(3.5-nlv_bar[j])*31.5)<16.&&
	     nlv_tcal0[j]>0.&&nlv_tcal1[j]>0.) {
	   
	   vdoubles[nlv_bar[j]+1]++;

	 }
       }
     }     
 
     Int_t cl_bar=0;

     totdetene=0.;

     //if (vmult<=1) continue;

     //if (doubleflag==0) continue;
     
     // ###########################   loop over clusters  #################################
     for (Int_t i=1;i<=clno;i++) {

       //if (!(z1->IsInside(clusters[i].etot,clusters[i].tf))) continue;

       if (fabs(clusters[i].xf)>115.||fabs(clusters[i].yf)>95.) continue;  // not shadowed by VETO

       if(((int)(clusters[i].zf/5))%2==0) {
	 cl_bar = (int)(clusters[i].zf/5)*50+(int)((clusters[i].yf+130.)/5);
       } else {
	 cl_bar = (int)(clusters[i].zf/5)*50+(int)((clusters[i].xf+130.)/5);
       }
            
       if(1==1) {  // condition on clusters, e.g. Tof vs Ehit cuts
       //if(cl_bar>350&&cl_bar<401) {  // condition on clusters, e.g. Tof vs Ehit cuts
       //if (clusters[i].dim>7) {

	 hTofvsEhitClust->Fill(clusters[i].etot,clusters[i].tf);
	 hTofvsEhitClustp[(cl_bar-1)/50+1]->Fill(clusters[i].etot,clusters[i].tf);
	 hPlane[(cl_bar-1)/50+1]->Fill(clusters[i].xf,clusters[i].yf);

	 int in_veto = 0;
	 int in_veto1 = 0;
	 int in_veto2 = 0;
	 for(int j=0;j<nlv_nbar;j++) {
	   if (nlv_bar[j]==8) continue;
	   
	   if (nlv_tcal0[j]>0.&&nlv_tcal1[j]>0.) {

	     if (fabs((Distance+clusters[i].zf)/clusters[i].tf-30.)<1.) { // gammas in VETO and NeuLAND
	       if (fabs((nlv_tcal0[j]+nlv_tcal1[j])/2.-37.)<5.) {
		 hVTofvsXgamma->Fill(clusters[i].xf,clusters[i].tf-(nlv_tcal0[j]+nlv_tcal1[j])/2.+70.);
		 hVTofvsBargamma->Fill(nlv_bar[j]+1,(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 hVTofvsEvntgamma->Fill(jentry,(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	       }
	     }

	     if (0.91*clusters[i].tf+2.-(nlv_tcal0[j]+nlv_tcal1[j])/2.<0.) {
	       hTofvsEhitClustBack->Fill(clusters[i].etot,clusters[i].tf);
	       hPlaneBack[(cl_bar-1)/50+1]->Fill(clusters[i].xf,clusters[i].yf);
	     }

	     hVEhitvsEhit->Fill(clusters[i].etot,sqrt(nlv_qcal0[j]*nlv_qcal1[j]));

	     hDYvsDX->Fill(clusters[i].xf*vDistance/(Distance+clusters[i].zf)-(3.5-nlv_bar[j])*31.5,
			   clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
	     
	     hVXvsX[(cl_bar-1)/50+1]->Fill(clusters[i].xf*vDistance/(Distance+clusters[i].zf),(3.5-nlv_bar[j])*31.5);
	     hVXvsXax[(cl_bar-1)/50+1]->Fill(clusters[i].ax*(vDistance-Distance-clusters[i].zf)+clusters[i].xf,(3.5-nlv_bar[j])*31.5);

	     hVTdiffvsY[(cl_bar-1)/50+1]->Fill(clusters[i].xf,clusters[i].yf+1*(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
	     
	     hVTofvsX->Fill(clusters[i].xf,clusters[i].tf-(nlv_tcal0[j]+nlv_tcal1[j])/2.+70.);

	     //check overlapping !!!!!!!
	     //if (vmult==1) {
	     if (vneighbour[nlv_bar[j]+1]==1) {
	       hVTdiffvsYneigh[(cl_bar-1)/50+1]->Fill(clusters[i].xf,clusters[i].yf+1*(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
	       hVTofvsXneigh->Fill(clusters[i].xf,clusters[i].tf-(nlv_tcal0[j]+nlv_tcal1[j])/2.+70.);
	       hVT1vsVT2->Fill(clusters[i].yf,-1*(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
	     }		 
	     
	     haycorr->Fill(clusters[i].yf,clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
	     haycorr2->Fill(clusters[i].yf/(Distance+clusters[i].zf)*vDistance,
			    clusters[i].yf/(Distance+clusters[i].zf)*vDistance+(nlv_tcal0[j]-nlv_tcal1[j])*0.9714/0.130);
	     if (clusters[i].dim==1)
	       haycorr1->Fill(clusters[i].ay*(vDistance-Distance-clusters[i].zf)+clusters[i].yf,(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
	     
	     
	     // ###################################### VETO CONDITIONS #######################################
	     if (fabs(clusters[i].xf*vDistance/(Distance+clusters[i].zf)-(3.5-nlv_bar[j])*31.5-1.1)<18.&& // <22 shift 2
		 fabs(clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130)<20.) {    //30
	       
	       in_veto=1;
	     }
	     if (fabs(clusters[i].xf*vDistance/(Distance+clusters[i].zf)-(3.5-nlv_bar[j])*31.5)<16.&&
		 ((vdoubles[nlv_bar[j]+1]==1&&fabs(clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130)<20.)
		  ||vdoubles[nlv_bar[j]+1]>1)) {
	       
	       in_veto1=1;
	     }
	     if (fabs(clusters[i].ax*(vDistance-Distance-clusters[i].zf)+clusters[i].xf-(3.5-nlv_bar[j])*31.5)<16.&&
		 ((vdoubles[nlv_bar[j]+1]==1&&fabs(clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130)<20.)
		  ||vdoubles[nlv_bar[j]+1]>1)) {
	       
	       in_veto2=1;
	     }
	     // ###################################### VETO CONDITIONS #######################################

	   } // condition for VETO hits
	 }  // loop over VETO paddles
	 
	 if (in_veto==0) { // ############## NEUTRAL CLUSTERS ####################
	   //if (in_veto1==0) {
	   nclno++;

	   if (vetoq==0) {
	     hTofvsEhitClustvstrong[(cl_bar-1)/50+1]->Fill(clusters[i].etot,clusters[i].tf);
	     hTofvsEhitClustVetoStrong->Fill(clusters[i].etot,clusters[i].tf);
	   }
	   hTofvsEhitClustVeto->Fill(clusters[i].etot,clusters[i].tf);
	   hTofvsBarClustVeto->Fill(cl_bar,clusters[i].tf);
	   hEvsBarClustVeto->Fill(cl_bar,clusters[i].etot);
	   hFirstZ->Fill(clusters[i].zf,clusters[i].tf);

	   hPlanecut[(cl_bar-1)/50+1]->Fill(clusters[i].xf,clusters[i].yf);
	   hTofvsEhitClustv[(cl_bar-1)/50+1]->Fill(clusters[i].etot,clusters[i].tf);

	   totdetene = totdetene + clusters[i].etot;
	   
	   for(int j=0;j<nlv_nbar;j++) {
	     if (nlv_bar[j]==8) continue;
	     
	     if (nlv_tcal0[j]>0.&&nlv_tcal1[j]>0.) {
	       //if (fabs(clusters[i].xf*vDistance/(Distance+clusters[i].zf)-(3.5-nlv_bar[j])*31.5-20.)<10.&& // remaining background
	       //   fabs(clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130)<20.) {
	       
		 
		 hVXvsXcut[(cl_bar-1)/50+1]->Fill(clusters[i].xf*vDistance/(Distance+clusters[i].zf),(3.5-nlv_bar[j])*31.5);
		 hVTdiffvsYcut[(cl_bar-1)/50+1]->Fill(clusters[i].xf,clusters[i].yf+1*(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
		 
		 hDYvsDXcut->Fill(clusters[i].xf*vDistance/(Distance+clusters[i].zf)-(3.5-nlv_bar[j])*31.5,
				  clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130);

		 hDXvsEvntcut->Fill(jentry,clusters[i].xf*vDistance/(Distance+clusters[i].zf)-(3.5-nlv_bar[j])*31.5);

		 hDYvsEvntcut->Fill(jentry,clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130);

		 hDYvsDXcorr->Fill(clusters[i].ax*(vDistance-Distance-clusters[i].zf)+clusters[i].xf-(3.5-nlv_bar[j])*31.5,
				   clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
		 		 
		 hDYvsVdoubles->Fill(vdoubles[nlv_bar[j]+1],
				     clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130);

		 hVetoEraw->Fill(nlv_bar[j]+1,sqrt(nlv_qraw0[j]*nlv_qraw1[j]));
		 hVetoEcal->Fill(nlv_bar[j]+1,sqrt(nlv_qcal0[j]*nlv_qcal1[j]));
		 hVetoEcal1->Fill(nlv_bar[j]+1,nlv_qcal0[j]);
		 hVetoEcal2->Fill(nlv_bar[j]+1,nlv_qcal1[j]);
		 
		 hDXvsX->Fill(clusters[i].xf,clusters[i].xf-(3.5-nlv_bar[j])*31.5);
		 hDXvsZ->Fill(clusters[i].zf,clusters[i].xf-(3.5-nlv_bar[j])*31.5);
		 hDYvsY->Fill(-1.*(nlv_tcal0[j]-nlv_tcal1[j])/0.130,clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
		 hDYvsZ->Fill(clusters[i].zf,clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
		 
		 hVTofvsTof->Fill(clusters[i].tf,(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 hVVelvsVel->Fill((Distance+clusters[i].zf)/clusters[i].tf,
				    vDistance/((nlv_tcal0[j]+nlv_tcal1[j])/2.));

		 hTofdiffvsEcal->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),clusters[i].tf-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 hTofdiffvsEhit->Fill(clusters[i].etot,clusters[i].tf-(nlv_tcal0[j]+nlv_tcal1[j])/2.);

		 hTofdiffvsVbar->Fill(nlv_bar[j]+1,0.91*clusters[i].tf+2.-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 hTofdiffvsBar->Fill(cl_bar,0.91*clusters[i].tf+2.-(nlv_tcal0[j]+nlv_tcal1[j])/2.);

		 hTofdiff->Fill(clusters[i].tf-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 hVTdiffEcal1[nlv_bar[j]+1]->Fill(nlv_tcal0[j]-nlv_tcal1[j],nlv_qcal0[j]);
		 hVTdiffEcal2[nlv_bar[j]+1]->Fill(nlv_tcal0[j]-nlv_tcal1[j],nlv_qcal1[j]);

		 hVTofvsEhitNoMatch->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	     }
	   }
	   //}
	 } else { // in_veto!=0  ################ CHARGED CLUSTERS #######################
	   cclno++;

	   hPlaneacut[(cl_bar-1)/50+1]->Fill(clusters[i].xf,clusters[i].yf);
	   hBaracut->Fill(cl_bar,clusters[i].yf);
	   hEvsBarClustNoVeto->Fill(cl_bar,clusters[i].etot);
	   hTofvsBarClustNoVeto->Fill(cl_bar,clusters[i].tf);
	   hTofvsEhitClustNoVeto->Fill(clusters[i].etot,clusters[i].tf);
	   
	   hXZclust->Fill(clusters[i].zf,clusters[i].xf);
	   hYZclust->Fill(clusters[i].zf,clusters[i].yf);

	   if (clusters[i].dim!=1) {
	     hclxz->Fill(clusters[i].xf,clusters[i].xl-clusters[i].xf);
	     hclyz->Fill(clusters[i].tl-clusters[i].tf,clusters[i].zl-clusters[i].zf);
	   }
	   
	   hDax->Fill(clusters[i].dim,
		      clusters[i].ax*(vDistance-Distance-clusters[i].zf)+clusters[i].xf-
		      clusters[i].xf*vDistance/(Distance+clusters[i].zf));
	   
	   hDay->Fill(clusters[i].dim,
		      clusters[i].ay*(vDistance-Distance-clusters[i].zf)+clusters[i].yf-
		      clusters[i].yf*vDistance/(Distance+clusters[i].zf));
	   
	   haxvsax->Fill(clusters[i].xf*vDistance/(Distance+clusters[i].zf),
			 clusters[i].ax*(vDistance-Distance-clusters[i].zf)+clusters[i].xf);
	   hayvsay->Fill(clusters[i].yf*vDistance/(Distance+clusters[i].zf),
			 clusters[i].ay*(vDistance-Distance-clusters[i].zf)+clusters[i].yf);
	   
	   hayvsax->Fill(clusters[i].ax*(vDistance-Distance-clusters[i].zf)+clusters[i].xf,
			 clusters[i].ay*(vDistance-Distance-clusters[i].zf)+clusters[i].yf);
	   
	   for(int j=0;j<nlv_nbar;j++) {
	     if (nlv_bar[j]==8) continue;
	     if (nlv_tcal0[j]>0.&&nlv_tcal1[j]>0.) {
	       
	       if (fabs(clusters[i].xf*vDistance/(Distance+clusters[i].zf)-(3.5-nlv_bar[j])*31.5-1.1)<18.&& // <22 shift 2
		   fabs(clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130)<20.) {    //30
		 

		 hTofvsEhitClustnv[(cl_bar-1)/50+1]->Fill(clusters[i].etot,clusters[i].tf);
		 hTofvsEhitClustNoVeto2->Fill(clusters[i].etot,clusters[i].tf);
		 hVTofvsEhitMatch2->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 
		 hVEhitvsEhitch->Fill(clusters[i].etot,sqrt(nlv_qcal0[j]*nlv_qcal1[j]));

		 hVTofvsEhitMatch->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 
		 hVetoTcal1->Fill(nlv_bar[j],nlv_tcal0[j]);
		 hVetoTcal2->Fill(nlv_bar[j],nlv_tcal1[j]);
		 
		 hVVelvsVelch->Fill((Distance+clusters[i].zf)/clusters[i].tf,
				    vDistance/((nlv_tcal0[j]+nlv_tcal1[j])/2.));
		 
		 hDYvsDXacut->Fill(clusters[i].xf*vDistance/(Distance+clusters[i].zf)-(3.5-nlv_bar[j])*31.5,
				   clusters[i].yf+(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
		 
		 hVXvsXch[(cl_bar-1)/50+1]->Fill(clusters[i].xf*vDistance/(Distance+clusters[i].zf),(3.5-nlv_bar[j])*31.5);
		 
		 hVTdiffvsYch[(cl_bar-1)/50+1]->Fill(clusters[i].xf,clusters[i].yf+1*(nlv_tcal0[j]-nlv_tcal1[j])/0.130);
		 
		 hVTofvsTofch->Fill(clusters[i].tf,(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 hTofdiffvsEcalch->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),clusters[i].tf-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 hTofdiffvsEhitch->Fill(clusters[i].etot,clusters[i].tf-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 
		 hTofdiffcorrvsEvnt->Fill(jentry,0.91*clusters[i].tf+2.-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 hTofdiffvsEvnt->Fill(jentry,clusters[i].tf-(nlv_tcal0[j]+nlv_tcal1[j])/2.);	       
		 
		 hTofdiffcorrvsEvntBar[nlv_bar[j]+1]->Fill(jentry,0.91*clusters[i].tf+2.-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 
		 hTofdiffvsVbarch->Fill(nlv_bar[j]+1,0.91*clusters[i].tf+2.-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 hTofdiffvsBarch->Fill(cl_bar,0.91*clusters[i].tf+2.-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		 
		 hTofdiffch->Fill(clusters[i].tf-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	       }
	     }
	   }

	   while (jentry+1 > current_event   && jentry < nentries) {
	     tpc_cnt++;
	     fChainx->GetEvent(tpc_cnt);
	     if (run != current_run) {
	       current_run = run;
	       printf("number of events  %d\n", current_event-event_offset);
	       event_offset = current_event;
	       printf("events do not match  run %d  offset %d\n", current_run, event_offset);
	     }
	     current_event = eventid + event_offset;
	   }	     

	   tpc_track_cnt = 0;

	   while (jentry+1 == current_event) {

	     if (mom!=0) {
	       
	       hdxvsdist->Fill(tpcdist, -1*dx*mom);
	       hdyvsdist->Fill(tpcdist, dy*mom);

	       if (sigma20z&&vz<-9.49569&&vz>-12.80121&&vx>-15&&vx<15&&vy<-206.06&&vy>-246.06
		   &&parentvid==vid&&tpcdist<5.) { //&&ndf>30&&abs(ea)<100) {
		 
		 hdyvsdx->Fill(-1*dx, dy);
		 hdyvsdz->Fill(dz, dy);
		 hdxvsdz->Fill(dz, -1*dx);
		 hYvsXTpc->Fill(-1.*dx*mom,dy*mom);
		 hYvsZTpc->Fill(mom*dz,mom*dy);
		 hXvsZTpc->Fill(mom*dz,-1.*mom*dx);
		 hDirTpc->Fill(atan2(-1.*dx,dz)*57.3);
		 
		 hDedxvsMom->Fill(mom/charge, dedx);
		 
		 hdxvsX->Fill(clusters[i].xf,-1.*mom*(-1*dx*cos(30./57.3)+dz*sin(30./57.3)));
		 hdyvsY->Fill(clusters[i].yf,dy*mom);

		 hTofvsMom->Fill(mom/charge,clusters[i].tf);
		 hDedxvsEhitClust->Fill(clusters[i].etot, dedx);
		 
		 hXYtpc->Fill(clusters[i].xf, clusters[i].yf);
		 
		 //if (dxvsXcut->IsInside(clusters[i].xf,(-1.*mom*(-1*dx*cos(30./57.3)+dz*sin(30./57.3))))
		 //  &&dyvsYcut->IsInside(clusters[i].yf,dy*mom/charge)) {
		 
		 
		 if (z1->IsInside(clusters[i].etot,clusters[i].tf)&&
		     protons_tpc->IsInside(mom/charge, dedx)) {

		   hDedxvsMomn->Fill(mom/charge, dedx);
		   hYvsXTpcn->Fill(-1.*mom*dx,mom*dy);
		   hYvsZTpcn->Fill(mom*dz,mom*dy);
		   hXvsZTpcn->Fill(mom*dz,-1.*mom*dx);
		   hdyvsdxn->Fill(-1*dx, dy);
		   hdyvsdzn->Fill(dz, dy);
		   hdxvsdzn->Fill(dz, -1*dx);
		   hDirTpcn->Fill(atan2(-1.*dx,dz)*57.3);
		   hdxvsXn->Fill(clusters[i].xf,-1.*mom*(-1*dx*cos(30./57.3)+dz*sin(30./57.3)));
		   hdyvsYn->Fill(clusters[i].yf,dy*mom);
		 }	 
	       
		 //if (!(fabs(sin(thetaL/57.3)*sin(phiL/57.3))<0.1&&
		 //    -1.*sin(thetaL/57.3)*cos(phiL/57.3)<-0.54&&
		 //    -1.*sin(thetaL/57.3)*cos(phiL/57.3)>-0.85)) {
		 //
		 
		 // if (fabs(dy)<0.1&&
		 //    -1.*dx<-0.54&&
		 //    -1.*dx>-0.85) hdyvsdxc->Fill(-1*dx*mom, dy*mom);
		 
		 //if (fabs(dy) < dx/5.&&-1*dx*mom<-250.
		 //   &&-1.*dx*mom<-0.5*dz*mom-210.&&-1.*dx*mom>-0.9*dz*mom-260.) {

		 if (YvsXtpccut->IsInside(-1.*dx*mom, dy*mom)&&
		     XvsZtpccut->IsInside(dz*mom, -1.*dx*mom)&&
		     dyvsdxtpccut->IsInside(-1.*dx, dy)&&
		     dyvsdztpccut->IsInside(dz, dy)) {
		   
		   //hdxvsdist->Fill(tpcdist, -1*dx*mom);
		   //hdyvsdist->Fill(tpcdist, dy*mom);
		   
		 hYvsXTpcc->Fill(-1.*dx*mom, dy*mom);
		 hYvsZTpcc->Fill(mom*dz,mom*dy);
		 hXvsZTpcc->Fill(mom*dz,-1.*mom*dx);

		 hdyvsdxc->Fill(-1*dx, dy);
		 hdyvsdzc->Fill(dz, dy);
		 hdxvsdzc->Fill(dz, -1*dx);
		   
		 hDirTpcc->Fill(atan2(-1.*dx,dz)*57.3);
		 
		 hdxvsMomc->Fill(mom, -1.*dx);
		 hdyvsMomc->Fill(mom, dy);
		 hdzvsMomc->Fill(mom, dz);

		 hDedxvsMomc->Fill(mom/charge, dedx);

		 if (protons_tpc->IsInside(mom/charge, dedx)) {
		 
		 if (z1->IsInside(clusters[i].etot,clusters[i].tf)) {
		   hXvsZTpccz1->Fill(mom*dz,-1.*mom*dx);
		   hDedxvsMomcz1->Fill(mom/charge, dedx);
		 }
		   
 		 if (protons->IsInside(clusters[i].etot,clusters[i].tf)) {
		   hXvsZTpccp->Fill(mom*dz,-1.*mom*dx);
		   hDedxvsMomcp->Fill(mom/charge, dedx);
		 }
		 
		 if (deuterons->IsInside(clusters[i].etot,clusters[i].tf)) {
		   hXvsZTpccd->Fill(mom*dz,-1.*mom*dx);
		   hDedxvsMomcd->Fill(mom/charge, dedx);
		 }
		 
 		 if (tritons->IsInside(clusters[i].etot,clusters[i].tf)) {
		   hXvsZTpcct->Fill(mom*dz,-1.*mom*dx);
		   hDedxvsMomct->Fill(mom/charge, dedx);
		 }
 
		 
		 hTofvsMomc->Fill(mom/charge,clusters[i].tf);
		 hDedxvsEhitClustc->Fill(clusters[i].etot, dedx);

		 
		 hXYtpcc->Fill(clusters[i].xf, clusters[i].yf);

		 hdxvsXc->Fill(clusters[i].xf,-1.*mom*(-1*dx*cos(30./57.3)+dz*sin(30./57.3)));
		 hdyvsYc->Fill(clusters[i].yf,dy*mom);

		 }
		 
		 hTofvsX->Fill(clusters[i].xf, clusters[i].tf);
		 hTofvsY->Fill(clusters[i].yf, clusters[i].tf);
		 
		 tpc_track_cnt++;
		 
		 for(int j=0;j<nlv_nbar;j++) {
		   if (nlv_bar[j]==8) continue;
		   if (nlv_tcal0[j]>0.&&nlv_tcal1[j]>0.) {
		     
		     if (helium->IsInside(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.)) {
		       hXvsZTpcche->Fill(mom*dz,-1.*mom*dx);
		       hDedxvsMomche->Fill(mom/charge, dedx);
		     }
		     hDedxvsVehit->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),dedx);
		     hVTofvsMomc->Fill(mom/charge,(nlv_tcal0[j]+nlv_tcal1[j])/2.);
		   }
		 }
	       }
	       }
	     
	     }
	    
	     
	     tpc_cnt++;
	     fChainx->GetEvent(tpc_cnt);
	     
	     if (run != current_run) {
	       current_run = run;
	       printf("number of events  %d\n", current_event-event_offset);
	       event_offset = current_event;
	       printf("events match  run %d  offset %d  eventid %d\n", current_run, event_offset, eventid);
	     }
	     
	     current_event = eventid + event_offset;
	     
	   } // while - loop over tpc tracks in the same event

	   hTrmult->Fill(tpc_track_cnt);

	   //tpc_tracks[jentry+1] = tpc_track_cnt;

		      }
       } // VETO geometry cut
     } // loop over clusters 
   } else {  // clno==0 ##############  nothing in NeuLAND ###################
     
     for(int j=0;j<nlv_nbar;j++) {
       if (nlv_bar[j]==8) continue;
       if (nlv_tcal0[j]>0.&&nlv_tcal1[j]>0.) {
	 
	 hVTofvsEhitStop->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	 
	 //hVTdiffvsVbar->Fill(nlv_bar[j]+1,nlv_tcal0[j]-nlv_tcal1[j]);
       }
     }
   }
   
   for(int j=0;j<nlv_nbar;j++) {
     if (nlv_bar[j]==8) continue;
     
     if (vdoublescut->IsInside(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.))
       hVdoublesc->Fill(nlv_bar[j]+1,vdoubles[nlv_bar[j]+1]);
     
     hVdoubles->Fill(nlv_bar[j]+1,vdoubles[nlv_bar[j]+1]);
     
     if (vdoubles[nlv_bar[j]+1]>1) 
       hVTofvsEhitMulti->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
   }
      
   for (iBar=1;iBar<=num_bars;iBar++) {
     if (fired[iBar]==2) {
       if(vetoq==0) {
	 //hPlanecut[(iBar-1)/50+1]->Fill(xhit[iBar],yhit[iBar]);
	 if(fabs(yhit[iBar])<95.&&fabs(xhit[iBar])<115.) {
	   hTofvsEhitallVetoq->Fill(ehit[iBar],tofc[iBar]);
	   hTofvsEhitvq[(iBar-1)/50+1]->Fill(ehit[iBar],tof[iBar]);
	   hEvsBarVetoq->Fill(iBar,ehit[iBar]);
	   hTofvsBarVetoq->Fill(iBar,tof[iBar]);
	   
	   //if (tof[iBar]>10000./(ehit[iBar]+60.)&&tof[iBar]<100.) {
	   hXYcut[(iBar-1)/50+1]->Fill(xhit[iBar],yhit[iBar]);
	   for(int j=0;j<nlv_nbar;j++) {
	     if (nlv_bar[j]==8) continue;
	     //hVEhitvsEhit->Fill(ehit[iBar],sqrt(nlv_qcal0[j]*nlv_qcal1[j]));
	   }
	 }
       } else {
	 if (plmultH[0]==0)
	   hXYacut[(iBar-1)/50+1]->Fill(xhit[iBar],yhit[iBar]);
       }
       
       if(veto==0) {
	 if(fabs(yhit[iBar])<95.&&fabs(xhit[iBar])<115.) {
	   hTofvsEhitallVeto->Fill(ehit[iBar],tof[iBar]);
	   hTofvsEhitv[(iBar-1)/50+1]->Fill(ehit[iBar],tof[iBar]);
	   hVelvsBarVeto->Fill(iBar,v[iBar]);
	   hTofvsBarVeto->Fill(iBar,tof[iBar]);
	 }
       } else {
	 hTofvsEhitnv[(iBar-1)/50+1]->Fill(ehit[iBar],tof[iBar]);
       }
     }
   }
   
   hMult->Fill(mult);     
   hMultPart->Fill(multpart);     
   hMultPmt->Fill(multpmt);     
   hMultV->Fill(multV);     
   hMultH->Fill(multH); 

   hMultv->Fill(vmult);
   hMultvq->Fill(vqmult);

      hClusterMultAfter->Fill(clno);
   hClusterMultNeutrons->Fill(nclno);
   hClusterMultCharged->Fill(cclno);
   hClNovsClNo->Fill(nclno,cclno);
   
   hClvsTotEne->Fill(totdetene,nclno);

   hNmultvsEvnt->Fill(jentry,nclno);
   hChmultvsEvnt->Fill(jentry,cclno);
   
   //hNmultvsTrmult->Fill(tpc_tracks[jentry+1],nclno);
   //hCmultvsTrmult->Fill(tpc_tracks[jentry+1],cclno);
   //hVmultvsTrmult->Fill(tpc_tracks[jentry+1],vmult);

  
   // if(mult==10) {
   //   for (iBar=1;iBar<=num_bars;iBar++) {
   //     if (fired[iBar]==2&&v[iBar]<29.&&v[iBar]>0.&&evtcnt<20000) {
   // 	 hEvent[evtcnt]->Fill(xhit[iBar]/5.+30.,yhit[iBar]/5.+30.,tof[iBar]);
   // 	 hEvent[evtcnt]->Fill(xhit[iBar]/5.+30.,zhit[iBar]/5.+60.,tof[iBar]);
   // 	 hEvent[evtcnt]->Fill(zhit[iBar]/5.+60.,yhit[iBar]/5.+30.,tof[iBar]);
   // 	 hEvntXY[evtcnt]->Fill(xhit[iBar],yhit[iBar],tof[iBar]);
   // 	 hEvntXZ[evtcnt]->Fill(xhit[iBar],zhit[iBar],tof[iBar]);
   // 	 hEvntYZ[evtcnt]->Fill(zhit[iBar],yhit[iBar],tof[iBar]);
   //     }
   //   }
   //   hClnovsEvnt->Fill(evtcnt,clno);
   //   evtcnt++;
   // }

 } 
}

Double_t QDC(Double_t channel)
{
  Double_t value=20*(9.e-6*channel*channel+0.0285*channel);
  return value;
}

Double_t wlk(Double_t x)
{
  Double_t y=0;

  Double_t par1=1500.; // +-0.2238
  Double_t par2=0.00075;//+-2.355e-05


  y=par1*pow(x,par2)-(par1*pow(400.,par2)); // Michael's
  //y=2.29083*log(x)-0.0870157*log(x)*log(x)-4.57824;  // mine

  //Int_t nseg;
  // xwalk[0]=0.;
  // for(nseg=1;nseg<=6;nseg++) {
 
  //   if (x>xwalk[nseg-1]&&x<=xwalk[nseg]) {
  //     y = -1*(awalk[nseg]*x+bwalk[nseg]-35.7); //35.7  // needs correction 0.07 -> 35.77 to match powlaw
  //   }
  // }

  return y;
}

void display_events() {

  TCanvas *C_histograms1;
  C_histograms1 = new TCanvas("NeuLAND histograms 1","NeuLAND histograms 1",0,0,1280,1024);
  C_histograms1->Clear();

  for(int i=0;i<20000;i++) {
    hEvent[i]->Draw("colz");
    pause();
  }
}

Double_t traw_to_tcal(struct TCalTable const *const a_tcal, Int_t const
    a_channel)
{
  if (a_channel < a_tcal->table[0].raw) {
    Int_t raw0 = a_tcal->table[0].raw;
    Int_t d_raw = a_tcal->table[1].raw - raw0;
    Double_t ns0 = a_tcal->table[0].ns;
    Double_t d_ns = a_tcal->table[1].ns - ns0;
    return d_ns / d_raw * (a_channel - raw0) + ns0;
  }
  int last = a_tcal->table_len - 1;
  if (a_channel > a_tcal->table[last].raw) {
    Int_t raw0 = a_tcal->table[last].raw;
    Int_t d_raw = a_tcal->table[last - 1].raw - raw0;
    Double_t ns0 = a_tcal->table[last].ns;
    Double_t d_ns = a_tcal->table[last - 1].ns - ns0;
    return d_ns / d_raw * (a_channel - raw0) + ns0;
  }
  for (uint i = 0; a_tcal->table_len - 1 > i; ++i) {
    if (a_tcal->table[i + 1].raw >= a_channel) {
      Int_t raw0 = a_tcal->table[i].raw;
      Int_t d_raw = a_tcal->table[i + 1].raw - raw0;
      Double_t ns0 = a_tcal->table[i].ns;
      Double_t d_ns = a_tcal->table[i + 1].ns - ns0;
      return d_ns / d_raw * (a_channel - raw0) + ns0;
    }
  }
  cerr << "Uhm, this should never happen," << endl;
  abort();
}

Double_t nnp_traw_to_tcal(Int_t const a_bar_id, Int_t const a_t_type, Int_t
    const a_channel)
{
  return traw_to_tcal(&tcal_nnp[a_bar_id][a_t_type], a_channel);
}

Double_t start_traw_to_tcal(Int_t const a_t_type, Int_t const a_channel)
{
  return traw_to_tcal(&tcal_start[a_t_type], a_channel);
}

void read_calibration()
{
  int pl,bar,pmt,nseg;
  float ped1,ped2;
  string cline;
  //  ifstream a_file ( "neuland_tcal.txt" );
  //  if(!a_file){
  //    cout << "Cannot open tcal file. Exit." << endl;
  //    exit(1);
  //  }
  //  getline (a_file,line);

  ifstream ped_file ( "ped0437.txt" );  
  if(!ped_file){
    cout << "Cannot open pedestal file. Exit." << endl;
    exit(1);
  }

  for(Int_t i=1;i<=num_bars;i++){
    ped_file>>pl>>bar>>pmt>>ped1;
    ped_file>>pl>>bar>>pmt>>ped2;
    pedestal1[(pl-1)*50+bar]=ped1;
    pedestal2[(pl-1)*50+bar]=ped2;
  }
  ped_file.close();	 

  ifstream tcal_file("nnp_tcal_r0437ridf.hh");
  if (!tcal_file) {
    cerr << "Cannot open tcal file. Exit." << endl;
    exit(1);
  }
  for (int line_no = 1;; ++line_no) {
    char line[256];
    tcal_file.getline(line, sizeof line);
    if (!tcal_file) {
      break;
    }
#define EXPECT(str) do {\
    if (0 != strcmp(token, str)) {\
      cerr << line_no << ": Missing " str", got \"" << token << "\".\n" << endl;\
      exit(1);\
    }\
  } while (0)
    char *p, *token;
    token = strtok_r(line, " (,)", &p);
    EXPECT("TIME_CALIB_POINT");
    token = strtok_r(p, " (,)", &p);
    EXPECT("SIGNAL_ID");
    token = strtok_r(p, " (,)", &p);
    bool is_start = false;
    int plane=0, bar_id;
    if (0 == strcmp(token, "START")) {
      is_start = true;
    } else if (0 == strcmp(token, "NNP")) {
      token = strtok_r(p, " (,)", &p);
      plane = strtol(token, NULL, 10);
      token = strtok_r(p, " (,)", &p);
      bar = strtol(token, NULL, 10);
      bar_id = (plane - 1) * 50 + (bar - 1) + 1;
    } else {
      cerr << line_no << ": Invalid tcal detector " << token << ".\n" << endl;
      exit(1);
    }
    token = strtok_r(p, " (,)", &p);
    int t_type = strtol(token, NULL, 10);
    token = strtok_r(p, " (,)", &p);
    int raw = strtol(token, NULL, 10);
    token = strtok_r(p, " (,)", &p);
    float ns = strtod(token, NULL);
    TCalTable *tcal;
    if (is_start) {
      tcal = &tcal_start[t_type];
    } else {
      tcal = &tcal_nnp[bar_id][t_type];
    }
    int entry = tcal->table_len;
    if (ENTRY_TABLE_LEN <= entry) {
      cerr << line_no << ": Collected " << entry << " entries." << endl;
      if (is_start) {
        cerr << "Start, t_type=" << t_type << endl;
      } else {
        cerr << "NNP, t_type=" << t_type << ", plane=" << plane << ", bar=" <<
            bar << "." <<endl;
      }
      exit(1);
    }
    tcal->table[entry].raw = raw;
    tcal->table[entry].ns  = ns;
    ++tcal->table_len;
  }
  tcal_file.close();
  // TCal tests.
//  cout << raw_to_tcal(1, 1, 1181) << endl;
//  cout << raw_to_tcal(1, 1, 1182) << endl;
//  cout << raw_to_tcal(1, 1, 1183) << endl;
//  cout << raw_to_tcal(1, 1, 1185) << endl;
//  cout << raw_to_tcal(1, 1, 1186) << endl;
//  cout << raw_to_tcal(1, 1, 1187) << endl;
//  cout << raw_to_tcal(1, 1, 3285) << endl;
//  cout << raw_to_tcal(1, 1, 3286) << endl;
//  cout << raw_to_tcal(1, 1, 3287) << endl;
//  exit(1);

  ifstream b_file ( "parameter_test/neuland_tdiff3.txt");
  //ifstream b_file ( "parameter_test/neuland_sync_gamma.txt");  // was _gamma OK!!! ig 30.05.2017
  //ifstream b_file ( "parameter_test/neuland_dtfine.txt");
  //ifstream b_file ( "parameter_test/neuland_dtdtfine.txt");
  //ifstream b_file ( "parameter_test/neuland_sync_gamma2.txt");
  if(!b_file){
    cout << "Cannot open tdiff/tsync file. Exit." << endl;
    exit(1);
  }
  getline (b_file,cline);
  ifstream e_file ( "parameter_test/neuland_esyncXXX.txt");
  if(!e_file){
    cout << "Cannot open ediff/esync file. Exit." << endl;
    exit(1);
  }
  getline (e_file,cline);

  for(Int_t i=1;i<=num_bars;i++){
    b_file>>bar>>tdiff[i]>>tsync[i]>>vscint[i];
    e_file>>bar>>ediff[i]>>esync[i]>>att[i];

    //tdiff[i]=0.;
    //tsync[i]=0.;
    //vscint[i]=1.;
    //ediff[i]=1.;
    //esync[i]=1.;
    //cout<< bar <<"     " << tdiff[i] << "     " << esync[i] << "     " << vscint[i] << "\n";
    //printf("tdiff %d %lf\n", i, tdiff[i]);
  }   
  b_file.close();	 
  e_file.close();

  ifstream w_file ( "parameter_test/neuland_walk.txt");
  if(!w_file){
    cout << "Cannot open walk file. Exit." << endl;
    exit(1);
  }

  for(int i=1;i<=6;i++) {
    w_file>>nseg>>xwalk[i]>>awalk[i]>>bwalk[i];
  }
  w_file.close();
}

void out_hist() {

  TFile *fhist = new TFile("nltpc.root","recreate");

  hdyvsdx->Write();
  hdyvsdz->Write();
  hdxvsdz->Write();
  hdyvsdxc->Write();
  hdyvsdzc->Write();
  hdxvsdzc->Write();
  hdyvsdxn->Write();
  hdyvsdzn->Write();
  hdxvsdzn->Write();
  
  hYvsXTpc->Write();
  hYvsZTpc->Write();
  hXvsZTpc->Write();
  
  hYvsXTpcn->Write();
  hYvsZTpcn->Write();
  hXvsZTpcn->Write();
  hYvsXTpcc->Write();
  hYvsZTpcc->Write();
  hXvsZTpcc->Write();
 
  hdzvsMomc->Write();
  hdxvsMomc->Write();
  hdyvsMomc->Write();

  hdxvsX->Write();
  hdyvsY->Write();
  
}

void loc2glo(float x, float y, float z, float *xglo, float *yglo, float *zglo) {

  float xg, yg, zg;
  
  xg = Distance*sin(th0) + x*cos(th0) + z*sin(th0);
  yg = y;
  zg = Distance*cos(th0) - x*sin(th0) + z*cos(th0);

  xglo = &xg;
  yglo = &yg;
  zglo = &zg;
}

void run()
{

  Int_t runs[] = {2841, 2843, 2844, 2845, 2846, 2848, 2849, 2850, 2851, 2852, 2855, 
                  2856, 2857, 2858, 2859, 2860, 2861, 2875, 2877, 2878, 2879, 2880, 
                  2881, 2882, 2883, 2884, 2887, 2888, 2889, 2890, 2891, 2892, 2893, 
                  2894, 2896, 2898, 2899, 2900, 2901, 2902, 2903, 2904, 2905, 2907, 
                  2914, 2916, 2917, 2919, 2920, 2921, 2922, 2924, 2925, 2926, 2927, 
                  2929, 2930, 2931, 2932, 2933, 2934, 2935, 2936, 2939, 2940, 2941, 
                  2942, 2943, 2944, 2945, 2946, 2948, 2955, 2956, 2958, 2959, 2960, 
                  2961, 2962, 2964, 2965, 2966, 2968, 2969, 2970, 2971, 2972, 2973, 
                  2975, 2976, 2977, 2978, 2979, 2980, 2981, 2982, 2983, 2984, 2985, 
                  2986, 2988, 2989, 2990, 2991, 2992, 2993, 2997, 2999, 3000, 3002, 
                  3003, 3007, 3039};

  Int_t numData = sizeof(runs)/sizeof(Int_t);

  char name[180]; 
  TChain *fChain = new TChain("nl");
  TChain *fChainx = new TChain("dedx");

  for(int i=0;i<2000;i++) {
    run2run[i]=0;
  }
  
  int cnt=-1;

  opencuts();
  
  RunNumber=999;
  
  read_calibration();
  
  for (int i = 0;i<20;i++) {
    //for (int i = 0;i<numData;i++) {
    
    sprintf(name,"/mnt/spirit/analysis/changj/SpiRITROOT.JWfix/macros/data/analysisCode-All-noLayerCut-GC-DS-GiordanoCommentOut-bShift-By-JWfix/dedxROOT/dedxSn132-All-noLayerCut_%d.root",i);
    fChainx->Add(name);
    printf(" file: %s\n",name);

    sprintf(name,"/mnt/spirit/analysis/gasparic/rootfiles/SMDAQ%04d.root",runs[i]);
    fChain->Add(name);
    cnt++;
    printf(" file: %s\n",name);
    run2run[cnt]=runs[i];

  }    

  for(int i=2841;i<=2841;i++) {  //  EOS2 all  2839 - 3166

    //for(int i=2839;i<=3166;i++) {  //  EOS2 all  2839 - 3166

  //for(int i=2839;i<=3055;i++) {  //  EOS2 132Sn on 124Sn  2839 - 3058

  //for(int i=3056;i<=3166;i++) {  //  EOS2 124Sn on 112Sn 3059 - 3166 

  //for(int i=3000;i<=3166;i++) {

    if (i==3017) continue;     //   very little counts
    if (i>=3018&&i<=3019) continue; // less counts than usual
    if (i>=3020&&i<=3023) continue; // very little counts
    if (i>=3024&&i<=3025) continue; // less counts than usual
    if (i>=3026&&i<=3030) continue;  // very little counts
    if (i>=3031&&i<=3036) continue;  // less counts than usual
    if (i==3037) continue;  // less counts than usual


    if (i>=3104&&i<=3112) continue; // a little less counts
    if (i>=3113&&i<=3126) continue; // almost no events    100k but no neutrons
    if (i>=3127&&i<=3133) continue; // much less counts
    if (i>=3134&&i<=3137) continue; // less counts

    //sprintf(name,"/mnt/spirit/analysis/gasparic/rootfiles/SMDAQ%04d.root",i);
    //fChain->Add(name);
    //cnt++;
    //printf(" file: %s\n",name);
    //run2run[cnt]=i;
  }


  // 2839 - 3009 #########################      132Sn on 124Sn

  // 2862 - 2874  no runs 
  // 2908 - 2915  almost nothing  shift in Tncorr - Tv after this break
  // 2950 - 2954  no runs         shift in Tncorr - Tv after this break 
  // 2972 - 2996  no runs

  // 3010 - 3012
  // 3013         Tof changes withing the run
  // 3014

  // 3015 - 3016  no runs
  // 3017 - 3025  some changes ?
  // 3026 - 3037  no runs
  // 3038 - 3039  changes

  // 3040 - 3055  no runs


  // 3056 - 3058  almost no events (12919 events)
  // 3059 - 3104 #########################      124Sn on 112Sn

  // 3105 - 3137  no runs          shift in Tncorr - Tv after this break

  // 3138 - 3166  OK
  
  loop(fChain, fChainx);
  //  out_hist();

  delete fChain;
  delete fChainx;
  
  printf("kraaaaj\n");

}

