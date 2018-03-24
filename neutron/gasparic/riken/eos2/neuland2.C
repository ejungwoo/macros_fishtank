#include <cmath>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <iostream>

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
#include        <algorithm>  // for sort 
#include        "NeulandCluster.h"

using namespace std;

//TFile *nfile= new TFile("nmult.root","RECREATE");
//TTree *ntree= new TTree("ntree","Tree with neutron multiplicities");
Int_t nrun,nevt,nevttot,nmult;

Int_t katmult[200000];
Int_t kyomult[200000];
Int_t ampmult[200000];

Int_t RunNumber=999;
Int_t fDisplayTrack=0;

Int_t run2run[3000];

int Nnekosan=0;

Double_t nnp_traw_to_tcal(Int_t,Int_t,Int_t);
Double_t start_traw_to_tcal(Int_t,Int_t);

//Double_t T0= 622.0;
//Double_t T0= 623.57-61.5+37.666;
Double_t T0 = 750.; // 745

const Double_t Distance = 844.8;  // EOS experiment
const Double_t vDistance = 812.3;  // EOS experiment distance 32.5 cm

//const Double_t c=29.979245; //cm/ns
const Double_t c=16.667; //cm/ns

const int num_planes=4;
const int num_bars=400;

const Int_t EhitRange=400; 

Double_t pedestal1[num_bars],pedestal2[num_bars];

Double_t vped1[8]={0.82,0.57,0.42,0.87,0.76,0.68,1.41,1.35};
Double_t vped2[8]={1.00,1.11,0.82,1.18,1.23,0.78,0.78,0.58};
Double_t vtoff1[8]={168.,164.,163.,168.,163.,165.,171.,117.};
Double_t vtoff2[8]={126.,126.,117.,121.,172.,172.,160.,168.};

Double_t vtdiff[8]={-0.16,1.33,-0.29,1.79,-0.02,4.78,-0.63,4.96};

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

Double_t xwalk[10],awalk[10],bwalk[10];

Double_t wlk(Double_t x);
Double_t QDC(Double_t channel);

TH1F* hMult = new TH1F("hMult","Multiplicity",100,-0.5,99.5);
TH1F* hMultPart = new TH1F("hMultPart","Multiplicity partial hits",100,-0.5,99.5);
TH1F* hMultPmt = new TH1F("hMultPmt","Multiplicity pmt",100,-0.5,99.5);
TH1F* hMultV = new TH1F("hMultV","Multiplicity vertical plane",100,-0.5,99.5);
TH1F* hMultH = new TH1F("hMultH","Multiplicity horizontal plane",100,-0.5,99.5);

TH1F* hMultv = new TH1F("hMultv","Veto time Multiplicity",10,-0.5,9.5);
TH1F* hMultvq = new TH1F("hMultvq","Veto charge Multiplicity",10,-0.5,9.5);
TH2F* hVdoubles = new TH2F("hVdoubles","Veto mult in bar vs bar",8,0.5,8.5,10,-0.5,9.5);

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

//TH2F* hT3vsEvnt[num_bars+1];
//TH2F* hT4vsEvnt[num_bars+1];
TH2F* hPed1vsEvnt[num_bars+1];
TH2F* hPed2vsEvnt[num_bars+1];

TH2F* hEhitvsPos[num_bars+1];
TH2F* hDT[num_bars+1];

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
TH2F* hTofvsRun = new TH2F("hTofvsRun","Tof vs run ",4000,0.5,4000.5,3000,0.,500.);

TH2F* hVTofvsEhit = new TH2F("hVTofvsEhit","Veto Tof vs Ehit",500,0.,150.,500,0.,400.);
TH2F* hVTofvsEhitStop = new TH2F("hVTofvsEhitStop","Veto Tof vs Ehit stopped",500,0.,150.,500,0.,400.);
TH2F* hVTofvsEhitMatch = new TH2F("hVTofvsEhitMatch","Veto Tof vs Ehit match",500,0.,150.,500,0.,400.);
TH2F* hVTofvsEhitNoMatch = new TH2F("hVTofvsEhitNoMatch","Veto Tof vs Ehit no match",500,0.,150.,500,0.,400.);
TH2F* hVTofvsTof = new TH2F("hVTofvsTof","Veto Tof vs Tof",500,0.,400.,500,0.,400.);
TH2F* hVEhitvsEhit = new TH2F("hVEhitvsEhit","Veto E vs E",500,0.,EhitRange,500,0.,50.);
TH2F* hVBarvsBar = new TH2F("hVBarvsBar","Veto Bars vs Bars",num_bars,0.5,num_bars+0.5,8,0.5,8.5);
TH2F* hVBarvsX = new TH2F("hVBarvsX","Veto Bars vs X",1000,-170.,170.,8,0.5,8.5);
TH2F* hVBarvsXcorr = new TH2F("hVBarvsXcorr","Veto Bars vs X corr",1000,-170.,170.,8,0.5,8.5);
TH2F* hVXvsX = new TH2F("hVVvsX","Veto X vs X",1000,-170.,170.,14,-175.,175.);

TH2F* hDYvsY = new TH2F("hDYvsY","Veto DY vs Y",1000,-170.,170.,500,-50.,50.);
TH2F* hDXvsX = new TH2F("hDXvsX","Veto DX vs X",1000,-170.,170.,50,-50.,50.);
TH2F* hDYvsZ = new TH2F("hDYvsZ","Veto DY vs Z",200,0.,40.,500,-50.,50.);
TH2F* hDXvsZ = new TH2F("hDXvsZ","Veto Dx vs Z",200,0.,40.,50,-50.,50.);

TH2F* hVTdiffvsYall = new TH2F("hVTdiffvsYall","Veto Tdiff vs Y",1000,-170.,170.,300,-30.,30.);
TH2F* hVetoEraw = new TH2F("hVetoEraw", "Veto Eraw",8,0.5,8.5,500,0,6000);
TH2F* hVetoEcal = new TH2F("hVetoEcal", "Veto Ecal",8,0.5,8.5,500,0,50);
TH2F* hVetoEcal1 = new TH2F("hVetoEcal1", "Veto Ecal1",8,0.5,8.5,500,0,50);
TH2F* hVetoEcal2 = new TH2F("hVetoEcal2", "Veto Ecal2",8,0.5,8.5,500,0,50);
TH2F* hVTdiffvsY[9];
TH2F* hVTdiffEcal1[9];
TH2F* hVTdiffEcal2[9];
TH2F* hEvsBarVeto = new TH2F("hEvsBarVeto","Energy vs Bars Veto",num_bars,0.5,num_bars+0.5,2000,0,200);
TH2F* hEvsBarVetoq = new TH2F("hEvsBarVetoq","Energy vs Bars Vetoq",num_bars,0.5,num_bars+0.5,2000,0,200);
TH2F* hVetovsRun = new TH2F("hVetovsRun","Veto Edep vs run ",4000,0.5,4000.5,3000,0.,50.);
TH2F* hVetoTofvsRun = new TH2F("hVetoTofvsRun","Veto Tof vs run ",4000,0.5,4000.5,3000,0.,300.);

TH2F* hTdiffvsBar = new TH2F("hTdiffvsBar","Tdiff vs Bars",num_bars,0.5,num_bars+0.5,8000,-400,400);
TH2F* hEdiffvsBar = new TH2F("hEdiffvsBar","Ediff vs Bars",num_bars,0.5,num_bars+0.5,1000,-10.,10.);
TH1F* hTdiff[num_bars+1];
TH1F* hEdiff[num_bars+1];
TH1F* hEcosm[num_bars+1];
TH1F* hE1[num_bars+1];
TH1F* hE2[num_bars+1];
TH1F* hE1cut[num_bars+1];
TH1F* hE2cut[num_bars+1];
TH1F* hTofbar[num_bars+1];
TH1F* hTofplane[9];

TH2F* htt1_tcal1 = new TH2F("htt1_tcal1","tt1 - tcal1",num_bars,0.5,num_bars+0.5,1000,-10,10);
TH2F* htt2_tcal2 = new TH2F("htt2_tcal2","tt2 - tcal2",num_bars,0.5,num_bars+0.5,1000,-10,10);

TH2F* hE1vsPos[num_bars+1];
TH2F* hE2vsPos[num_bars+1];

TH2F* hTdiffvsXBar[num_bars+1];
TProfile* hTdiffvsXBarprof[num_bars+1];
TH2F* hEdiffvsXBar[num_bars+1];
TProfile* hEdiffvsXBarprof[num_bars+1];

TH1F* hTsync[num_bars+1][num_bars+1];
TH1F* hTsync_corr[num_bars+1];
TH1F* hTsync_corrH[num_planes*2+1];
TH1F* hTsync_corrV[num_planes*2+1];
TH1F* hTsync_corr175_225 = new TH1F("hTsync_corr175_225","hTsync_corr175_225",2000,-40.,40.);

TH2F* hTsync2D[num_bars+1];
TH2F* hTsync2D_corr = new TH2F("hTsync2D_corr","Tsync corr vs Bars",num_bars,0.5,num_bars+0.5,500,-50,50);
TH2F* hTdiffvsEvnt[num_bars+1];
TH2F* hCosmall  = new TH2F("hCosmall","Thit vs Path",1000,0,500,5000,-150,150);
TProfile* hCosmallprof = new TProfile("hCosmallprof","Thit vs Path",250,0,500,-100,100);
TH2F* hBarvsThit[num_planes+1][num_planes+1];

TH2F* hAvsEdepV;//  = new TH2F("hAvsEdepV","A vs Edep V",500,0,1000,1000,0,10);
TH2F* hAvsEdepH;//  = new TH2F("hAvsEdepH","A vs Edep H",500,0,1000,1000,0,10);
TH2F* hEdepvsPathV;//  = new TH2F("hEdepvsPathV","Edep vs Path V",1000,0,10,500,0,100);
TH2F* hEdepvsPathH;//  = new TH2F("hEdepvsPathH","Edep vs Path H",1000,0,10,500,0,100);
TProfile* hEdepvsPathVprof;//  = new TProfile("hEdepvsPathVprof","Edep vs Path V prof",200,0,10,0,100);
TProfile* hEdepvsPathHprof;//  = new TProfile("hEdepvsPathHprof","Edep vs Path H prof",200,0,10,0,100);

TH1F* hClusterNo=new TH1F("hClusterNo","Cluster multiplicity",40,-0.5,39.5);
TH1F* hClusterDim=new TH1F("hClusterDim","Cluster dimension",50,0.,50.);
TH2F* hTofvsEhitClust = new TH2F("hTofvsEhitClust","Tof vs Ehit clusters",500,0.,1000.,3000,0.,300.);
TH2F* hTofvsEhitClustLast = new TH2F("hTofvsEhitClustLast","Tof vs Ehit clusters last",1500,0.,500.,3000,0.,300.);
TH2F* hFirstZ = new TH2F("hFirstZ","Tof vs Z clusters",500,0.,40.,3000,0.,300.);
TH2F* hTofvsEhitClustVeto = new TH2F("hTofvsEhitClustVeto","Tof vs Ehit neutrons",500,0.,1000.,3000,0.,300.);
TH2F* hTofvsEhitClustNoVeto = new TH2F("hTofvsEhitClustNoVeto","Tof vs Ehit charged",500,0.,1000.,3000,0.,300.);
TH2F* hTofvsBarClustVeto = new TH2F("hTofvsBarClustVeto","Tof vs Bar clusters",num_bars,0.5,num_bars+0.5,3000,0.,300.);
TH2F* hEvsBarClustVeto = new TH2F("hEvsBarClustVeto","E vs Bar neutrons",num_bars,0.5,num_bars+0.5,500,0.,1000.);
TH2F* hEvsBarClustNoVeto = new TH2F("hEvsBarClustNoVeto","E vs Bar charged",num_bars,0.5,num_bars+0.5,500,0.,1000.);
TH2F* hTofvsEhitClustv[num_planes*2+1];

TH1F* hTdiffClust = new TH1F("hTdiffClust","tdiff in cluster",1000,-50,50);

TH1F* hClusterNoNeutrons=new TH1F("hClusterNoNeutrons","Cluster multiplicity",40,-0.5,39.5);
TH1F* hClusterNoCharged=new TH1F("hClusterNoCharged","Cluster multiplicity",40,-0.5,39.5);
TH2F* hClNovsClNo=new TH2F("hClNovsClNo","Cluster mult veto vs neuland",40,-0.5,39.5,40,-0.5,39.5);

TH2F* hNLvsKA=new TH2F("hNLvsKA","Multiplicity vs KA mult",15,-0.5,14.5,40,-0.5,39.5);
TH2F* hNLvsKY=new TH2F("hNLvsKY","Multiplicity vs KY mult",35,-0.5,34.5,40,-0.5,39.5);
TH2F* hKYvsKA=new TH2F("hKYvsKA","KY vs KA mult",15,-0.5,14.5,35,-0.5,34.5);
TH2F* hNLvsKAmix=new TH2F("hNLvsKAmix","Multiplicity vs KA mult mix",15,-0.5,14.5,40,-0.5,39.5);
TH2F* hNLvsKYmix=new TH2F("hNLvsKYmix","Multiplicity vs KY mult mix",35,-0.5,34.5,40,-0.5,39.5);
TH2F* hKYvsKAmix=new TH2F("hKYvsKAmix","KY vs KA mult mix",15,-0.5,14.5,35,-0.5,34.5);
TH2F* hNLvsKAmix1=new TH2F("hNLvsKAmix1","Multiplicity vs KA mult mix",15,-0.5,14.5,40,-0.5,39.5);
TH2F* hNLvsKYmix1=new TH2F("hNLvsKYmix1","Multiplicity vs KY mult mix",35,-0.5,34.5,40,-0.5,39.5);
TH2F* hKYvsKAmix1=new TH2F("hKYvsKAmix1","KY vs KA mult mix",15,-0.5,14.5,35,-0.5,34.5);
TH2F* hNLvsKYdiff=new TH2F("hNLvsKYdiff","Multiplicity vs KY mult diff",35,-0.5,34.5,40,-0.5,39.5);
TH2F* hNLvsKAdiff=new TH2F("hNLvsKAdiff","Multiplicity vs KA mult diff",15,-0.5,14.5,40,-0.5,39.5);
TH2F* hKYvsKAdiff=new TH2F("hKYvsKAdiff","KY vs KA mult diff",15,-0.5,14.5,35,-0.5,34.5);
TH2F* hNLvsKYdiff1=new TH2F("hNLvsKYdiff1","Multiplicity vs KY mult diff",35,-0.5,34.5,40,-0.5,39.5);
TH2F* hNLvsKAdiff1=new TH2F("hNLvsKAdiff1","Multiplicity vs KA mult diff",15,-0.5,14.5,40,-0.5,39.5);
TH2F* hKYvsKAdiff1=new TH2F("hKYvsKAdiff1","KY vs KA mult diff",15,-0.5,14.5,35,-0.5,34.5);

TH1F* hTofdiff=new TH1F("hTofdiff","Tn - Tv",4000,-200,200);
TH1F* hTofdiffn=new TH1F("hTofdiffn","Tn - Tv",4000,-200,200);
TH2F* hTofdiffvsEcal=new TH2F("hTofdiffvsEcal","Tn - Tv",500,0,50,1000,-200,200);
TH2F* hTofdiffvsEhit=new TH2F("hTofdiffvsEhit","Tn - Tv",1000,0,500,1000,-200,200);

TH2F* hPlane[2*num_planes+1];
TH2F* hPlanecut[2*num_planes+1];
TH2F* hPlaneacut[2*num_planes+1];

TH2F* hXY = new TH2F("hXY","XY hits",170,-170,170,170,-170,170); 
TH2F* hXYclust = new TH2F("hXYclust","XY hits cluster",170,-170,170,170,-170,170); 
TH2F* hXZclust = new TH2F("hXZclust","XY hits cluster",8,0.,40.,170,-170,170); 
TH2F* hYZclust = new TH2F("hYZclust","XY hits cluster",8,0.,40.,170,-170,170); 

TH2F* hDTvsE[num_bars+1];
TH2F* hDTvsEall = new TH2F("hDTvsEall","Delta T vs E all",300,0.,500.,1000,-80,80);

TH1F* hNstart = new TH1F("hNstart","Multiplicity start",20,-0.5,19.5);
TH1F* hTstart = new TH1F("hTstart","Time start",2000,-42000,42000);

// cluster
TH1F* hnCluster = new TH1F("hnCluster","nCluster",40,-0.5,39.5);
TH2F* hCausality = new TH2F("hCausality", "Distance vs Delta T", 300,0.,15., 400,0.,200.);
TH2F* hCausalityGamma = new TH2F("hCausalityGamma", "Distance vs Delta T", 300,0.,15., 400,0.,200.);
TH1F* hVelCluster = new TH1F("hVelCluster","Velocity between clusters",400,0,40.);
TH1F* hVelClusterGamma = new TH1F("hVelClusterGamma","Velocity between clusters",400,0,40.);
TH2F* hCVelvsDist = new TH2F("hCVelvsDist", "ClusterVelocity vs Distance", 400,0.,200., 200,0.,200.);
TH2F* hCVelvsDistGamma = new TH2F("hCVelvsDistGamma", "ClusterVelocity vs Distance", 400,0.,200., 200,0.,200.);

static const Int_t nMaxClust=num_bars;
NeulandCluster clust[nMaxClust];

void loop(TChain *fChain)
{
  char name1[60]; 
  char name2[60];

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
  }
  for(int i=1;i<=num_planes*2;i++){
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
      sprintf(name1,"hTofvsEhitClustv%d",i);
      sprintf(name2,"Tof vs Ehit clusters %d",i);     
      hTofvsEhitClustv[i] = new TH2F(name1,name2, 200,0.,EhitRange,3000,0.,300.); 
  }
  for(int i=1;i<=2*num_planes;i++){
      sprintf(name1,"plane%i",i);
      sprintf(name2,"Xhit vs Yhit %i",i);     
      hPlane[i] = new TH2F(name1,name2,170,-170,170,170,-170,170); 
  }
  for(int i=1;i<=2*num_planes;i++){
      sprintf(name1,"plane%icut",i);
      sprintf(name2,"Xhit vs Yhit %i",i);     
      hPlanecut[i] = new TH2F(name1,name2,170,-170,170,170,-170,170); 
  }
  for(int i=1;i<=2*num_planes;i++){
      sprintf(name1,"plane%iacut",i);
      sprintf(name2,"Xhit vs Yhit %i",i);     
      hPlaneacut[i] = new TH2F(name1,name2,170,-170,170,170,-170,170); 
  }
  // for(int i=0;i<=num_bars;i++){
  //   sprintf(name1,"TdiffBar%i",i);
  //   sprintf(name2,"Time diff of bar %i",i);     
  //   hTdiff[i] = new TH1F(name1,name2,1000,-100.,100.); // was 2000 bins
  // }
  // for(int i=0;i<=num_bars;i++){
  //   sprintf(name1,"EdiffBar%i",i);
  //   sprintf(name2,"Energy diff of bar %i",i);     
  //   hEdiff[i] = new TH1F(name1,name2,1000,-10.,10.); // was 2000 bins
  // }
  for(int i=0;i<=num_bars;i++){
    sprintf(name1,"Ecosm%i",i);
    sprintf(name2,"Energy of bar %i cosm",i);     
    hEcosm[i] = new TH1F(name1,name2,1000,0.,2000.); // was 2000 bins
  }
  for(int i=0;i<=num_bars;i++){
    sprintf(name1,"E1cosm%i",i);
    sprintf(name2,"Energy1 of bar %i cosm",i);     
    hE1[i] = new TH1F(name1,name2,500,0.,1000.); // was 2000 bins
    sprintf(name1,"E1cut%i",i);
    sprintf(name2,"Energy1 of bar %i cosm cut",i);     
    hE1cut[i] = new TH1F(name1,name2,1000,0.,1000.); // was 2000 bins
  }
  for(int i=0;i<=num_bars;i++){
    sprintf(name1,"Tofbar%i",i);
    sprintf(name2,"Gamma tof of bar %i",i);     
    hTofbar[i] = new TH1F(name1,name2,100,26.,32.); // was 2000 bins
  }
  for(int i=0;i<=8;i++){
    sprintf(name1,"Tofplane%i",i);
    sprintf(name2,"Gamma tof of plane %i",i);     
    hTofplane[i] = new TH1F(name1,name2,100,26.,32.); // was 2000 bins
  }
  for(int i=0;i<=num_bars;i++){
    sprintf(name1,"E2cosm%i",i);
    sprintf(name2,"Energy2 of bar %i cosm",i);     
    hE2[i] = new TH1F(name1,name2,500,0.,1000.); // was 2000 bins
    sprintf(name1,"E2cut%i",i);
    sprintf(name2,"Energy2 of bar %i cosm cut",i);     
    hE2cut[i] = new TH1F(name1,name2,1000,0.,1000.); // was 2000 bins
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
    hVTdiffvsY[i] = new TH2F(name1,name2,500,-170.,170.,200,-30,30);
  }

  for(int i=0;i<=num_bars;i++){
    sprintf(name1,"EhitvsPos%i",i);
    sprintf(name2,"Ehit vs position %i",i);     
    hEhitvsPos[i] = new TH2F(name1,name2,200,-170,170,200,0,EhitRange); // was 2000 bins
  }
  for(int i=0;i<=num_bars;i++){
     sprintf(name1,"E1vsPos%i",i);
     sprintf(name2,"Energy1 of bar %i vs Position",i);     
     hE1vsPos[i] = new TH2F(name1,name2,150,-200,200,1000,0.,1000.);
   }
   for(int i=0;i<=num_bars;i++){
     sprintf(name1,"E2vsPos%i",i);
     sprintf(name2,"Energy2 of bar %i vs Position",i);     
     hE2vsPos[i] = new TH2F(name1,name2,150,-200,200,1000,0.,1000.);
   }
  // for(int i=0;i<=num_bars;i++){
  //     sprintf(name1,"E1vsPosMult%i",i);
  //     sprintf(name2,"Energy1 of bar %i vs Position Mult",i);     
  //     hE1vsPosMult[i] = new TH2F(name1,name2,150,-200,200,1000,0.,1000.);
  //   }
  //   for(int i=0;i<=num_bars;i++){
  //     sprintf(name1,"E2vsPosMult%i",i);
  //     sprintf(name2,"Energy2 of bar %i vs Position Mult",i);     
  //     hE2vsPosMult[i] = new TH2F(name1,name2,150,-200,200,1000,0.,1000.);
  //   }

  for(int i=0;i<=num_bars;i++){
    for(int j=0;j<=num_bars;j++){
      sprintf(name1,"Tsync%ito%i",i,j);
      sprintf(name2,"Time sync of bar %i to bar %i",i,j);     
      hTsync[i][j] = new TH1F(name1,name2,1000,-50.,50.); // was 2000 bins
    }
  }
  for(int i=0;i<=num_bars;i++){
    sprintf(name1,"Tsync_corr%i",i);
    sprintf(name2,"Time sync corrected of bar %i",i);     
    hTsync_corr[i] = new TH1F(name1,name2,500,-50.,50.); // was 2000 bins
  }
  for(int i=0;i<=num_planes*2;i++){
    sprintf(name1,"Tsync_corrH%i",i);
    sprintf(name2,"Time sync corrected of H plane %i",i);     
    hTsync_corrH[i] = new TH1F(name1,name2,2000,-40.,40.); // was 2000 bins
    sprintf(name1,"Tsync_corrV%i",i);
    sprintf(name2,"Time sync corrected of V plane %i",i);     
    hTsync_corrV[i] = new TH1F(name1,name2,2000,-40.,40.); // was 2000 bins
  }

  for(int i=0;i<=num_bars;i++){
    sprintf(name1,"TdiffvsXBar%i",i);
    sprintf(name2,"Tdiff vs X bar %i",i);     
    hTdiffvsXBar[i] = new TH2F(name1,name2,50,0.5,50.5,3000,-300,300);
    sprintf(name1,"TdiffvsXBarprof%i",i);
    sprintf(name2,"Tdiff vs X bar prof %i",i);     
    hTdiffvsXBarprof[i] = new TProfile(name1,name2,50,0.5,50.5,-300,300);
    sprintf(name1,"EdiffvsXBar%i",i);
    sprintf(name2,"Ediff vs X bar %i",i);     
    hEdiffvsXBar[i] = new TH2F(name1,name2,50,0.5,50.5,500,-5,5);
    sprintf(name1,"EdiffvsXBarprof%i",i);
    sprintf(name2,"Ediff vs X bar prof %i",i);     
    hEdiffvsXBarprof[i] = new TProfile(name1,name2,50,0.5,50.5,-5,5);
  }
  for(int i=0;i<=num_bars;i++){
    sprintf(name1,"Tsync2DBar%i",i);
    sprintf(name2,"Time sync of bar %i",i);     
    hTsync2D[i] = new TH2F(name1,name2,50,0.5,50.5,500,-50,50); // was 2000 bins
  }
  for(int i=0;i<=num_bars;i++){
    sprintf(name1,"DT%i",i);
    sprintf(name2,"DT vs bar %i",i);     
    hDT[i] = new TH2F(name1,name2,400,0.5,400.5,500,-10.,10.); // was 2000 bins
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
 
 int cosm_up, cosm_down;

 Double_t dist = -1.e20;
 Double_t tstart;
 Int_t nstart;
 Double_t s, v[num_bars+1], beta;
 Double_t dt;

 // char name[20];

 TRandom3 rnd;
 
 Double_t a,b,path,sx, sy, sz, sx2, sy2, sz2, sxy, syz;
 int n;
 Double_t ah,bh,pathh,sxh, syh, sx2h, sy2h, sxyh;
 int nh;
 
 Long64_t nentries = fChain->GetEntries();
 //nentries = 1000000;
 
 //TH2F* hBarvsEvnt = new TH2F("hBarvsEvnt","Bar vs event",1000,0,nentries,num_bars,0.5,num_bars+0.5);
 //TH2F* hBarvsEvnt2 = new TH2F("hBarvsEvnt2","Bar vs event 2",1000,0,nentries,num_bars,0.5,num_bars+0.5);
 TH2F* hTofvsEvnt = new TH2F("hTofvsEvnt","Tof vs event ",1000,0,nentries,5000,0.,1000.);
 TH2F* hEvntvsRun = new TH2F("hEvntvsRun","Evnt vs Run ",4000,0.5,4000.5,1000,0,nentries);
 TH2F* hVetovsEvnt = new TH2F("hVetovsEvnt","Veto Edep vs event ",1000,0,nentries,3000,0.,50.);
 TH2F* hVetoTofvsEvnt = new TH2F("hVetoTofvsEvnt","Veto Tof vs event ",1000,0,nentries,3000,0.,300.);
 TH2F* hTstartvsEvnt = new TH2F("hTstartvsEvnt","Tstart vs event ",1000,0,nentries,2000,0.,2000.);

 TH2F* hClmultvsEvnt = new TH2F("hClmultvsEvnt","Clmult vs event ",1000,0,nentries,50,-0.5,49.5);
 TH2F* hNmultvsEvnt = new TH2F("hNmultvsEvnt","Nmult vs event ",1000,0,nentries,50,-0.5,49.5);
 TH2F* hChmultvsEvnt = new TH2F("hChmultvsEvnt","Chmult vs event ",1000,0,nentries,50,-0.5,49.5);

 // for(int i=0;i<num_bars+1;i++){
 //   sprintf(name1,"TdiffvsEvntc%i",i);
 //   sprintf(name2,"Tdiff of tac channel %i",i);     
 //   hTdiffvsEvnt[i] = new TH2F(name1,name2,100,0,nentries,200,-300.,300.);
 // }
 
 // for(int i=0;i<num_bars+1;i++){
 //   sprintf(name1,"T3vsEvntc%i",i);
 //   sprintf(name2,"T3 of tac channel %i",i);     
 //   hT3vsEvnt[i] = new TH2F(name1,name2,100,0,nentries,500,0.,2000.);
 //   sprintf(name1,"T4vsEvntc%i",i);
 //   sprintf(name2,"T4 of tac channel %i",i);     
 //   hT4vsEvnt[i] = new TH2F(name1,name2,1000,0,nentries,100,0.,2000.);
 // }
 for(int i=0;i<num_bars+1;i++){
   sprintf(name1,"hPed1vsEvntc%i",i);
   sprintf(name2,"Ped1 of tac channel %i",i);     
   hPed1vsEvnt[i] = new TH2F(name1,name2,100,0,nentries,200,-10.,200.);
 }
 for(int i=0;i<num_bars+1;i++){
   sprintf(name1,"hPed2vsEvntc%i",i);
   sprintf(name2,"Ped2 of tac channel %i",i);     
   hPed2vsEvnt[i] = new TH2F(name1,name2,100,0,nentries,200,-10.,200.);
 }

 Int_t n2n=0;

 Int_t runold=0;
 Int_t evt_cnt=0;

 char katname[100];
 char kyoname[100];
 char ampname[100];

 for (Long64_t jentry=0; jentry<nentries;jentry++) {

   fChain->GetEvent(jentry);
   //printf("----------------\n");	  
   if ((float(jentry)/20000.)==int(jentry/20000)) {
     cout << "event: " << jentry << " of " << nentries << endl;
   }

   RunNumber = run2run[fChain->GetTreeNumber()];

   evt_cnt++;
   if(RunNumber!=runold) {
     evt_cnt=1;
     runold=RunNumber;

//      sprintf(katname,"lukasik/run_%04d.trigbit.root",RunNumber);
//      for(Long64_t ikat=0;ikat<200000;ikat++) katmult[ikat]=0;
//      TChain *tkat= new TChain("kat");
//      tkat->Add(katname);
//      Int_t katm;
//      TBranch *b_katm;
//      tkat->SetBranchAddress("nhit",&katm,&b_katm);
//      Int_t kate;
//      TBranch *b_kate;
//      tkat->SetBranchAddress("evt",&kate,&b_kate);
//      Long64_t nkat = tkat->GetEntries();
//      for(Long64_t ikat=0;ikat<nkat;ikat++) {
//        tkat->GetEvent(ikat);
//        katmult[ikat+1]=katm;
//        if(ikat!=kate) printf("mismatch katana\n");
//      }

//      sprintf(kyoname,"lukasik/run_%04d.kyoto.root",RunNumber);
//      for(Long64_t ikyo=0;ikyo<200000;ikyo++) kyomult[ikyo]=0;
//      TChain *tkyo= new TChain("tka");
//      tkyo->Add(kyoname);
//      Int_t kyom;
//      TBranch *b_kyom;
//      tkyo->SetBranchAddress("nhit",&kyom,&b_kyom);
//      Int_t kyoe;
//      TBranch *b_kyoe;
//      tkyo->SetBranchAddress("eveid",&kyoe,&b_kyoe);
//      Long64_t nkyo = tkyo->GetEntries();
//      for(Long64_t ikyo=0;ikyo<nkyo;ikyo++) {
//        tkyo->GetEvent(ikyo);
//        kyomult[ikyo+1]=kyom;
//        if(ikyo!=kyoe) printf("mismatch kyoto\n");
//      }

     //sprintf(ampname,"lukasik/run%04d.katamp.root",RunNumber);
     // for(Long64_t iamp=0;iamp<200000;iamp++) ampmult[iamp]=0;
     // TChain *tamp= new TChain("ntree");
     // tamp->Add(ampname);
     // Int_t ampe;
     // TBranch *b_ampe;
     // tamp->SetBranchAddress("nevt",&ampe,&b_ampe);
     // Float_t ampm;
     // TBranch *b_ampm;
     // tamp->SetBranchAddress("namax",&ampm,&b_ampm);
     // Long64_t namp = tamp->GetEntries();
     // for(Long64_t iamp=0;iamp<namp;iamp++) {
     //   tamp->GetEvent(iamp);
     //   ampmult[iamp+1]=(Int_t)(ampm*0.025);
     //   //if(iamp!=ampe) 
     //   printf("mismatch amp  %d    %d    %f\n",iamp, ampe, ampm);
     // }
   }

   if (RunNumber<10) {
     T0 = 756.5-7.;
   } else if (RunNumber<1736) {
     T0 = 756.5;
   } else if (RunNumber<1745) {
     T0 = 756.5+40.-185.;
   } else if (RunNumber<1747) {
     T0 = 756.5;
   } else if (RunNumber<1781) {
     T0 = 756.5+40.-185.;
   } else if (RunNumber<1785) {
     T0 = 756.5;
   } else if (RunNumber<1787) {
     T0 = 756.5-5.;
   } else if (RunNumber<1801) {
     T0 = 756.5;
   } else if (RunNumber<1809) {
     T0 = 756.5+40.-185.;
   } else if (RunNumber<1811) {
     T0 = 756.5-5.;
   } else if (RunNumber<1818) {
     T0 = 756.5+40.-185.;
   } else if (RunNumber<1847) {
     T0 = 756.5-5.;
   } else if (RunNumber<1882) {
     T0 = 756.5+40.-270.;
   } else if (RunNumber<2180) {
     T0 = 756.5+40.-340.;
   } else if (RunNumber<2210) {
     T0 = 360.2;
   } else if (RunNumber<2839) {
     T0 = 560.2;
   } else if (RunNumber<3116) {
     T0 = 584.5;
   } else if (RunNumber<3137) {
     T0 = 607.5;
   } else {
     T0 = 582.7;
   }

   //T0 = 756.5;
 
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
   nstart = 0;

   hTstartvsEvnt->Fill(jentry,tstart);

   Int_t vdoubles[9];

   for (iHit=0;iHit<nlv_nbar;iHit++) {
     if (nlv_bar[iHit]<8) {
       if (RunNumber<1906) {
	 nlv_tcal0[iHit] = 4.*nlv_tcal0[iHit] - 120.;
	 nlv_tcal1[iHit] = 4.*nlv_tcal1[iHit] - 120.;
       }
       nlv_qcal0[iHit] = nlv_qcal0[iHit] - vped1[nlv_bar[iHit]];
       nlv_qcal1[iHit] = nlv_qcal1[iHit] - vped2[nlv_bar[iHit]];

       nlv_tcal0[iHit] = nlv_tcal0[iHit] - vtoff1[nlv_bar[iHit]]+66.-vtdiff[nlv_bar[iHit]]/2.;
       nlv_tcal1[iHit] = nlv_tcal1[iHit] - vtoff2[nlv_bar[iHit]]+66.+vtdiff[nlv_bar[iHit]]/2.;

       hVTofvsEhit->Fill(sqrt(nlv_qcal0[iHit]*nlv_qcal1[iHit]),(nlv_tcal0[iHit]+nlv_tcal1[iHit])/2.);
       hVetovsRun->Fill(RunNumber,sqrt(nlv_qcal0[iHit]*nlv_qcal1[iHit]));
       hVetovsEvnt->Fill(jentry,sqrt(nlv_qcal0[iHit]*nlv_qcal1[iHit]));
       hVetoTofvsRun->Fill(RunNumber,(nlv_tcal0[iHit]+nlv_tcal1[iHit])/2.);
       hVetoTofvsEvnt->Fill(jentry,(nlv_tcal0[iHit]+nlv_tcal1[iHit])/2.);
       
       vdoubles[nlv_bar[iHit]+1]=0;

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

   cosm_up=0;
   cosm_down=0;

   veto = 0;
   vetoq= 0;

   for(int j=0;j<nlv_nbar;j++) {
     if (nlv_bar[j]==8) continue;
     
     if((nlv_qcal0[j]>0.1||nlv_qcal0[j]!=nlv_qcal0[j])||
	(nlv_qcal1[j]>0.1||nlv_qcal1[j]!=nlv_qcal1[j])) {
       vetoq=1;
       vqmult++;
     }

     if(nlv_traw0[j]>0||nlv_traw1[j]>0) {
       veto=1;
       vmult++;
     }
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
       hPed1vsEvnt[iBar]->Fill(jentry,e1[iBar]);
       hPed2vsEvnt[iBar]->Fill(jentry,e2[iBar]);

       // sync level

       ee1[iBar]=e1[iBar]*ediff[iBar]*esync[iBar];
       ee2[iBar]=e2[iBar]/ediff[iBar]*esync[iBar];
       
       tt1[iBar]=t1[iBar]-tdiff[iBar]/2.-tsync[iBar]+wlk(e1[iBar]);
       tt2[iBar]=t2[iBar]+tdiff[iBar]/2.-tsync[iBar]+wlk(e2[iBar]);
       tt1nowalk[iBar]=t1[iBar]-tdiff[iBar]/2.-tsync[iBar];
       tt2nowalk[iBar]=t2[iBar]+tdiff[iBar]/2.-tsync[iBar];

       htt1_tcal1->Fill(iBar,tt1[iBar]-tstart-tcal1[iBar]);
       htt2_tcal2->Fill(iBar,tt2[iBar]-tstart-tcal2[iBar]);

       // hit level
       pos[iBar]=(tt1[iBar]-tt2[iBar])*vscint[iBar];
       thit[iBar]=(tt1[iBar]+tt2[iBar])/2.;
 
       // PMT saturation
       ee1[iBar] = ee1[iBar]/(1-0.0096*ee1[iBar]);
       ee2[iBar] = ee2[iBar]/(1-0.0096*ee2[iBar]);
       ehit[iBar]=sqrt(ee1[iBar]*ee2[iBar]);

       if (ehit[iBar]<0) ehit[iBar]=10000.;

       if (fabs(pos[iBar])>500.) continue;

       hEvsBar->Fill(iBar,ehit[iBar]);

       hE1vsPos[iBar]->Fill(pos[iBar],e1[iBar]);
       hE2vsPos[iBar]->Fill(pos[iBar],e2[iBar]);
       hEhitvsPos[iBar]->Fill(pos[iBar],ehit[iBar]);

       if (((iBar-1)/50)%2==1) {
	 xhit[iBar] = (iBar-76-plane*100)*5.+ rnd.Uniform(5.);
	 yhit[iBar] = pos[iBar];
	 zhit[iBar] = 5. + plane*10.+ rnd.Uniform(5.);
       } else {
	 xhit[iBar] = pos[iBar];
	 yhit[iBar] = (iBar-26-plane*100)*5.+ rnd.Uniform(5.);
	 zhit[iBar] = plane*10.+ rnd.Uniform(5.);
       }

       if(yhit[iBar]>50.) cosm_up=1;
       if(yhit[iBar]<-50.) cosm_down=1;

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

       hTofvsX->Fill(xhit[iBar],tof[iBar]);
       hTofvsY->Fill(yhit[iBar],tof[iBar]);
       
       hTofvsEhit[(iBar-1)/50+1]->Fill(ehit[iBar],tof[iBar]);
       hTofvsEhitall->Fill(ehit[iBar],tofc[iBar]);

       hTofvsEvnt->Fill(jentry,tof[iBar]);
       hTofvsRun->Fill(RunNumber,tof[iBar]);
       hEvntvsRun->Fill(RunNumber,jentry);

       hTofvsPath->Fill(s,tof[iBar]);
       
       if(ehit[iBar]>0.){
	 hTof[(iBar-1)/50+1]->Fill(tof[iBar]);
	 hTofc[(iBar-1)/50+1]->Fill(tofc[iBar]);
	 if(ehit[iBar]>10.0) hTofch[(iBar-1)/50+1]->Fill(tofc[iBar]);
	 hEdiffvsBar->Fill(iBar,log(ee1[iBar]/ee2[iBar]));
	 hTofvsZall->Fill(zhit[iBar],tof[iBar]);
	 //hTdiffvsEvnt[iBar]->Fill(jentry,tt1[iBar]-tt2[iBar]);     
       }

       // gamma tof calibration
       if (ehit[iBar]>5.) {
	 hTofbar[iBar]->Fill(tofc[iBar]);
	 hTofplane[(iBar-1)/50+1]->Fill(tofc[iBar]);
       }

       // walk correction calibration
       hTofvsErawall->Fill(e1[iBar],tofnowalk[iBar]);
       hTofvsErawall->Fill(e2[iBar],tofnowalk[iBar]);
       hTofvsErawallcorr->Fill(e1[iBar],tofc[iBar]);
       hTofvsErawallcorr->Fill(e2[iBar],tofc[iBar]);

       //printf("%d %lf\n",iBar,ehit[iBar],tof[iBar]);

       // Christiaan's method
       // for(int j=0;j<nlv_nbar;j++) {
       // 	 if (nlv_bar[j]==8) continue;

       // 	 if (//cl_t-(nlv_tcal0[j]+nlv_tcal1[j])/2.<0.&&
       // 	     fabs(xhit[iBar]*vDistance/Distance-(3.5-nlv_bar[j])*31.5)<18.&&
       // 	     fabs(yhit[iBar]+(nlv_tcal0[j]-nlv_tcal1[j])/0.131)<7.5&&
       // 	     nlv_tcal0[j]>0.&&nlv_tcal1[j]>0.&&
       // 	     ((nlv_qcal0[j]>0.1)||(nlv_qcal0[j]!=nlv_qcal0[j])||
       // 	      (nlv_qcal1[j]>0.1)||(nlv_qcal1[j]!=nlv_qcal1[j]))) {
       // 	   v[iBar]=100.;
       // 	 }
       //}


       if (vetoq>=0) { // ******************* CLUSTERS ********************
	 
	 //if(fabs(yhit[iBar])<95.&&fabs(xhit[iBar])<115.&&v[iBar]<28.) {
	 if(v[iBar]<28.) {
	   
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
	       
	       //if (sqrt((xhit[iBar]-xcl)*(xhit[iBar]-xcl)+
	       //	(yhit[iBar]-ycl)*(yhit[iBar]-ycl)+
	       //	(zhit[iBar]-zcl)*(zhit[iBar]-zcl))<10.) {
	       if (fabs(xhit[iBar]-xcl)<7.5&&
		   fabs(yhit[iBar]-ycl)<7.5&&
		   fabs(zhit[iBar]-zcl)<7.5
		   &&fabs(thit[iBar]-tcl-0.4)<0.6) {

		 //hTdiffClust->Fill(sqrt((xhit[iBar]-xcl)*(xhit[iBar]-xcl)+
		 //		   (yhit[iBar]-ycl)*(yhit[iBar]-ycl)+
		 //			(zhit[iBar]-zcl)*(zhit[iBar]-zcl)));
		 hTdiffClust->Fill(thit[iBar]-tcl);
		 
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
   
   Double_t cl_e, cl_t, cl_x, cl_y, cl_z, cl_ehit, cl_elast, cl_tlast;
   Int_t cl_bar=0;

   hClusterNo->Fill(clno);
   hClmultvsEvnt->Fill(jentry,clno);
   int nclno=0, cclno=0;

   if (clno>0) {
     for (Int_t i=1;i<=clno;i++) {
       hClusterDim->Fill(cl[i][0]);
       cl_e = 0.;
       cl_elast = 0.;
       cl_ehit = 0.;
       cl_t = 10000.;
       cl_tlast = 0.;
       cl_x = -1000.;
       cl_y = -1000.;
       cl_z = -1000.;
       for (Int_t j=1;j<=cl[i][0];j++) {
	 //if (ehit[cl[i][j]]>5.)
	 cl_e = cl_e + ehit[cl[i][j]];
	 if (tof[cl[i][j]]<cl_t) {
	   cl_t = tof[cl[i][j]];
	   cl_x = xhit[cl[i][j]];
	   cl_y = yhit[cl[i][j]];
	   cl_z = zhit[cl[i][j]];
	   cl_ehit = ehit[cl[i][j]];
	   cl_bar = cl[i][j];
	 }
	 if (tof[cl[i][j]]>cl_tlast) {
	   cl_tlast = tof[cl[i][j]];
	   cl_elast = ehit[cl[i][j]];
	 }	    
       }

 
       if(fabs(cl_x)<115.&&fabs(cl_y)<90.) {
	 hTofvsEhitClust->Fill(cl_e,cl_t);
	 hTofvsEhitClustLast->Fill(cl_elast,cl_t);
	 
	 hXYclust->Fill(cl_x,cl_y);
	 hXZclust->Fill(cl_z,cl_x);
	 hYZclust->Fill(cl_z,cl_y);
	 
	 int in_veto = 0;
	 for(int j=0;j<nlv_nbar;j++) {
	   if (nlv_bar[j]==8) continue;
	   
	   if (//cl_t-(nlv_tcal0[j]+nlv_tcal1[j])/2.<0.&&
	       fabs(cl_x*vDistance/Distance-(3.5-nlv_bar[j])*31.5)<18.&&
	       //fabs(cl_y+(nlv_tcal0[j]-nlv_tcal1[j])/0.131)<7.5&&
	       //nlv_tcal0[j]>0.&&nlv_tcal1[j]>0.&&
	       ((nlv_qcal0[j]>0.1)||(nlv_qcal0[j]!=nlv_qcal0[j])||
		(nlv_qcal1[j]>0.1)||(nlv_qcal1[j]!=nlv_qcal1[j]))) {
	     
	     in_veto=1;
	     
	     vdoubles[nlv_bar[j]+1]++;
	     
	     hTofdiffn->Fill(cl_t-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	     
	     hVXvsX->Fill(cl_x,(3.5-nlv_bar[j])*25.);
	     
	     hVBarvsX->Fill(cl_x,nlv_bar[j]+1);
	     hVBarvsXcorr->Fill(cl_x*vDistance/(Distance+cl_z),nlv_bar[j]+1);
	     hVTdiffvsY[nlv_bar[j]+1]->Fill(cl_y,nlv_tcal0[j]-nlv_tcal1[j]);
	     hVTdiffvsYall->Fill(cl_y,nlv_tcal0[j]-nlv_tcal1[j]);
	     
	     hVetoEraw->Fill(nlv_bar[j]+1,sqrt(nlv_qraw0[j]*nlv_qraw1[j]));
	     hVetoEcal->Fill(nlv_bar[j]+1,sqrt(nlv_qcal0[j]*nlv_qcal1[j]));
	     hVetoEcal1->Fill(nlv_bar[j]+1,nlv_qcal0[j]);
	     hVetoEcal2->Fill(nlv_bar[j]+1,nlv_qcal1[j]);
	   }
	 }
	 
	 if (in_veto==0) {
	   nclno++;
	   hEvsBarClustVeto->Fill(cl_bar,cl_e);
	   hTofvsEhitClustv[(cl_bar-1)/50+1]->Fill(ehit[cl_bar],tof[cl_bar]);
	   hTofvsBarClustVeto->Fill(cl_bar,cl_t);
	   
	   for(int j=0;j<nlv_nbar;j++) {
	     if (nlv_bar[j]==8) continue;
	     
	     if ((nlv_qcal0[j]>0.1||nlv_qcal0[j]!=nlv_qcal0[j])||
		 (nlv_qcal1[j]>0.1||nlv_qcal1[j]!=nlv_qcal1[j])) {
	       
	       hVTofvsEhitMatch->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	       
	       hTofdiffvsEcal->Fill(nlv_qcal0[j],cl_t-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	       hTofdiffvsEcal->Fill(nlv_qcal1[j],cl_t-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	       hTofdiffvsEhit->Fill(cl_e,cl_t-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	       
	       hTofdiff->Fill(cl_t-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	       hVTdiffEcal1[nlv_bar[j]+1]->Fill(nlv_tcal0[j]-nlv_tcal1[j],nlv_qcal0[j]);
	       hVTdiffEcal2[nlv_bar[j]+1]->Fill(nlv_tcal0[j]-nlv_tcal1[j],nlv_qcal1[j]);
	       
	       hTofdiffvsEcal->Fill(nlv_qcal0[j],cl_t-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	       hTofdiffvsEcal->Fill(nlv_qcal1[j],cl_t-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	       hTofdiffvsEhit->Fill(cl_e,cl_t-(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	       
	       hDXvsX->Fill(cl_x*vDistance/Distance,cl_x-(3.5-nlv_bar[j])*31.5);
	       hDXvsZ->Fill(cl_z,cl_x-(3.5-nlv_bar[j])*31.5);
	       hDYvsY->Fill(-1.*(nlv_tcal0[j]-nlv_tcal1[j])/0.131,cl_y+(nlv_tcal0[j]-nlv_tcal1[j])/0.131);
	       hDYvsZ->Fill(cl_z,cl_y+(nlv_tcal0[j]-nlv_tcal1[j])/0.131);
	       
	       hTofvsEhitClustVeto->Fill(cl_e,cl_t);
	       hFirstZ->Fill(cl_z,cl_t);
	     }
	   }
	 } else {
	   cclno++;
 	   hEvsBarClustNoVeto->Fill(cl_bar,cl_e);
	   hTofvsEhitClustNoVeto->Fill(cl_e,cl_t);

	   hPlaneacut[(cl_bar-1)/50+1]->Fill(cl_x,cl_y);
	   
	   for(int j=0;j<nlv_nbar;j++) {
	     if (nlv_bar[j]==8) continue;
	     
	     if ((nlv_qcal0[j]>0.1||nlv_qcal0[j]!=nlv_qcal0[j])||
		 (nlv_qcal1[j]>0.1||nlv_qcal1[j]!=nlv_qcal1[j])) {
	       
	       hVTofvsEhitNoMatch->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
	     }
	   }
	 }
       }
     }
   } else {
     // for(int j=0;j<nlv_nbar;j++) {
     //   if (nlv_bar[j]==8) continue;
     //   if ((nlv_qcal0[j]>0.1||nlv_qcal0[j]!=nlv_qcal0[j])||
     // 	   (nlv_qcal1[j]>0.1||nlv_qcal1[j]!=nlv_qcal1[j])) 
     // 	 hVTofvsEhitStop->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
     // }
   }
 
   for(int j=0;j<nlv_nbar;j++) {
     if (nlv_bar[j]==8) continue;
     
     hVdoubles->Fill(nlv_bar[j]+1,vdoubles[nlv_bar[j]+1]);
     
     if (vdoubles[nlv_bar[j]+1]>1) 
       hVTofvsEhitStop->Fill(sqrt(nlv_qcal0[j]*nlv_qcal1[j]),(nlv_tcal0[j]+nlv_tcal1[j])/2.);
   }


   hClusterNoNeutrons->Fill(nclno);
   hClusterNoCharged->Fill(cclno);
   hClNovsClNo->Fill(nclno,cclno);

   hNmultvsEvnt->Fill(jentry,nclno);
   hChmultvsEvnt->Fill(jentry,cclno);

   sx=0.;
   sy=0.;
   sx2=0.;
   sy2=0.;
   sxy=0.;
   n=0;
   sz=0.;
   sz2=0.;
   syz=0.;
   nh=0;

   for (iBar=1;iBar<=num_bars;iBar++) {
     if (fired[iBar]==2) {
       n++;
       sx = sx + xhit[iBar];
       sy = sy + yhit[iBar]; 
       sx2 = sx2 + xhit[iBar]*xhit[iBar];
       sy2 = sy2 + yhit[iBar]*yhit[iBar];
       sxy = sxy + xhit[iBar]*yhit[iBar];
       
       sz = sz + zhit[iBar];
       sz2 = sz2 + zhit[iBar]*zhit[iBar];
       syz = syz + zhit[iBar]*yhit[iBar];
     }
   }

   a = (n*sxy-sx*sy)/(n*sx2-sx*sx);  //  X-Y plane
   b = (sy-a*sx)/n;
   
   ah = (n*syz-sz*sy)/(n*sz2-sz*sz); //  Z-Y plane
   bh = (sy-ah*sz)/n;
   
   for (iBar=1;iBar<=num_bars;iBar++) {
     if (fired[iBar]==2) {

       if(mult>0) { //&&(fabs(a)>1.0||fabs(ah)>1.0)) {  //&&cosm_up+cosm_down==2) {
	 if (ehit[iBar]>0.) {
	   hEcosm[iBar]->Fill(ehit[iBar]);
	   hEvsBarcosm->Fill(iBar,ehit[iBar]);
	   hEvsBarcosmprof->Fill(iBar,ehit[iBar]);
	   hE1[iBar]->Fill(e1[iBar]);
	   hE2[iBar]->Fill(e2[iBar]);
	   hTdiffvsBar->Fill(iBar,tt1[iBar]-tt2[iBar]);
	 }
	 
	 for (jBar=1;jBar<=num_bars;jBar++) {
	   if ((fired[jBar]==2)&&(jBar!=iBar)) {
	     if (ehit[iBar]>0.&&ehit[jBar]>0.) {
	       
	       // if((abs(iBar-jBar)%100)==0&&((iBar-1)%100)<50) {
	       if((iBar-1)/50==(jBar-1)/50) {
		 
		 dist = sqrt((xhit[iBar]-xhit[jBar])*(xhit[iBar]-xhit[jBar])+
			     (yhit[iBar]-yhit[jBar])*(yhit[iBar]-yhit[jBar])+
			     (zhit[iBar]-zhit[jBar])*(zhit[iBar]-zhit[jBar]));
	       } else {
		 //dist = (yhit[iBar]-yhit[jBar])/fabs(yhit[iBar]-yhit[jBar])*
		 dist = -1.*(zhit[iBar]-zhit[jBar])/fabs(zhit[iBar]-zhit[jBar])*
		   sqrt((xhit[iBar]-xhit[jBar])*(xhit[iBar]-xhit[jBar])+
			(yhit[iBar]-yhit[jBar])*(yhit[iBar]-yhit[jBar])+
			(zhit[iBar]-zhit[jBar])*(zhit[iBar]-zhit[jBar]));
	       }
	       //hDT[iBar]->Fill(jBar,thit[jBar]-thit[iBar]-dist/c);
	       hDT[iBar]->Fill(jBar,thit[jBar]-thit[iBar]);
	       hTsync[iBar][jBar]->Fill(thit[jBar]-thit[iBar]-dist/c);
	       //hTsync[iBar][jBar]->Fill(thit[jBar]-thit[iBar]);

	       //if (fabs(a)>1.0||fabs(ah)>1.0) {
		 hCosmall->Fill(dist,thit[jBar]-thit[iBar]);
		 hCosmallprof->Fill(dist,thit[jBar]-thit[iBar]);
		 hDTvsEall->Fill(ehit[jBar],thit[jBar]-thit[iBar]);
		 //} 
	     }
	   }
	 } 
       }

       hPlane[(iBar-1)/50+1]->Fill(xhit[iBar],yhit[iBar]);

       if(vetoq==0) {
	 //hPlanecut[(iBar-1)/50+1]->Fill(xhit[iBar],yhit[iBar]);
	 if(fabs(yhit[iBar])<95.&&fabs(xhit[iBar])<115.) {
	   hTofvsEhitallVetoq->Fill(ehit[iBar],tofc[iBar]);
	   hTofvsEhitvq[(iBar-1)/50+1]->Fill(ehit[iBar],tof[iBar]);
	   hEvsBarVetoq->Fill(iBar,ehit[iBar]);
	   hTofvsBarVetoq->Fill(iBar,tof[iBar]);
	   
	   //if (tof[iBar]>10000./(ehit[iBar]+60.)&&tof[iBar]<100.) {
	   hXY->Fill(xhit[iBar],yhit[iBar]);
	   for(int j=0;j<nlv_nbar;j++) {
	     if (nlv_bar[j]==8) continue;
	     hVEhitvsEhit->Fill(ehit[iBar],sqrt(nlv_qcal0[j]*nlv_qcal1[j]));
	     // hVetoEraw->Fill(nlv_bar[j]+1,sqrt(nlv_qraw0[j]*nlv_qraw1[j]));
	     // hVetoEcal->Fill(nlv_bar[j]+1,sqrt(nlv_qcal0[j]*nlv_qcal1[j]));
	     // hVetoEcal1->Fill(nlv_bar[j]+1,nlv_qcal0[j]);
	     // hVetoEcal2->Fill(nlv_bar[j]+1,nlv_qcal1[j]);
	     //hVTdiffEcal1[nlv_bar[j]+1]->Fill(nlv_tcal0[j]-nlv_tcal1[j],nlv_qcal0[j]);
	     //hVTdiffEcal2[nlv_bar[j]+1]->Fill(nlv_tcal0[j]-nlv_tcal1[j],nlv_qcal1[j]);
	   }
	 }
       } else {
	 hPlaneacut[(iBar-1)/50+1]->Fill(xhit[iBar],yhit[iBar]);
       }
	 
       if(veto==0) {
	 if(fabs(yhit[iBar])<95.&&fabs(xhit[iBar])<115.) {
	   hTofvsEhitallVeto->Fill(ehit[iBar],tof[iBar]);
	   hTofvsEhitv[(iBar-1)/50+1]->Fill(ehit[iBar],tof[iBar]);
	   hVelvsBarVeto->Fill(iBar,v[iBar]);
	   hTofvsBarVeto->Fill(iBar,tof[iBar]);
	 }
       } else {
	 if(fabs(yhit[iBar])<95.&&fabs(xhit[iBar])<115.)
	   hTofvsEhitnv[(iBar-1)/50+1]->Fill(ehit[iBar],tof[iBar]);
       }

       //int in_veto = 0;
       //for(int j=0;j<nlv_nbar;j++) {
       //if (nlv_bar[j]==8) continue;
       //if(fabs(yhit[iBar])>95.||fabs(xhit[iBar])>115.) continue;
       //if (in_veto==1) continue;
       //if (fabs(xhit[iBar]*vDistance/(Distance+zhit[iBar])-(3.5-nlv_bar[j])*25.)<12.5&&
       //  ((nlv_qcal0[j]>0.1||nlv_qcal0[j]!=nlv_qcal0[j])||
       //     (nlv_qcal1[j]>0.1||nlv_qcal1[j]!=nlv_qcal1[j]))) continue;
       //in_veto=1;
       //}
       
     }
   }     

   sx=0.;
   sy=0.;
   sx2=0.;
   sy2=0.;
   sxy=0.;
   n=0;
   sxh=0.;
   syh=0.;
   sx2h=0.;
   sy2h=0.;
   sxyh=0.;
   nh=0;
 
   for (int i=1;i<=50;i++) {
     for (int p=0;p<num_planes;p++) {
       iBar=p*100+i;
       for (int j=51;j<=100;j++) {
	 jBar=p*100+j;
	 if (fired[iBar]==2&&fired[jBar]==2){	
	   //hTsync2D[iBar]->Fill((j-50),thit[jBar]-thit[iBar]); 
	   //hTsync2D[jBar]->Fill(i,thit[jBar]-thit[iBar]); 
	   hTsync2D[iBar]->Fill((j-50),tofc[jBar]-tofc[iBar]); 
	   hTsync2D[jBar]->Fill(i,tofc[jBar]-tofc[iBar]); 
	   if (mult>10&&(fabs(a)>2.0||fabs(ah)>2.0)) { // mult>4 for cosmics
	     //if (plmultH[p]==1&&plmultV[p]==1) {
	     hTdiffvsXBar[iBar]->Fill(j-50,tt1[iBar]-tt2[iBar]);
	     hTdiffvsXBarprof[iBar]->Fill(j-50,tt1[iBar]-tt2[iBar]);
	     hEdiffvsXBar[iBar]->Fill(j-50,log(e1[iBar]/e2[iBar]));
	     hEdiffvsXBarprof[iBar]->Fill(j-50,log(e1[iBar]/e2[iBar]));
	     
	     hTdiffvsXBar[jBar]->Fill(i,tt1[jBar]-tt2[jBar]);
	     hTdiffvsXBarprof[jBar]->Fill(i,tt1[jBar]-tt2[jBar]);
	     hEdiffvsXBar[jBar]->Fill(i,log(e1[jBar]/e2[jBar]));
	     hEdiffvsXBarprof[jBar]->Fill(i,log(e1[jBar]/e2[jBar]));
	     //}
	   }
	 }
       }
       
       if (plmultH[p]>10&&plmultV[p]==0) {
	 if (fired[iBar]==2) {
	   sxh = sxh + xhit[iBar];
	   syh = syh + yhit[iBar]; 
	   sx2h = sx2h + xhit[iBar]*xhit[iBar];
	   sy2h = sy2h + yhit[iBar]*yhit[iBar];
	   sxyh = sxyh + xhit[iBar]*yhit[iBar];
	   nh++;
	   hE1cut[iBar]->Fill(e1[iBar]);
	   hE2cut[iBar]->Fill(e2[iBar]);
	 }
       }

       if (plmultV[p]>10&&plmultH[p]==0) {
	 if (fired[iBar+50]==2) {
	   sx = sx + xhit[iBar+50];
	   sy = sy + yhit[iBar+50]; 
	   sx2 = sx2 + xhit[iBar+50]*xhit[iBar+50];
	   sy2 = sy2 + yhit[iBar+50]*yhit[iBar+50];
	   sxy = sxy + xhit[iBar+50]*yhit[iBar+50];
	   n++;
	   hE1cut[iBar+50]->Fill(e1[iBar+50]);
	   hE2cut[iBar+50]->Fill(e2[iBar+50]);
	 }
       }
     }
   }
   
   a = (n*sxy-sx*sy)/(n*sx2-sx*sx);
   b = (sy-a*sx)/n;
   path = sqrt(1.0+a*a);
   
   ah = (nh*sxyh-sxh*syh)/(nh*sy2h-syh*syh);
   bh = (sxh-ah*syh)/nh;
   pathh = sqrt(1.0+ah*ah);
   
   for (int p=0;p<num_planes;p++) {  
     for (int i=1;i<=50;i++) {
       iBar=p*100+i;

       if (plmultH[p]>10&&plmultV[p]==0) {
	 if (fired[iBar]==2) {
	   // hAvsEdepH->Fill(ehit[iBar],fabs(ah));
	   // hEdepvsPathH->Fill(pathh,ehit[iBar]);
	   // hEdepvsPathHprof->Fill(pathh,ehit[iBar]);
	 }
       }
       
       if (plmultV[p]>10&&plmultH[p]==0) {
	 if (fired[iBar+50]==2) {
	   // hAvsEdepV->Fill(ehit[iBar+50],fabs(a));
	   // hEdepvsPathV->Fill(path,ehit[iBar+50]);
	   // hEdepvsPathVprof->Fill(path,ehit[iBar+50]);
	 }
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

   hNLvsKA->Fill(katmult[evt_cnt],mult);
   hNLvsKY->Fill(kyomult[evt_cnt],mult);
   hKYvsKA->Fill(katmult[evt_cnt],kyomult[evt_cnt]);
   hNLvsKAmix->Fill(katmult[evt_cnt+1],mult);
   hNLvsKYmix->Fill(kyomult[evt_cnt+1],mult);
   hKYvsKAmix->Fill(katmult[evt_cnt],kyomult[evt_cnt+1]);
   hNLvsKAmix1->Fill(katmult[evt_cnt+2],mult);
   hNLvsKYmix1->Fill(kyomult[evt_cnt+2],mult);
   hKYvsKAmix1->Fill(katmult[evt_cnt+2],kyomult[evt_cnt+1]);

   nrun = RunNumber;
   nevt = evt_cnt;
   nevttot = jentry;
   nmult = clno;

   //ntree->Fill();

   //if (vetoq==0)
   continue;

   Double_t ethr=0.; //7.   // tracking / clustering threshold (delta E).
   Double_t velthr=28.; // tracking / clustering threshold (velocity). : gamma rejection

   // cluster finding
   Int_t nClust = 0;
   for(Int_t iClust=0; iClust<nMaxClust; iClust++){
     clust[iClust].SetnHit(0);
   }
   NeulandHit hit;

   for (iBar=1; iBar<=num_bars; iBar++) {

     if(fired[iBar]==2 && ehit[iBar]>ethr && v[iBar]<velthr){  // Cut low-energy component / target gamma ray
       //       printf("%d %lf\n",iBar,tof[iBar]);
       bool first = true;
       Int_t iClustFirst= -1;
       hit.setT(tof[iBar]);
       hit.setX(xhit[iBar]);
       hit.setY(yhit[iBar]);
       hit.setZ(zhit[iBar]);

       for(Int_t iClust=0; iClust<nClust; iClust++){
	 if(clust[iClust].IsInCluster(hit)){
	   //if this Hit belongs to clust[iClust]

	   if(first){ 
	     // if clust[iClust] is the fist cluster
	     // -> Add hit to this cluster
	     first = false;
	     clust[iClust].AddHit(hit);
	     iClustFirst=iClust;
	   } else { 
	     // if clust[iClust] is not the fist cluster
	     // -> Add clust[iClust] into clust[iClustFirst]
	     clust[iClustFirst].AddCluster(clust[iClust]);
	     clust[iClust].SetnHit(0);
	   }
	 }
       }
       if(first == true){
	 // hit does not belong to any existing cluster
	 // -> put in a new cluster
	 clust[nClust].AddHit(hit);
	 nClust++;
       }

     }
   }  // for (iBar=1...  

   // sort hits in each cluster with the hit time
   for(Int_t iClust=0; iClust<nClust; iClust++){
     clust[iClust].TimeSort();
   }
   // sort clusters with the cluster size to reject size 0
   sort(clust, clust + nClust,GreaterNeulandCluster());
   for(Int_t iClust=0; iClust<nClust; iClust++){
     if(clust[iClust].GetnHit()==0) {
       nClust=iClust;
       continue;
     }
   }

   hnCluster->Fill(nClust);
   if(nClust==2) {
     dist = sqrt( pow((clust[0].GetHit(0).getX()-clust[1].GetHit(0).getX()),2.0)
                 +pow((clust[0].GetHit(0).getY()-clust[1].GetHit(0).getY()),2.0)
                 +pow((clust[0].GetHit(0).getZ()-clust[1].GetHit(0).getZ()),2.0));
     dt = abs(clust[0].GetHit(0).getT()-clust[1].GetHit(0).getT());
     hCausality->Fill(dt,dist);
     hVelCluster->Fill(dist/dt);
     hCVelvsDist->Fill(dist,dist/dt);

     if( clust[0].GetnHit()==1 || clust[1].GetnHit()==1 ) {
       hCausalityGamma->Fill(dt,dist);
       hVelClusterGamma->Fill(dist/dt);
       hCVelvsDistGamma->Fill(dist,dist/dt);
     }
   }
   
   // printf("nClust %d\n",nClust);
   // for(Int_t iClust=0; iClust<nClust; iClust++){
   //   printf("  iClust: %d, size: %d\n", iClust, clust[iClust].GetnHit());
   //   for(iHit=0; iHit<clust[iClust].GetnHit(); iHit++){
   //     printf("    iHit: %d %lf %lf %lf %lf\n", 
   // 	      iHit, clust[iClust].GetHit(iHit).getT(), clust[iClust].GetHit(iHit).getX(), clust[iClust].GetHit(iHit).getY(), clust[iClust].GetHit(iHit).getZ());
   //   }
   // }


   //if(nClust!=2) continue;  // only nClust=2
   //if(!causality(clust[0],clust[1])) continue;   // causality cut
   //if(clust[0].GetnHit()<2 || clust[1].GetnHit()<2) continue; // only cluster_size>=2

   n2n++;
   //   if(n2n%10!=1) continue;

   if(! fDisplayTrack) continue;

   // cluster display by ascii output
   char HitXZ[9][51];
   char HitYZ[9][51];
   char char1[1];
   //   if ((float(jentry)/10000.)==int(jentry/10000)) {
   if (1) {
     //     printf("\r event: %ud of %ud\n", jentry,nentries);
     //     cout << "event: " << jentry << " of " << nentries << endl;

     dist = sqrt( pow((clust[0].GetHit(0).getX()-clust[1].GetHit(0).getX()),2.0)
                 +pow((clust[0].GetHit(0).getY()-clust[1].GetHit(0).getY()),2.0)
                 +pow((clust[0].GetHit(0).getZ()-clust[1].GetHit(0).getZ()),2.0));
     dt = abs(clust[0].GetHit(0).getT()-clust[1].GetHit(0).getT());

     printf("nCluster: %d\n", nClust);
     printf(" [FirstHits]\n");
     printf("  Clust0: T = %8.3lf\n", clust[0].GetHit(0).getT());
     printf("  Clust1: T = %8.3lf\n", clust[1].GetHit(0).getT());
     printf("        dT  = %8.3lf\n", dt);
     printf("       dist = %8.3lf\n", dist);
     printf("    dist/dT = %8.3lf\n", dist/dt);

     printf("       ZX        ||         ZY     \n");
     for (int i=1;i<=8;i++) {
       for (int j=1;j<=50;j++) {
   	 HitXZ[i][j]=' ';
   	 HitYZ[i][j]=' ';
       }
     }

     for(Int_t iClust=0; iClust<nClust; iClust++){
       for(iHit=0; iHit<clust[iClust].GetnHit(); iHit++){
   	 int ix = (clust[iClust].GetHit(iHit).getX() + 125)/5 + 1;
   	 int iy = (clust[iClust].GetHit(iHit).getY() + 125)/5 + 1;
   	 int iz = (clust[iClust].GetHit(iHit).getZ()      )/5 + 1;

   	 sprintf(char1,"%1d",iClust);
	 HitXZ[iz][51-ix]=char1[0];
	 HitYZ[iz][51-iy]=char1[0];
       }
     }

     for (int j=1;j<=50;j++) {
       for (int i=1;i<=8;i++) {
   	 printf("%c ",HitXZ[i][j]);
       }
       printf(" %02d  ",51-j);
       for (int i=1;i<=8;i++) {
   	 printf("%c ",HitYZ[i][j]);
       }
       printf("\n");
     }
   }
   printf("\n");

 } 
 printf("Number of 2n candidates: %d\n", n2n);
}

void calibrate_tdiff()
{
   cout<<"performing calibration of bars!"<<endl;

  Int_t bins, Ntotal, l, r;
  Double_t xcenter, ave, left, right, center;
  for(Int_t i=1; i<=num_bars;i++){
 // Finding time btcal. It should be zero for a fully illuminated bar.
    Int_t nbins=hTdiff[i]->GetNbinsX(); 
    Ntotal=0;
    center=0;
    ave=0;
    bins=0;
    l=0;
    r=0;
    
    for(Int_t k=1; k<nbins; k++){ 
      Double_t binContent = hTdiff[i]->GetBinContent(k);  
      Ntotal += binContent;
      xcenter = hTdiff[i]->GetXaxis()->GetBinCenter(k);
      center+=xcenter*binContent;
      if(binContent>0) bins+=1;
      ave+=binContent;
    }
    if(bins>0){
      ave=ave/bins;
    }
    else{
      ave=0;
    }
    cout <<"ave " << ave << endl;
    for(Int_t k=1; k<nbins; k++){ 
      Double_t binContent = hTdiff[i]->GetBinContent(k);  
      if(binContent>ave/2. && l==0.) l=k;
      if(binContent<ave/2. && l!=0. && r==0. && k>l+10) r=k;
    }
    
    left=hTdiff[i]->GetXaxis()->GetBinCenter(l);
    right=hTdiff[i]->GetXaxis()->GetBinCenter(r);
    if(abs(right-left)!=0){
      vscint[i]=270./abs(right-left);
    }
    else{
      vscint[i]=0.;
    }
    
    tdiff[i]=(left+right)/2.;
    
    cout <<i<< "  left, right :"<< left << "  " << right << endl;  
    cout << "tdiff, v: " <<tdiff[i]<<"  "<< vscint[i] << "\n";
  }

  string line;
  ofstream a_file ( "./parameter_test/neuland_sync_1.txt" );
  a_file<<"#bar  tdiff_offset[ns]  sync_offset[ns]  vscint[cm/ns] \n";
    
  for(Int_t i=1;i<=num_bars;i++){
    a_file<< i <<"     " << tdiff[i] <<"     " << tsync[i]<< "     "<< vscint[i] << "\n";
  }   
  a_file.close();
}


void calibrate_tdiff2()
{

  Double_t offset,slope;
  ofstream a_file ( "./parameter_test/neuland_syncXX_2.txt" );
  a_file<<"#bar  tdiff_offset[ns]  sync_offset[ns]  vscint[cm/ns] \n";

  for(int i=1;i<=num_bars;i++) {

    TF1 *f1 = new TF1("f1", "pol1");
    hTdiffvsXBarprof[i]->Fit("f1");
    offset=f1->GetParameter(0);
    slope=f1->GetParameter(1);

    tdiff[i] = offset + slope*25.5;
    vscint[i] = 5./fabs(slope);
    
    a_file<< i <<"     " << tdiff[i] <<"     " << tsync[i]<< "     "<< vscint[i] << "\n";
  }
  a_file.close();
}

void calibrate_tdiff3()
{

  Int_t ntot, peakbin, nbins;
  ofstream a_file ( "./parameter_test/neuland_tdiff3.txt" );
  a_file<<"#bar  ediff_offset[ns]  esync_offset[ns]  att[cm/ns] \n";

  Double_t x,y,a,b,sx, sy, sx2, sy2, sxy;
  int n;
  
  for(int j=1;j<=num_bars;j++) {

    nbins=hTdiffvsXBar[j]->GetNbinsX(); 
 
    sx=0.;
    sy=0.;
    sx2=0.;
    sy2=0.;
    sxy=0.;
    n=0;

    for(int i=1;i<nbins;i++) {
      TH1F *cch;
      cch = (TH1F*)hTdiffvsXBar[j]->ProjectionY("cch",i,i+1);
      ntot = cch->GetEntries();
    
      if (ntot>0) {
	peakbin=cch->GetMaximumBin();
	
	x=hTdiffvsXBar[j]->GetXaxis()->GetBinCenter(i);
	y=cch->GetXaxis()->GetBinCenter(peakbin);
	
	sx = sx + x;
	sy = sy + y; 
	sx2 = sx2 + x*x;
	sy2 = sy2 + y*y;
	sxy = sxy + x*y;
	n++;
      }
    }
    a = (n*sxy-sx*sy)/(n*sx2-sx*sx);
    b = (sy-a*sx)/n;
  
    tdiff[j] = b + a*25.5;
    vscint[j] = 5./fabs(a);

    a_file<< j <<"     " << tdiff[j] <<"     " << tsync[j]<< "     "<< vscint[j] << "\n";
  }
  
  a_file.close();
}

void calibrate_tsync0()
{
  Double_t mean[num_bars+1][num_bars+1];
  Double_t hmean;
  Double_t hmeanH[4], hmeanV[4], hmean175_225;
  Int_t nh, nv, nall;
  
  Int_t peakbin;
  Double_t peakh, left, right;
  
  for (int i=1;i<=num_bars;i++) {
    for (int j=1;j<=num_bars;j++) {

      mean[i][j] = 0.;
      if(i==j) continue;
  
      peakbin=hTsync[i][j]->GetMaximumBin();
      peakh=hTsync[i][j]->GetXaxis()->GetBinCenter(peakbin);
      
      left = peakh - 10.;
      right = peakh + 10.;
      
      //TF1 *f1 = new TF1("f1", "gaus", left, right);
      //hTsync[i][j]->Fit("f1", "R");
      //mean[i][j] = f1->GetParameter(1); //value of 1st parameter	
      
      mean[i][j] = peakh;	
    }
  }
  
  for (int p=0;p<num_planes;p++) {
    for (int i=1;i<=50;i++) {
      for (int j=51;j<=100;j++) {
   	  hTsync_corr[p*100+j]->Fill(mean[p*100+i][p*100+j]-mean[p*100+i][p*100+75]);
   	  hTsync_corr[p*100+i]->Fill(mean[p*100+i][p*100+j]-mean[p*100+25][p*100+j]);
 
  	  hTsync2D_corr->Fill(p*100+j,mean[p*100+i][p*100+j]-mean[p*100+i][p*100+75]);
   	  hTsync2D_corr->Fill(p*100+i,mean[p*100+i][p*100+j]-mean[p*100+25][p*100+j]);
      }
    }
  }

  for (int p=0;p<num_planes;p++) {
    hmeanV[p]=0.;
    hmeanH[p]=0.;
    hmean175_225=0.;
    nh=0;
    nv=0;
    nall=0;
    for (int i=1;i<=50;i++) {
      hmeanH[p] = hmeanH[p] + (mean[p*100+i][p*100+75]-mean[p*100+i][175]);
      nh = nh + 1;
      hTsync_corrH[p+1]->Fill(mean[p*100+i][p*100+75]-mean[p*100+i][175]);
      hmeanV[p] = hmeanV[p] +  (mean[p*100+25][p*100+i+50]-mean[225][p*100+i+50]);
      nv = nv + 1;
      hTsync_corrV[p+1]->Fill(mean[p*100+i][p*100+75]-mean[p*100+i][175]);
      hmean175_225 = hmean175_225 + (mean[p*100+i][175]-mean[p*100+i][225]);
      hmean175_225 = hmean175_225 + (mean[p*100+i+50][175]-mean[p*100+i+50][225]);
      nall = nall + 2;
      hTsync_corr175_225->Fill(mean[p*100+i][175]-mean[p*100+i][225]);
      hTsync_corr175_225->Fill(mean[p*100+i+50][175]-mean[p*100+i+50][225]);
    }
    hmeanH[p]=hmeanH[p]/nh;
    hmeanV[p]=hmeanV[p]/nv;
    hmean175_225=hmean175_225/nall;
  }

  for (int p=0;p<num_planes;p++) {
    for (int i=1;i<=100;i++) {
      hmean = hTsync_corr[p*100+i]->GetMean();
      
      printf("before  %d     %f     %f\n",p*100+i,tsync[p*100+i],hmean);
      if (i<=50) {
	//tsync[p*100+i]=tsync[p*100+i]-hmean + mean[25][p*100+25];//+mean[p*100+25][p*100+75];
	tsync[p*100+i]=tsync[p*100+i] - hmean - hmeanV[p];//+mean[p*100+25][p*100+75];
	//printf("vs 25     %f\n",mean[25][p*100+25]); 	
      } else {
	//tsync[p*100+i]=tsync[p*100+i]+hmean + mean[25][p*100+75];// +mean[p*100+25][p*100+75];
	tsync[p*100+i]=tsync[p*100+i] + hmean + hmeanH[p];// + hmean175_225;// +mean[p*100+25][p*100+75];
	//printf("vs 25     %f\n",mean[25][p*100+75]);
      }
      printf("after           %f\n",tsync[p*100+i]);      
    }
  }


  ofstream a_file ( "./parameter_test/neuland_sync_2x0.txt" );
  a_file<<"#bar  tdiff_offset[ns]  sync_offset[ns]  vscint[cm/ns] \n";
  for(Int_t i=1;i<=num_bars;i++){
    a_file<< i <<"     " << tdiff[i] <<"     " << tsync[i]<< "     "<< vscint[i] << "\n";
  }   
  a_file.close();	 

}


void calibrate_tsync()
{
  Double_t mean[num_bars+1][num_bars+1];
  Double_t hmean;
  Double_t hmeanH[4], hmeanV[4], hmean175_225;
  Int_t nh, nv, nall;
 
  Int_t peakbin;
  Double_t peakh, left, right;
  
  for (int i=1;i<=num_bars;i++) {
    for (int j=1;j<=num_bars;j++) {

      mean[i][j] = -10000.;
  
      if (i!=j&&(hTsync[i][j]->GetEntries())>0.) {
	peakbin=hTsync[i][j]->GetMaximumBin();
	peakh=hTsync[i][j]->GetXaxis()->GetBinCenter(peakbin);
	
	left = peakh - 10.;
	right = peakh + 10.;
	
	//TF1 *f1 = new TF1("f1", "gaus", left, right);
	//hTsync[i][j]->Fit("f1", "R");
	//mean[i][j] = f1->GetParameter(1); //value of 1st parameter	
	
	mean[i][j] = peakh;	
      }
    }
  }

   // for (int p=0;p<num_planes;p++) {
   //   for (int i=1;i<=50;i++) {
   //     for (int j=51;j<=100;j++) {
   // 	if(1) {
   // 	  hTsync_corr[p*100+j]->Fill(mean[p*100+i][p*100+j]-mean[p*100+i][p*100+75]);
   // 	  hTsync_corr[p*100+i]->Fill(mean[p*100+i][p*100+j]-mean[p*100+25][p*100+j]);
   // 	  hTsync2D_corr->Fill(p*100+j,mean[p*100+i][p*100+j]-mean[p*100+i][p*100+75]);
   // 	  hTsync2D_corr->Fill(p*100+i,mean[p*100+i][p*100+j]-mean[p*100+25][p*100+j]);
   // 	}
   //     }
   //   }
   // }
  hmean175_225=0.;
  nall=0;
  for (int i=1;i<=400;i++) {
    if(mean[i][175]!=-10000&&mean[i][225]!=-10000.) {
      hmean175_225 = hmean175_225 + (mean[i][175]-mean[i][225]);
      nall = nall + 1;
      hTsync_corr175_225->Fill(mean[i][175]-mean[i][225]);
    }
    for (int j=1;j<=400;j++) {
    if(i!=j&&mean[i][j]!=-10000.&&mean[i][((j-1)/50)*50+25]!=-10000.) {
    	hTsync_corr[j]->Fill(mean[i][j]-mean[i][((j-1)/50)*50+25]);
    	hTsync2D_corr->Fill(j,mean[i][j]-mean[i][((j-1)/50)*50+25]);
     }
    }
  }
  hmean175_225=hmean175_225/nall;

  for (int p=0;p<num_planes;p++) {
    hmeanV[p]=0.;
    hmeanH[p]=0.;
    nh=0;
    nv=0;
    for (int i=1;i<=400;i++) {
      if(mean[i][p*100+75]!=-10000&&mean[i][175]!=-10000.) {
	hmeanH[p] = hmeanH[p] + (mean[i][p*100+75]-mean[i][175]);
	nh = nh + 1;
	hTsync_corrH[p+1]->Fill(mean[i][p*100+75]-mean[i][175]);
      }
     if(mean[p*100+25][i]!=-10000&&mean[225][i]!=-10000.) {
       hmeanV[p] = hmeanV[p] +  (mean[p*100+25][i]-mean[225][i]);
       nv = nv + 1;
       hTsync_corrV[p+1]->Fill(mean[p*100+25][i]-mean[225][i]);
     }
    }
    hmeanH[p]=hmeanH[p]/nh;
    hmeanV[p]=hmeanV[p]/nv;
  }


  for (int p=0;p<num_planes;p++) {
    for (int i=1;i<=100;i++) {
      //hmean = hTsync_corr[p*100+i]->GetMean();
      
      peakbin=hTsync_corr[p*100+i]->GetMaximumBin();
      peakh=hTsync_corr[p*100+i]->GetXaxis()->GetBinCenter(peakbin);
      
      left = peakh - 0.5;
      right = peakh + 0.5;
      //TF1 *f1 = new TF1("f1", "gaus", left, right);
      //hTsync_corr[p*100+i]->Fit("f1", "R");
      //hmean = f1->GetParameter(1); //value of 1st parameter	

      hmean = peakh;

      printf("before  %d     %f     %f\n",p*100+i,tsync[p*100+i],hmean);
      if (i<=50) {
	//tsync[p*100+i]=tsync[p*100+i]-hmean + mean[25][p*100+25];//+mean[p*100+25][p*100+75];
	tsync[p*100+i]=tsync[p*100+i]+hmean - hmeanV[p];//+mean[p*100+25][p*100+75];
	//printf("vs 25     %f\n",mean[25][p*100+25]); 	
      } else {
	//tsync[p*100+i]=tsync[p*100+i]+hmean + mean[25][p*100+75];// +mean[p*100+25][p*100+75];
	tsync[p*100+i]=tsync[p*100+i]+hmean + hmeanH[p] + hmean175_225;// +mean[p*100+25][p*100+75];
	//printf("vs 25     %f\n",mean[25][p*100+75]);
      }
      printf("after           %f\n",tsync[p*100+i]);      
    }
  }


  ofstream a_file ( "./parameter_test/neuland_sync_2x.txt" );
  a_file<<"#bar  tdiff_offset[ns]  sync_offset[ns]  vscint[cm/ns] \n";
  for(Int_t i=1;i<=num_bars;i++){
    a_file<< i <<"     " << tdiff[i] <<"     " << tsync[i]<< "     "<< vscint[i] << "\n";
  }   
  a_file.close();	 

}


void calibrate_tsync_gamma()
{
  Int_t peakbin;
  Double_t peakh, left, right;

  Double_t offset;
  ofstream a_file ( "./parameter_test/neuland_sync_gamma.txt" );
  a_file<<"#bar  tdiff_offset[ns]  sync_offset[ns]  vscint[cm/ns] \n";

  for(int i=1;i<=num_bars;i++) {

    offset=0.;

    if((hTofbar[i]->Integral())>100) {
      
      peakbin=hTofbar[i]->GetMaximumBin();
      peakh=hTofbar[i]->GetXaxis()->GetBinCenter(peakbin);
      
      left = peakh - 1.5;
      right = peakh + 1.5;
      //TF1 *f1 = new TF1("f1", "gaus", left, right);
      TF1 *f1 = new TF1("f1", "gaus", 26., 32.);
      
      hTofbar[i]->Fit("f1","R","",left, right);
      offset=f1->GetParameter(1);
 
      left = offset - 1.;
      right = offset + 1.;
      hTofbar[i]->Fit("f1","R","",left, right);
      offset=f1->GetParameter(1);
      
      offset = offset-28.5;

    }
    
    if (fabs(offset)<5.) {
      tsync[i] = tsync[i]+offset;
    } else {
      tsync[i] = tsync[i];
    }
    a_file<< i <<"     " << tdiff[i] <<"     " << tsync[i]<< "     "<< vscint[i] << "\n";   
  }
  a_file.close();
}

void calibrate_tsync_gamma2()
{

  Int_t peakbin, bar;
  Double_t peakh, left, right;

  Double_t offset;
  ofstream a_file ( "./parameter_test/neuland_sync_gamma2.txt" );
  a_file<<"#bar  tdiff_offset[ns]  sync_offset[ns]  vscint[cm/ns] \n";

  for(int i=1;i<=8;i++) {

    offset=0.;

    if((hTofplane[i]->Integral())>100) {
      
      peakbin=hTofplane[i]->GetMaximumBin();
      peakh=hTofplane[i]->GetXaxis()->GetBinCenter(peakbin);
      
      left = peakh - 0.8;
      right = peakh + 0.5;

      TF1 *f1 = new TF1("f1", "gaus", left, right);
      
      hTofplane[i]->Fit("f1","R");
      offset=f1->GetParameter(1);
      //offset = offset-37.666;
      
      offset = peakh - 61.5;  // maximum

    }
    
    for(int j=1;j<=50;j++) {
      bar=(i-1)*50+j;
      if (fabs(offset)<10.) {
	tsync[bar] = tsync[bar]+offset;
      } else {
	tsync[bar] = tsync[bar];
      }
      a_file<< bar <<"     " << tdiff[bar] <<"     " << tsync[bar]<< "     "<< vscint[bar] << "\n";   
    }
  }
  a_file.close();
}

void calibrate_tsync_fine()
{
  Int_t peakbin;
  Double_t peakh, offset, dtmust;

  Int_t plane;

  ofstream a_file ( "./parameter_test/neuland_sync_2xx.txt" );
  a_file<<"#bar  tdiff_offset[ns]  sync_offset[ns]  vscint[cm/ns] \n";

  for(int i=1;i<=num_bars;i++) {

    offset = 0.;  

    plane = (i-1)/50 + 1;

    if(((i-1)%50+1<11||(i-1)%50+1>40)&&(plane%2)==1) {

      peakbin=hTsync[(plane-1)*50+25][i]->GetMaximumBin();
      peakh=hTsync[(plane-1)*50+25][i]->GetXaxis()->GetBinCenter(peakbin);

      dtmust = (((plane-1)*50+25)-i)*5./c;

      //offset = peakh - dtmust;
      offset = peakh; // already -dist/c in hTsync[][]

      printf("%d      %f      %f      %f\n",i,peakh,dtmust,offset);

    }
    
    tsync[i] = tsync[i] + offset;

    a_file<< i <<"     " << tdiff[i] <<"     " << tsync[i]<< "     "<< vscint[i] << "\n";
  }   
  a_file.close();	 
}

void calibrate_ediff()
{
   cout<<"performing calibration of bars!"<<endl;

  Int_t bins, Ntotal, l, r;
  Double_t xcenter, ave, left, right, center;
  for(Int_t i=1; i<=num_bars;i++){
    // Finding time btcal. It should be zero for a fully illuminated bar.
    Int_t nbins=hEdiff[i]->GetNbinsX(); 
    Ntotal=0;
    center=0;
    ave=0;
    bins=0;
    l=0;
    r=0;
    
    for(Int_t k=1; k<nbins; k++){ 
      Double_t binContent = hEdiff[i]->GetBinContent(k);  
      Ntotal += binContent;
      xcenter = hEdiff[i]->GetXaxis()->GetBinCenter(k);
      center+=xcenter*binContent;
      if(binContent>0) bins+=1;
      ave+=binContent;
    }
    if(bins>0){
      ave=ave/bins;
    }
    else{
      ave=0;
    }
    cout <<"ave " << ave << endl;
    for(Int_t k=1; k<nbins; k++){ 
      Double_t binContent = hEdiff[i]->GetBinContent(k);  
      if(binContent>ave/2. && l==0.) l=k;
      if(binContent<ave/2. && l!=0. && r==0. && k>l+10) r=k;
    }
    
    left=hEdiff[i]->GetXaxis()->GetBinCenter(l);
    right=hEdiff[i]->GetXaxis()->GetBinCenter(r);
    if(abs(right-left)!=0){
      att[i]=2*270./abs(right-left);
    }
    else{
      att[i]=0.;
    }
    
    ediff[i]=exp(-1*(left+right)/4.);
    
    cout << "left, right :"<< left << "  " << right << endl;  
    cout << "ediff, att: " <<ediff[i]<<"  "<< att[i] << "\n";
  }

  ofstream a_file ( "./parameter_test/neuland_esync.txt" );
  a_file<<"#bar  ediff_offset[ns]  esync_offset[ns]  att[cm/ns] \n";
    
  for(Int_t i=1;i<=num_bars;i++){
    a_file<< i <<"     " << ediff[i] <<"     " << esync[i]<< "     "<< att[i] << "\n";
  }
  a_file.close();	 
}

void calibrate_ediff2()
{

  Double_t offset,slope;
  ofstream a_file ( "./parameter_test/neuland_esyncXX.txt" );
  a_file<<"#bar  ediff_offset[ns]  esync_offset[ns]  att[cm/ns] \n";

  for(int i=1;i<=num_bars;i++) {

    TF1 *f1 = new TF1("f1", "pol1");
    hEdiffvsXBarprof[i]->Fit("f1");
    offset=f1->GetParameter(0);
    slope=f1->GetParameter(1);

    ediff[i] = exp(-(offset + slope*25.5)/2.);
    att[i] = 2.*5./fabs(slope);
    
    a_file<< i <<"     " << ediff[i] <<"     " << esync[i]<< "     "<< att[i] << "\n";
  }
  a_file.close();
}

void calibrate_ediff3()
{

  Int_t ntot, peakbin, nbins;
  ofstream a_file ( "./parameter_test/neuland_ediff3.txt" );
  a_file<<"#bar  ediff_offset[ns]  esync_offset[ns]  att[cm/ns] \n";

  Double_t x,y,a,b,sx, sy, sx2, sy2, sxy;
  int n;
  
  for(int j=1;j<=num_bars;j++) {

    nbins=hEdiffvsXBar[j]->GetNbinsX(); 
 
    sx=0.;
    sy=0.;
    sx2=0.;
    sy2=0.;
    sxy=0.;
    n=0;

    for(int i=1;i<nbins;i++) {
      TH1F *cch;
      cch = (TH1F*)hEdiffvsXBar[j]->ProjectionY("cch",i,i+1);
      ntot = cch->GetEntries();
    
      if (ntot>0) {
	peakbin=cch->GetMaximumBin();
	
	x=hEdiffvsXBar[j]->GetXaxis()->GetBinCenter(i);
	y=cch->GetXaxis()->GetBinCenter(peakbin);
	
	sx = sx + x;
	sy = sy + y; 
	sx2 = sx2 + x*x;
	sy2 = sy2 + y*y;
	sxy = sxy + x*y;
	n++;
      }
    }
    a = (n*sxy-sx*sy)/(n*sx2-sx*sx);
    b = (sy-a*sx)/n;
  
    ediff[j] = exp(-(b + b*25.5)/2.);
    att[j] = 2.*5./fabs(a);
    
    a_file<< j <<"     " << ediff[j] <<"     " << esync[j]<< "     "<< att[j] << "\n";
  }
  
  a_file.close();
}
void calibrate_walk()
{

  Int_t ntot, peakbin, nbins;
  ofstream a_file ( "./parameter_test/neuland_walk.txt" );

  Double_t x,y,a,b,sx, sy, sx2, sy2, sxy;
  int n, nseg;
  
  nbins=hTofvsErawall->GetNbinsX(); 

  nseg=0;
  x=0.;
  sx=0.;
  sy=0.;
  sx2=0.;
  sy2=0.;
  sxy=0.;
  n=0;

  for(int i=1;i<nbins;i++) {
    TH1F *cch;
    cch = (TH1F*)hTofvsErawall->ProjectionY("cch",i,i+1);
    ntot = cch->GetEntries();
    
    if (ntot>100) {
      
      peakbin=cch->GetMaximumBin();
      
      x=hTofvsErawall->GetXaxis()->GetBinCenter(i);
      y=cch->GetXaxis()->GetBinCenter(peakbin);
      
      sx = sx + x;
      sy = sy + y; 
      sx2 = sx2 + x*x;
      sy2 = sy2 + y*y;
      sxy = sxy + x*y;
      n++;
    }
    
    if (n==10) {
      
      a = (n*sxy-sx*sy)/(n*sx2-sx*sx);
      b = (sy-a*sx)/n;
      
      nseg = nseg + 1;
      
      a_file<< nseg << "     " << x <<"     " << a <<"     " << b << "\n";
      
      sx=0.;
      sy=0.;
      sx2=0.;
      sy2=0.;
      sxy=0.;
      n=0;
    }
  }
  
  a_file.close();
}

void calibrate_esync()
{

  Double_t peak;
  ofstream a_file ( "./parameter_test/neuland_esyncXY.txt" );
  a_file<<"#bar  ediff_offset[ns]  esync_offset[ns]  att[cm/ns] \n";

  for(int i=1;i<=num_bars;i++) {
    peak = hEvsBarcosmprof->GetBinContent(i);
    
    esync[i]=11.5/peak*esync[i];
    
    a_file<< i <<"     " << ediff[i] <<"     " << esync[i]<< "     "<< att[i] << "\n";
  }
  a_file.close();
}

void calibrate_esync2()
{

  Double_t peak, peakh=0;
  Int_t ntot, ntoth, peakbin, nbins, nbinsh;
  bool firstbin=true;
  ofstream a_file ( "./parameter_test/neuland_esyncXXX.txt" );
  a_file<<"#bar  ediff_offset[ns]  esync_offset[ns]  att[cm/ns] \n";

  Double_t x,y,a,b,sx, sy, sx2, sy2, sxy;
  int n;
  
  sx=0.;
  sy=0.;
  sx2=0.;
  sy2=0.;
  sxy=0.;
  n=0;

  nbins=hEdepvsPathV->GetNbinsX(); 
  nbinsh=hEdepvsPathH->GetNbinsX(); 


  for(int i=1;i<nbinsh;i++) {
    TH1F *cch;
    cch = (TH1F*)hEdepvsPathH->ProjectionY("cch",i,i+1);
    ntoth = cch->GetEntries();
    
    if (ntoth>0&&firstbin){
      peakbin=cch->GetMaximumBin();
      peakh=cch->GetXaxis()->GetBinCenter(peakbin);
      firstbin =false;
    }
  }
  printf("horizontal: %f\n",peakh);

  for(int i=1;i<nbins;i++) {
    TH1F *ccc;
    
    ccc = (TH1F*)hEdepvsPathV->ProjectionY("ccc",i,i+1);
    ntot = ccc->GetEntries();

    if (ntot>0&&n<10) {
      peakbin=ccc->GetMaximumBin();
      
      x=hEdepvsPathV->GetXaxis()->GetBinCenter(i);
      y=ccc->GetXaxis()->GetBinCenter(peakbin);
      
      sx = sx + x;
      sy = sy + y; 
      sx2 = sx2 + x*x;
      sy2 = sy2 + y*y;
      sxy = sxy + x*y;
      n++;
    }
  }

  a = (n*sxy-sx*sy)/(n*sx2-sx*sx);
  b = (sy-a*sx)/n;

  peak = a+b;
  printf("verticals: %f\n",peak);
 
  for(int p=0;p<num_planes;p++) {
    for(int i=1;i<=100;i++) {
      if (i>50) {
	esync[p*100+i]=11.5/peak*esync[p*100+i];
      } else {
	esync[p*100+i]=11.5/peakh*esync[p*100+i];
      }
      cout<< p*100+i <<"     " << ediff[p*100+i] <<"     " << esync[p*100+i]<< "     "<< att[p*100+i] << endl;
    a_file<< p*100+i <<"     " << ediff[p*100+i] <<"     " << esync[p*100+i]<< "     "<< att[p*100+i] << "\n";
    }
  }
  a_file.close();
}

void find_thr()
{

  char name[60]; 
  Int_t nbin, norm;
  Int_t binContent;
  Bool_t done;
  Float_t thr, xmin, xmax;

  sprintf(name,"thr/neuland_thr_curr.txt");

  FILE *thrfile;
  thrfile=fopen(name,"w");

  for(Int_t i=1;i<=num_bars;i++){

    printf("bar %i\n",i);

    xmin=hE1[i]->GetXaxis()->GetXmin();
    xmax=hE1[i]->GetXaxis()->GetXmax();
    nbin=hE1[i]->GetNbinsX(); 
    norm=hE1[i]->Integral(150,nbin);
    hE1[i]->Scale(100000./norm);
    done = false;
    for(Int_t k=1; k<=nbin; k++) {
      if (done) continue;
      binContent = hE1[i]->GetBinContent(k);
      if (binContent>1000) {  // was 200
	thr = k*(xmax-xmin)/nbin-(xmax-xmin)/nbin/2.;
	done = true;
	fprintf(thrfile,"%d     1     %f\n",i,thr);
	printf("1 %i    %i   %f\n",i,k,thr);
      }
      if (k==nbin) {
	printf("reci nekaj\n");
	fprintf(thrfile,"%d     1     %f\n",i,0.);
      }
    }
    
    xmin=hE2[i]->GetXaxis()->GetXmin();
    xmax=hE2[i]->GetXaxis()->GetXmax();
    nbin=hE2[i]->GetNbinsX(); 
    norm=hE2[i]->Integral(150,nbin);
    hE2[i]->Scale(100000./norm);
    done = false;
    for(Int_t k=1; k<=nbin; k++) {
      if (done) continue;
      binContent = hE2[i]->GetBinContent(k);
      if (binContent>1000) { // was 200
	thr = k*(xmax-xmin)/nbin-(xmax-xmin)/nbin/2.;
	done = true;
	fprintf(thrfile,"%d     2     %f\n",i,thr);
	printf("2 %i    %i   %f\n",i,k,thr);
      }
      if (k==nbin) {
	printf("reci nekaj\n");
	fprintf(thrfile,"%d     2     %f\n",i,0.);
      }
    }
  }
  fclose(thrfile);
}

void get_hv()
{

  int pl,bar,pm;
  TF1 *gfit;
  ofstream curr_file ( "hv/neuland_curr.txt" );

  FILE *readout;
  char var[80];
  char var2[80];
  float hv;
  
  Double_t mean, sigma, xmin, xmax;
  int nbin;
  
  for(int i=1;i<=400;i++) {
    
    pl = (i-1)/50+1;
    bar = (i-1)%50+1;
    pm = 1;
    sprintf(var,"caget NEULAND:HV:NNP%02d_%02d_%01d:vmon",pl,bar,pm);
    readout = popen(var,"r");
    fscanf(readout,"%s  %f",var2,&hv);
    printf("NNP%02d_%02d_%d HV = %f \n",pl,bar,pm,hv);
    pclose(readout);

    xmin=hE1[i]->GetXaxis()->GetXmin();
    xmax=hE1[i]->GetXaxis()->GetXmax();
    nbin=hE1[i]->GetNbinsX();

    mean = hE1[i]->GetMaximumBin();
    mean = mean*(xmax-xmin)/nbin-(xmax-xmin)/nbin/2.;
    sigma = mean/2.;

    hE1[i]->SetAxisRange(mean-sigma,mean+sigma);
    gfit = new TF1("Gaussian","gaus");
    hE1[i]->Fit("Gaussian");
    mean = gfit->GetParameter(1);
    sigma = gfit->GetParameter(2);

    curr_file<< i <<"     " << pm <<"     " << mean << "     "<< hv << "\n";

    pm = 2;
    sprintf(var,"caget NEULAND:HV:NNP%02d_%02d_%01d:vmon",pl,bar,pm);
    readout = popen(var,"r");
    fscanf(readout,"%s  %f",var2,&hv);
    printf("NNP%02d_%02d_%d HV = %f \n",pl,bar,pm,hv);
    pclose(readout);

    xmin=hE2[i]->GetXaxis()->GetXmin();
    xmax=hE2[i]->GetXaxis()->GetXmax();
    nbin=hE2[i]->GetNbinsX();

    mean = hE2[i]->GetMaximumBin();
    mean = mean*(xmax-xmin)/nbin-(xmax-xmin)/nbin/2.;
    sigma = mean/2.;
    hE2[i]->SetAxisRange(mean-sigma,mean+sigma);
    gfit = new TF1("Gaussian","gaus",mean-sigma/2.,mean+sigma/2.);
    hE2[i]->Fit("Gaussian");
    mean = gfit->GetParameter(1);
    sigma = gfit->GetParameter(2);
 
    curr_file<< i <<"     " << pm <<"     " << mean << "     "<< hv << "\n";
        
  }
  curr_file.close();
}

void fine_tune_hv()
{
  ifstream fl;
  fl.open("hv/neuland_hv_233.txt");
  ifstream fh;
  fh.open("hv/neuland_hv_234.txt");
  ofstream fnew("hv/hv_new.txt");

  int barl, barh;
  int pmtl, pmth;
  float peakl[800];
  float vl[800];
  float peakh[800];
  float vh[800];
  float hv_new[800];

  float a[800];
  float b[800];

  for(int i=0;i<800;i++) {
    fl>>barl>>pmtl>>peakl[i]>>vl[i];
    fh>>barh>>pmth>>peakh[i]>>vh[i];
    
    if ((barl!=barh)||(pmtl!=pmth)){
      printf("something wrong in gainmatching hv\n");
      printf("%d    %d   %d   %d   %d\n",i,barl,pmtl, barh,pmth);
    }
    
    if (fabs(peakh[i]-peakl[i])>5.) {
      a[i]=(vh[i]-vl[i])/(peakh[i]-peakl[i]);
    } else {
      a[i]=0.;
    }
    
    b[i]=vl[i]-a[i]*peakl[i];

     if ((barl-1)%100+1<=50) { // horizontal
      hv_new[i] = a[i]*200.+b[i];   // 200 is reference for horizontal

      if (hv_new[i]>1400.) {
	printf("HV too high %d  %d  %f\n", barl, pmtl, hv_new[i]);
	hv_new[i]=1400.;
      }
      fnew<<barl<<"    "<<pmtl<<"    "<<hv_new[i]<<endl;
    } else {                  // vertical
      hv_new[i] = a[i]*232.+b[i];   // 232 is reference for vertical

      if (hv_new[i]>1400.) {
	printf("HV too high %d  %d  %f\n", barl, pmtl, hv_new[i]);
	hv_new[i]=1400.;
      }
      fnew<<barl<<"    "<<pmtl<<"    "<<hv_new[i]<<endl;
    }  
    
  }
  fnew.close();
  fl.close();
  fh.close();
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

  Int_t nseg;

  y=par1*pow(x,par2)-(par1*pow(400.,par2)); // Michael's
  //y=2.29083*log(x)-0.0870157*log(x)*log(x)-4.57824;  // mine

  // xwalk[0]=0.;
  // for(nseg=1;nseg<=6;nseg++) {
 
  //   if (x>xwalk[nseg-1]&&x<=xwalk[nseg]) {
  //     y = -1*(awalk[nseg]*x+bwalk[nseg]-35.7); //35.7  // needs correction 0.07 -> 35.77 to match powlaw
  //   }
  // }

  return y;
}


void plot_histograms()
{
  
  char name[60]; 

  TCanvas *C_histograms1;
  C_histograms1 = new TCanvas("NeuLAND histograms 1","NeuLAND histograms 1",0,0,1280,1024);
  C_histograms1->Clear();
  C_histograms1->Divide(2,2);
  C_histograms1->cd(1);
  hPlane[1]->Draw("colz");
  C_histograms1->cd(2);
  hPlane[2]->Draw("colz");
  C_histograms1->cd(3);
  hPlane[3]->Draw("colz");
  C_histograms1->cd(4);
  hPlane[4]->Draw("colz");
  sprintf(name,"fig/NeuLANDhistograms1_run%03d.png",RunNumber);
  C_histograms1->Print(Form(name), "png");
  
  TCanvas *C_histograms2;
  C_histograms2 = new TCanvas("NeuLAND histograms 2","NeuLAND histograms 2",0,0,1280,1024);
  C_histograms2->Clear();
  C_histograms2->Divide(2,2);
  C_histograms2->cd(1);
  hPlane[5]->Draw("colz");
  C_histograms2->cd(2);
  hPlane[6]->Draw("colz");
  C_histograms2->cd(3);
  hPlane[7]->Draw("colz");
  C_histograms2->cd(4);
  hPlane[8]->Draw("colz");
  sprintf(name,"fig/NeuLANDhistograms2_run%03d.png",RunNumber);
  C_histograms2->Print(Form(name), "png");
  
  TCanvas *C_histograms3;
  C_histograms3 = new TCanvas("NeuLAND histograms 3","NeuLAND histograms 3",0,0,1200,1000);
  C_histograms3->Clear();
  gStyle->SetOptLogz();
  hTdiffvsBar->Draw("colz");
  sprintf(name,"fig/NeuLANDhistograms3_run%03d.png",RunNumber);
  C_histograms3->Print(Form(name), "png");

  TCanvas *C_histograms4;
  C_histograms4 = new TCanvas("NeuLAND histograms 4","NeuLAND histograms 4",0,0,1200,1000);
  C_histograms4->Clear();
  gStyle->SetOptLogz();
  C_histograms4->SetLogz();
  C_histograms4->Divide(1,2);

  C_histograms4->cd(1);
  hE1vsBar->GetYaxis()->SetRangeUser(0,1000);
  hE1vsBar->Draw("colz");
  C_histograms4->cd(2);
  hE2vsBar->GetYaxis()->SetRangeUser(0,1000);
  hE2vsBar->Draw("colz");
  sprintf(name,"fig/NeuLANDhistograms4_run%03d.png",RunNumber);
  C_histograms4->Print(Form(name), "png");
  
  TCanvas *C_histograms5;
  C_histograms5 = new TCanvas("NeuLAND histograms 5","NeuLAND histograms 5",0,0,1200,1000);
  C_histograms5->Clear();
  gStyle->SetOptLogz();
  C_histograms5->SetLogz();
  C_histograms5->Divide(2,2);

  C_histograms5->cd(1);
  hPlane[2]->Draw("colz");
  C_histograms5->cd(2);
  C_histograms5->cd(3);
  C_histograms5->cd(4);
  hTofvsEhit[2]->Draw("colz");

  sprintf(name,"fig/NeuLANDhistograms5_run%03d.png",RunNumber);
  C_histograms5->Print(Form(name), "png");
  

  TCanvas *C_histograms6;
  C_histograms6 = new TCanvas("NeuLAND histograms 6","NeuLAND histograms 6",0,0,1200,1000);
  C_histograms6->Clear();
  gStyle->SetOptLogz();
  C_histograms6->SetLogz();
  C_histograms6->Divide(2,2);

  C_histograms6->cd(1);
  hTofvsEhit[1]->Draw("colz");
  C_histograms6->cd(2);
  hTofvsEhit[2]->Draw("colz");
  C_histograms6->cd(3);
  hTofvsEhit[3]->Draw("colz");
  C_histograms6->cd(4);
  hTofvsEhit[4]->Draw("colz");

  sprintf(name,"fig/NeuLANDhistograms6_run%03d.png",RunNumber);
  C_histograms6->Print(Form(name), "png");


  TCanvas *C_histograms7;
  C_histograms7 = new TCanvas("NeuLAND histograms 7","NeuLAND histograms 7",0,0,1200,1000);
  C_histograms7->Clear();
  gStyle->SetOptLogz();
  C_histograms7->SetLogz();
  C_histograms7->Divide(2,2);

  C_histograms7->cd(1);
  hTofvsEhit[5]->Draw("colz");
  C_histograms7->cd(2);
  hTofvsEhit[6]->Draw("colz");
  C_histograms7->cd(3);
  hTofvsEhit[7]->Draw("colz");
  C_histograms7->cd(4);
  hTofvsEhit[8]->Draw("colz");

  sprintf(name,"fig/NeuLANDhistograms7_run%03d.png",RunNumber);
  C_histograms7->Print(Form(name), "png");

  return;
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
      cerr << line_no << ": Missing "str", got \"" << token << "\".\n" << endl;\
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

  ifstream b_file ( "parameter_test/neuland_sync_gamma.txt");  // was _gamma
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

void out_hist()
{
  char name[60];
  sprintf(name,"hst/hist_run%03d.root",RunNumber);
  printf("Saving histogram: %s\n",name);

  TFile* fhist = new TFile(name,"recreate");
  hMult->Write();
  hMultPmt->Write();
  hMultPart->Write();
  hMultV->Write();
  hMultH->Write();
  hEvsBar->Write();
  hEvsBarcosm->Write();
  hEvsBarcosmprof->Write();
  hBars->Write();
  hT1vsBar->Write();
  hT2vsBar->Write();
  hE1vsBar->Write();
  hE2vsBar->Write();
  hT3vsBar->Write();
  hT4vsBar->Write();
  hT5vsBar->Write();
  hT6vsBar->Write();
  hTdiffvsBar->Write();
  hEdiffvsBar->Write();
  hTsync2D_corr->Write();
  hDTvsEall->Write();

  for(int i=0;i<=num_planes;i++){
    for(int j=0;j<=num_planes;j++){
      //hCosmEdepX[i][j]->Write();
    }
  }
  for(int i=0;i<=num_planes;i++){
    for(int j=0;j<=num_planes;j++){
      //hCosmEdepY[i][j]->Write();
    }
  }
  for(int i=0;i<=2*num_planes;i++){
      hPlane[i]->Write();
  }
  for(int i=0;i<=num_bars;i++){
    hEcosm[i]->Write();
  }
  for(int i=0;i<=num_bars;i++){
    hE1[i]->Write();
    hE1cut[i]->Write();
  }
  for(int i=0;i<=num_bars;i++){
    hE2[i]->Write();
    hE2cut[i]->Write();
  }

  for(int p=0;p<num_planes;p++) {
    for(int i=1;i<=50;i++){
      for(int j=51;j<=100;j++){
	hTsync[p*100+i][p*100+j]->Write();
      }
    }
  }
  for(int i=0;i<=num_bars;i++){
    hTsync_corr[i]->Write();
  }
  for(int i=0;i<=num_bars;i++){
    hTdiffvsXBar[i]->Write();
    hTdiffvsXBarprof[i]->Write();
    hEdiffvsXBar[i]->Write();
    hEdiffvsXBarprof[i]->Write();
  }
  for(int i=0;i<=num_bars;i++){
    hTsync2D[i]->Write();
  }
  for(int i=0;i<=num_bars;i++){
    hDTvsE[i]->Write();
  }
  for(int i=0;i<=num_bars;i++){
    hEhitvsPos[i]->Write();
  }
  for(int i=1;i<=num_planes*2;i++){
    hTof[i]->Write();
    hTofc[i]->Write();
    hTofch[i]->Write();
    hTofvsEhit[i]->Write();
    hTofvsEhitv[i]->Write();
  }
  hTofvsBar->Write();
  hTofcvsBar->Write();
  
  fhist->Close();
}

void fit_gamma(TH1F* h1, Double_t sigma, Double_t center, Double_t lower, Double_t upper)
{
  Double_t height = h1->GetBinContent(h1->GetXaxis()->FindBin(center));
  Double_t ylower = h1->GetBinContent(h1->GetXaxis()->FindBin(lower));
  Double_t yupper = h1->GetBinContent(h1->GetXaxis()->FindBin(upper));

  TF1 *f1=new TF1("f1","gaus(0)+pol1(3)");
  f1->SetParameter(0, height);
  f1->SetParameter(1, center);
  f1->SetParameter(2, sigma);
  f1->SetParameter(3, ylower-(yupper-ylower)/(upper-lower)*lower);
  f1->SetParameter(4, (yupper-ylower)/(upper-lower));

  h1->Fit(f1, "", "", lower, upper);
}

void subhisto()
{
  int nx,ny;
  double c1,c2;

  nx = hNLvsKY->GetXaxis()->GetNbins(); 
  ny = hNLvsKYmix->GetYaxis()->GetNbins(); 
  for (int i=1;i<=nx;i++) { 
    for (int j=1;j<=ny;j++) { 
      c1 =hNLvsKY->GetBinContent(i,j); 
      if (c1 > 0) { 
	c2 =hNLvsKYmix->GetBinContent(i,j); 
	//if (c1-c2 > 0) 
	  hNLvsKYdiff->SetBinContent(i,j,c1-1.0*c2); 
      }
    }
  }
  nx = hNLvsKA->GetXaxis()->GetNbins(); 
  ny = hNLvsKAmix->GetYaxis()->GetNbins(); 
  for (int i=1;i<=nx;i++) { 
    for (int j=1;j<=ny;j++) { 
      c1 =hNLvsKA->GetBinContent(i,j); 
      if (c1 > 0) { 
	c2 =hNLvsKAmix->GetBinContent(i,j); 
	//if (c1-c2 > 0) 
	  hNLvsKAdiff->SetBinContent(i,j,c1-1.0*c2); 
      }
    }
  }
  nx = hKYvsKA->GetXaxis()->GetNbins(); 
  ny = hKYvsKAmix->GetYaxis()->GetNbins(); 
  for (int i=1;i<=nx;i++) { 
    for (int j=1;j<=ny;j++) { 
      c1 =hKYvsKA->GetBinContent(i,j); 
      if (c1 > 0) { 
	c2 =hKYvsKAmix->GetBinContent(i,j); 
	//if (c1-c2 > 0) 
	  hKYvsKAdiff->SetBinContent(i,j,c1-1.0*c2); 
      }
    }
  }
  nx = hNLvsKYmix1->GetXaxis()->GetNbins(); 
  ny = hNLvsKYmix->GetYaxis()->GetNbins(); 
  for (int i=1;i<=nx;i++) { 
    for (int j=1;j<=ny;j++) { 
      c1 =hNLvsKYmix1->GetBinContent(i,j); 
      if (c1 > 0) { 
	c2 =hNLvsKYmix->GetBinContent(i,j); 
	//if (c1-c2 > 0) 
	  hNLvsKYdiff1->SetBinContent(i,j,c1-1.0*c2); 
      }
    }
  }
  nx = hNLvsKAmix1->GetXaxis()->GetNbins(); 
  ny = hNLvsKAmix->GetYaxis()->GetNbins(); 
  for (int i=1;i<=nx;i++) { 
    for (int j=1;j<=ny;j++) { 
      c1 =hNLvsKAmix1->GetBinContent(i,j); 
      if (c1 > 0) { 
	c2 =hNLvsKAmix->GetBinContent(i,j); 
	//if (c1-c2 > 0) 
	  hNLvsKAdiff1->SetBinContent(i,j,c1-1.0*c2); 
      }
    }
  }
  nx = hKYvsKAmix1->GetXaxis()->GetNbins(); 
  ny = hKYvsKAmix->GetYaxis()->GetNbins(); 
  for (int i=1;i<=nx;i++) { 
    for (int j=1;j<=ny;j++) { 
      c1 =hKYvsKAmix1->GetBinContent(i,j); 
      if (c1 > 0) { 
	c2 =hKYvsKAmix->GetBinContent(i,j); 
	//if (c1-c2 > 0) 
	  hKYvsKAdiff1->SetBinContent(i,j,c1-1.0*c2); 
      }
    }
  }
}


void run(Int_t run=999)
{
  char name[100]; 
  TChain *fChain = new TChain("nl");

  for(int i=0;i<2000;i++) {
    run2run[i]=0;
  }

  int cnt=-1;

  RunNumber=run;
 
  read_calibration();

  // for(int i=2614;i<2634;i++) {
  //   sprintf(name,"../eos/SMDAQ%04d.root",i);
  //   fChain->Add(name);
  //   cnt++;
  //   printf(" file: %s\n",name);
  //   run2run[cnt]=i;
  // }
  for(int i=2839;i<3167;i++) {  //  EOS2 2839-3166
    sprintf(name,"rootfiles/SMDAQ%04d.root",i);
    fChain->Add(name);
    cnt++;
    printf(" file: %s\n",name);
    run2run[cnt]=i;
  }

  //ntree->Branch("nrun", &nrun , "nrun/I");
  //ntree->Branch("nevt", &nevt , "nevt/I");
  //ntree->Branch("nevttot", &nevttot , "nevttot/I");
  //ntree->Branch("nmult", &nmult , "nmult/I");
  
  loop(fChain);
  //  out_hist();

  subhisto();

  //ntree->AutoSave();
  //nfile->Close();
}

