#include <iostream>
#include "TChain.h"
#include "TCanvas.h"
#include "cmath"

Double_t PI=TMath::Pi();
Double_t NLt=29.597;
void JonNeuland(){
  Double_t maxNLr=100;
  gStyle->SetOptFit();
  Int_t runs[]={2841};
  /*
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
  */
  Int_t numData = sizeof(runs)/sizeof(Int_t);
  Int_t events[numData];
  Int_t tracks[numData];
  char name[180];

  TChain *fChainNL = new TChain("NL");
  TChain *fChainSTPC = new TChain("dedx");
  //define histograms
  auto *histNLXY=new TH2D("histNLXY","TPC projection to NeuLAND;x'(mm);y'(mm)",100,-2000,2000,100,-2000,2000);
  auto *histNLXZ=new TH2D("histNLXZ","TPC projection to NeuLAND;x'(mm);yz'(mm)",1000,-2000,2000,1000,-2000,2000);
  auto *histnlxy=new TH2D("histnlxy","NeuLAND X/Y distribution;nlx(mm);nly(mm)",100,-2000,2000,100,-2000,2000);
  auto *histAOE=new TH2D("histAOE","TPC angle of emission;#theta (deg);#phi(deg)",90,0,90,360,0,360);
  auto *histDeDx= new TH2D("histDeDx","dE/dx vs. p/Z;p/Z (MeV/c);dE/dx (ADC/mm)", 300, -2500, 2500, 300, 0, 900);
  
  auto *histNLz=new TH1D("histNLz","z distribution TPC projection;z'",1000,-500,500);  
  auto *histNLr=new TH1D("histNLr","distance between points",1000,0,10000);

  auto *histNLXvZ=new TH2D("histNLXvZ","Projected Z vs NeuLAND X;NL x(mm);projected Z (mm)",240,-1200,1200,200,-1000,1000);
  auto *histXvX=new TH2D("histXvX","NeuLAND X vs x';x'(mm);X (NeuLAND) (mm)",100,-3000,3000,100,-3000,3000);
  auto *histYvY=new TH2D("histYvY","NeuLAND Y vs Y';y'(mm);Y (NeuLAND) (mm)",100,-2000,2000,100,-2000,2000);
  auto *histDYvDX=new TH2D("histDYvDX","dy v dx;dx;dy",100,-1,1,100,-1,1);
  
  
  //define NL branches
  Double_t nlx; fChainNL->SetBranchAddress("NLx", &nlx);
  Double_t nly; fChainNL->SetBranchAddress("NLy", &nly);
  Int_t nlrun; fChainNL->SetBranchAddress("NLrun", &nlrun);
  Int_t nlevent; fChainNL->SetBranchAddress("NLevent", &nlevent);
  
  
  //define branches
  Int_t run; fChainSTPC->SetBranchAddress("run", &run);
  Int_t eventid; fChainSTPC->SetBranchAddress("eventid", &eventid);
  Double_t NLx; fChainSTPC->SetBranchAddress("NLx", &NLx);
  Double_t NLy; fChainSTPC->SetBranchAddress("NLy", &NLy);
  Double_t NLz; fChainSTPC->SetBranchAddress("NLz", &NLz);
  Double_t mom; fChainSTPC->SetBranchAddress("mom", &mom);
  Int_t charge; fChainSTPC->SetBranchAddress("charge", &charge);
  Double_t dedx; fChainSTPC->SetBranchAddress("dedx", &dedx);
  Double_t dx; fChainSTPC->SetBranchAddress("dx", &dx);
  Double_t dy; fChainSTPC->SetBranchAddress("dy", &dy);
  Double_t dz; fChainSTPC->SetBranchAddress("dz", &dz);
  Double_t vx; fChainSTPC->SetBranchAddress("vx", &vx);
  Double_t vy; fChainSTPC->SetBranchAddress("vy", &vy);
  Double_t vz; fChainSTPC->SetBranchAddress("vz", &vz);
  Double_t dist; fChainSTPC->SetBranchAddress("dist", &dist);
  Double_t ndf; fChainSTPC->SetBranchAddress("ndf", &ndf);
  
  Int_t  TotalEvents=0;
  Int_t LocalEvents=0;
  Int_t  TotalTracks=0;
  Int_t LocalTracks=0;
  Double_t theta, phi;
  Int_t inNL=0;
  Double_t NLr;
  Double_t NLxp, NLyp, NLzp;


  for (int i = 0;i<numData;i++) {

    sprintf(name,"/mnt/spirit/analysis/changj/SpiRITROOT.JWfix/macros/data/analysisCode-All-noLayerCut-GC-DS-GiordanoCommentOut-bShift-By-JWfix/dedxROOT/dedxSn132-All-noLayerCut_%d.root",i);
    //sprintf(name,"/mnt/spirit/analysis/barneyj/SpiRITROOT.latest/macros/data/analysisCode-Sn132-noLayerCut-GC-DS-GiordanoCommentOut-bShift-By-anodeCalib/dedxROOT2/dedxSn132-All-noLayerCut_%d.root",i);
    //sprintf(name,"/mnt/spirit/analysis/barneyj/SpiRITROOT.latest/macros/data/analysisCode-Sn132-noLayerCut-GC-DS-GiordanoCommentOut-bshift-By-anodeCalib/dedxROOT/dedxSn132-All-noLayerCut_%d.root",i);
    
    fChainSTPC->Add(name);
    //printf(" file: %s\n",name);
    LocalTracks=fChainSTPC->GetEntries()-TotalTracks;
    tracks[i]=LocalTracks;
    TotalTracks=fChainSTPC->GetEntries();

    
    //sprintf(name,"/mnt/spirit/analysis/gasparic/rootfiles/SMDAQ%04d.root",runs[i]);
    fChainNL->Add("NLout.root");
    LocalEvents=fChainNL->GetEntries()-TotalEvents;
    events[i]=LocalEvents;
    TotalEvents=fChainNL->GetEntries();
    //cnt++;
    //printf(" file: %s\n",name);
    //run2run[cnt]=runs[i];

  }

  Int_t cnt=0;
  Int_t cnt2=0;
  fChainSTPC->GetEntry(cnt2);//get first event for STPC here to set "eventid" variable
  cout << eventid << endl;
  //note that eventid starts from 1
  for(int iRun=0;iRun<2/*numData*/;iRun++){
    LocalEvents=events[iRun];
    Int_t cnt3=0;
    for(int iEvent=1;iEvent<=LocalEvents;iEvent++){
      if(iEvent%10000==0) cout << "Current: " << iEvent << "/" << LocalEvents << endl;
      
      fChainNL->GetEntry(cnt);
      //should do NeuLAND analysis here
      //printf("NeuLAND: %i\n",iEvent);
      //if( nlx>-9998) histnlxy->Fill(nlx*10,nly*10);
      
      while(eventid==iEvent&& cnt2<tracks[iRun] ){
	fChainSTPC->GetEntry(cnt2);
	cnt2++;
	//NLr=sqrt((NLx-nlx*10)*(NLx-nlx*10)+(NLy-nly*10)*(NLy-nly*10)+NLz*NLz);
	if(eventid!=iEvent) continue;//if events are not matching, skip to end of loop 
	//if(eventid>LocalEvents-1 || eventid<2) printf("%i,%i\n",eventid,iEvent);//this is just a check that the first and last events are analyzed
	//track by track analysis here
	//if(NLz>-9998)
	phi=atan2(dy,dx)*180/PI+180*(1-abs(atan2(dy,dx))/atan2(dy,dx));
	theta=acos(dz/sqrt(dz*dz+dy*dy+dx*dx))*180/PI;
	if( nlx>-9998) histnlxy->Fill(nlx*10,nly*10);
	if(nlx>-9998 && NLz>0){
	//if(NLr<1500){
	  //if((phi>300 || phi<60)){
	  if( vz<-9.49569&&vz>-12.80121&&vx>-15&&vx<15&&vy<-206.06&& vy>-246.06 /*&& dist<5 && ndf>30*/){	
	    histNLXZ->Fill(NLx,NLz);
	    NLyp=NLy+213.3;
	    NLz=NLz-8025.7;
	    NLx=NLx-4225.9;
	    histNLXZ->Fill(NLx,NLz);
	    //histDzvx->Fill(NLx,NLz);	    
	    NLxp=NLx*TMath::Cos(NLt*PI/180.)-NLz*TMath::Sin(NLt*PI/180.);
	    NLzp=NLx*TMath::Sin(NLt*PI/180.)+NLz*TMath::Cos(NLt*PI/180.);

	    //histNLXZ->Fill(NLxp,NLzp);
	    NLr=sqrt((NLxp-nlx*10)*(NLxp-nlx*10)+(NLyp-nly*10)*(NLyp-nly*10)+NLzp*NLzp);
	    if(NLr<maxNLr){
	      histNLXY->Fill(NLx,NLy);
	      histNLXZ->Fill(NLxp,NLzp);
	      histNLXvZ->Fill(nlx*10,NLzp);

	      histNLz->Fill(NLzp);
	      histNLr->Fill(NLr);
	      histAOE->Fill(theta,phi);
	      histDeDx->Fill(mom/charge,dedx);
	    
	      histXvX->Fill(NLxp,nlx*10);
	      histYvY->Fill(NLyp,nly*10);
	      histDYvDX->Fill(dx,dy);	    
	      if(NLxp>-1500 && NLxp<1500 /*&& NLy>-1000 && NLy<1000*/) inNL++;
	    }
	    }
	}
      }//while events are matching

      cnt++;
    }//loop thru events
  }//loop thru runs
  //histNLXY->Draw("colz");

  
  TLine *NLtop=new TLine(-1250,1000,1250,1000);
  TLine *NLleft=new TLine(-1250,1000,-1250,-1000);
  TLine *NLright=new TLine(1250,1000,1250,-1000);
  TLine *NLbottom=new TLine(-1250,-1000,1250,-1000);

  auto cvs1=new TCanvas("cvs1","",0,0,1600,1000);
  cvs1->Divide(2,2);
  cvs1->cd(1);
  histNLXY->Draw("colz");
  
  NLtop->Draw("same");
  NLbottom->Draw("same");
  NLright->Draw("same");
  NLleft->Draw("same");
  cvs1->cd(3);
  histnlxy->Draw("colz");
  
  cvs1->cd(2);
  histDeDx->Draw("colz");
  //histNLz->Draw();
  //histXvX->Draw("colz");
  cvs1->cd(4);
  
  histNLr->Draw();
  //histYvY->Draw("colz");
  cvs1->SaveAs(Form("distance%f.png",maxNLr));  
  auto cvs2=new TCanvas("cvs2","",0,0,1600,1000);
  cvs2->Divide(2,2);
  cvs2->cd(1);
  histXvX->Draw("colz");
  cvs2->cd(2);
  histYvY->Draw("colz");
  cvs2->cd(3);
  histNLXvZ->Draw();
  //histNLXvZ->Fit("pol1");
  cvs2->cd(4);
  //histNLXZ->Draw("colz");
  histDYvDX->Draw("colz");
  cout << "in Neuland: " << inNL << endl;  
  cvs2->SaveAs(Form("corr%f.png",maxNLr));
  

}


