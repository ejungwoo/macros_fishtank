#include "style.h"
#include <cmath>
using namespace style;

Double_t EvalBichsel(TString species, Double_t value)
{
  if (species == "pi")
    return 1.52844-0.00533277*value+1.95981e-05*pow(value,2)-3.28283e-08*pow(value,3)+2.07723e-11*pow(value,4);
  else if(species == "p")
    return 1.68615-0.00124295*value+8.34557e-07*pow(value,2)-2.447e-10*pow(value,3)+2.57546e-14*pow(value,4);
  else if(species == "d")
    return 1.74617-0.0007082*value+2.50968e-07*pow(value,2)-3.93875e-11*pow(value,3)+2.24589e-15*pow(value,4);
  else if(species == "t")
    return 1.78826-0.000504648*value+1.18876e-07*pow(value,2)-1.21566e-11*pow(value,3)+4.4748e-16*pow(value,4);
  else{
    std::cout << "It needs implementation, sorry!" << std::endl;
    exit(0);}
}

ROOT::Math::Interpolator* GetInterpolator(TString species)
{
  const volatile Int_t n = 75;
  int &nn = const_cast <int &> (n);
  if (species == "pi") nn = 69; 

  Double_t betagamma, eloss, pmin, pmax;
  Double_t p, avgdEdx, avgdEdx2, dEdx, dEdx2;
  Double_t pB[n], dEdxB[n], pG[n], dEdxG[n], ratio[n];

  ifstream file;
  if (species == "pi") file.open("./bgamma_pid_pions.dat");
  else file.open("./bgamma_pid.dat");

  ifstream file2; 

  // Bichsel file  
  for(Int_t i=0; i<n; i++){
    file>>betagamma>>eloss;
    if (species == "p")
      pB[i]=betagamma*938.27;
    else if(species == "pi")
      pB[i]=betagamma*139.57;
    else if(species == "d")
      pB[i]=betagamma*(938.27+939.57);
    else if(species == "t")
      pB[i]=betagamma*(938.27+2*939.57);
    else{
      std::cout << "It needs implementation, sorry!" << std::endl;
      exit(0);}

    dEdxB[i]=eloss/10000; // keV/cm to MeV/mm    

  }

  fPmin = pB[0];
  fPmax = pB[n-1];

  // support variables for pion case
  Int_t j=0;
  Double_t pG2[63], ratio2[63];      

  // load geant4 file for calculating the ration
  if (species == "p")
    file2.open("./geant4protons_75_events.txt");
  else if(species == "pi")
    file2.open("./geant4piminus_75_events.txt");  
  else if(species == "d")
    file2.open("./geant4deuterons_75_events.txt");      
  else if(species == "t")
    file2.open("./geant4tritons_75_events.txt");
  else{
    std::cout << "It needs implementation, sorry!" << std::endl;
    exit(0);}

  for(Int_t i=0; i<n; i++){
    file2>>p>>avgdEdx>>avgdEdx2>>dEdx>>dEdx2;
    pG[i]=p;
    dEdxG[i]=dEdx;
    ratio[i]=dEdxB[i]/dEdxG[i];
    if (species == "pi"){
      if(ratio[i]<1.5){
        pG2[j]=pG[i];
        ratio2[j]=ratio[i];
        j++;                    
      }
    }
  }

  // interpolation
  ROOT::Math::Interpolator* myinter = new ROOT::Math::Interpolator(n, ROOT::Math::Interpolation::kCSPLINE );
  if (species == "pi") myinter->SetData(63, pG2, ratio2);
  else myinter->SetData(n, pG, ratio);

  return myinter;

} 

void pidline()
{
  auto file = new TFile("justpid.root","read");
  file -> Get("canvas-0") -> Draw();
  gstat(0);

  for (auto useSpline : {1,0})
  {
    auto g_pi = new TGraph(); g_pi -> SetLineColor(kRed);
    auto g_p  = new TGraph(); g_p  -> SetLineColor(kBlack);
    auto g_d  = new TGraph(); g_d  -> SetLineColor(kBlue);
    auto g_t  = new TGraph(); g_t  -> SetLineColor(kGreen);

    TH2D *frame = nullptr;
    if (useSpline) {
      auto in_pi = GetInterpolator("pi");
      auto in_p  = GetInterpolator("p");
      auto in_d  = GetInterpolator("d");
      auto in_t  = GetInterpolator("t");
      for (auto x=10.; x<1500.; x+=1) { double val=in_pi->Eval(x); if (!std::isnan(val)) g_pi->SetPoint(g_pi->GetN(),x,val); }
      for (auto x=10.; x<1500.; x+=1) { double val=in_p ->Eval(x); if (!std::isnan(val)) g_p ->SetPoint(g_p ->GetN(),x,val); }
      for (auto x=10.; x<1500.; x+=1) { double val=in_d ->Eval(x); if (!std::isnan(val)) g_d ->SetPoint(g_d ->GetN(),x,val); }
      for (auto x=10.; x<1500.; x+=1) { double val=in_t ->Eval(x); if (!std::isnan(val)) g_t ->SetPoint(g_t ->GetN(),x,val); }
      cc("c");
      frame = new TH2D("frame","interpolation;p;dE/dx",100,0,1500,100,0.8,2);
    } else {
      for (auto x=10.; x<1500.; x+=1) g_pi->SetPoint(g_pi->GetN(),x,EvalBichsel("pi",x));
      for (auto x=10.; x<1500.; x+=1) g_p ->SetPoint(g_p ->GetN(),x,EvalBichsel("p" ,x));
      for (auto x=10.; x<1500.; x+=1) g_d ->SetPoint(g_d ->GetN(),x,EvalBichsel("d" ,x));
      for (auto x=10.; x<1500.; x+=1) g_t ->SetPoint(g_t ->GetN(),x,EvalBichsel("t" ,x));
      cc("c2");
      frame = new TH2D("frame2","fit function;p;dE/dx",100,0,1500,100,0.8,2);
    }

    frame -> Draw();
    g_pi  -> Draw("lsame");
    g_p   -> Draw("lsame");
    g_d   -> Draw("lsame");
    g_t   -> Draw("lsame");
  }
}
