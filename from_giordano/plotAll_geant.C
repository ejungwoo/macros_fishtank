
TGraph* GetMCP(TString filename="", Int_t type = 0);
TGraph* GetBichsel(TString species= "", TString filename="", Int_t type = 0);
TGraph* correction(TString filename="", TString filename2="");
TGraph* correction(TString species = "", TString filename="", TString filename2="");  
TGraph* correctionSpline(TString filename, TString filename2, TString species, const Int_t ni);
  
void plotAll_geant()
{

  gStyle -> SetOptStat(0);
  gStyle -> SetTitleOffset(2.4,"y");
  gStyle -> SetTitleOffset(1.2,"x");
  gStyle -> SetPadLeftMargin(0.17);
  gStyle -> SetLabelSize(0.02,"X");
  
  // dedx without secondary electrons
  //  TGraph* graphMCPE = GetMCP("geant4piminus_75_events.txt",0);
  //  TGraph* graphMCPE = GetMCP("geant4protons_75_events.txt",0);
  //  TGraph* graphMCPE = GetMCP("geant4deuterons_75_events.txt",0);
  TGraph* graphMCPE = GetMCP("geant4tritons_75_events.txt",0);
  graphMCPE -> SetMarkerStyle(24);
  graphMCPE -> SetMarkerColor(1);
  // dedx with secondary electrons  
  /*
  //  TGraph* graphMCPE2 = GetMCP("geant4piminus_75_events.txt",1);
  TGraph* graphMCPE2 = GetMCP("geant4protons_75_events.txt",1);
  graphMCPE2 -> SetMarkerStyle(24);
  graphMCPE2 -> SetMarkerColor(3);
  */
  // dedx from bichsel new
  //  TGraph* graphMCBP = GetBichsel("pi", "bgamma_pid_pions.dat", 0);
  //  TGraph* graphMCBP = GetBichsel("p", "bgamma_pid.dat", 0);
  //  TGraph* graphMCBP = GetBichsel("d", "bgamma_pid.dat", 0);
  TGraph* graphMCBP = GetBichsel("t", "bgamma_pid.dat", 0);
  graphMCBP -> SetMarkerStyle(29);
  graphMCBP -> SetMarkerColor(2);
  // dedx from bichsel old
  /*
  TGraph* graphMCBP2 = GetBichsel("", "bichsel_pion.txt", 1);
  //  TGraph* graphMCBP2 = GetBichsel("", "bichsel_p.txt", 1);  
  graphMCBP2 -> SetMarkerStyle(22);
  graphMCBP2 -> SetMarkerColor(4); 
  */
  
  TLegend *legendE = new TLegend(0.55,0.80,0.9,0.9);
  /*
  legendE -> AddEntry(graphMCPE, "dedx pion (MC)","PL");
  //  legendE -> AddEntry(graphMCPE2, "dedx pion (MC) with secondary","PL");
  legendE -> AddEntry(graphMCBP, "dedx pion Bichsel new","PL");
  //  legendE -> AddEntry(graphMCBP2, "dedx pion MC Bichsel","PL");    
  */
  /*
  legendE -> AddEntry(graphMCPE, "dedx proton (MC)","PL");
  //  legendE -> AddEntry(graphMCPE2, "dedx proton (MC) with secondary","PL");
  legendE -> AddEntry(graphMCBP, "dedx proton Bichsel new","PL");
  //  legendE -> AddEntry(graphMCBP2, "dedx proton (MC) Bichsel old","PL");  
  legendE -> SetFillColor(0);
  */
  /*
  legendE -> AddEntry(graphMCPE, "dedx deuteron (MC)","PL");
  //  legendE -> AddEntry(graphMCPE2, "dedx proton (MC) with secondary","PL");
  legendE -> AddEntry(graphMCBP, "dedx deuteron Bichsel new","PL");
  //  legendE -> AddEntry(graphMCBP2, "dedx proton (MC) Bichsel old","PL");  
  legendE -> SetFillColor(0);
  */
  legendE -> AddEntry(graphMCPE, "dedx triton (MC)","PL");
  //  legendE -> AddEntry(graphMCPE2, "dedx proton (MC) with secondary","PL");
  legendE -> AddEntry(graphMCBP, "dedx triton Bichsel new","PL");
  //  legendE -> AddEntry(graphMCBP2, "dedx proton (MC) Bichsel old","PL");  
  legendE -> SetFillColor(0);  
  
  TCanvas* cvsP = new TCanvas("cvsP","cvsP",700,700);
  cvsP -> SetGrid();
  //  graphMCPE  -> SetMinimum(0.0001);
  //  graphMCPE2  -> SetMaximum(0.011);
  graphMCPE  -> GetXaxis() -> SetTitle("Momentum (MeV/c)");
  graphMCPE  -> GetYaxis() -> SetTitle("Energy Loss (MeV/mm)");
  graphMCPE  -> Draw("AP");
  //  graphMCPE2   -> Draw("P SAME");
  graphMCBP   -> Draw("P SAME");
  //  graphMCBP2  -> Draw("P SAME");  
  legendE     -> Draw("SAME");

  ///////////////////////////////////////////////////////////
  // plot for correction
  ///////////////////////////////////////////////////////////
  
  TGraph* graphCorrection = correction("pi","bgamma_pid_pions.dat","geant4piminus_75_events.txt");  
  //  TGraph* graphCorrection = correction("p","bgamma_pid.dat","geant4protons_75_events.txt");
  //  TGraph* graphCorrection = correction("d","bgamma_pid.dat","geant4deuterons_75_events.txt");
  //  TGraph* graphCorrection = correction("t","bgamma_pid.dat","geant4tritons_75_events.txt");
  graphCorrection -> SetMinimum(0.8);
  //  graphCorrection -> SetMaximum(1.4);
  graphCorrection -> SetMarkerStyle(24);
  graphCorrection -> SetMarkerColor(1);

  
  TGraph* graphCorrection2 = correctionSpline("bgamma_pid_pions.dat","geant4piminus_75_events.txt","pi",69);
  //  TGraph* graphCorrection2 = correctionSpline("bgamma_pid.dat","geant4protons_75_events.txt","p",75);
  //  TGraph* graphCorrection2 = correctionSpline("bgamma_pid.dat","geant4deuterons_75_events.txt","d",75);
  //  TGraph* graphCorrection2 = correctionSpline("bgamma_pid.dat","geant4tritons_75_events.txt","t",75);
  graphCorrection2 -> SetMarkerStyle(29);
  graphCorrection2 -> SetMarkerColor(2);
  
  TCanvas* cvsP2 = new TCanvas("cvsP2","cvsP2",700,700);
  cvsP2 -> SetGrid();
  graphCorrection  -> Draw("AP");
  graphCorrection  -> GetXaxis() -> SetTitle("Momentum (MeV/c)");
  graphCorrection  -> GetYaxis() -> SetTitle("Bichels/GEANT4 Energy Loss");  
  graphCorrection2  -> Draw("P SAME");
  
}

TGraph* correction(TString species = "", TString filename="", TString filename2=""){
  Int_t i=0;
  const Int_t n=69;
  //  const Int_t n=75;
  Double_t betagamma, eloss, dEdxB[n], pB[n], dEdxG[n], pG[n], ratio[n];
  TGraph* graph = new TGraph();
  ifstream file(filename);
  ifstream file2(filename2);
  // bichsel file
  for(Int_t i=0; i<n; i++){    
    file>>betagamma>>eloss;
    if (species == "p")
      pB[i]=betagamma*938.27;
    else if(species == "pi")
      pB[i]=betagamma*139.57;
    else if(species == "d")
      p=betagamma*(938.27+939.57);
    else if(species == "t")
      p=betagamma*(938.27+2*939.57);
    else{
      std::cout << "It needs implementation, sorry!" << std::endl;
      exit(0);
    }
    
    dEdxB[i]=eloss/10000; // keV/cm to MeV/mm 
    //    cout << pB[i] << " " << dEdxB[i] << endl;
  }
  // geant file
  Double_t p, avgdEdx, avgdEdx2, dEdx, dEdx2;
  for(Int_t i=0; i<n; i++){
    file2>>p>>avgdEdx>>avgdEdx2>>dEdx>>dEdx2;
    pG[i]=p;
    dEdxG[i]=dEdx;
    ratio[i]=dEdxB[i]/dEdxG[i];
    //    cout << pG[i] << " " << dEdxG[i] << " " << ratio[i] << endl;
    if (ratio[i]<1.5)
      graph->SetPoint(i,pG[i],ratio[i]);
  }
  
  return graph;
}

/*
TGraph* correction(TString filename="", TString filename2="")
{
  //  const Int_t n = 69;
  const Int_t n = 75;  
  Int_t i=0;
  Double_t pb, eloss, dEdxB[n], pB[n], dEdxG[n], pG[n], ratio[n];
  TGraph* graph = new TGraph();
  ifstream file(filename);
  ifstream file2(filename2);
  // bichsel file
  for(Int_t i=0; i<n; i++){
    file>>pb>>eloss;
    pB[i]=pb;
    dEdxB[i]=eloss;
    //    cout << pB[i] << " " << dEdxB[i] << endl;
  }
  // geant file
  Double_t p, avgdEdx, avgdEdx2, dEdx, dEdx2;
  for(Int_t i=0; i<n; i++){    
    file2>>p>>avgdEdx>>avgdEdx2>>dEdx>>dEdx2;
    pG[i]=p;
    dEdxG[i]=dEdx;
    ratio[i]=dEdxB[i]/dEdxG[i];
    //    cout << pG[i] << " " << dEdxG[i] << " " << ratio[i] << endl;
    if (ratio[i]<2.)
      graph->SetPoint(i,pG[i],ratio[i]);
  }

  return graph;
}
*/

TGraph* GetBichsel(TString species = "", TString filename="", Int_t type = 0) // Bichsel dE/dx
{
  Int_t i=0;
  Double_t betagamma, dEdx, p;
  TGraph* graph = new TGraph();
  ifstream file(filename);

  if(type == 0){
    while(file>>betagamma>>dEdx){
      if (species == "p")
	p=betagamma*938.27;
      else if(species == "pi")
	p=betagamma*139.57;
      else if(species == "d")
	p=betagamma*(938.27+939.57);
      else if(species == "t")
	p=betagamma*(938.27+2*939.57);
      else{
	std::cout << "It needs implementation, sorry!" << std::endl;
	exit(0);
      }
      
      dEdx=dEdx/10000; // keV/cm to MeV/mm
      graph->SetPoint(i++,p,dEdx);
    }
  }
  else if (type == 1){
    while(file>>p>>dEdx){
      graph->SetPoint(i++,p,dEdx);      
    }
  }
  
  return graph;
}

TGraph* GetMCP(TString filename="", Int_t type = 0) // MC Proton dE/dx
{
  Int_t i=0;
  Double_t p, avgdEdx, avgdEdx2, dEdx, dEdx2;
  TGraph* graph = new TGraph();
  ifstream file(filename);

  // dEdx, secondaries not included
  if(type==0) 
    while(file>>p>>avgdEdx>>avgdEdx2>>dEdx>>dEdx2) 
      graph->SetPoint(i++,p,dEdx);

  // dEdx, secondaries include
  else if(type==1)
    while(file>>p>>avgdEdx>>avgdEdx2>>dEdx>>dEdx2)     
      graph->SetPoint(i++,p,dEdx+dEdx2);

  // average dEdx, secondaries not included
  else if(type==2)
    while(file>>p>>avgdEdx>>avgdEdx2>>dEdx>>dEdx2)  
      graph->SetPoint(i++,p,avgdEdx);

  // average dEdx, secondaries included
  else if(type==3)
    while(file>>p>>avgdEdx>>avgdEdx2>>dEdx>>dEdx2)    
      graph->SetPoint(i++,p,avgdEdx+avgdEdx2);

  return graph;
}

TGraph* correctionSpline(TString filename, TString filename2, TString species, const Int_t ni){
  Int_t i=0;
  Double_t betagamma, eloss, dEdxB[ni], pB[ni], dEdxG[ni], pG[ni], ratio[ni], pmin, pmax;
  TGraph* graph = new TGraph();
  ifstream file(filename);
  ifstream file2(filename2);
  // bichsel file
  for(Int_t i=0; i<ni; i++){
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
      exit(0);
    }
    
    dEdxB[i]=eloss/10000; // keV/cm to MeV/mm

  }

  pmin = pB[0];
  pmax = pB[ni-1];

  cout << pmin << " " << pmax << endl;
  
  // geant file
  Int_t j=0;
  Double_t pG2[63], ratio2[63];
  Double_t p, avgdEdx, avgdEdx2, dEdx, dEdx2;
  for(Int_t i=0; i<ni; i++){
    file2>>p>>avgdEdx>>avgdEdx2>>dEdx>>dEdx2;
    pG[i]=p;
    dEdxG[i]=dEdx;
    ratio[i]=dEdxB[i]/dEdxG[i];

    //    cout << pG[i] << " " << ratio[i] << endl;
    if(ratio[i]<1.5){
      pG2[j]=pG[i];
      ratio2[j]=ratio[i];
      j++;
    }
    //    graph->SetPoint(i,pG[i],ratio[i]);
    //    j++;
  }
  //  cout << j << endl;
  
  // interpolation
  const Int_t nb = 1000;
  Double_t ix[nb], iy[nb];

  ROOT::Math::Interpolator myinter(ni, ROOT::Math::Interpolation::kCSPLINE );
  //  myinter.SetData(ni, pG, ratio);
  myinter.SetData(63, pG2, ratio2);
  
  for (Int_t i=0; i< nb; i++){
    ix[i]  = (Double_t) i*(pmax-pmin)/(nb-1) + pmin;
    iy[i] = myinter.Eval(ix[i]);
    graph->SetPoint(i,ix[i],iy[i]);
  }  
  
  return graph;

}
