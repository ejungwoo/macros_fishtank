auto mpi = 139.57018;
auto mp = 938.2720813;
auto mn = 939.565346;
auto md = 1875.612762;
auto mt = 2808.921112;
auto mhe3 = 2808.39132;
auto mal = 3727.379378;
auto NLtheta=29.579;
auto NLr=8561.05;
auto TPC_S_x=-2.;
auto TPC_S_y=213.3;
auto TPC_S_z=-580.4;

void dedxSn132_all(Int_t iRun = 0) {
  /*
  Int_t runs[] = {2843, 2844, 2845, 2846, 2848, 2849, 2850, 2851, 2852, // 9
                  2855, 2856, 2857, 2858, 2859, 2860, 2861,             // 7
                  2877, 2878, 2879, 2880, 2881, 2882, 2883, 2884,       // 8
                  2887, 2888, 2889, 2890, 2891, 2892, 2893, 2894};      // 8
  Int_t numFiles[] = {42, 10, 38, 38, 5, 37, 36, 38, 42,
                      38, 6, 43, 45, 38, 43, 43,
                      46, 44, 37, 41, 38, 44, 46, 31,
                      46, 46, 45, 44, 38, 7, 37, 48}; 
  */
  Int_t runs[]={2841,2843};
  /*  Int_t runs[] = {2841, 2843, 2844, 2845, 2846, 2848, 2849, 2850, 2851, 2852, 2855, 
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
  Int_t numFiles[]={43,42};
  /*
 Int_t numFiles[] = {43, 42, 10, 38, 38, 5, 37, 36, 38, 42, 38, 6, 43, 45, 38, 43, 43, 47, 
                      46, 44, 37, 41, 38, 44, 46, 31, 46, 46, 45, 44, 38, 7, 37, 48, 33, 58, 
                      36, 44, 45, 43, 52, 42, 46, 35, 10, 42, 43, 47, 46, 44, 39, 30, 43, 45, 
                      47, 21, 41, 42, 44, 44, 46, 45, 45, 45, 42, 45, 44, 46, 45, 47, 46, 46, 
                      44, 45, 46, 45, 45, 43, 44, 49, 50, 31, 15, 42, 44, 49, 45, 45, 44, 45, 
                      11, 45, 47, 45, 24, 43, 44, 3, 48, 49, 45, 37, 48, 47, 44, 17, 37, 47, 
                      26, 41, 41, 40, 23};
 
  */
 Int_t numData = sizeof(runs)/sizeof(Int_t);
  TChain **chain = new TChain*[numData];
  TChain **beam = new TChain*[numData];
  for (Int_t i = 0; i < numData; i++) {
    chain[i] = new TChain("cbmsim");
    beam[i] = new TChain("beam");
  }

/*
  auto pCutFile = new TFile("analysisCode/pdFuncCut.root");
  auto pFuncCut = (TCutG *) pCutFile -> FindObjectAny("pFuncCut");

  auto piCutFile = new TFile("analysisCode/pionCut.root");
  auto pimCut = (TCutG *) piCutFile -> FindObjectAny("pim");
  auto pipCut = (TCutG *) piCutFile -> FindObjectAny("pip");
  */

//  auto pipBgCutFile = new TFile("analysisCode/pipBg.root");
//  auto pipBgCut = (TCutG *) pipBgCutFile -> FindObjectAny("pipBg");
  auto *histXZ=new TH2D("histXZ","",1000,-10000,2000,1000,-10000,2000);
  TClonesArray *trackArray = nullptr;
  TClonesArray *vertexArray = nullptr;
  STEventHeader *header = nullptr;
  Bool_t sigma10, sigma15, sigma20, sigma20z, proton, pip, pim, pipbg;
  Double_t vx, vy, vz, projmevu, projpx, projpy, projpz, proja, projb, projx, projy, projz;

  for (Int_t iRun = 0; iRun < numData; iRun++) {
    for (Int_t iFile = 0; iFile < numFiles[iRun]; iFile++) {
      //chain[iRun] -> AddFile(Form("Sn132-noLayerCut-GC-DS-GiordanoCommentOut-bShift-By-anodeCalib/run%d_s%d.reco.v1.04.root", runs[iRun], iFile));
//      chain[iRun] -> AddFile(Form("Sn132-noLayerCut-GC-DS-GiordanoCommentOut-bShift-By/run%d_s%d.reco.v1.04.root", runs[iRun], iFile));
      chain[iRun] -> AddFile(Form("run%d_s%d.reco.v1.04.root", runs[iRun], iFile));
      
    }
    beam[iRun] -> AddFile(Form("/mnt/spirit/analysis/changj/BeamAnalysis/macros/refined/Sn132_all/run%d.refined.root", runs[iRun]));
  }

  for (Int_t i = 0; i < numData; i++) {
    chain[i] -> SetBranchAddress("STEventHeader", &header);
    chain[i] -> SetBranchAddress("STRecoTrack", &trackArray);
    chain[i] -> SetBranchAddress("STVertex", &vertexArray);
    beam[i] -> SetBranchAddress("sigma10", &sigma10);
    beam[i] -> SetBranchAddress("sigma15", &sigma15);
    beam[i] -> SetBranchAddress("sigma20", &sigma20);
    beam[i] -> SetBranchAddress("sigma20z", &sigma20z);
    beam[i] -> SetBranchAddress("projmevu", &projmevu);
    beam[i] -> SetBranchAddress("projpx", &projpx);
    beam[i] -> SetBranchAddress("projpy", &projpy);
    beam[i] -> SetBranchAddress("projpz", &projpz);
    beam[i] -> SetBranchAddress("proja", &proja);
    beam[i] -> SetBranchAddress("projb", &projb);
    beam[i] -> SetBranchAddress("projx", &projx);
    beam[i] -> SetBranchAddress("projy", &projy);
    beam[i] -> SetBranchAddress("projz", &projz);
  }

//  auto cvs = new TCanvas("cvs", "", 800, 500);
//  auto hist = new TH2D("hist", "", 500, -500, 2000, 100, 0, 500);

//  auto file = new TFile("dedxSn132-LayerCut10to90.root", "recreate");
  //auto file = new TFile(Form("dedxSn132-All-noLayerCut_%d.root", iRun), "recreate");

  for (auto iRun = 0; iRun < numData; iRun++) {
    
  auto file = new TFile(Form("dedxROOT/dedxSn132-All-noLayerCut_%d.root", iRun), "recreate");
  
  Double_t dedx, mom, ndf, dx, dy, dz, dist, pt, pzCM, rapL, rapCM, rapCMNorm, KECM, phiL, thetaL, phiCM, thetaCM, pocavx, pocavy, pocavz, ea, NLx, NLy, NLz;
  Int_t charge, pid, stpid, run, eventid, vid, parentvid, mult;
 
  auto tree = new TTree("dedx", "");
//  tree -> SetMaxTreeSize(100000000000LL/2);
  tree -> Branch("run", &run);
  tree -> Branch("eventid", &eventid);
  tree -> Branch("dedx", &dedx);
  tree -> Branch("mom", &mom);
  tree -> Branch("charge", &charge);
  tree -> Branch("pid", &pid);
  tree -> Branch("stpid", &stpid);
  tree -> Branch("ndf", &ndf);
  tree -> Branch("dx", &dx); // mom direction
  tree -> Branch("dy", &dy);
  tree -> Branch("dz", &dz);
  tree -> Branch("vid", &vid); // vertex id 
  tree -> Branch("parentvid", &parentvid); // track's parent vertex id 
  tree -> Branch("vx", &vx); // vertex position
  tree -> Branch("vy", &vy);
  tree -> Branch("vz", &vz);
  tree -> Branch("pocavx", &pocavx); // POCA vertex position
  tree -> Branch("pocavy", &pocavy);
  tree -> Branch("pocavz", &pocavz);
  tree -> Branch("dist", &dist); // poca and vertex distance
  tree -> Branch("sigma10", &sigma10);
  tree -> Branch("sigma15", &sigma15);
  tree -> Branch("sigma20", &sigma20);
  tree -> Branch("sigma20z", &sigma20z);
  tree -> Branch("mult", &mult);
  tree -> Branch("pt", &pt);
  tree -> Branch("pzCM", &pzCM);
  tree -> Branch("KECM", &KECM);
  tree -> Branch("rapL", &rapL);
  tree -> Branch("rapCM", &rapCM);
  tree -> Branch("rapCMNorm", &rapCMNorm);
  tree -> Branch("projx", &projx);
  tree -> Branch("projy", &projy);
  tree -> Branch("projz", &projz);
  tree -> Branch("phiL", &phiL);
  tree -> Branch("thetaL", &thetaL);
  tree -> Branch("phiCM", &phiCM);
  tree -> Branch("thetaCM", &thetaCM);
  tree -> Branch("proton", &proton);
  tree -> Branch("pip", &pip);
  tree -> Branch("pim", &pim);
  tree -> Branch("pipbg", &pipbg);
  tree -> Branch("NLx", &NLx);
  tree -> Branch("NLy", &NLy);
  tree -> Branch("NLz", &NLz);
  //tree -> Branch("ea", &ea);

  Double_t amu2mev = 931.494028;
  Double_t beamZ = 50;
  Double_t beamN = 82;
  Double_t beamA = beamZ + beamN;
  Double_t targetZ = 50;
  Double_t targetN = 74;
  Double_t targetA = targetZ + targetN;
  Double_t aMass = 0.5*(mp + mn) - 8; // MeV


    run = runs[iRun];
    auto numEntries = chain[iRun] -> GetEntries();
    cout << "Run: " << runs[iRun] << " Entries: " << numEntries << endl;
    for (auto iEntry = 0; iEntry < numEntries; iEntry++) {
      proton = 0;
      pim = 0;
      pip = 0;
      pipbg = 0;

      beam[iRun] -> GetEntry(iEntry);
      chain[iRun] -> GetEntry(iEntry);
      eventid = header -> GetEventID();

      Double_t beamKE = projmevu*0.99917; // MeV
      Double_t beamM = beamA*aMass;
      Double_t beamE = beamKE*beamA + beamM;
      Double_t beamP = TMath::Sqrt(beamE*beamE - beamM*beamM); // MeV/c
      Double_t beamRapidity = 0.5*log((beamE + beamP)/(beamE - beamP));
      Double_t E_CM = beamE + targetA*aMass; 
      Double_t P_CM = beamP;
      Double_t rapidity_CM = 0.5*log((E_CM + P_CM)/(E_CM - P_CM));
      Double_t beta = -(P_CM/E_CM);
      TVector3 betaVec(0., 0., beta);

      if (iEntry%10000 == 0)
        cout << "Current: " << iEntry << "/" << numEntries << endl;

  //    if (!sigma20)
  //      continue;

//        if (!(posVertex.Z() < -9.49569 && posVertex.Z() > -12.80121))
//          continue;

//        if (!(posVertex.X() > -15 && posVertex.X() < 15 && posVertex.Y() < -206.06 && posVertex.Y() > -246.06))
//          continue;

      auto numTracks = trackArray -> GetEntries();
      mult = numTracks;
      for (auto iTrack = 0; iTrack < numTracks; iTrack++) {
        auto track = (STRecoTrack *) trackArray -> At(iTrack);

        parentvid = track -> GetVertexID();
        if (parentvid < 0) {
          vid = -99998;
          vx = -99998;
          vy = -99998;
          vz = -99998;
          dist = -99998;
        } else {
          vid = parentvid;
          auto vertex = (STVertex *) vertexArray -> At(parentvid);
          auto posVertex = vertex -> GetPos();

          vx = posVertex.X();
          vy = posVertex.Y();
          vz = posVertex.Z();

          auto pocaVertex = track -> GetPOCAVertex();
          pocavx = pocaVertex.X();
          pocavy = pocaVertex.Y();
          pocavz = pocaVertex.Z();
          dist = (pocaVertex - posVertex).Mag();
        }

//          if (track -> GetParentID() != iVert)
//            continue;

//          if ((pocaVertex - posVertex).Mag() > 5.)
//            continue;

  //      if (track -> GetNDF() < 30)
  //        continue;

        dedx = track -> GetdEdxWithCut(0, 0.7);
        charge = track -> GetCharge();
        mom = track -> GetMomentum().Mag();
        pid = STPID::GetPDG(track -> GetPID());
        stpid = track -> GetPID();

        Double_t mass = 0;
        switch (pid) {
          case 2212:
            mass = mp;
            break;

          case 211:
          case -211:
            mass = mpi;
            break;

          case 1000010020:
            mass = md;
            break;

          case 1000010030:
            mass = mt;
            break;

          case 1000020030:
            charge = 2;
            mass = mhe3;
            break;

          case 1000020040:
            charge = 2;
            mass = mal;
            break;

          default:
            mass = 0;
            break;
        }

        auto momVec = track -> GetMomentum();

        if (momVec.Z() < 0)
          momVec = -momVec;

        momVec.RotateY(-proja/1000.);
        momVec.RotateX(projb/1000.);
        pt = momVec.Pt();
        Double_t e = TMath::Sqrt(mass*mass + mom*mom);
        TLorentzVector vec(momVec, e);
        rapL = vec.Rapidity();
        vec.Boost(betaVec);
        pzCM = vec.Pz();
        KECM = vec.E() - vec.M();
        rapCM = vec.Rapidity();
        rapCMNorm = vec.Rapidity()/(beamRapidity - rapidity_CM);
        phiCM = vec.Phi()*180./TMath::Pi();
        phiCM = (phiCM < 0 ? phiCM + 360 : phiCM);
        thetaCM = vec.Theta()*180./TMath::Pi();

        momVec = momVec.Unit();
        ndf = track -> GetClusterIDArray() -> size();
        dx = momVec.X();
        dy = momVec.Y();
        dz = momVec.Z();
        phiL = momVec.Phi()*180./TMath::Pi();
        phiL = (phiL < 0 ? phiL + 360 : phiL);
        thetaL = momVec.Theta()*180./TMath::Pi();

	TVector3 neuland = track->GetPosNeuland();
	TVector3 katana = track->GetPosKatana();
	
	NLx=-9999; NLy=-9999;NLz=-9999;
	histXZ->Fill(neuland.X(),neuland.Z());
	//if(neuland.Z()>1 && neuland.Y()>-1500 && neuland.Y()<1500){
	  //cout << neuland.X() << "," << neuland.Y() << "," << neuland.Z() << endl;
	
	NLx=neuland.X();//+TPC_S_x;
	NLy=neuland.Y();//+TPC_S_y;
	NLz=neuland.Z();//+TPC_S_z;
	
	//NLx=NLx*(.87)+NLz*.494;//rotate to NL frame
	//NLz=-NLx*.494+NLz*.87;//rotate to NL frame
	
	/*
	NLx=neuland.X()+2;//translate to SAMURAI frame
	  NLy=neuland.Y()+213.3;//translate to SAMURAI frame
	  NLz=neuland.Z()-580.4;//translate to SAMURAI frame
	  NLx=NLx-NLr*TMath::Sin(NLtheta*TMath::Pi()/180.);//translate to NL origin
	  NLz=NLz-NLr*TMath::Cos(NLtheta*TMath::Pi()/180.);//translate to NL origin
	*/
	//NLx=NLx*TMath::Cos(NLtheta*TMath::Pi()/180.)-NLz*TMath::Sin(NLtheta*TMath::Pi()/180.);//rotate to NL frame
	//NLz=NLx*TMath::Sin(NLtheta*TMath::Pi()/180.)+NLz*TMath::Cos(NLtheta*TMath::Pi()/180.);//rotate to NL frame

	//NLx=neuland.X()-4225.9;
	//NLy=neuland.Y()+213.3;
	//NLz=NLz-NLr;
	
	histXZ->Fill(NLx,NLz);	
	  

//}//if approx in right area for NeuLAND 
        //ea = track -> GetEffectiveArea();

//        proton = pFuncCut -> IsInside(mom/charge, dedx);
//        pim = pimCut -> IsInside(mom/charge, dedx);
//        pip = pipCut -> IsInside(mom/charge, dedx);
//        pipbg = pipBgCut -> IsInside(mom/charge, dedx);

  //      hist -> Fill(mom/charge, dedx);
        tree -> Fill();
      }

    }


    file -> cd();
    tree -> Write();
    tree->Delete();
  }


  
  histXZ->Draw("colz");
}
