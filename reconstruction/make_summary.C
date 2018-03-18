void make_summary(int iRun = 0)
{
  auto input_tree = new TChain("cbmsim");
  auto input_beam = new TChain("beam");
  /*
  TString tag = "trackChargeFix.1499.d50474b";  // before matching helix and reco track index
  TString tag = "trackChargeFix.1500.8a7d00d";  // before format version fix
  TString tag = "trackChargeFix.1501.7e5c17a";  // JungWoo - cut gg
  TString tag = "trackChargeFix.1502.9030344";  // PhysicsRun - use gg
  */

  //TString tag = "trackChargeFix.1503.4a82220";  // with parentID
  //TString tag = "trackChargeFix.1508.bf6c33b";
  //TString tag = "trackChargeFix.1504.9a9e4cc";  // cluster error on=
  TString tag = "develop.1522.ac6f82e";

  Int_t runs[] = {2900,2901,2902,2903,2904,2905};
  Int_t numEventsInSplit = 2000;
  Int_t numEventsInRun[] = {86916, 88750, 85788, 102451, 82746, 91334};
  Int_t numSplits[] = {44, 45, 43, 52, 42, 46};

  for (auto split = 0; split < numSplits[iRun]; ++split)
    input_tree -> Add(Form("data/run%d_s%d.reco.",runs[iRun],split)+tag+".root");

  input_beam -> Add(Form("/mnt/spirit/analysis/changj/BeamAnalysis/macros/refined/Sn132_all/run%d.refined.root",iRun));

  TClonesArray *vertexArray = nullptr;
  TClonesArray *recoTrackArray = nullptr;
  TClonesArray *vaddTrackArray = nullptr;

  input_tree -> SetBranchAddress("STVertex", &vertexArray);
  input_tree -> SetBranchAddress("STRecoTrack", &recoTrackArray);
  input_tree -> SetBranchAddress("VATracks", &vaddTrackArray);

  bool sigma20;
  input_beam -> SetBranchAddress("sigma20", &sigma20);

  auto fileName = TString("data/summary.")+tag+".merged2.root";
  TFile *output_file = new TFile(fileName,"recreate");
  TTree *output_tree = new TTree("data","");

  bool goodVertex;
  Int_t runID = run, splitID, eventID, entry, numTracks, vertexID, numdEdx;
  TVector3 posV;

  //reco
  Int_t charge1, pdg1, trackID1;
  Double_t distV1, distV21, dedx1, length1;
  TVector3 p1, pocav1;

  //va
  Int_t charge2, pdg2, trackID2;
  Double_t distV2, distV22, dedx2, length2;
  TVector3 p2, pocaV2, pocaV22;

  output_tree -> Branch("runid",      &runID);
  output_tree -> Branch("splitID",    &splitID);
  output_tree -> Branch("eventID",    &eventID);
  output_tree -> Branch("numdEdx",    &numdEdx);
  output_tree -> Branch("goodVertex", &goodVertex);
  output_tree -> Branch("posV",       "TVector3", &posV);

  output_tree -> Branch("charge",   &charge1);

  output_tree -> Branch("pdg1",     &pdg1);
  output_tree -> Branch("dedx1",    &dedx1);
  output_tree -> Branch("length1",  &length1);
  output_tree -> Branch("p1",       "TVector3", &p1);
  output_tree -> Branch("p1",       "TVector3", &p1);
  output_tree -> Branch("p1",       "TVector3", &p1);
  output_tree -> Branch("pocav1",   "TVector3", &pocav1);
  output_tree -> Branch("pocav1",   "TVector3", &pocav1);

  output_tree -> Branch("pdg2",     &pdg2);
  output_tree -> Branch("dedx2",    &dedx2);
  output_tree -> Branch("distV2",   &distV2);
  output_tree -> Branch("length2",  &length2);
  output_tree -> Branch("p2",       "TVector3", &p2);

  //output_tree -> Branch("entry",     &entry);
  //output_tree -> Branch("numTracks", &numTracks);
  //output_tree -> Branch("vertexID",  &vertexID);
  //output_tree -> Branch("trackID1", &trackID1);
  //output_tree -> Branch("trackID2", &trackID2);
  //output_tree -> Branch("distV1",   &distV1);
  //output_tree -> Branch("distV21",  &distV21);
  //output_tree -> Branch("distV22",  &distV22);
  //output_tree -> Branch("charge2",  &charge2);

  Int_t numEvents = input_tree -> GetEntries();
  for (entry = 0; entry < numEvents; ++entry)
  //for (entry = 0; entry < 100; ++entry)
  {
    if (entry % 1000 == 0)
      cout << "Event " << entry << " / " << numEvents << " -> " << recoTrackArray -> GetEntries() << " " << vaddTrackArray -> GetEntries() << endl;

    input_tree -> GetEntry(entry);

    auto idxRun = 0;
    for (;idxRun < 6; ++idxRun) {
      auto tempID = eventID - numEventsInRun[idxRun];
      if (tempID < 0)
        break;
      else
        eventID = tempID;
    }
    runID = runs[idxRun];
    splitID = eventID/numEventsInSplit;
    eventID = eventID - numEventsInSplit*splitID;

    if (vertexArray -> GetEntries() != 1) {
      continue;
    }

    auto vertex = (STVertex *) vertexArray -> At(0);
    auto posV = vertex -> GetPos();

    goodVertex = false;
    if (posV.Z()>-13.8&&posV.Z()<-10.3)
      goodVertex = true;

    vertexID = 0; //XXX
    //vertexID = ((STRecoTrack *) vaddTrackArray -> At(0)) -> GetVertexID();

    numTracks = vaddTrackArray -> GetEntries();

    //if (numTracks != vaddTrackArray -> GetEntries()) { continue; }

    for (trackID2 = 0; trackID2 < numTracks; ++trackID2)
    {
      auto track2 = (STRecoTrack *) vaddTrackArray -> At(trackID2);
      if (track2 -> GetVertexID() != vertexID)
        continue;

      trackID1 = track2->GetParentID();
      auto track1 = (STRecoTrack *) recoTrackArray -> At(trackID1);

      auto dedxArray0 = track1 -> GetdEdxPointArray();
      numdEdx = dedxArray0 -> size();
      if (numdEdx < 4) {
        continue;
      }

      pdg1 = STPID::GetPDG(track1 -> GetPID());
      pdg2  = STPID::GetPDG(track2  -> GetPID());

      pocav1 = track1 -> GetPOCAVertex() - posV;
      pocav2 = track2 -> GetPOCAVertex() - posV;
      dist1 = pocav1.Mag();
      dist2 = pocav2.Mag();

      auto dedxFirst = dedxArray0 -> at(0);
      auto dedxLast = dedxArray0 -> at(numdEdx-1);
      length1 = dedxLast.fLength - dedxFirst.fLength;
      length2 = dedxLast.fLength;

      charge1 = track1 -> GetCharge();
      charge2  = track2  -> GetCharge();

      p1 = track1 -> GetMomentum();
      p2 = track2 -> GetMomentum();
      dedx1 = track1 -> GetdEdxWithCut(0,0.7);
      dedx2 = track2 -> GetdEdxWithCut(0,0.7);

      output_tree -> Fill();
    }
  }

  output_file -> cd();
  output_tree -> Write();
  cout << ">> " << fileName << endl;
}
