#include "style.h"
#include "../reconstruction/config_draw.h"
using namespace style;

void make_hitcluster_summary()
{
  auto input_tree = new TChain("cbmsim");
  //input_tree -> Add("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s200.reco.trackChargeFix.1505.c48e1f8.test.root");
  input_tree -> Add("/mnt/spirit/analysis/user/leej/macros/hitcluster/data/run2900_s200.reco.develop.1521.dce79d3.test.root");

  TClonesArray *clusterArray = nullptr;
  TClonesArray *helixArray = nullptr;
  input_tree -> SetBranchAddress("STHitCluster", &clusterArray);
  input_tree -> SetBranchAddress("STHelixTrack", &helixArray);

  auto fileName = TString("data/summary.root");
  TFile *output_file = new TFile(fileName,"recreate");
  TTree *output_tree = new TTree("data","");

  TVector3 pos, cov;
  Double_t charge;
  Int_t n, clusterID;
  Long64_t eventID;
  bool layerCluster;
  //bool lastCluster;

  output_tree -> Branch("eventID",&eventID);
  output_tree -> Branch("clusterID",&clusterID);
  output_tree -> Branch("pos","TVector3",&pos);
  output_tree -> Branch("cov","TVector3",&cov);
  output_tree -> Branch("layerCluster",&layerCluster);
  //output_tree -> Branch("lastCluster",&lastCluster);
  output_tree -> Branch("charge",&charge);
  output_tree -> Branch("n",&n);

  auto numEvents = input_tree -> GetEntries();
  for (eventID = 0; eventID < numEvents; ++eventID)
  //for (auto e: {3,4,7})
  {
    //eventID = e;
    input_tree -> GetEntry(eventID);
    cout << eventID << " / " << numEvents << endl;

    auto numHelix = helixArray -> GetEntries();
    for (auto iHelix = 0; iHelix < numHelix; ++iHelix)
    {
      auto helix = (STHelixTrack *) helixArray -> At(iHelix);
      auto ids = helix -> GetClusterIDArray();
      for (auto id : *ids) {
        clusterID = id;
        auto cluster = (STHitCluster *) clusterArray -> At(id);

        layerCluster = true;
        if (cluster -> GetRow() != -1)
          layerCluster = false;


        pos = cluster -> GetPosition();
        auto mat = cluster -> GetCovMatrix();
        cov.SetXYZ(mat(0,0),mat(1,1),mat(2,2));
        charge = cluster -> GetCharge();
        n = cluster -> GetNumHits();

        output_tree -> Fill();
      }
      //lastCluster = true;
      //output_tree -> Fill();
      //lastCluster = false;
    }
  }

  output_file -> cd();
  output_tree -> Write();
  cout << ">> " << fileName << endl;

  gApplication -> Terminate();
}
