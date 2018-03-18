#include "style.h"
using namespace style;

void new_clustering_test()
{
  gstat(0);

  gStyle -> SetPalette(kBrownCyan);

  auto file = new TFile("/mnt/spirit/analysis/user/leej/SpiRITROOT.trackChargeFix/macros/data/mc0_s0.reco.v1.04.root");
  auto tree = (TTree *) file -> Get("cbmsim");
  TClonesArray *hitArray = nullptr;
  TClonesArray *helixArray = nullptr;
  tree -> SetBranchAddress("STHit",&hitArray);
  tree -> SetBranchAddress("STHelixTrack",&helixArray);

  tree -> GetEntry(0);
  auto helix = (STHelixTrack *) helixArray -> At(0);
  auto ids = helix -> GetHitIDArray();
  helix -> FinalizeHits(); //XXX

  TObjArray hits;
  for (auto id : *ids) {
    auto hit = (STHit *) hitArray -> At(id);
    hits.Add(hit);
  }

  if (hits.GetEntriesFast() < 10)
    return;

  TVector3 q;
  Double_t alpha;

  bool isRow = 0;
  auto section = 0;
  bool currentBuildIsRow = 1; // should be different from isRow so it can init as son as loop start.
  auto currentBuildSection = -999; // should be different from section so it can init as son as loop start.

  TObjArray buildingClusters[2];
  Int_t buildingSection[2] = {-1,-1}; // 0: layer, 1: row
  Int_t numBuildingClusters[2] = {0}; // 0: layer, 1: row
  Int_t buildBoundary[2][2];
  buildBoundary[0][0] =  999;
  buildBoundary[0][1] = -999;
  buildBoundary[1][0] =  999;
  buildBoundary[1][1] = -999;

  TClonesArray *fClusterArray = new TClonesArray("STHitClusterRich",100);
  fClusterArray -> Clear();
  auto numClusters = 0;

  hits.Sort();
  TIter itHits(&hits);
  STHit *hit;

  vector<Int_t> sectionCollection;

  //auto hist = new TH2D("pp","",112,0,1344,108,-432,432); //XXX
  //auto hist = new TH2D("Pad Plane","Hit-Clustering; z (mm); x (mm)",112,0,1344,54,0,432); //XXX
  auto hist = new TH2D("Pad Plane","Hit-Clustering; z (mm); x (mm)",56,0,672,54,0,432); //XXX
  while ((hit = (STHit *) itHits.Next()))
  {
    helix -> ExtrapolateToPointAlpha(hit -> GetPosition(), q, alpha);

    section = (Int_t)(alpha/(TMath::Pi()/4));
    isRow = abs(section)%2;

    if (section != currentBuildSection) // init
    {
      for (auto section0 : sectionCollection)
        if (section0 == section)
          continue;

      sectionCollection.push_back(section);

      cout << "section: " << currentBuildSection << " -> " << section << endl;

      currentBuildIsRow = isRow;
      currentBuildSection = section;
      vector<Int_t> initList; // 1 for row, 0 for layer

      if (currentBuildSection-section > 1) // if jumped section
      {
        if (buildingSection[1] > buildingSection[0]) {
          initList.push_back(0);
          initList.push_back(1);
        } else {
          initList.push_back(1);
          initList.push_back(0);
        }
      }
      else if (isRow) initList.push_back(1);
      else            initList.push_back(0);

      for (auto rl : initList) {
        buildingClusters[rl].Clear();
        numBuildingClusters[rl] = 0;
        buildBoundary[abs(rl-1)][0] = 999;
        buildBoundary[abs(rl-1)][1] = -999;
      }

      //TODO
      continue;
    }
    else // continue build with same option
    {
      if (isRow) {
        auto row = hit -> GetRow();
        auto layer = hit -> GetLayer();

        //shoud check opposite boundary because next build will use shade area
        if (layer < buildBoundary[0][0]) buildBoundary[0][0] = layer;
        if (layer > buildBoundary[0][1]) buildBoundary[0][1] = layer;

        // check before build
        bool foundCluster = false;
        for (auto iCluster = numBuildingClusters[0]-1; iCluster >= 0; --iCluster) {
          auto cluster = (STHitClusterRich *) buildingClusters[0].At(iCluster);
          if (cluster -> GetLayer() == layer) {
            cluster -> SetClusterID(-1);
            foundCluster = true;
            break;
          }
        }
        if(foundCluster)
          continue;

        // check this build
        if (row >= buildBoundary[1][0] && row <= buildBoundary[1][1])
          continue;

        foundCluster = false;
        for (auto iCluster = numBuildingClusters[1]-1; iCluster >= 0; --iCluster) {
          auto cluster = (STHitClusterRich *) buildingClusters[1].At(iCluster);
          if (cluster -> GetRow() == row) {
            foundCluster = true;
            cluster -> AddHit(hit);
            break;
          }
        }
        if (!foundCluster) {
          auto cluster = (STHitClusterRich *) fClusterArray -> ConstructedAt(numClusters++);
          cluster -> SetClusterID(1);
          cluster -> AddHit(hit);
          cluster -> SetRow(row);
          cluster -> SetLayer(-1);
          buildingClusters[1].Add(cluster);
          ++numBuildingClusters[1];
        }
      } ////////////////////////////////////////////////////// was for row build
      else {
        auto row = hit -> GetRow();
        auto layer = hit -> GetLayer();

        //shoud check opposite boundary because next build will use shade area
        if (row < buildBoundary[1][0]) buildBoundary[1][0] = row;
        if (row > buildBoundary[1][1]) buildBoundary[1][1] = row;

        // check before build
        bool foundCluster = false;
        for (auto iCluster = numBuildingClusters[1]-1; iCluster >= 0; --iCluster) {
          auto cluster = (STHitClusterRich *) buildingClusters[1].At(iCluster);
          if (cluster -> GetRow() == row) {
            cluster -> SetClusterID(-1);
            foundCluster = true;
            break;
          }
        }
        if(foundCluster)
          continue;

        // check this build
        if (layer >= buildBoundary[0][0] && layer <= buildBoundary[0][1])
          continue;

        foundCluster = false;
        for (auto iCluster = numBuildingClusters[0]-1; iCluster >= 0; --iCluster) {
          auto cluster = (STHitClusterRich *) buildingClusters[0].At(iCluster);
          if (cluster -> GetLayer() == layer) {
            foundCluster = true;
            cluster -> AddHit(hit);
            if (cluster -> GetZ() == hit -> GetZ())
              cout << cluster -> GetZ() << " " << hit -> GetZ() << endl;
            break;
          }
        }
        if (!foundCluster) {
          auto cluster = (STHitClusterRich *) fClusterArray -> ConstructedAt(numClusters++);
          cluster -> SetClusterID(1);
          cluster -> AddHit(hit);
          cluster -> SetRow(-1);
          cluster -> SetLayer(layer);
          buildingClusters[0].Add(cluster);
          ++numBuildingClusters[0];
        }
      }
    } ////////////////////////////////////////////////////// was for layer build

    hist -> Fill(hit -> GetZ(), hit -> GetX(), hit -> GetCharge());
    //if (isRow) hist -> Fill(hit -> GetZ(), hit -> GetX(),2); //XXX
    //else       hist -> Fill(hit -> GetZ(), hit -> GetX()); //XXX
  }

  auto gr = new TGraphAsymmErrors();
  gr -> SetMarkerStyle(20);
  gr -> SetMarkerSize(0.4);
  auto gl = new TGraphAsymmErrors();
  gl -> SetMarkerStyle(20);
  gl -> SetMarkerSize(0.4);
  for (auto iCluster = 0; iCluster < numClusters; ++iCluster) {
    auto cluster = (STHitClusterRich *) fClusterArray -> At(iCluster);
    if (cluster -> GetClusterID() == -1) {
      fClusterArray -> Remove(cluster);
    }
    if (cluster -> IsRowCluster()) {
      auto z = cluster -> GetZ();
      auto max = cluster -> GetZMax();
      auto min = cluster -> GetZMin();
      auto zh = max - z + 6;
      auto zl = z - min + 6;
      cout << "z: " << zl << " " << zh << endl;
      gr -> SetPoint(gr->GetN(),cluster->GetZ(),cluster->GetX());
      gr -> SetPointError(gr->GetN()-1,zl,zh,0.,0.);
    }
    else if (cluster -> IsLayerCluster()) {
      auto x = cluster -> GetX();
      auto max = cluster -> GetXMax();
      auto min = cluster -> GetXMin();
      auto xh = max - x + 4;
      auto xl = x - min + 4;
      cout << "x: " << xl << " " << xh << endl;
      gl -> SetPoint(gl->GetN(),cluster->GetZ(),cluster->GetX());
      gl -> SetPointError(gl->GetN()-1,0.,0.,xl,xh);
    }
    else
      cout << "asd;lfkjasd;fkjawe;oivjscnlvkjahwef;oaiwjef awef" << endl;
  }
  fClusterArray -> Compress();

  auto graph = new TGraph();
  for (auto i = 0.; i < 1; ++i)
    helix ->  InterpolateByRatio(Double_t r) const;

  auto cvs = cc();//new TCanvas("cvs","",1150,700); //XXX
  cvs -> SetLogz();
  make(hist) -> Draw("colz"); //XXX
  gr -> Draw("psame");
  gl -> Draw("psame");

  save(cvs);
}
