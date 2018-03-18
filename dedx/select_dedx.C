void select_dedx()
{
  TString cut = "ndf>50&&dist<5&&dedx>50&&dedx<110&&mom>500&&mom<500.1&&charge>0";

  auto tree = new TChain("dedx");
  //for (auto i = 0; i < 56; ++i)
  for (auto i = 0; i < 56; ++i)
    tree -> AddFile(Form("dedxROOT/dedxSn132-noLayerCut_%d.root",i));

  auto cvs = new TCanvas("cvs","",700,500);
  tree -> Draw("dedx:mom>>hist(400,-800,3000,400,0,1000)",cut,"colz");
  cvs -> SetLogz();

  auto file1 = new TFile("data/selected.root","recreate");
  auto tree1 = (TTree *) tree -> CopyTree(cut);
  tree1 -> Write();
}
