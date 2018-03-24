//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jan  9 17:10:20 2018 by ROOT version 6.08/00
// from TTree dedx/
// found on file: /mnt/spirit/analysis/changj/SpiRITROOT.JWfix/macros/data/analysisCode-All-noLayerCut-GC-DS-GiordanoCommentOut-bShift-By-JWfix/dedxROOT/dedxSn132-All-noLayerCut_0.root
//////////////////////////////////////////////////////////

#ifndef tpc_h
#define tpc_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class tpc {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
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
   Double_t        dist;
   Bool_t          sigma10;
   Bool_t          sigma15;
   Bool_t          sigma20;
   Bool_t          sigma20z;
   Int_t           mult;
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

   // List of branches
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
   TBranch        *b_dist;   //!
   TBranch        *b_sigma10;   //!
   TBranch        *b_sigma15;   //!
   TBranch        *b_sigma20;   //!
   TBranch        *b_sigma20z;   //!
   TBranch        *b_mult;   //!
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

   tpc(TTree *tree=0);
   virtual ~tpc();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef tpc_cxx
tpc::tpc(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/mnt/spirit/analysis/changj/SpiRITROOT.JWfix/macros/data/analysisCode-All-noLayerCut-GC-DS-GiordanoCommentOut-bShift-By-JWfix/dedxROOT/dedxSn132-All-noLayerCut_0.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/mnt/spirit/analysis/changj/SpiRITROOT.JWfix/macros/data/analysisCode-All-noLayerCut-GC-DS-GiordanoCommentOut-bShift-By-JWfix/dedxROOT/dedxSn132-All-noLayerCut_0.root");
      }
      f->GetObject("dedx",tree);

   }
   Init(tree);
}

tpc::~tpc()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t tpc::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t tpc::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void tpc::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("eventid", &eventid, &b_eventid);
   fChain->SetBranchAddress("dedx", &dedx, &b_dedx);
   fChain->SetBranchAddress("mom", &mom, &b_mom);
   fChain->SetBranchAddress("charge", &charge, &b_charge);
   fChain->SetBranchAddress("pid", &pid, &b_pid);
   fChain->SetBranchAddress("stpid", &stpid, &b_stpid);
   fChain->SetBranchAddress("ndf", &ndf, &b_ndf);
   fChain->SetBranchAddress("dx", &dx, &b_dx);
   fChain->SetBranchAddress("dy", &dy, &b_dy);
   fChain->SetBranchAddress("dz", &dz, &b_dz);
   fChain->SetBranchAddress("vid", &vid, &b_vid);
   fChain->SetBranchAddress("parentvid", &parentvid, &b_parentvid);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("pocavx", &pocavx, &b_pocavx);
   fChain->SetBranchAddress("pocavy", &pocavy, &b_pocavy);
   fChain->SetBranchAddress("pocavz", &pocavz, &b_pocavz);
   fChain->SetBranchAddress("dist", &dist, &b_dist);
   fChain->SetBranchAddress("sigma10", &sigma10, &b_sigma10);
   fChain->SetBranchAddress("sigma15", &sigma15, &b_sigma15);
   fChain->SetBranchAddress("sigma20", &sigma20, &b_sigma20);
   fChain->SetBranchAddress("sigma20z", &sigma20z, &b_sigma20z);
   fChain->SetBranchAddress("mult", &mult, &b_mult);
   fChain->SetBranchAddress("pt", &pt, &b_pt);
   fChain->SetBranchAddress("pzCM", &pzCM, &b_pzCM);
   fChain->SetBranchAddress("KECM", &KECM, &b_KECM);
   fChain->SetBranchAddress("rapL", &rapL, &b_rapL);
   fChain->SetBranchAddress("rapCM", &rapCM, &b_rapCM);
   fChain->SetBranchAddress("rapCMNorm", &rapCMNorm, &b_rapCMNorm);
   fChain->SetBranchAddress("projx", &projx, &b_projx);
   fChain->SetBranchAddress("projy", &projy, &b_projy);
   fChain->SetBranchAddress("projz", &projz, &b_projz);
   fChain->SetBranchAddress("phiL", &phiL, &b_phiL);
   fChain->SetBranchAddress("thetaL", &thetaL, &b_thetaL);
   fChain->SetBranchAddress("phiCM", &phiCM, &b_phiCM);
   fChain->SetBranchAddress("thetaCM", &thetaCM, &b_thetaCM);
   fChain->SetBranchAddress("proton", &proton, &b_proton);
   fChain->SetBranchAddress("pip", &pip, &b_pip);
   fChain->SetBranchAddress("pim", &pim, &b_pim);
   fChain->SetBranchAddress("pipbg", &pipbg, &b_pipbg);
   fChain->SetBranchAddress("ea", &ea, &b_ea);
   Notify();
}

Bool_t tpc::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void tpc::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t tpc::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef tpc_cxx
