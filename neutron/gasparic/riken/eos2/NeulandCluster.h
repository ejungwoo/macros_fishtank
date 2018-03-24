#include <algorithm>  // for sort 

static const Int_t nMaxSize = 400;
static const Int_t tTolerance = 10.0; // [ns]
static const Int_t dTolerance = 10.0; // [cm]
static const Int_t cc = 29.979;       // [cm/ns]

static const Double_t tReso = 0.3;   // timing resolution [ns]
static const Double_t dReso = 10.0;  // [cm]

class NeulandHit{
 private:
  Double_t t, x, y, z;
 public:
  NeulandHit();
  Double_t getT() const {return t;}
  Double_t getX(){return x;}
  Double_t getY(){return y;}
  Double_t getZ(){return z;}
  void setT(Double_t tin){t = tin;}
  void setX(Double_t xin){x = xin;}
  void setY(Double_t yin){y = yin;}
  void setZ(Double_t zin){z = zin;}
};

NeulandHit::NeulandHit(){
  t = 1.e20;
  x = 1.e20;
  y = 1.e20;
  z = 1.e20;
}

bool operator<(const NeulandHit &nh1, const NeulandHit &nh2){  
  return nh1.getT() < nh2.getT();
}



class NeulandCluster{

private:
  Int_t nHit;
  NeulandHit nharray[nMaxSize];
  
public:
  NeulandCluster();
  void AddHit(NeulandHit nh);
  void AddCluster(NeulandCluster cluster);
  bool IsInCluster(NeulandHit hit);
  void TimeSort();
  Int_t GetnHit() const {return nHit;};
  void SetnHit(Int_t n){nHit = n;};
  NeulandHit GetHit(Int_t iHit){return nharray[iHit];};
};

NeulandCluster::NeulandCluster(){
  nHit = 0;
}

void NeulandCluster::AddHit(NeulandHit nh){
  nharray[nHit].setT(nh.getT());
  nharray[nHit].setX(nh.getX());
  nharray[nHit].setY(nh.getY());
  nharray[nHit].setZ(nh.getZ());
  nHit++;
}

void NeulandCluster::AddCluster(NeulandCluster cluster){
  //  printf("AddCluster: before nhit %d\n",nHit);
  //  printf("AddCluster: add %d hit\n", cluster.GetnHit());
  for(Int_t iHit=0; iHit<cluster.GetnHit(); iHit++){
    this->AddHit(cluster.GetHit(iHit));
  }
  //  printf("AddCluster: after nhit %d\n",nHit);
}

bool NeulandCluster::IsInCluster(NeulandHit hit){
  for(Int_t iHit=0; iHit<nHit; iHit++){
    Double_t dist;
    dist =  pow((hit.getX()-nharray[iHit].getX()),2.0)
          + pow((hit.getY()-nharray[iHit].getY()),2.0)
          + pow((hit.getZ()-nharray[iHit].getZ()),2.0);
    dist = sqrt(dist);
    if(dist<dTolerance && abs(hit.getT()-nharray[iHit].getT())<tTolerance){
      return true;
    }
  }
  return false;
}

void NeulandCluster::TimeSort(){
  sort(nharray, nharray + nHit);
}

class GreaterNeulandCluster{
 public:
  bool operator()(const NeulandCluster &nc1, const NeulandCluster &nc2) const {
    return nc1.GetnHit() > nc2.GetnHit();
  }
};


bool causality(NeulandCluster nc1, NeulandCluster nc2){
  Double_t dist;
  dist = sqrt( pow((nc1.GetHit(0).getX()-nc2.GetHit(0).getX()),2.0)
              +pow((nc1.GetHit(0).getY()-nc2.GetHit(0).getY()),2.0)
              +pow((nc1.GetHit(0).getZ()-nc2.GetHit(0).getZ()),2.0));
  return dist-dReso > cc * (abs(nc1.GetHit(0).getT()-nc2.GetHit(0).getT())+tReso*sqrt(2.0));
}
