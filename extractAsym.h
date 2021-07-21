#ifndef EXTRACT_ASYM_H_
#define EXTRACT_ASYM_H_
#include "TTree.h"


void setBranchAddresses(TTree* mTree,diHadTreeFields& fields)
{
  mTree->SetBranchAddress("Q2",&fields.Q2);
  mTree->SetBranchAddress("trueQ2",&fields.trueQ2);
  mTree->SetBranchAddress("x",&fields.x);
  mTree->SetBranchAddress("trueX",&fields.trueX);
  mTree->SetBranchAddress("y",&fields.y);
  mTree->SetBranchAddress("W",&fields.W);
  mTree->SetBranchAddress("Mx",&fields.Mx);
  mTree->SetBranchAddress("evtNr",&fields.evtNr);
  mTree->SetBranchAddress("polarization",&fields.polarization);
  mTree->SetBranchAddress("numHadronPairs",&fields.numHadronPairs);
  mTree->SetBranchAddress("phiR",fields.phiR);
  mTree->SetBranchAddress("phiS",fields.phiS);
  mTree->SetBranchAddress("truePhiR",fields.truePhiR);
  mTree->SetBranchAddress("truePhiS",fields.truePhiS);
  mTree->SetBranchAddress("trueM",fields.trueM);
  mTree->SetBranchAddress("trueZ",fields.trueZ);

  mTree->SetBranchAddress("phiH",fields.phiH);
  mTree->SetBranchAddress("z",fields.z);
  mTree->SetBranchAddress("M",fields.M);
  mTree->SetBranchAddress("theta",fields.theta);
  mTree->SetBranchAddress("xF",fields.xF);
  mTree->SetBranchAddress("pT",fields.pT);
  mTree->SetBranchAddress("weight",fields.weight);
  mTree->SetBranchAddress("weightUpperLimit",fields.weightUpperLimit);
  mTree->SetBranchAddress("weightLowerLimit",fields.weightLowerLimit);
  mTree->SetBranchAddress("rawWeight",fields.rawWeight);
  mTree->SetBranchAddress("rawWeightUnc",fields.rawWeightUnc);
  mTree->SetBranchAddress("pairType",fields.pairType);
}


pair<double,double> getA(double** vals, vector<float>& phiBins,int kinBin, const char* recTypeName,const char* boundName, int binning,char* graphName)
{
  double y[20];
  double ey[20];
  double x[20];
  double ex[20];

  char buffer[200];
  
  int numPhiBins=phiBins.size();
  cout <<"numphi bins: " << numPhiBins <<endl;
  for(int phiBin=0;phiBin<phiBins.size();phiBin++)
    {
      cout <<" getting nup for bin " << phiBin <<endl;
      double Nup=vals[phiBin][1];
      double Ndown=vals[phiBin][0];
      cout <<"nup+ndown: "<< Nup+Ndown <<endl;
      double A=(Nup-Ndown)/(Nup+Ndown);
      
      cout <<"phiBin: "<< phiBin << " Nup: "<< Nup <<" Ndown: " << Ndown <<" A: " << A <<endl;
      y[phiBin]=A;
      x[phiBin]=(phiBin+0.5)*2*TMath::Pi()/numPhiBins;
      ex[phiBin]=0.0;
      
      double eU=sqrt(Nup);
      double eD=sqrt(Ndown);
      
      double uDeriv=2*Ndown/((Nup+Ndown)*(Nup+Ndown));
      double dDeriv=2*Nup/((Nup+Ndown)*(Nup+Ndown));
      
      ey[phiBin]=sqrt(uDeriv*uDeriv*eU*eU+dDeriv*dDeriv*eD*eD);
      
    }
  TGraphErrors g(phiBins.size(),x,y,ex,ey);
  
  gStyle->SetOptFit(111);
  cout <<"saving as " << buffer <<endl;
  sprintf(buffer,"graphFor_rec_%s_binning%d_kinBin%d_bound_%s_%s.png",recTypeName,binning,kinBin,boundName,graphName);
  //sprintf(buffer,"graphFor_rec_%s_binning%d_kinBin%d_bound_%s.png",recTypeNames[i].c_str(),binning,kinBin,boundNames[iBound].c_str());
  TCanvas c1;
  TF1 f1("f1","[0]*sin(x)",0,2*M_PI);
  f1.SetParameters(0,0.0);
  g.Fit(&f1);
  g.Draw("AP");
    c1.SaveAs(buffer);
  double A=f1.GetParameter(0);
  double err=f1.GetParError(0);
  return pair<double,double>(A,err);
}


template<class T> T** allocateArray(int dim1, int dim2);
template<class T> T*** allocateArray(int dim1, int dim2, int dim3);
template<class T> T**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template<class T> T***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template<class T> T****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template<class T> T******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);

#endif
