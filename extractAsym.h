#ifndef EXTRACT_ASYM_H_
#define EXTRACT_ASYM_H_
#include "TTree.h"


void setBranchAddresses(TTree* mTree,diHadTreeFields& fields)
{
  mTree->SetBranchAddress("Q2",&fields.Q2);
  mTree->SetBranchAddress("x",&fields.x);
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

  mTree->SetBranchAddress("phiH",fields.phiH);
  mTree->SetBranchAddress("z",fields.z);
  mTree->SetBranchAddress("M",fields.M);
  mTree->SetBranchAddress("theta",fields.theta);
  mTree->SetBranchAddress("xF",fields.xF);
  mTree->SetBranchAddress("pT",fields.pT);
  mTree->SetBranchAddress("weight",fields.weight);
  mTree->SetBranchAddress("weightUpperLimit",fields.weightUpperLimit);
  mTree->SetBranchAddress("weightLowerLimit",fields.weightLowerLimit);
  mTree->SetBranchAddress("pairType",fields.pairType);
}


template<class T> T** allocateArray(int dim1, int dim2);
template<class T> T*** allocateArray(int dim1, int dim2, int dim3);
template<class T> T**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template<class T> T***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template<class T> T****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template<class T> T******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);

#endif
