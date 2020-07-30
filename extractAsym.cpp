#include "TLorentzVector.h"
#include "TF1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TLegend.h"
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <set>
using namespace std;
#include "AUTweight.h"
//#include "studyPIDs.h"

#include "TGraph.h"
#include "TFile.h"
#include "TStyle.h"
#include "TLine.h"
#include "TClonesArray.h"

#include <stdio.h>      /* printf, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"

#include "studyXQ2RecDelph.h"

#include "extractAsym.h"

int main(int argc, char** argv)
{
  int colors[]={kRed,kBlue,kGreen,kBlack,kCyan,kMagenta,kOrange,kYellow};
  int markerStyles[]={20,21,22,23,43,33,34,47};

  string recTypeNames[]={"elec","hadronic","da","mixed","mostlyLepton","trueBoost","generatedWAcc","generatedWOAcc"};
  if(argc<2)
    {
      cout <<"need to supply input file " <<endl;
      exit(0);
    }
  TFile* mFile(argv[0]);

  TTree* mTrees[8];
  diHadTreeFields treeFields[8];

  int treeIndex=0;
  //8 trees, 3 binnings (x,z,m), up to 10 kin bins, up to 16 azimuthal bins and 2 polarization bins
  //  double counts[8][4][10][16][2];
  //  double countsUpper[8][4][10][16][2];
  //  double countsLower[8][4][10][16][2];

  //this also sets them to zero and can be dynamic
  double***** counts=allocateArray<double>(8,4,10,16,2);
  double***** countsUpper=allocateArray<double>(8,4,10,16,2);
  double***** countsLower=allocateArray<double>(8,4,10,16,2);


  double*** kinMeans=allocateArray<double>(8,4,10);
  int*** kinCounts=allocateArray<double>(8,4,10);

  double*** amps=allocateArray(8,4,10);
  double*** ampErrs=allocateArray(8,4,10);
  
  const int xBinning=0;
  const int zBinning=1;
  const int mBinning=2;
  
  vector<float> xBins;
  vector<float> zBins;
  vector<float> mBins;
  vector<float> phiBins;

  int numPhiBins=16;
  for(int i=0;i<numPhiBins;i++)
    {
      float bin=(i+1)*2*TMath::Pi()/numPhiBins;
      phiBins.push_back(bin);
    }
  
  xBins.push_back(0.1);
  xBins.push_back(0.2);
  xBins.push_back(0.3);
  xBins.push_back(50.0);

  zBins.push_back(0.3);
  zBins.push_back(0.4);
  zBins.push_back(0.5);
  zBins.push_back(100.0);

  mBins.push_back(0.5);
  mBins.push_back(0.7);
  mBins.push_back(0.9);
  mBins.push_back(200.0);

  
  
  for(int i=0;i<recTypeEnd;i++)
    {
      sprintf(buffer,"tree_%s",recTypeNames[i].c_str());
      setBranchAddresses(mTrees[treeIndex],treeFields[treeIndex]);


      long numEntries=mTrees[treeIndex].GetEntries();
      for(int i=0;i<numEntries;i++)
	{
	  mTrees[treeIndex]->GetEntry(i);
	  //model not valid below 0.2
	  if(treeFields[treeIndex].z<0.2)
	    continue;

	  int xBin=getBin(xBins,treeFields[treeIndex].x);
	  int zBin=getBin(zBins,treeFields[treeIndex].z);
	  int mBin=getBin(mBins,treeFields[treeIndex].M);
	  float ang=treeFields[treeIndex].phiR+treeFields[treeIndex].phiS;
	  while(ang>2*TMath::Pi())
	    {
	      ang-=(2*TMath::Pi());
	    }
	  while(ang<0)
	    {
	      ang+=(2*TMath::Pi());
	    }
	  int phiBin=getBin(phiBins,ang);
	  int polIndex=0;
	  if(treeFields[treeIndex].polarization>0)
	    {
	      polIndex=1;
	    }

	  kinMeans[i][xBinning][xBin]+=(treeFields[treeIndex].weight*treeFields[treeIndex].x);
	  kinMeans[i][mBinning][xBin]+=(treeFields[treeIndex].weight*treeFields[treeIndex].M);
	  kinMeans[i][zBinning][xBin]+=(treeFields[treeIndex].weight*treeFields[treeIndex].z);

	  kinCounts[i][xBinning][xBin]+=(treeFields[treeIndex].weight);
	  kinCounts[i][mBinning][xBin]+=(treeFields[treeIndex].weight);
	  kinCounts[i][zBinning][xBin]+=(treeFields[treeIndex].weight);

	  
	  counts[i][xBinning][xBin][phiBin][polIndex]+=treeFields[treeIndex].weight;
	  countsUpper[i][xBinning][xBin][phiBin][polIndex]+=treeFields[treeIndex].weightUpperLimit;
	  countsLower[i][xBinning][xBin][phiBin][polIndex]+=treeFields[treeIndex].weightLowerLimit;
	  
	  counts[i][mBinning][mBin][phiBin][polIndex]+=treeFields[treeIndex].weight;
	  countsUpper[i][mBinning][mBin][phiBin][polIndex]+=treeFields[treeIndex].weightUpperLimit;
	  countsLower[i][mBinning][mBin][phiBin][polIndex]+=treeFields[treeIndex].weightLowerLimit;


	  counts[i][zBinning][zBin][phiBin][polIndex]+=treeFields[treeIndex].weight;
	  countsUpper[i][zBinning][zBin][phiBin][polIndex]+=treeFields[treeIndex].weightUpperLimit;
	  countsLower[i][zBinning][zBin][phiBin][polIndex]+=treeFields[treeIndex].weightLowerLimit;

	  
	}
    }
  
  for(int i=0;i<recTypeEnd;i++)
    {
      for(int binning=0;binning<3;binning++)
	{
	  int numKinBins=xBins.size();
	  for(int kinBin=0;kinBin<numKinBins;kinBin++)
	    {

	      double y[20];
	      double ey[20];
	      double x[20];
	      double ex[20];
	      
	      for(int phiBin=0;phiBin<phiBins.size();phiBins++)
		{
		  Nup=counts[i][binning][kinBin][phiBin][0];
		  Ndown=counts[i][binning][kinBin][phiBin][0];
		  double A=(Nup-Ndown)/(Nup+Ndown);
		  y[phiBin]=A;
		  x[phiBin]=(phiBin+1)*2*TMath::Pi()/numPhiBins;
		  ex[phiBin]=0.0;

		  float eU=sqrt(Nup);
		  float eD=sqrt(Ndown);

		  float uDeriv=2*Ndown/((Nup+Ndown)*(Nup+Ndown));
		  float dDeriv=2*Nup/((Nup+Ndown)*(Nup+Ndown));

		  ey[phiBin]=sqrt(uDeriv*uDeriv*eU*eU+dDeriv*dDeriv*eD*eD);
		  
		}

	      TGraphErrors g(phiBins.size(),x,y,ex,ey);
	      sprintf(buffer,"graphFor_rec_%s_binning%d_kinBin%d.png",recTypesNames[i],binning,kinBin);
	      TCanvas c1;
	      TF2 f1("f1","[0]*sin(x)",0,2*M_PI);
	      f1.SetParameters(0.0);
	      g.Fit(&f1);

	      g.Draw(AP);
	      c1.SaveAs(buffer);
	      amps[i][binning][kinBin]=f1.GetParameter(0);
	      ampErrs[i][binning][kinBin]=f1.GetParError(0);
	    }

	}

    }

  int numKinBins=xBins.size();
  for(int binning=0;binning<3;binning++)
    {
      for(int i=0;i<recTypeEnd;i++)
	{
	  double y[10];
	  double x[10];
	  double ey[10];
	  double ex[10];
	  for(int iKin=0;iKin<numKinBins;iKin++)
	    {
	      y[iKin]=amps[i][binning][iKin];
	      ey[iKin]=ampErrs[i][binning][iKin];
	      x[iKin]=kinMeans[i][binning][iKin]/kinCounts[i][binning][iKin];
	      ex[iKin]=0.0;
	    }
	  TGraphErrors g(numKinBins,x,y,ex,ey);
	  sprintf(buffer,"amp_binning%d_recType_%s");
	  g.SetName(buffer);
	  

	}

}



template<class T> T** allocateArray(int dim1, int dim2)
{
  T** ret=new T*[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=new T[dim2];
    }

  for(int i=0;i<dim1;i++)
    {
      for(int j=0;j<dim2;j++)
	{
	  ret[i][j]=0;
	}
    }
  return ret;
};

template<class T> T*** allocateArray(int dim1, int dim2, int dim3)
{
  T*** ret=new T**[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3);
    }
  return ret;
};

template<class T> T**** allocateArray(int dim1, int dim2, int dim3, int dim4)
{
  T**** ret=new T***[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4);
    }
  return ret;
};


template<class T> T***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5)
{
 T***** ret=new T****[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4,dim5);
    }
  return ret;
};


template<class T> T****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6)
{
 T****** ret=new T*****[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4,dim5,dim6);
    }
  return ret;
};

template<class T> T******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7)
{
 T******* ret=new T******[dim1];
  for(int i=0;i<dim1;i++)
    {
      ret[i]=allocateArray<T>(dim2,dim3,dim4,dim5,dim6,dim7);
    }
  return ret;
};

template int** allocateArray(int dim1, int dim2);
template int*** allocateArray(int dim1, int dim2, int dim3);
template int**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template int***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template int****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template int******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);

template unsigned long** allocateArray(int dim1, int dim2);
template unsigned long*** allocateArray(int dim1, int dim2, int dim3);
template unsigned long**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template unsigned long***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template unsigned long****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template unsigned long******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);


template float** allocateArray(int dim1, int dim2);
template float*** allocateArray(int dim1, int dim2, int dim3);
template float**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template float***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template float****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template float******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);

template double** allocateArray(int dim1, int dim2);
template double*** allocateArray(int dim1, int dim2, int dim3);
template double**** allocateArray(int dim1, int dim2, int dim3, int dim4);
template double***** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5);
template double****** allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6);
template double******* allocateArray(int dim1, int dim2, int dim3, int dim4, int dim5,int dim6, int dim7);
