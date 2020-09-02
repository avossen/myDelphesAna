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
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TF1.h"
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
enum RecType{elec, hadronic, da, mixed, mostlyLepton, truth,gen, genNoAcc,recTypeEnd};
enum bound{mean,upper,lower,endBound};
int main(int argc, char** argv)
{
  long downCount=0;
  long upCount=0;
  int colors[]={kRed,kBlue,kGreen,kBlack,kCyan,kMagenta,kOrange,kYellow};
  int markerStyles[]={20,21,22,23,43,33,34,47};
  string binningNames[]={"x","z","M"};
  cout <<"??" <<endl;
  char buffer[1000];

  string recTypeNames[]={"elec","hadronic","da","mixed","mostlyLepton","trueBoost","generatedWAcc","generatedWOAcc"};
  string boundNames[]={"mean","upper","lower"};
  if(argc<2)
    {
      cout <<"need to supply input file " <<endl;
      exit(0);
    }
  TFile mFile(argv[1]);
  TTree* mTrees[8];


  for(int i=0;i<recTypeEnd;i++)
    {
      sprintf(buffer,"tree_%s",recTypeNames[i].c_str());
      mTrees[i]=(TTree*)mFile.Get(buffer);
      if(mTrees[i]==0)
	cout <<" didn't load " << buffer <<endl;
    }
  
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
  double*** kinCounts=allocateArray<double>(8,4,10);

  //first dimention is bound (mean, upper, lower)
  double**** amps=allocateArray<double>(3,8,4,10);
  double**** ampErrs=allocateArray<double>(3,8,4,10);
  
  const int xBinning=0;
  const int zBinning=1;
  const int mBinning=2;
  
  vector<float> xBins;
  vector<float> zBins;
  vector<float> mBins;
  vector<float> phiBins;

  int numPhiBins=8;
  for(int i=0;i<numPhiBins;i++)
    {
      float bin=(i+1)*2*TMath::Pi()/numPhiBins;
      phiBins.push_back(bin);
    }
  
  xBins.push_back(0.15);
  xBins.push_back(0.2);
    xBins.push_back(0.25);
        xBins.push_back(0.3);
  xBins.push_back(0.35);
  xBins.push_back(0.45);
    xBins.push_back(0.6);

  xBins.push_back(50.0);

  zBins.push_back(0.3);
    zBins.push_back(0.35);
  zBins.push_back(0.4);
    zBins.push_back(0.45);
  zBins.push_back(0.5);
    zBins.push_back(0.6);
  zBins.push_back(100.0);
  
  mBins.push_back(0.3);
  mBins.push_back(0.5);
    mBins.push_back(0.6);
  mBins.push_back(0.7);
    mBins.push_back(0.8);
  mBins.push_back(0.9);
    mBins.push_back(1.2);
  mBins.push_back(200.0);


  
  for(int i=0;i<recTypeEnd;i++)
    {
      treeIndex=i;
      sprintf(buffer,"tree_%s",recTypeNames[i].c_str());
      //      cout <<"setting branching addresse" <<endl;
      setBranchAddresses(mTrees[treeIndex],treeFields[treeIndex]);
      //      cout <<"done " <<endl;


      long numEntries=mTrees[treeIndex]->GetEntries();
      //      cout <<"we have " << numEntries << " entries " <<endl;
      cout <<"getting " << numEntries <<" for recTypeName " << recTypeNames[i]<<endl;
      for(int ie=0;ie<numEntries;ie++)
	{
	  mTrees[treeIndex]->GetEntry(ie);


	  int xBin=getBin(xBins,treeFields[treeIndex].x);
	  if(treeFields[treeIndex].x<0.0)
	    {
	      continue;
	    }
	  if(treeFields[treeIndex].x>1.0)
	    {
	      continue;
	   }
	  



	  //some temporary event cuts
	  if(treeFields[treeIndex].x<0.05)
	    {
	      continue;
	    }
	  if(treeFields[treeIndex].y>0.1)
	      {
		continue;
	      }
	

	  //	  cout <<"x: "<< treeFields[treeIndex].x <<" bin: " << xBin<<endl;
	  //	  cout <<"looking at " << treeFields[treeIndex].numHadronPairs<<" pairs " <<endl;

	  for(int iPair=0;iPair<treeFields[treeIndex].numHadronPairs;iPair++)
	    {
	      int polIndex=0;
	      if(treeFields[treeIndex].polarization[iPair]>0)
		{
		  polIndex=1;
		}

	      diHadTreeFields* fields=&(treeFields[treeIndex]);
	      if(treeIndex==elec)
		{
		  //		  cout <<" elec true phiR: "<< fields->truePhiR[iPair] <<" rec: "<< fields->phiR[iPair] << " true phiS: "<< fields->truePhiS[iPair] <<" rec: " << fields->phiS[iPair] <<endl;
		}


	      if(isnan(treeFields[treeIndex].weight[iPair]))
		continue;
	      if(isnan(treeFields[treeIndex].weightUpperLimit[iPair]))
		continue;
	      //model not valid below 0.2
	      if(treeFields[treeIndex].z[iPair]<0.2)
		continue;
	      
	      if(treeFields[treeIndex].z[iPair]>1.0)
		continue;
	      if(treeFields[treeIndex].M[iPair]>3.0)
		continue;
	      if(treeFields[treeIndex].M[iPair]<0.1)
		continue;
	      
	      int zBin=getBin(zBins,treeFields[treeIndex].z[iPair]);
	      int mBin=getBin(mBins,treeFields[treeIndex].M[iPair]);
	      float ang=treeFields[treeIndex].phiR[iPair]+treeFields[treeIndex].phiS[iPair];
	      //	      float ang=treeFields[treeIndex].phiR[iPair];
	      //	      	      	      float ang=treeFields[treeIndex].phiS[iPair];
	      //	      cout <<"phiS: "<< ang <<endl;

	      //	      cout <<" z: " << treeFields[treeIndex].z[iPair] <<endl;
	      //	      cout <<" M: " << treeFields[treeIndex].M[iPair] <<endl;
	      
	      //	      cout <<"1" <<endl;
	      while(ang>2*TMath::Pi())
		{
		  ang-=(2*TMath::Pi());
		}
	      while(ang<0)
		{
		  ang+=(2*TMath::Pi());
		}
	      //	      cout <<" weight: "<< treeFields[treeIndex].weight[iPair] <<" upper: " << treeFields[treeIndex].weightUpperLimit[iPair] <<" lower: "<< treeFields[treeIndex].weightLowerLimit[iPair]<<endl;

	      if(isnan(ang))
		continue;
	      int phiBin=getBin(phiBins,ang);
	      
	      float weight=treeFields[treeIndex].weight[iPair];
	      
	      //	      weight=1.0;
	      //	      cout <<"2" <<endl;
	      kinMeans[i][xBinning][xBin]+=(weight*treeFields[treeIndex].x);
	      kinMeans[i][mBinning][mBin]+=(weight*treeFields[treeIndex].M[iPair]);
	      kinMeans[i][zBinning][zBin]+=(weight*treeFields[treeIndex].z[iPair]);

	      if(isnan(treeFields[treeIndex].x) || isnan(treeFields[treeIndex].M[iPair] || treeFields[treeIndex].z[iPair]))
		cout <<"some pair nan" <<endl;
		 
	      if(isnan(weight))
		cout <<"weight nan" <<endl;

	      if(isnan(treeFields[treeIndex].weightUpperLimit[iPair]))
		cout <<"weight nan upper" <<endl;

	      if(isnan(treeFields[treeIndex].weightLowerLimit[iPair]))
		cout <<"weight nan lower" <<endl;

	      kinCounts[i][xBinning][xBin]+=(weight);
	      kinCounts[i][mBinning][mBin]+=(weight);
	      kinCounts[i][zBinning][zBin]+=(weight);
	      

	      counts[i][xBinning][xBin][phiBin][polIndex]+=weight;

	      countsUpper[i][xBinning][xBin][phiBin][polIndex]+=treeFields[treeIndex].weightUpperLimit[iPair];

	      countsLower[i][xBinning][xBin][phiBin][polIndex]+=treeFields[treeIndex].weightLowerLimit[iPair];

	      counts[i][mBinning][mBin][phiBin][polIndex]+=weight;

	      countsUpper[i][mBinning][mBin][phiBin][polIndex]+=treeFields[treeIndex].weightUpperLimit[iPair];

	      countsLower[i][mBinning][mBin][phiBin][polIndex]+=treeFields[treeIndex].weightLowerLimit[iPair];

	      //	      cout <<"3" <<endl;
	      //	      cout <<" i: "<< i <<" binning: "<< zBinning <<" zbin: "<< zBin <<" phiBin: "<< phiBin <<" pol: "<< polIndex <<endl;
	      counts[i][zBinning][zBin][phiBin][polIndex]+=weight;
	      countsUpper[i][zBinning][zBin][phiBin][polIndex]+=treeFields[treeIndex].weightUpperLimit[iPair];
	      countsLower[i][zBinning][zBin][phiBin][polIndex]+=treeFields[treeIndex].weightLowerLimit[iPair];
	      //	      cout <<"3+" <<endl;
	    }
	}
    }

  //  	      cout <<"4" <<endl;
  for(int i=0;i<recTypeEnd;i++)
    {
      cout <<" looking at " << recTypeNames[i] <<endl;
      for(int binning=0;binning<3;binning++)
	{
	  int numKinBins=xBins.size();
	  for(int kinBin=0;kinBin<numKinBins;kinBin++)
	    {

	      double y[3][20];
	      double ey[3][20];
	      double x[3][20];
	      double ex[3][20];

	      double** locCounts[3];
	      locCounts[mean]=counts[i][binning][kinBin];
	      locCounts[upper]=countsUpper[i][binning][kinBin];
	      locCounts[lower]=countsLower[i][binning][kinBin];
	      
	      for(int iBound=mean;iBound<endBound;iBound++)
		{
		  for(int phiBin=0;phiBin<phiBins.size();phiBin++)
		    {
		      double Nup=locCounts[iBound][phiBin][1];
		      double Ndown=locCounts[iBound][phiBin][0];
		      double A=(Nup-Ndown)/(Nup+Ndown);
		      
		      cout <<"phiBin: "<< phiBin << " Nup: "<< Nup <<" Ndown: " << Ndown <<" A: " << A <<endl;
		      y[iBound][phiBin]=A;
		      x[iBound][phiBin]=(phiBin+0.5)*2*TMath::Pi()/numPhiBins;
		      ex[iBound][phiBin]=0.0;

		      double eU=sqrt(Nup);
		      double eD=sqrt(Ndown);
		      
		      double uDeriv=2*Ndown/((Nup+Ndown)*(Nup+Ndown));
		      double dDeriv=2*Nup/((Nup+Ndown)*(Nup+Ndown));
		      
		      ey[iBound][phiBin]=sqrt(uDeriv*uDeriv*eU*eU+dDeriv*dDeriv*eD*eD);
		  
		    }
		  TGraphErrors g(phiBins.size(),x[iBound],y[iBound],ex[iBound],ey[iBound]);
		  
		  gStyle->SetOptFit(111);
		  sprintf(buffer,"graphFor_rec_%s_binning%d_kinBin%d_bound_%s.png",recTypeNames[i].c_str(),binning,kinBin,boundNames[iBound].c_str());
		  TCanvas c1;
		  TF1 f1("f1","[0]*sin(x)",0,2*M_PI);
		  f1.SetParameters(0,0.0);
		  g.Fit(&f1);
		  g.Draw("AP");
		  c1.SaveAs(buffer);
		  amps[iBound][i][binning][kinBin]=f1.GetParameter(0);
		  ampErrs[iBound][i][binning][kinBin]=f1.GetParError(0);
		  cout << recTypeNames[i] <<", binning: "<< binning<< " kinBin: "<< kinBin <<endl;
		  cout <<"amp for bound " << iBound<<" " << amps[iBound][i][binning][kinBin] <<endl;


		  
		}

	    }
	}
    }
  int numKinBins=xBins.size();
  TH1D* ampDiffHistos[8];
  TH1D* ampDiffHistosTotal[8];
  for(int i=0;i<recTypeEnd;i++)
    {
      sprintf(buffer,"diffHisto_%s",recTypeNames[i].c_str());
      ampDiffHistos[i]=new TH1D(buffer,buffer,100,-3,3);
      sprintf(buffer,"diffHistoTotal_%s",recTypeNames[i].c_str());
      ampDiffHistosTotal[i]=new TH1D(buffer,buffer,100,-0.03,0.03);
    }
  
  for(int binning=0;binning<3;binning++)
    {
      
      TCanvas c;
      float maxX=0.8;

      
      for(int i=0;i<recTypeEnd;i++)
	{
	  double y[10];
	  double x[10];
	  double ey[10];
	  double ex[10];

	  //the rest for the bounds are the same
	  double boundEYUpper[10];
	  double boundEYLower[10];

	  for(int iKin=0;iKin<numKinBins;iKin++)
	    {
	      y[iKin]=amps[mean][i][binning][iKin];
	      boundEYUpper[iKin]=abs(amps[upper][i][binning][iKin]-y[iKin]);
	      boundEYLower[iKin]=abs(amps[lower][i][binning][iKin]-y[iKin]);
	      cout <<"ikin: " <<iKin <<" bound upper: "<< boundEYUpper[iKin] <<" lower: "<< boundEYLower[iKin]<<endl;
	      ey[iKin]=ampErrs[mean][i][binning][iKin];
	      x[iKin]=kinMeans[i][binning][iKin]/kinCounts[i][binning][iKin]+i*0.01;
	      cout <<"getting x for binning " << binning <<" : means: "<< kinMeans[i][binning][iKin] << " counts: "<< kinCounts[i][binning][iKin] <<endl;
	      ex[iKin]=0.0;
	      cout <<"iKin: "<< iKin << " y: "<< y[iKin] <<" x:" << x[iKin] << " ey: "<< ey[iKin]<<endl;


	      if(i!=gen)
		{
		  double diff=(y[iKin]-amps[mean][gen][binning][iKin])/sqrt(ey[iKin]*ey[iKin]+ampErrs[mean][i][binning][iKin]*ampErrs[mean][i][binning][iKin]);

		  ampDiffHistos[i]->Fill(diff);
		  //no normalization to uncertainty
		  diff=(y[iKin]-amps[mean][gen][binning][iKin]);
		  ampDiffHistosTotal[i]->Fill(diff);
		}
	      
	    }

	  sprintf(buffer,"ampGraphBounds_binning%s_recType_%s",binningNames[binning].c_str(),recTypeNames[i].c_str());
	  TGraphAsymmErrors* gBounds=new TGraphAsymmErrors(numKinBins,x,y,0,0,boundEYLower,boundEYUpper);
	  gBounds->SetName(buffer);
	  gBounds->SetLineWidth(4);
	  gBounds->SetLineColor(kGreen);

	  TGraphErrors* g=new TGraphErrors(numKinBins,x,y,ex,ey);
	  sprintf(buffer,"ampGraph_binning%s_recType_%s",binningNames[binning].c_str(),recTypeNames[i].c_str());
	  
	  g->SetName(buffer);
	  g->SetMarkerStyle(markerStyles[i]);
	  g->SetMarkerColor(colors[i]);

	  if(binning==0)
	    {
	      g->GetYaxis()->SetTitle("A");
	      g->GetXaxis()->SetTitle(binningNames[binning].c_str());
	      g->GetYaxis()->SetRangeUser(-0.3,0.3);
	      if(i==0)
		{
		gPad->DrawFrame(0.0,-0.3, 1.0, 0.3);
		g->GetXaxis()->SetLimits(0.0,0.8);
		}
	      
	      g->GetXaxis()->SetRangeUser(0.0,0.8);
	    }
	  else
	    {
	      if(i==0)
		{

		}
	      g->GetYaxis()->SetRangeUser(-0.3,0.3);
	    }

	  if(binning==1)
	    {
	      if(i==0)
		{
		  		  		maxX=0.85;
		  g->GetXaxis()->SetLimits(0.0,maxX);

		gPad->DrawFrame(0.0,-0.3, maxX, 0.3);

		}
	      g->GetXaxis()->SetRangeUser(0.0,0.85);
	      g->GetYaxis()->SetRangeUser(-0.1,0.15);
	    }
	  if(binning==2)
	    {
	      if(i==0)
		{
		  maxX=1.5;
		gPad->DrawFrame(0.0,-0.3, maxX, 0.3);
		g->GetXaxis()->SetLimits(0.0,maxX);
		}
	      g->GetXaxis()->SetRangeUser(0.0,maxX);
	      g->GetYaxis()->SetRangeUser(-0.1,0.2);
	    }

	  if(i==0)   
	    {
	      cout <<"drawing " <<endl;

	      g->Draw("AP");
	    }
	  //	  else
	  //	    {
	  //	      g->Draw("SAME P");
	  //	    }
	  if(i==gen)
	    {
	      gBounds->Draw("SAME E");
	    }
	      //redraw on top
	      {
		
		g->Draw("SAME P");
	      }
	    
	  //	        sprintf(buffer,"amps_binning_%s_recType_%s.png",binningNames[binning].c_str(),recTypeNames[i].c_str());
	}
      TLine *line = new TLine(0,0,maxX,0);
      line->Draw();
      sprintf(buffer,"amps_binning_%s.png",binningNames[binning].c_str());
      c.SaveAs(buffer);
    }

  TCanvas c;
  for(int i=0;i<recTypeEnd;i++)
    {
      sprintf(buffer,"ampDiffHisto_%s.png",recTypeNames[i].c_str());
      ampDiffHistos[i]->Draw();
      c.SaveAs(buffer);
      sprintf(buffer,"ampDiffHistoTotal_%s.png",recTypeNames[i].c_str());
      ampDiffHistosTotal[i]->Draw();
      c.SaveAs(buffer);

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
