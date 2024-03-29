#include "TText.h"
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
#include "AUTweight.h"


//enum RecType{elec, hadronic, da, mixed, mostlyLepton, truth,gen, genNoAcc,recTypeEnd};
enum RecType{elec, gen,truth, recTypeEnd, hadronic, da, mixed, mostlyLepton, genNoAcc};
enum bound{mean,upper,lower,endBound};
int main(int argc, char** argv)
{
  //polarization

  bool debug=false;
  double minZ=0.2;
  double minM=0.1;
  ofstream outFile;
  outFile.open ("extractionOutput.txt");
  bool useRealAsym=false;
  //plot the inner pad of the projections vs m (otherwise z)
  bool vsM=true;
  int minCounts=10;
  gStyle->SetOptStat(0);
  AUTweight*  m_weights;
  m_weights=new AUTweight();
  
  ////----for the fancy plotting
  //  double minXCut=0.001;
  //following Marco's mesage
  
  double minXCut=0.0001;
  double minQ2Cut=1.0;
  double minQ2_fancy=0.1;
  double maxQ2_fancy=10000;
  double logMinQ2=log10(minQ2_fancy);
  double logMaxQ2=log10(maxQ2_fancy);
  //should be smaller then the actual cuts, so one can still see the axis
  double minX_fancy=0.0001;
  double maxX_fancy=1.0;
  double logMinX=log10(minX_fancy);
  double logMaxX=log10(maxX_fancy);
  double maxQ2Fancy=1000000;
  double maxXFancy=1.0;
  ////-----
  long downCount=0;
  long upCount=0;
  int colors[]={kRed,kBlue,kGreen,kBlack,kCyan,kMagenta,kOrange,kYellow};
  int markerStyles[]={20,21,22,23,43,33,34,47};
  string binningNames[]={"Q2","x","z","M"};
  cout <<"??" <<endl;
  char buffer[1000];
  char buffer1[1000];
string recTypeNames[]={"elec","generatedWAcc","truth","hadronic","da","mixed","mostlyLepton","generatedWOAcc"};
//  string recTypeNames[]={"elec", "generatedWAcc","hadronic","da","mixed","mostlyLepton","generatedWOAcc"};
  string boundNames[]={"mean","upper","lower"};
  if(argc<2)
    {
      cout <<"need to supply input file " <<endl;
      exit(0);
    }



  char* listfile=argv[1];
  ifstream fListFile(listfile);
  string line;


  diHadTreeFields treeFields[8];

  int treeIndex=0;
  //8 trees, 3 binnings (x,z,m), up to 10 kin bins, up to 16 azimuthal bins and 2 polarization bins
  //  double counts[8][4][10][16][2];
  //  double countsUpper[8][4][10][16][2];
  //  double countsLower[8][4][10][16][2];
  vector<float> Q2Bins;
  vector<float> xBins;
  vector<float> zBins;
  vector<float> mBins;
  vector<float> phiBins;


  double xBinsPerDecade=3;

  double q2BinsPerDecade=2;
	 // double q2BinsPerDecade=1;
    for(int i=1;i<9;i++)
  //  for(int i=1;i<2;i++)
    {
      cout <<" adding q2Bin: "<< pow(10,i/q2BinsPerDecade)<<endl;
      Q2Bins.push_back(pow(10,+i/q2BinsPerDecade));
      }
	 //    Q2Bins.push_back(1000);
	 
	 //	 Q2Bins.push_back(10.0);
	 Q2Bins.push_back(100000);

  
	 //   Q2Bins.push_back(100);
	 //   Q2Bins.push_back(100000);


    for(int i=1;i<12;i++)
//	//       for(int i=1;i<3;i++)
     {
       cout <<" adding xbin: "<< pow(10,-4+i/xBinsPerDecade) <<endl;
       xBins.push_back(pow(10,-4.0+i/xBinsPerDecade));
     }
   


    //        xBins.push_back(0.1);
    //	   xBins.push_back(0.3);
    //	   xBins.push_back(0.4);
       xBins.push_back(100.0);
//
//just for the plots
////  xBins.push_back(0.2);
////  xBins.push_back(0.3);
////  xBins.push_back(0.4);
////  xBins.push_back(0.5);
////  xBins.push_back(100.00);
////  
  
    zBins.push_back(0.3);
  //            zBins.push_back(0.35);
    zBins.push_back(0.4);
	//                  zBins.push_back(0.45);
        zBins.push_back(0.5);
        zBins.push_back(0.6);
        zBins.push_back(0.7);
        zBins.push_back(0.8);
        zBins.push_back(0.9);
//       zBins.push_back(0.4);
       zBins.push_back(100.0);
       
    	mBins.push_back(0.3);
	mBins.push_back(0.5);
	mBins.push_back(0.6);
	mBins.push_back(0.7);
	mBins.push_back(0.8);
	mBins.push_back(0.9);
	mBins.push_back(1.2);
	mBins.push_back(1.5);
  mBins.push_back(2000.0);
  cout <<"mbins size: "<< mBins.size() <<endl;


  int numXBins=xBins.size();
  int numQ2Bins=Q2Bins.size();
  int numZBins=zBins.size();
  int numMBins=mBins.size();


  int maxBins=numXBins;
  if(numZBins> maxBins)
    maxBins=numZBins;
  if(numMBins>maxBins)
    maxBins=numMBins;
  if(numQ2Bins>maxBins)
    maxBins=numQ2Bins;
  
  //this also sets them to zero and can be dynamic
  //		  counts[i][xBinning][xBin][phiBin][polIndex]+=weight*weightFile;
  double***** counts=allocateArray<double>(8,4,maxBins,16,2);

  double***** countsUpper=allocateArray<double>(8,4,maxBins,16,2);
  double***** countsLower=allocateArray<double>(8,4,maxBins,16,2);

  //let's do this just for the electrons
  double****** countsAll=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins,16,2);
  double****** countsAllTrueAng=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins,16,2);
  double****** countsAllTrueAngAndKin=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins,16,2);
  double****** countsUpperAll=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins,16,2);
  double****** countsLowerAll=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins,16,2);

  double **** asymAll=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  double **** asymAllTrueAng=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  double **** asymAllTrueAngAndKin=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  double **** asymErrAll=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  double **** asymErrAllTrueAng=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  double **** asymErrAllTrueAngAndKin=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  double **** asymUpperAll=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  double **** asymLowerAll=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);

  enum kinBin{kinBinQ2,kinBinX,kinBinZ,kinBinM,kinBinEnd};
    //try to just take the mean weight, instead of the value of the actual fit
  double **** meanWeightAll=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
    double **** meanWeightAllTrueFromElec=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
    double **** meanWeightAllTrue=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
    double **** meanWeightUncAll=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  double ***** kinMeansAll=allocateArray<double>(4,numQ2Bins,numXBins,numZBins,numMBins);
  double **** yMeansAll=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  double ***** kinCountsAll=allocateArray<double>(4,numQ2Bins,numXBins,numZBins,numMBins);
  double **** kinCountsAllNoFileWeight=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  double **** kinCountsAllTrue=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
    double **** kinCountsAllTrueNoFileWeight=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
  
  //need to keep track of unweighted counts to get the uncertainties of the estimates




  double*** kinMeans=allocateArray<double>(8,4,10);
  double*** kinCounts=allocateArray<double>(8,4,10);
  double*** kinCountsNoFileWeight=allocateArray<double>(8,4,10);
  //first dimention is bound (mean, upper, lower)
  double**** amps=allocateArray<double>(3,8,4,10);
  double**** ampErrs=allocateArray<double>(3,8,4,10);

  //legacy variables
  const int xBinning=kinBinX;
  const int zBinning=kinBinZ;
  const int mBinning=kinBinM;

  //temporarily..
  int numPhiBins=16;
  for(int i=0;i<numPhiBins;i++)
    {
      float bin=(i+1)*2*TMath::Pi()/numPhiBins;
      phiBins.push_back(bin);
    }


  ///---
  //    vector<string> vNoISRFilenames;
  //   char* listfile=argv[1];
  //  ifstream fListFile(listfile);
  //  string line;
  
  double polarization=0.7;
  double dilution=1.0;

  while(getline(fListFile,line))
    {
      //get filename, upper and lower Q2 cut, weight
      auto iss=istringstream(line);
      auto str=std::string();
      int minQ2=0;
      int maxQ2=0;
      double weightFile=0;
      //not too clean, since polarization and dilution are aplied to all...
      iss>> str >> minQ2 >> maxQ2 >> weightFile >> polarization >> dilution;
      cout <<"reading: " << str << " minQ2: " << minQ2 <<" maxQ2 " << maxQ2 <<" weight: " << weightFile << " polarization: "<< polarization <<" dilution: " << dilution<<endl;
      if(str.length()<10)
	{
	  cout <<"bail..." <<endl;
	  break;
	}
      
  
      //      TFile mFile(argv[1]);
      cout <<"opening root file: "<< str <<endl;
      TFile mFile(str.c_str());
      TTree* mTrees[8];
      
      for(int i=0;i<recTypeEnd;i++)
	{
	  sprintf(buffer,"tree_%s",recTypeNames[i].c_str());
	  mTrees[i]=(TTree*)mFile.Get(buffer);
	  if(mTrees[i]==0)
	    cout <<" didn't load " << buffer <<endl;
	}

      
      for(int i=0;i<recTypeEnd;i++)
	{

	  if(debug)
	    cout <<" get rec type " << i <<endl;
	  bool isTrueWithAcc=false;
	  bool isElecMethod=false;
	  if(i==0)
	    isElecMethod=true;
	  if(i==truth)
	    isTrueWithAcc=true;
	  
	  treeIndex=i;
	  if(debug)

	  sprintf(buffer,"tree_%s",recTypeNames[i].c_str());
	  if(debug)
	    cout << "getting tree: "<< buffer <<endl;
	  //      cout <<"setting branching addresse" <<endl;
	  setBranchAddresses(mTrees[treeIndex],treeFields[treeIndex]);
	  //      cout <<"done " <<endl;
	  
	  
	  long numEntries=mTrees[treeIndex]->GetEntries();
	  if(debug)
	    {
	      cout <<"we have " << numEntries << " entries " <<endl;
	      cout <<"getting " << numEntries <<" for recTypeName " << recTypeNames[i]<<endl;
	    }
	  //	  	  	  	  numEntries=100;
	  for(int ie=0;ie<numEntries;ie++)
	    {
	      if(debug)
		cout <<" entry : "<< ie <<" nr: " << treeFields[treeIndex].evtNr <<endl;
	      
	      mTrees[treeIndex]->GetEntry(ie);
	      double Q2=treeFields[treeIndex].Q2;
	      double trueQ2=treeFields[treeIndex].trueQ2;
	      double x=treeFields[treeIndex].x;
	      double trueX=treeFields[treeIndex].trueX;
	      //	      cout <<" x: "<< x <<" Q2: "<< Q2 <<endl;
	      //needed to incorporate files with different Q2 ranges
	      if(Q2< minQ2 || Q2> maxQ2)
		{
		  //~	  		  cout <<"q2 cut: "<< Q2 <<endl;
		  continue;
		}
	      else
		{
		  //		  cout <<" after Q2 : " << Q2 <<endl;
		}
	      if(debug)
		cout <<" Q2: "<< Q2 <<" x: "<< x <<endl;
	      int Q2Bin=getBin(Q2Bins,Q2);
	      int trueQ2Bin=getBin(Q2Bins,trueQ2);
	      int xBin=getBin(xBins,x);
	      int trueXBin=getBin(xBins,trueX);
	      if(xBin==0 && Q2Bin==0)
		{
		  //		  cout <<"we have xbin 0 and q2 ==0 " <<endl;
		}
	      //	      cout <<"1" <<endl;
	      if(treeFields[treeIndex].x<0.0)
		{
		  continue;
		}
	      if(treeFields[treeIndex].x>1.0)
		{
		  continue;
		}
	      //	      	      cout <<"2" <<endl;
	      if(treeFields[treeIndex].x>minXCut && treeFields[treeIndex].x<0.000160)
		{
				cout << "found x ; " << treeFields[treeIndex].x<<endl;
		}
	      //some temporary event cuts
	      if(treeFields[treeIndex].x<minXCut)
		{
		  	      continue;
		}
	      if(treeFields[treeIndex].Q2<minQ2Cut)
		{
		  	      continue;
		}
	      //	      cout <<"3" <<endl;
	      //conservative...
	      //	      cout <<"y: " << treeFields[treeIndex].y <<endl;
	      	      if(treeFields[treeIndex].y>0.85)
			{
			  continue;
			}

	      if(treeFields[treeIndex].y<0.1)
		{
		  //		  cout <<"cont due to 0.1 " << endl;
		  continue;
		}
	      if(treeFields[treeIndex].y<0.01)
		{
		  continue;
		}
	      
	      	  if(debug)
		    {
	            cout <<"x: "<< treeFields[treeIndex].x <<" bin: " << xBin<<endl;
		    cout <<"looking at " << treeFields[treeIndex].numHadronPairs<<" pairs " <<endl;
		    }

	      for(int iPair=0;iPair<treeFields[treeIndex].numHadronPairs;iPair++)
		{
		  if(debug)
		    cout <<" looking at pair " << iPair <<" of " << treeFields[treeIndex].numHadronPairs <<endl;
		  int polIndex=0;
		  double z=treeFields[treeIndex].z[iPair];
		  double M=treeFields[treeIndex].M[iPair];

		  if(treeFields[treeIndex].polarization[iPair]>0)
		    {
		      polIndex=1;
		    }

		  diHadTreeFields* fields=&(treeFields[treeIndex]);
		  
		  if(isnan(treeFields[treeIndex].weight[iPair]))
		    continue;
		  if(isnan(treeFields[treeIndex].weightUpperLimit[iPair]))
		    continue;
		  //model not valid below 0.2
		  if(treeFields[treeIndex].z[iPair]<minZ)
		    continue;
		  
		  if(treeFields[treeIndex].z[iPair]>1.0)
		    continue;
		  if(treeFields[treeIndex].M[iPair]>3.0)
		    continue;
		  if(treeFields[treeIndex].M[iPair]<minM)
		    continue;


		  if(debug)
		    {
		    cout <<" looking at pair with z: " << treeFields[treeIndex].z[iPair]<<" M: " << treeFields[treeIndex].M[iPair]<< " weight: " << treeFields[treeIndex].rawWeight[iPair]<<endl;
		    cout <<" x; " << x <<" Q2: "<< Q2 <<endl;
		    }

		  float trueM=treeFields[treeIndex].trueM[iPair];
		  float trueZ=treeFields[treeIndex].trueZ[iPair];
		  float trueQ2=treeFields[treeIndex].trueQ2;
		  float trueX=treeFields[treeIndex].trueX;
		  
		  
		  int zBin=getBin(zBins,treeFields[treeIndex].z[iPair]);
		  int trueZBin=getBin(zBins,treeFields[treeIndex].trueZ[iPair]);
		  int mBin=getBin(mBins,treeFields[treeIndex].M[iPair]);
		  int trueMBin=getBin(mBins,treeFields[treeIndex].trueM[iPair]);
		  float ang=treeFields[treeIndex].phiR[iPair]+treeFields[treeIndex].phiS[iPair];
		  float trueAng=treeFields[treeIndex].truePhiR[iPair]+treeFields[treeIndex].truePhiS[iPair];
		  //	      float ang=treeFields[treeIndex].phiR[iPair];
	      //	      	      	      float ang=treeFields[treeIndex].phiS[iPair];
		  //	            cout <<"phiS: "<< ang <<endl;

		  //	      	      cout <<" z: " << treeFields[treeIndex].z[iPair] <<endl;
		  //	      	      cout <<" M: " << treeFields[treeIndex].M[iPair] <<endl;
	      
	      //	      cout <<"1" <<endl;
		  while(ang>2*TMath::Pi())
		    {
		      ang-=(2*TMath::Pi());
		    }
		  while(ang<0)
		    {
		      ang+=(2*TMath::Pi());
		    }
		  while(trueAng>2*TMath::Pi())
		    {
		      trueAng-=(2*TMath::Pi());
		    }
		  while(trueAng<0)
		    {
		      trueAng+=(2*TMath::Pi());
		    }

		  if(debug)
		    cout <<" weight: "<< treeFields[treeIndex].weight[iPair] <<" upper: " << treeFields[treeIndex].weightUpperLimit[iPair] <<" lower: "<< treeFields[treeIndex].weightLowerLimit[iPair]<<endl;

		  if(isnan(ang) || isnan(trueAng))
		    {
		      //		      cout <<" ang is nan" <<endl;
		      continue;
		    }
		  int phiBin=getBin(phiBins,ang);
		  int truePhiBin=getBin(phiBins,trueAng);
		  
		  float weight=treeFields[treeIndex].weight[iPair];

		  double meanWeight=0;
		  double unc=0;


		  //		  cout <<"true Q2: " << sqrt(trueQ2) <<" trueM: " << trueM<<" trueX: "<< trueX <<" tureZ: "<< trueZ <<endl;
		  m_weights->getWeight(sqrt(trueQ2),trueM,trueX,trueZ,meanWeight,unc);

		  //		  cout <<"rawWeight " << treeFields[treeIndex].rawWeight[iPair] <<" rec weight: "  <<meanWeight <<endl;
		  //		  cout <<"weight " <<  treeFields[treeIndex].weight[iPair] <<" recweight: "<< 1+treeFields[treeIndex].polarization[iPair]*meanWeight*sin(trueAng)<<endl;
		  
	      //	      weight=1.0;
		  	  if(debug)
			    cout <<"2" <<endl;
		  kinMeans[i][xBinning][xBin]+=(weight*weightFile*treeFields[treeIndex].x);
		  kinMeans[i][mBinning][mBin]+=(weight*weightFile*treeFields[treeIndex].M[iPair]);
		  kinMeans[i][zBinning][zBin]+=(weight*weightFile*treeFields[treeIndex].z[iPair]);
		  //		  cout <<"2>" <<endl;

		  if(isnan(treeFields[treeIndex].x) || isnan(treeFields[treeIndex].M[iPair] || treeFields[treeIndex].z[iPair]))
		    {
		      cout <<"some pair nan" <<endl;
		    }
		  
		  if(isnan(weight))
		    {
		      		      cout <<" weight is nan" <<endl;
		    cout <<"weight nan" <<endl;
		    }
		  
		  if(isnan(treeFields[treeIndex].weightUpperLimit[iPair]))
		    cout <<"weight nan upper" <<endl;
		  
		  if(isnan(treeFields[treeIndex].weightLowerLimit[iPair]))
		    cout <<"weight nan lower" <<endl;


		  if(isTrueWithAcc)
		    {
		      double rawWeight=treeFields[treeIndex].rawWeight[iPair];
		      if(debug)
			{
			cout <<"mean weight true before " << meanWeightAllTrue[Q2Bin][xBin][zBin][mBin]<<endl;
			cout <<" kin counts before: "<<  kinCountsAllTrue[Q2Bin][xBin][zBin][mBin]<<endl;
			cout <<" adding : " <<  rawWeight*weightFile << " and count; "<< weight*weightFile<<endl;
			}
		      //		      if(rawWeight>0)
			{
			  meanWeightAllTrue[Q2Bin][xBin][zBin][mBin]+=(rawWeight*weightFile);
			  kinCountsAllTrue[Q2Bin][xBin][zBin][mBin]+=weight*weightFile;
			  kinCountsAllTrueNoFileWeight[Q2Bin][xBin][zBin][mBin]+=weight;
			}
		      if(debug)
			{
			  cout <<" meman weight now: "  <<  meanWeightAllTrue[Q2Bin][xBin][zBin][mBin]<<endl;
			  cout <<" kin counts now: "<<  kinCountsAllTrue[Q2Bin][xBin][zBin][mBin]<<endl;
			    
}
			 
		      //		      cout <<"adding true pair weight: "<< rawWeight <<" q2bin: "<< Q2Bin<< " xbin: "<< xBin <<" zBin: "<< zBin <<" mBin: " << mBin <<endl;
		    }
		
		  
		  if(isElecMethod)
		    {
		      double rawWeight=treeFields[treeIndex].rawWeight[iPair];
		      double rawWeightUnc=treeFields[treeIndex].rawWeightUnc[iPair];

		      if(xBin==0)
			{
			  //			  cout <<"getting raw Weight :" << rawWeight <<" for x " << x <<endl;
			}
		      if(debug)
			{
			  cout <<"2>>" <<endl;
			  cout <<"xbin: "<< xBin <<" zBin: " << zBin << " mBin : " << mBin <<endl;
			
			}
		      if(debug)
			{
			  cout <<"fdsafds q2bin: "<< Q2Bin <<" xbin: "<< xBin <<" zbin: "<< zBin <<" mbin: "<< mBin <<endl;
			  cout <<"mean weight elec before " << meanWeightAll[Q2Bin][xBin][zBin][mBin]<<endl;
			  cout <<" kin counts before: "<<  kinCountsAll[Q2Bin][xBin][zBin][mBin]<<endl;
			  cout <<" adding : " <<  rawWeight*weightFile << " and count; "<< weight*weightFile<<endl;
			}

		      //		      if(rawWeight>0)
		      //			{
			  meanWeightAll[Q2Bin][xBin][zBin][mBin]+=(rawWeight*weightFile);
		      //		      cout <<"adding rec pair weight: "<< rawWeight <<" q2bin: "<< Q2Bin<< " xbin: "<< xBin <<" zBin: "<< zBin <<" mBin: " << mBin <<endl;

		      
			  meanWeightUncAll[Q2Bin][xBin][zBin][mBin]+=(rawWeightUnc*weightFile);


			  yMeansAll[Q2Bin][xBin][zBin][mBin]+=(weight*weightFile*treeFields[treeIndex].y);
		      kinMeansAll[kinBinQ2][Q2Bin][xBin][zBin][mBin]+=(weight*weightFile*Q2);
		      
		      kinMeansAll[kinBinX][Q2Bin][xBin][zBin][mBin]+=(weight*weightFile*x);
		      kinMeansAll[kinBinZ][Q2Bin][xBin][zBin][mBin]+=(weight*weightFile*z);
		      if(debug)
			cout <<"2>>0-" <<endl;

		      kinMeansAll[kinBinM][Q2Bin][xBin][zBin][mBin]+=(weight*weightFile*M);
		      kinCountsAll[kinBinQ2][Q2Bin][xBin][zBin][mBin]+=weight*weightFile;
		      kinCountsAllNoFileWeight[Q2Bin][xBin][zBin][mBin]+=1;
		      if(debug)
			{
			cout <<"mean weight elec now " << meanWeightAll[Q2Bin][xBin][zBin][mBin]<<endl;
			cout <<" kin counts now: "<<  kinCountsAll[kinBinQ2][Q2Bin][xBin][zBin][mBin]<<endl;
			}


		      if(debug)
			cout <<"2>>>" <<endl;
		      
		      countsAll[Q2Bin][xBin][zBin][mBin][phiBin][polIndex]+=weight;
		      countsAllTrueAng[Q2Bin][xBin][zBin][mBin][truePhiBin][polIndex]+=weight;
		      //don't need to recalculate the weight, just put it at the correct kinematics...
		      if(Q2Bin!=trueQ2Bin || xBin!=trueXBin || mBin!=trueMBin || zBin!=trueZBin)
			{
			  //			  cout <<"Q2Bin: "<< Q2Bin <<" trueBin: "<< trueQ2Bin <<" xbin : " << xBin <<" trueXBin: "<< trueXBin <<" zbin: "<< zBin <<" true z bin : "<< trueZBin <<" mbin: " << mBin <<" trueMBin " << trueMBin <<endl;
			}
		      else
			{
			  //			  cout <<"same bins " <<endl;
			}


		      countsAllTrueAngAndKin[trueQ2Bin][trueXBin][trueZBin][trueMBin][truePhiBin][polIndex]+=weight;
		      countsUpperAll[Q2Bin][xBin][zBin][mBin][phiBin][polIndex]+=treeFields[treeIndex].weightUpperLimit[iPair];
		      countsLowerAll[Q2Bin][xBin][zBin][mBin][phiBin][polIndex]+=treeFields[treeIndex].weightLowerLimit[iPair];
		      //		      kinCountsNoFileWeight[Q2Bin][xBin][zBin][mBin]+=1;
		      //		      cout <<"2>>>>" <<endl;
		    }


		  //		  cout <<"2-" <<endl;
		  kinCounts[i][xBinning][xBin]+=(weight*weightFile);
		  kinCounts[i][mBinning][mBin]+=(weight*weightFile);
		  kinCounts[i][zBinning][zBin]+=(weight*weightFile);
		  kinCountsNoFileWeight[i][xBinning][xBin]+=(weight);
		  kinCountsNoFileWeight[i][mBinning][mBin]+=(weight);
		  kinCountsNoFileWeight[i][zBinning][zBin]+=(weight);
		  
		  //		  cout <<"2-0" <<endl;
		  //		  cout <<"2-0-" <<endl;
		  //		  cout <<" fields: "<< treeFields[treeIndex].x <<endl;
		  //		  cout <<"ipair: " << iPair <<endl;
		  //		  cout << " treeF: " << treeFields[treeIndex].weightUpperLimit[iPair]*weightFile <<endl;
		  //		  cout <<" we want to fill " << i <<" xBinning " << xBinning <<" xbin: "<< xBin <<" phiBin: "<< phiBin <<" polIndex: "<< polIndex <<endl;
		  counts[i][xBinning][xBin][phiBin][polIndex]+=weight;
		  countsUpper[i][xBinning][xBin][phiBin][polIndex]+=treeFields[treeIndex].weightUpperLimit[iPair];
		  //		  		  cout <<"2-0-0" <<endl;
		  countsLower[i][xBinning][xBin][phiBin][polIndex]+=treeFields[treeIndex].weightLowerLimit[iPair];

		  
		  //		  		  cout <<"2-0-0-" <<endl;
		  counts[i][mBinning][mBin][phiBin][polIndex]+=weight;
		  //		  		  cout <<"2-00" <<endl;
		  countsUpper[i][mBinning][mBin][phiBin][polIndex]+=treeFields[treeIndex].weightUpperLimit[iPair];
		  
		  countsLower[i][mBinning][mBin][phiBin][polIndex]+=treeFields[treeIndex].weightLowerLimit[iPair];
		  
		  //	      	      cout <<"3" <<endl;
	      //	      cout <<" i: "<< i <<" binning: "<< zBinning <<" zbin: "<< zBin <<" phiBin: "<< phiBin <<" pol: "<< polIndex <<endl;
		  counts[i][zBinning][zBin][phiBin][polIndex]+=weight;
		  countsUpper[i][zBinning][zBin][phiBin][polIndex]+=treeFields[treeIndex].weightUpperLimit[iPair];
		  countsLower[i][zBinning][zBin][phiBin][polIndex]+=treeFields[treeIndex].weightLowerLimit[iPair];
		  //	      cout <<"3+" <<endl;
		}
	    }
	}
    }
  cout <<"after all the events... " << endl;
  //  	      cout <<"4" <<endl;
  for(int i=0;i<recTypeEnd;i++)
    {
      
      bool isElecMethod=false;
      if(i==0)
	isElecMethod=true;	

      cout <<" looking at " << recTypeNames[i] <<endl;
      //  enum kinBin{kinBinQ2,kinBinX,kinBinZ,kinBinM,kinBinEnd};
      cout <<"what???" <<endl;
      for(int binning=kinBinQ2;binning<kinBinEnd;binning++)
	{
	  int numKinBins=xBins.size();
	  switch(binning)
	    {
	    case kinBinQ2:
	      if(debug)
		cout <<" q2 binss " << endl;
	      numKinBins=Q2Bins.size();
	      break;
	    case kinBinX:
	      if(debug)
		cout <<" x binss " << endl;

	      numKinBins=xBins.size();
	      break;
	    case kinBinZ:
	      if(debug)
		cout <<" z binss " << endl;

	      numKinBins=zBins.size();
	      break;
	    case kinBinM:
	      if(debug)
		cout <<" m binss " << endl;

	      numKinBins=mBins.size();
	      break;
	    }
	  /////-----

	  cout << "we have " << numKinBins <<" bins" <<endl;
	  //////---
	  
	  for(int kinBin=0;kinBin<numKinBins;kinBin++)
	    {
//	      double y[3][20];
//	      double ey[3][20];
//	      double x[3][20];
//	      double ex[3][20];

	      
	      double** locCounts[3];
	      cout <<" looking at i: "<< i << " binning: "<< binning <<" kinBin : " << kinBin <<endl;
	      locCounts[mean]=counts[i][binning][kinBin];
	      locCounts[upper]=countsUpper[i][binning][kinBin];
	      locCounts[lower]=countsLower[i][binning][kinBin];
	      double fileWeightFactor=sqrt(kinCountsNoFileWeight[i][binning][kinBin]/kinCounts[i][binning][kinBin]);

	      for(int iBound=mean;iBound<endBound;iBound++)
		{

		  cout <<"looking at bound " << iBound <<" i: "<< i <<" binning: "<< binning <<" kinBin: "<< kinBin << endl;
		  pair<double,double> fitRes= getA(locCounts[iBound], phiBins,kinBin,recTypeNames[i].c_str(),boundNames[iBound].c_str(),binning,"normal");
		  cout <<"after fr " <<endl;
		  cout <<"setting amps to " << fitRes.first <<endl;
		  amps[iBound][i][binning][kinBin]=fitRes.first;
		  ampErrs[iBound][i][binning][kinBin]=fitRes.second*fileWeightFactor;
		  cout <<"amps and amp err " << iBound <<endl;
/////		  for(int phiBin=0;phiBin<phiBins.size();phiBin++)
/////		    {
/////		      double Nup=locCounts[iBound][phiBin][1];
/////		      double Ndown=locCounts[iBound][phiBin][0];
/////		      double A=(Nup-Ndown)/(Nup+Ndown);
/////		      
/////		      cout <<"phiBin: "<< phiBin << " Nup: "<< Nup <<" Ndown: " << Ndown <<" A: " << A <<endl;
/////		      y[iBound][phiBin]=A;
/////		      x[iBound][phiBin]=(phiBin+0.5)*2*TMath::Pi()/numPhiBins;
/////		      ex[iBound][phiBin]=0.0;
/////
/////		      double eU=sqrt(Nup);
/////		      double eD=sqrt(Ndown);
/////		      
/////		      double uDeriv=2*Ndown/((Nup+Ndown)*(Nup+Ndown));
/////		      double dDeriv=2*Nup/((Nup+Ndown)*(Nup+Ndown));
/////		      
/////		      ey[iBound][phiBin]=sqrt(uDeriv*uDeriv*eU*eU+dDeriv*dDeriv*eD*eD);
/////		  
/////		    }
/////		  TGraphErrors g(phiBins.size(),x[iBound],y[iBound],ex[iBound],ey[iBound]);
/////		  
/////		  gStyle->SetOptFit(111);
/////		  sprintf(buffer,"graphFor_rec_%s_binning%d_kinBin%d_bound_%s.png",recTypeNames[i].c_str(),binning,kinBin,boundNames[iBound].c_str());
/////		  TCanvas c1;
/////		  TF1 f1("f1","[0]*sin(x)",0,2*M_PI);
/////		  f1.SetParameters(0,0.0);
/////		  g.Fit(&f1);
/////		  g.Draw("AP");
/////		  c1.SaveAs(buffer);
/////		  amps[iBound][i][binning][kinBin]=f1.GetParameter(0);
/////		  ampErrs[iBound][i][binning][kinBin]=f1.GetParError(0);
/////		  cout << recTypeNames[i] <<", binning: "<< binning<< " kinBin: "<< kinBin <<endl;
/////		  cout <<"amp for bound " << iBound<<" " << amps[iBound][i][binning][kinBin] <<endl;
/////		}
		}
	      
	    }
	}
    }
  for(int q2Bin=0;q2Bin<Q2Bins.size();q2Bin++)
    {
      for(int xBin=0;xBin<xBins.size();xBin++)
	{
	  for(int zBin=0;zBin<zBins.size();zBin++)
	    {
	      for(int mBin=0;mBin<mBins.size();mBin++)
		{

		  sprintf(buffer,"q2Bin_%i_xbin_%i_zbin_%i_mbin_%i",q2Bin,xBin,zBin,mBin);
		  sprintf(buffer1,"%s_all",buffer);
		  double fileWeightFactor=sqrt(kinCountsAllNoFileWeight[q2Bin][xBin][zBin][mBin]/kinCountsAll[kinBinQ2][q2Bin][xBin][zBin][mBin]);

		  
		  //	      pair<double,double> fitRes= getA(locCounts[mean], phiBins,kinBin,recTypeNames[i].c_str(),boundNames[iBound].c_str(),binning);
		  //	      		  countsAll[Q2Bin][xBin][zBin][mBin][phiBin][polIndex];
		  pair<double,double> fitRes=getA(countsAll[q2Bin][xBin][zBin][mBin],phiBins,0,"_all","_mean",0,buffer1);
		  sprintf(buffer1,"%s_trueAng",buffer);
		  pair<double,double> fitResTrueAng=getA(countsAllTrueAng[q2Bin][xBin][zBin][mBin],phiBins,0,"_all","_mean",0,buffer1);
		  sprintf(buffer1,"%s_trueAngAndKin",buffer);
		  pair<double,double> fitResTrueAngAndKin=getA(countsAllTrueAngAndKin[q2Bin][xBin][zBin][mBin],phiBins,0,"_all","_mean",0,buffer1);
		  sprintf(buffer1,"%s_normal",buffer);

		  pair<double,double> fitResUpper=getA(countsUpperAll[q2Bin][xBin][zBin][mBin],phiBins,0,"_all","_upper",0,buffer1);
		  
		  pair<double,double> fitResLower=getA(countsLowerAll[q2Bin][xBin][zBin][mBin],phiBins,0,"_all","_lower",0,buffer1);

		  cout <<"fit for all q2bin: "<< q2Bin <<" xBin: "<< xBin <<" zBin " << zBin <<" mBin: "<< mBin << " A: " <<fitRes.first<<endl;
		  
		  asymAll[q2Bin][xBin][zBin][mBin]=fitRes.first;
		  asymAllTrueAng[q2Bin][xBin][zBin][mBin]=fitResTrueAng.first;
		  asymAllTrueAngAndKin[q2Bin][xBin][zBin][mBin]=fitResTrueAngAndKin.first;
		  asymUpperAll[q2Bin][xBin][zBin][mBin]=fitResUpper.first;
		  asymLowerAll[q2Bin][xBin][zBin][mBin]=fitResLower.first;
		  //only care for uncertainties of the mean asym
		  asymErrAll[q2Bin][xBin][zBin][mBin]=fitRes.second*fileWeightFactor;
		  asymErrAllTrueAng[q2Bin][xBin][zBin][mBin]=fitResTrueAng.second/fileWeightFactor;
		  asymErrAllTrueAngAndKin[q2Bin][xBin][zBin][mBin]=fitResTrueAngAndKin.second*fileWeightFactor;;
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
  
  for(int binning=kinBinQ2;binning<kinBinEnd;binning++)
    {
      TCanvas c;
      float maxX=0.8;
      for(int i=0;i<recTypeEnd;i++)
	{
	  double y[30];
	  double x[30];
	  double ey[30];
	  double ex[30];

	  //the rest for the bounds are the same
	  double boundEYUpper[30];
	  double boundEYLower[30];

	  switch(binning)
	    {
	    case kinBinQ2:
	      numKinBins=Q2Bins.size();
	      break;
	    case kinBinX:
	      numKinBins=xBins.size();
	      break;
	    case kinBinZ:
	      numKinBins=zBins.size();
	      break;
	    case kinBinM:
	      numKinBins=mBins.size();
	      break;
	    }

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

	  cout <<" moving on.. " <<endl;
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

	  if(binning==kinBinX)
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

	  if(binning==kinBinZ)
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
	  if(binning==kinBinM)
	    {
	      if(i==0)
		{
		  maxX=2.0;
		gPad->DrawFrame(0.0,-0.3, maxX, 0.3);
		g->GetXaxis()->SetLimits(0.0,maxX);
		}
	      g->GetXaxis()->SetRangeUser(0.0,maxX);
	      g->GetYaxis()->SetRangeUser(-0.02,0.15);
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
	    
	  //	        spxbrintf(buffer,"amps_binning_%s_recType_%s.png",binningNames[binning].c_str(),recTypeNames[i].c_str());
	}
      TLine *line = new TLine(0,0,maxX,0);
      line->Draw();
      sprintf(buffer,"amps_binning_%s.png",binningNames[binning].c_str());
      c.SaveAs(buffer);
    }

  cout <<"about to start fancy plots " <<endl;
  //////---------
  bool isLogFancyX=false;
  TCanvas fancyPlot("fancy","fancy",10,10,6000,4000);
  ///make the axis
  TH2D axis("myAx","",100,minX_fancy,maxX_fancy,100,minQ2_fancy,maxQ2_fancy);
  axis.GetXaxis()->SetTitle("x");
  axis.GetYaxis()->SetTitle("Q^{2} [GeV]");

  axis.Draw();
  //needed to get the updated coordinates
  fancyPlot.Update();
  //  (xbins, Q2 bins
  

  if(isLogFancyX)
    {
      fancyPlot.SetLogx();
    }
  fancyPlot.SetLogy();
  fancyPlot.Update();
  
  cout <<"xBins size: "<< xBins.size() <<endl;
  for(int iQ=0;iQ<Q2Bins.size();iQ++)
	{
	    for(int iX=0;iX<xBins.size();iX++)
    {
      cout <<"bin " << xBins[iX] <<endl;

	  sprintf(buffer,"pad_xbin_%d_q2bin_%d",iX,iQ);

	  //	  logMinQ2, logMinX
	  double xlow=log10(minXCut);
	  if(!isLogFancyX)
	    {
	      xlow=minXCut;
	    }
	  double ylow=log10(minQ2Cut);

	  double xlowUser=minXCut;
	  double ylowUser=minQ2Cut;
	  double xupUser=xBins[iX];
	  double yupUser=Q2Bins[iQ];
	  if(iX>0)
	    {
	      xlow=xBins[iX-1];
	      if(isLogFancyX)
		{
		  xlow=log10(xlow);
		}
	      xlowUser=xBins[iX-1];
	    }
	  if(iQ>0)
	    {
	      ylow=log10(Q2Bins[iQ-1]);
	      ylowUser=Q2Bins[iQ-1];
	    }
	  cout << " iX: "<< iX <<" iQ: "<< iQ <<" xlow: "<< xlow <<" ylow: "<< ylow <<endl;
	  ///this might not be necessary, since the user coordinate system of the canvas can be used
	  double xup=xBins[iX];
	  if(isLogFancyX)
	    {
	      xup=log10(xup);
	    }
	  double yup=log10(Q2Bins[iQ]);
	  if(xup>maxX_fancy)
	    xup=maxX_fancy;
	  if(isLogFancyX)
	    {
	      if(xup>logMaxX)
		xup=logMaxX;
	    }
	  
	  if(yup>logMaxQ2)
	    yup=logMaxQ2;
	  //if we updated after the log setting, the getX1 is actually giving the decade.
	  double plotXRange=fancyPlot.GetX2()-fancyPlot.GetX1();
	  double plotQ2Range=fancyPlot.GetY2()-fancyPlot.GetY1();
	  //logmaxQ2
	  
	  cout <<" logmaxX: "<< logMaxX <<" logminx: "<< logMinX << " plotxrange: "<< plotXRange <<endl;
	  cout <<" xlow: " << xlow <<endl;
	  xlow=(xlow-fancyPlot.GetX1())/plotXRange;
	  xup=(xup-fancyPlot.GetX1())/plotXRange;

	  ylow=(ylow-fancyPlot.GetY1())/plotQ2Range;
	  yup=(yup-fancyPlot.GetY1())/plotQ2Range;
	  //////
	  
	  //might happen with the last bin, since it has a 'catch-all'bound
	  if(yup > 1.0)
	    yup=1.0;
	  if(xup>1.0)
	    xup=1.0;
	  if(yupUser>maxQ2Fancy)
	    yupUser=maxQ2Fancy;
	  if(xupUser>maxXFancy)
	    xupUser=maxXFancy;

	  cout <<"xupUser: " << xupUser <<" xlowUser: "<< xlowUser <<" yupUser: "<< yupUser <<" ylowUser: " << ylowUser <<endl;
	  cout <<"x1: "<< fancyPlot.GetX1() <<" x2: "<< fancyPlot.GetX2() <<" y1: " << fancyPlot.GetY1() <<" y2: " << fancyPlot.GetY2()<<endl;
	  ////we have the co-ordinates in user space, so let's convert to NDC
	  //	  xup=(log10(xupUser)-log10(fancyPlot.GetX1()))/(log10(fancyPlot.GetX2())-log10(fancyPlot.GetX1()));
	  //	  xlow=(log10(xlowUser)-log10(fancyPlot.GetX1()))/(log10(fancyPlot.GetX2())-log10(fancyPlot.GetX1()));
	  //	  yup=(log10(yupUser)-log10(fancyPlot.GetY1()))/(log10(fancyPlot.GetY2())-log10(fancyPlot.GetY1()));
	  //	  ylow=(log10(ylowUser)-log10(fancyPlot.GetY1()))/(log10(fancyPlot.GetY2())-log10(fancyPlot.GetY1()));
	  cout <<"pad with xlow: "<< xlow <<" ylow: " << ylow <<" xup: "<< xup <<" yup: " << yup <<endl;
	  TPad* myPad=new TPad(buffer,buffer,xlow, ylow, xup, yup,kBlue,2,1);//,kBlue,2,1);
	  myPad->Draw();
	  myPad->cd();

	  ////////------------------divide again for M-----------
	  ///num m bins
	  //////--------
	  vector<float> outerBinning=mBins;
	  vector<float> innerBinning=zBins;

	  cout <<"innerBinning size: "<< innerBinning.size()<<endl;
	  if(vsM)
	    {
	      cout <<" vs M switch " << endl;
	      outerBinning=zBins;
	      innerBinning=mBins;
	      cout <<"innerBinning size now : "<< innerBinning.size()<<endl;
	      cout <<"mBins size now : "<< mBins.size()<<endl;
	    }
	  	      cout <<"innerBinning size and now : "<< innerBinning.size()<<endl;
	    for(int io=0;io<outerBinning.size();io++)
	      {
		xlow=0.0;
		ylow=io*1.0/outerBinning.size();
		xup=1.0;
		yup=(io+1)*1.0/outerBinning.size();
		cout <<"put pad in pad with xlow: "<< xlow <<" ylow: "<< ylow <<" xup: " << xup <<" yup: "<< yup <<endl;
		
		sprintf(buffer, "innerPad_iq%d_ix%d_io%d",iQ,iX,io);
		TPad* myPadInner=new TPad(buffer,buffer,xlow, ylow, xup, yup);//,kBlue,2,1);
		//
		if(io!=outerBinning.size()-1)
		  {
		    myPadInner->SetTopMargin(0.0);
		  }
		if(io!=0)
		  {
		    myPadInner->SetBottomMargin(0.0);
		    double xmin=0;
		    double ymin=0;
		    double xmax=0;
		    double ymax=0;
		    myPadInner->GetRangeAxis(xmin, ymin, xmax ,ymax);
		    cout <<" got axis ranges: xmin: "<< xmin << " ymin: " << ymin <<" xmax: "<< xmax <<" ymax: " << ymax <<endl;
		    //no point of this
		    //		    myPadInner->RangeAxis(xmin,ymin,xmax,ymax);
		    myPadInner->Update();
		  }
		myPadInner->Update();
		myPadInner->Draw();
		//never know with root...just do it twice to be sure
		myPadInner->Update();
		myPadInner->Draw();
		myPadInner->cd();
		//	  myPad->SetFillColor(kBlue);
		///draw on the pad...
		double x[40];
		double ex[40];
		double y[40];
		double yTrue[40];
		double yTrueKins[40];
		double yTrueAng[40];
		double yTrueAngAndKin[40];
		double ySys[40];
		double ey[40];
		double eyTrueAng[40];
		double eyTrueAngAndKin[40];
	     
		double meanWeight=0;
		double unc=0;

		int filledInnerBins=0;
		cout <<"inner binning size: "<<innerBinning.size()<<endl;
		for(int ii=0;ii<innerBinning.size();ii++)
		  {
		    int iz=ii;
		    int im=io;
		    if(vsM)
		      {
			iz=io;
			im=ii;
		      }
		    //iz is the ex inner
		    cout <<" ii: "<< ii <<" filledInnerBin: "<< filledInnerBins <<endl;
		    double counts=0;
		    counts=kinCountsAll[kinBinQ2][iQ][iX][iz][im];
		    
		    if(counts<minCounts)
		      {
			cout <<"mincounts iQ: "<< iQ <<" iX :"<< iX <<" iz: " << iz <<" im: "<< im <<" counts: "<< counts <<endl;
			continue;
		      }
		    //	      	      asymAll[q2Bin][xBin][zBin][mBin]=fitRes.first;
		    //		      kinMeansAll[kinBinM][Q2Bin][xBin][zBin][mBin]+=(weight*weightFile*M);
		    //counts All for Q2 binning, since the counts are all the same
		    double locQ2=0;
		    double locM=0;
		    double locZ=0;
		    double locX=0;
		    double locY=0;
		    double locMeanWeightAll;
		    double locMeanWeightAllTrue;
		    double locMeanWeightUncAll;


		    locY=yMeansAll[iQ][iX][iz][im]/kinCountsAll[kinBinQ2][iQ][iX][iz][im];
		    locQ2=kinMeansAll[kinBinQ2][iQ][iX][iz][im]/kinCountsAll[kinBinQ2][iQ][iX][iz][im];
		    locM=kinMeansAll[kinBinM][iQ][iX][iz][im]/kinCountsAll[kinBinQ2][iQ][iX][iz][im];		    
		    locZ=kinMeansAll[kinBinZ][iQ][iX][iz][im]/kinCountsAll[kinBinQ2][iQ][iX][iz][im];	      
		    locX=kinMeansAll[kinBinX][iQ][iX][iz][im]/kinCountsAll[kinBinQ2][iQ][iX][iz][im];
		      
		    locMeanWeightAll=meanWeightAll[iQ][iX][iz][im]/kinCountsAll[kinBinQ2][iQ][iX][iz][im];
		    locMeanWeightAllTrue=meanWeightAllTrue[iQ][iX][iz][im]/kinCountsAllTrue[iQ][iX][iz][im];
		    cout <<"locMeanWeightAllTrue: "<< locMeanWeightAllTrue <<" meanWeight: "<< meanWeightAllTrue[iQ][iX][iz][im] <<" kincounts: "<< kinCountsAllTrue[iQ][iX][iz][im]<<endl;		      
		    locMeanWeightUncAll=meanWeightUncAll[iQ][iX][iz][im]/kinCountsAll[kinBinQ2][iQ][iX][iz][im];
		    if(iX==0)
		      {
			cout <<"mean weight all:  "<< meanWeightAll[iQ][iX][iz][im] <<" counts: "<< kinCountsAll[kinBinQ2][iQ][iX][iz][im]<<endl;
			cout << iQ <<" im: "<< im <<" iz : "<< iz <<endl;
		      }
		    
		    x[filledInnerBins]=locZ;
		    if(vsM)
		      x[filledInnerBins]=locM;
		    
		
		    //
		    //
		    //	      cout <<"iz: "<< iz <<" x[iz]: " << x[iz] <<endl;
		    ex[filledInnerBins]=0.0;
		    //just put the model expections on there
		    //
		    cout <<"getting weight for point " << sqrt(locQ2) << " m: "<< locM <<" x: "<< locX <<" z: "<< locZ <<endl;
		    
		    m_weights->getWeight(sqrt(locQ2),locM,locX,locZ,meanWeight,unc);
		    if(useRealAsym)
		      {cout <<"setting real asym for q2: "<< iQ <<" ix: "<< iX <<" iz: " << iz <<" im: "<< im <<" to " << asymAll[iQ][iX][iz][im]<<endl;
			y[filledInnerBins]=asymAll[iQ][iX][iz][im];
			yTrue[filledInnerBins]=locMeanWeightAllTrue;
			yTrueAng[filledInnerBins]=asymAllTrueAng[iQ][iX][iz][im];
			yTrueAngAndKin[filledInnerBins]=asymAllTrueAngAndKin[iQ][iX][iz][im];
			
		      }
		    else
		      {
			//			y[filledZBins]=meanWeight;
			cout <<"ix: "<< iX <<" iQ : "<< iQ << " iz: " << iz <<" im: "<< im << " mean weight: "<< locMeanWeightAll<<endl;
			y[filledInnerBins]=locMeanWeightAll;
			yTrue[filledInnerBins]=locMeanWeightAllTrue;
		      }
		    //		    ySys[filledZBins]=unc;
		    ySys[filledInnerBins]=locMeanWeightUncAll;
		    cout <<"using weight: "<< meanWeight <<" unc: "<< unc<<endl;
		    if(useRealAsym)
		      {
			  ey[filledInnerBins]=asymErrAll[iQ][iX][iz][im];
			  eyTrueAng[filledInnerBins]=asymErrAllTrueAng[iQ][iX][iz][im];
			  eyTrueAngAndKin[filledInnerBins]=asymErrAllTrueAngAndKin[iQ][iX][iz][im];
		      }
		    else
		      {
			cout <<"doing ey" <<endl;
		    //factor 1.5 is empirical
			  ey[filledInnerBins]=1.5/sqrt(kinCountsAll[kinBinQ2][iQ][iX][iz][im]);

			cout <<"done with ey " << endl;
		      }
		    //account for polarization
		    ey[filledInnerBins]/=(polarization*dilution);
		    cout <<"dada " <<endl;

		    //		      double **** kinCountsNoFileWeight=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
		    cout <<" hey hey " <<endl;
		    cout <<"wanting to get iQ: "<< iQ <<" iX: "<< iX <<" iz: "<< iz <<" im:" << im <<endl;
		    //		      double **** kinCountsNoFileWeight=allocateArray<double>(numQ2Bins,numXBins,numZBins,numMBins);
		    //		    cout <<"we want to divide by " <<sqrt(kinCountsNoFileWeight[iQ][iX][iz][im]) <<endl;
		    
		    double relSimuUncert=0;
		    relSimuUncert=1/sqrt(kinCountsAllNoFileWeight[iQ][iX][iz][im]);
		    cout <<"writing to file .." <<endl;
		    double rel=0;
		    if(yTrue[filledInnerBins]>0)
		      rel=(y[filledInnerBins]-yTrue[filledInnerBins])/yTrue[filledInnerBins];
		    
		    //		    outFile << iQ <<" " << iX << " " <<im <<" " << iz <<  " " << locQ2 <<" " << locX << " " << locM <<" " <<locZ << " " << locY<<" "   y[filledInnerBins] <<" +- "<< ey[filledInnerBins] << " +- " << y[filledInneinc
		    outFile << iQ <<" " << iX << " " <<iz <<" " << im <<  " " << locQ2 <<" " << locX << " " << locZ <<" " <<locM << " " <<locY << " " <<  y[filledInnerBins] <<" +- "<< ey[filledInnerBins] << " " <<relSimuUncert  <<  " "<<kinCountsAllNoFileWeight[iQ][iX][iz][im]<<endl;
		    //		    cout <<"asymErr: "<< asymErrAll[iQ][iX][iz][im] <<" naive: "<< 1.0/sqrt(kinCountsAll[kinBinQ2][iQ][iX][iz][im]) <<endl;
		    cout << "iQ:  " <<iQ <<" ix: "<< iX <<" im: "<< im <<" iz: "<< iz <<" y val: : " <<  y[filledInnerBins]  <<" uncertainty: "<<  ey[filledInnerBins] <<endl;
		    filledInnerBins++;
		  }
	     
		cout <<"we have " << filledInnerBins << " filled inner bins " << endl;
	      	    
		sprintf(buffer,"fancyGraph_xbin_%d_q2bin_%d",iX,iQ);
		TGraphErrors* gr=new TGraphErrors(filledInnerBins,x,y,ex,ey);
		TGraphErrors* grTruth=new TGraphErrors(filledInnerBins,x,yTrue,ex,ey);
		TGraphErrors* grTrueAng=new TGraphErrors(filledInnerBins,x,yTrueAng,ex,ey);
		TGraphErrors* grTrueAngAndKin=new TGraphErrors(filledInnerBins,x,yTrueAngAndKin,ex,ey);
		TGraphErrors* grSys=new TGraphErrors(filledInnerBins,x,y,ex,ySys);
		gr->SetTitle("");
		gr->SetMarkerStyle(markerStyles[0]);
		gr->SetMarkerColor(colors[0]);

		grTruth->SetTitle("");
		grTrueAng->SetTitle("");
		grTrueAngAndKin->SetTitle("");
		grTruth->SetMarkerStyle(markerStyles[2]);
		grTrueAng->SetMarkerStyle(markerStyles[2]);
		grTrueAngAndKin->SetMarkerStyle(markerStyles[2]);
		grTruth->SetMarkerColor(kBlack);
		grTrueAng->SetMarkerColor(kBlue);
		grTrueAngAndKin->SetMarkerColor(kMagenta);


		grSys->SetTitle("sys");
		grSys->SetMarkerStyle(markerStyles[0]);
		grSys->SetMarkerColor(kGreen);
		grSys->SetLineWidth(0);
		grSys->SetLineColor(kGreen);
		grSys->SetFillColor(kGreen);
		
		myPadInner->cd();
		gr->Draw("AP");
		if(io!=0)
		  {
		    gr->GetXaxis()->SetDrawOption("B");
		  }
		double maxX=0.85;
		if(vsM)
		  maxX=1.5;
		gr->GetXaxis()->SetLimits(0.0,maxX);
		
		TLine *line = new TLine(0,0,maxX,0);
		line->Draw();
		gr->GetXaxis()->SetRangeUser(0.0,maxX);
		gr->GetYaxis()->SetRangeUser(-0.02,0.2);
		gr->GetXaxis()->SetTitle("z");
		gr->GetXaxis()->SetTitleSize(0.05);
		gr->GetXaxis()->SetTitleOffset(0.32);
		gr->GetXaxis()->SetLabelSize(0.05);
		if(vsM)
		  gr->GetXaxis()->SetTitle("M");
		gr->GetYaxis()->SetTitle("A");
		gr->SetName(buffer);
		gr->Draw("AP");
		grSys->Draw("SAME 3");
		gr->Draw("SAME P");
		///		grTruth->Draw("SAME P");
		///		grTrueAng->Draw("SAME P");
		//		grTrueAngAndKin->Draw("SAME P");
		line->Draw();
		if(iX==0)
		  {
		    float lowerM=0.0;
		    float lowerZ=0.0;
		    if(vsM)
		      {
			if(io==0)
			  lowerZ=minZ;
			else
			  lowerZ=zBins[io-1];
		      }
		    else
		      {
			if(io==0)
			  lowerM=minM;
			else
			  lowerM=mBins[io-1];
		      }
		    
		    sprintf(buffer,"%.2f < M < %.2f",lowerM,outerBinning[io]);
		    if(vsM)
		      {
			sprintf(buffer,"%.2f < z < %.2f",lowerZ,outerBinning[io]);
		      }
		   TText *t = new TText(.3,.8,buffer);
		   t->SetNDC(true);
		   t->SetTextAlign(22);
		   t->SetTextColor(kRed+2);
		   t->SetTextFont(43);
		   t->SetTextSize(15);
		   //		   t->SetTextAngle(45);
		   t->Draw();
		  }
		//go back to mother pad to place the new one.
		myPad->cd();
	      }
	  ///---
	  //otherwise we'll do canvas in canvas...
	  fancyPlot.cd();
	        cout <<"ix: " << iX << " xbins size: " << xBins.size() <<endl;
	}

    }
  outFile.flush();
  outFile.close();
  char vs[10];
  if(vsM)
    sprintf(vs,"M");
  else
    sprintf(vs,"z");
  
  sprintf(buffer,"fancyPlot_Vs%s.png",vs);
  fancyPlot.SaveAs(buffer);
  sprintf(buffer,"fancyPlot_Vs%s.pdf",vs);
  fancyPlot.SaveAs(buffer);
  sprintf(buffer,"fancyPlot_Vs%s.root",vs);
  fancyPlot.SaveAs(buffer);

  //////--------


  
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
