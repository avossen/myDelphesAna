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
using namespace std;


enum RecType{elec, hadronic, da, mixed, mostlyLepton, truth,gen, genNoAcc,recTypeEnd};
//  int colors[]={kRed,kBlue,kGreen,kBlack,kCyan,kMagenta,kOrange,kYellow};
string recTypeNames[]={"elec","hadronic","da","mixed","mostlyLepton","trueBoost","generatedWAcc","generatedWOAcc"};
class ExRootTreeReader;
class ExRootResult;
const int corrBinsX=20;
const int corrBinsQ2=10;
struct TestPlots
{
  TH2* angCorr;
  TH1* angRes;
  TH1* angTrue;
  TH1* angRec;
  TH1* yRes;
  TH2* yCorr;
  TH1* yRec;
  TH1* yReal;
  TH2* Q2VsXSmear;
  TH2* Q2VsXSmearDA;
  TH2* Q2VsXSmearJB;
  TH2* Q2VsXSmearMixed;
  TH2* Q2VsXSmearNorm;
  TH2* corrConventionalQ2;
  TH2* corrConventionalX;


  TH2* Q2CorrElec;
  TH2* Q2CorrJB;
  TH2* Q2CorrDA;
  TH2* XCorrElec;
  TH2* XCorrJB;
  TH2* XCorrDA;
  TH2* XCorrMixed;

  TH2* corrJBQ2;
  TH2* corrJBQ2_2;
  TH2* corrJBX;

  TH2* corrDAQ2;
  TH2* corrDAX;

  TH2* corrMixedQ2;
  TH2* corrMixedX;


  TH2* zMeans[8];
  TH2* zDiffs[8];
  
  TH2* ptMeans[8];
  //only need counter for one hadron pair related variable...
  TH2* ptCounts[8];
  TH2* ptDiffs[8];

  TH2* phiRMeans[8];
  TH2* phiRDiffs[8];

  TH1*** phiRHistos[8];
  
  vector<float> zBins;
  vector<float> yBins;
};

void BookHistograms(ExRootResult *result, TestPlots *plots,int beamEnergyI, int  hadronBeamEnergyI)
{
  TLegend *legend;
  TPaveText *comment;
  char buffer[300];
  sprintf(buffer,"Fraction of events staying in bin (%dx%d)",beamEnergyI,hadronBeamEnergyI);
  plots->yRes=result->AddHist1D("yRes","yRes","#Delta y","counts",200,-2.0,2.0);
  plots->yReal=result->AddHist1D("yReal","yReal","corr y","counts",200,-2.0,2.0);
  plots->yRec=result->AddHist1D("yRec","yRec","Rec y","counts",500,0.0,50.0);
  plots->yCorr=result->AddHist2D("yCorr","yCorr","yRec","yReal",50,0,1,50,0,1,0,0);
  plots->angCorr=result->AddHist2D("angCorr","angCorr","angRec","angReal",50,0,6,50,0,6,0,0);
  plots->angRes=result->AddHist1D("angRes","angRes","#Delta ang","counts",200,-2.0,7.0);
  plots->angTrue=result->AddHist1D("angTrue","angTrue"," ang true","counts",200,0,6.5);
  plots->angRec=result->AddHist1D("angRec","angRec"," ang rec","counts",200,0,6.5);
  plots->Q2VsXSmear=result->AddHist2D("Q2VsXSmear",buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);
  plots->Q2VsXSmearDA=result->AddHist2D("Q2VsXSmearDA",buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);
  plots->Q2VsXSmearJB=result->AddHist2D("Q2VsXSmearJB",buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);
  plots->Q2VsXSmearMixed=result->AddHist2D("Q2VsXSmearMixed",buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);
  plots->Q2VsXSmearNorm=result->AddHist2D("Q2VsXSmearNorm",buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);

  for(int i=0;i<recTypeEnd;i++)
    {
      plots->phiRHistos[i]=new TH1**[plots->zBins.size()];
      for(int j=0;j<plots->zBins.size();j++)
	{
	  plots->phiRHistos[i][j]=new TH1*[plots->yBins.size()];
	  for(int k=0;k<plots->yBins.size();k++)
	    {
	      sprintf(buffer,"phi_R_%s_zBin_%d_yBin_%d",recTypeNames[i].c_str(),j,k);
	      plots->phiRHistos[i][j][k]=result->AddHist1D(buffer,buffer,"phi_R","counts",30,0,2*TMath::Pi());
	    }
	}
    }
  
  for(int i=0;i<recTypeEnd;i++)
    {
      sprintf(buffer,"pt_mean_%s",recTypeNames[i].c_str());
      plots->ptMeans[i]=result->AddHist2D(buffer,buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);
      sprintf(buffer,"pt_counts_%s",recTypeNames[i].c_str());
       plots->ptCounts[i]=result->AddHist2D(buffer,buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);
      sprintf(buffer,"pt_diffs_%s",recTypeNames[i].c_str());
       plots->ptDiffs[i]=result->AddHist2D(buffer,buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);


      sprintf(buffer,"z_mean_%s",recTypeNames[i].c_str());
      plots->zMeans[i]=result->AddHist2D(buffer,buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);
      sprintf(buffer,"z_diffs_%s",recTypeNames[i].c_str());
       plots->zDiffs[i]=result->AddHist2D(buffer,buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);

       
       sprintf(buffer,"phiR_mean_%s",recTypeNames[i].c_str());
       plots->phiRMeans[i]=result->AddHist2D(buffer,buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);

       sprintf(buffer,"phiR_diffs_%s",recTypeNames[i].c_str());
       plots->phiRDiffs[i]=result->AddHist2D(buffer,buffer,"x","Q^{2}",corrBinsX,0.0001,1,corrBinsQ2,0.5,10000,1,1);
    }

  cout <<"done with my histoos " << endl;

  plots->Q2CorrElec=result->AddHist2D("q2corrElec","q2corrElec","Q2 true", "Q2 rec",100,1,200,100,1,200,1,1);
  plots->Q2CorrJB=result->AddHist2D("q2corrJB","q2corrJB","Q2 true", "Q2 rec",100,1,200,100,1,200,1,1);
  plots->Q2CorrDA=result->AddHist2D("q2corrDA","q2corrDA","Q2 true", "Q2 rec",100,1,200,100,1,200,1,1);

  plots->XCorrElec=result->AddHist2D("XcorrElec","XcorrElec","x true", "x false",100,0.0001,1.0,200,0.0001,1.0,1,1);
  plots->XCorrJB=result->AddHist2D("XcorrJB","XcorrJB","x true", "x false",100,0.0001,1.0,200,0.0001,1.0,1,1);
  plots->XCorrDA=result->AddHist2D("XcorrDA","XcorrDA","x true", "x false",100,0.0001,1.0,200,0.0001,1.0,1,1);
  plots->XCorrMixed=result->AddHist2D("XcorrMixed","XcorrMixed","x true", "x false",100,0.0001,1.0,200,0.0001,1.0,1,1);

  plots->corrConventionalQ2=result->AddHist2D("convQ2","convQ2","Q2 true", "Q2 rec",100 ,1,200,100 ,1,200,1,1);
  plots->corrConventionalX=result->AddHist2D("convX","convX","x true", "x rec", 100 ,0.0001,1.0,100 ,-1.0,1.0,1,1);

  plots->corrJBQ2=result->AddHist2D("jbQ2","jbQ2","Q2 true", " Q2 rec", 100 ,1,200,100 ,-1,1,1,1);
  plots->corrJBQ2_2=result->AddHist2D("jbQ2_2","jbQ2_2","Q2 true", "Q2 rec", 100 ,1,200,100 ,1,200,1,1);
  plots->corrJBX=result->AddHist2D("jbX","jbX","x true", "x rec",100 ,0.0001,1.0,100 ,-1.0,1.0,1,1);

  plots->corrDAQ2=result->AddHist2D("daQ2","daQ2","Q2 true", "Q2 rec", 100 ,1,200,100 ,-1,1,1,1);
  plots->corrDAX=result->AddHist2D("daX","daX","x true", "x rec",100 ,0.0001,1.0,100 ,-1.0,1.0,1,1);

  plots->corrMixedQ2=result->AddHist2D("mixedQ2","mixedQ2","Q2 true", "Q2 rec", 100 ,1,200,100 ,-1,1,1,1);
  plots->corrMixedX=result->AddHist2D("mixedX","mixedX","x true"," x rec",100 ,0.0001,1.0,100 ,-1.0,1.0,1,1);


  //  BinLog(plots->Q2VsXSmearNormEt.GetXaxis());
  BinLog(plots->Q2VsXSmear->GetXaxis());
  BinLog(plots->Q2VsXSmearDA->GetXaxis());
  BinLog(plots->Q2VsXSmearJB->GetXaxis());
  BinLog(plots->Q2VsXSmearMixed->GetXaxis());
  BinLog(plots->Q2VsXSmearNorm->GetXaxis());


  BinLog(plots->Q2VsXSmearMixed->GetYaxis());
  BinLog(plots->Q2VsXSmear->GetYaxis());
  BinLog(plots->Q2VsXSmearDA->GetYaxis());
  BinLog(plots->Q2VsXSmearJB->GetYaxis());
  BinLog(plots->Q2VsXSmearNorm->GetYaxis());
  //  BinLog(plots->Q2VsXSmearEt->GetYaxis());
  //  BinLog(plots->Q2VsXSmearNormEt->GetYaxis());


  BinLog(plots->Q2CorrElec->GetXaxis());
  BinLog(plots->Q2CorrElec->GetYaxis());

  BinLog(plots->Q2CorrJB->GetXaxis());
  BinLog(plots->Q2CorrJB->GetYaxis());

  BinLog(plots->Q2CorrDA->GetXaxis());
  BinLog(plots->Q2CorrDA->GetYaxis());

  BinLog(plots->XCorrElec->GetXaxis());
  BinLog(plots->XCorrElec->GetYaxis());

  BinLog(plots->XCorrJB->GetXaxis());
  BinLog(plots->XCorrJB->GetYaxis());

  BinLog(plots->XCorrDA->GetXaxis());
  BinLog(plots->XCorrDA->GetYaxis());
  BinLog(plots->XCorrMixed->GetXaxis());
  BinLog(plots->XCorrMixed->GetYaxis());


  BinLog(plots->corrConventionalQ2->GetXaxis());
  BinLog(plots->corrConventionalX->GetXaxis());
  BinLog(plots->corrJBQ2_2->GetXaxis());
  BinLog(plots->corrJBQ2_2->GetYaxis());
  BinLog(plots->corrJBQ2->GetXaxis());
  BinLog(plots->corrJBX->GetXaxis());

  BinLog(plots->corrDAQ2->GetXaxis());
  BinLog(plots->corrDAX->GetXaxis());

  BinLog(plots->corrMixedX->GetXaxis());

  cout <<"setting logs " <<endl;

  cout <<endl;
    for(int i=0;i<recTypeEnd;i++)
    {
      cout <<" doing i " << i <<endl;
      BinLog(plots->ptMeans[i]->GetXaxis());
      BinLog(plots->ptMeans[i]->GetYaxis());
      
      BinLog(plots->ptDiffs[i]->GetXaxis());
      BinLog(plots->ptDiffs[i]->GetYaxis());
      
      BinLog(plots->ptCounts[i]->GetXaxis());
      BinLog(plots->ptCounts[i]->GetYaxis());
      cout <<" done pt " << endl;
      BinLog(plots->zMeans[i]->GetXaxis());
      BinLog(plots->zMeans[i]->GetYaxis());
      BinLog(plots->zDiffs[i]->GetXaxis());
      BinLog(plots->zDiffs[i]->GetYaxis());
      cout <<" done z " << endl;
      BinLog(plots->phiRMeans[i]->GetXaxis());
      BinLog(plots->phiRMeans[i]->GetYaxis());
      
      BinLog(plots->phiRDiffs[i]->GetXaxis());
      BinLog(plots->phiRDiffs[i]->GetYaxis());
      cout <<" done with this iteration " <<endl;
    }
  cout <<"done .." <<endl;
  
  //  BinLog(plots->corrConventionalQ2SmallY->GetXaxis());
  //  BinLog(plots->corrConventionalXSmallY->GetXaxis());

  //  BinLog(plots->corrJBQ2SmallY->GetXaxis());
  //  BinLog(plots->corrJBXSmallY->GetXaxis());

  //  BinLog(plots->corrDAQ2SmallY->GetXaxis());
  //  BinLog(plots->corrDAXSmallY->GetXaxis());

  //  BinLog(plots->corrMixedXSmallY->GetXaxis());


  //  BinLog(plots->plots->Q2VsXSmear->GetXaxis());
  //  BinLog(plots->Q2VsXSmear->GetYaxis());
				      
  //
  //  plots->fElectronDeltaPT = result->AddHist1D(
  //    "electron_delta_pt", "(p_{T}^{particle} - p_{T}^{electron})/p_{T}^{particle}",
  //    "(p_{T}^{particle} - p_{T}^{electron})/p_{T}^{particle}", "number of electrons",
  //    100, -0.1, 0.1);
  //
  //  plots->fElectronDeltaEta = result->AddHist1D(
  //    "electron_delta_eta", "(#eta^{particle} - #eta^{electron})/#eta^{particle}",
  //    "(#eta^{particle} - #eta^{electron})/#eta^{particle}", "number of electrons",
  //    100, -0.1, 0.1);
  //
  //  plots->fJetDeltaPT = result->AddHist1D(
  //    "jet_delta_pt", "(p_{T}^{jet} - p_{T}^{constituents})/p_{T}^{jet}",
  //    "(p_{T}^{jet} - p_{T}^{constituents})/p_{T}^{jet}", "number of jets",
  //100, -1.0e-1, 1.0e-1);
  //
  //   plots->res_dphi_z_10_20 = result->AddHist2D("dphi_res_z_10_20", "Collins Angle resolution, 10 <E_{jet}< 20\
  // GeV", "#phi_{C}^{gen}-#phi_{C}^{reco} [rad]", "generated z", 51, -mindphi,mindphi, 10,0.0,1.0);
  //
}
void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
}
void AnalyzeEvents(ExRootTreeReader *treeReader, TestPlots *plots, double beamEnergy, double hadronBeamEnergy,double sqrtS)
{
  cout <<"in ana events " <<endl;
    bool useTruth=false;
  //       bool useTruth=true;
  float logQ2Range=4.3; //for 0.5 to 10k
  float logXRange=4;

  TClonesArray *branchEvent=treeReader->UseBranch("Event");
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
  TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");

  cout <<"after use branch.." <<endl;
  Long64_t allEntries = treeReader->GetEntries();

  GenParticle *particle;
  Electron *electron;

  Track *track;
  Tower *tower;

  Jet *jet;
  Jet *genjet;
  Jet *matchedjet;
  TObject *object;

  TLorentzVector momentum;

  Float_t Eem, Ehad;
  Bool_t skip;
  Long64_t entry;

  Int_t i, j, pdgCode;
  //      cout <<"original e px: " << ePxOrig <<" py: " << ePyOrig<< " pz: "<< ePzOrig<<endl;


  for(entry=0;entry < allEntries;++entry)
    {
      if(entry%5000==0)
	{
	  cout <<"looking at event " << entry <<" of " << allEntries <<" : " <<100*entry/allEntries <<" % " <<endl;
	}

      if(entry/(double)allEntries>0.1)
	{
	  //	  break;
	}
      treeReader->ReadEntry(entry);

      //	     # four-momenta of proton, electron, virtual photon/Z^0/W^+-.                                                                       //this is the truth                                         
      //      cout <<"going through entry " << entry <<endl;
      //      cout <<" is there a pointer here? " << branchParticle->At(0) <<endl;

      //      cout <<"branch event has " << branchEvent->GetEntries() <<endl;

    
      if(branchParticle->GetEntries()>10)
 	{	
 	  for(int i=0;i<10;i++)
	    {
	      GenParticle* g=((GenParticle*)branchParticle->At(i));
	      //	      cout <<"particle at index " <<i   <<" pid: " << g->PID << " energy: "<< g->E <<endl;
	    }
	}
    
      //      TLorentzVector  pProton    =  ((GenParticle*)branchParticle->At(0))->P4();// #these numbers 0 , 3, 5 are hardcoded in Pythia8
      TLorentzVector  pProton    =  ((GenParticle*)branchParticle->At(5))->P4();// #these numbers 0 , 3, 5 are hardcoded in Pythia8
	     
      //      TLorentzVector	  pleptonIn    = ((GenParticle*)branchParticle->At(3))->P4();
      //seems to be different in dire as well...
      TLorentzVector	  pleptonIn    = ((GenParticle*)branchParticle->At(2))->P4();
      //      cout <<"proton... has energy: " << pProton.E() <<" lepton in: " << pleptonIn.E() <<endl;
      //      cout <<"lepton..." <<endl;
      //      TLorentzVector	  pleptonOut   = ((GenParticle*)branchParticle->At(5))->P4();
      //in the dire output this is at 4 (probably because coordinate system was switched originally)
      TLorentzVector	  pleptonOut   = ((GenParticle*)branchParticle->At(4))->P4();
      //      cout <<"lep out..." <<endl;
      TLorentzVector	  pPhoton      = pleptonIn - pleptonOut;
      //      cout <<"photon..." <<endl;
      //Q2, W2, Bjorken x, y, nu.                                                                                                                                                  
      double  Q2 = -pPhoton.M2();
      double  W2 = (pProton + pPhoton).M2();
      double  x = Q2 / (2. * pProton.Dot(pPhoton));
      double  y = (pProton.Dot(pPhoton)) / (pProton.Dot(pleptonIn));
      //       cout <<"true q2: " << Q2 <<" x: " << x <<" y: " << y <<endl;
    
      TLorentzVector e;
      Electron* electron;
      if(branchElectron->GetEntries()>0)
 	{
 	  electron = ((Electron*)branchElectron->At(0));
 	  e =electron->P4();
 	}
      else
 	{
 	  continue;
 	}
      //       cout <<"study rec .." <<endl;
      studyReconstruction(branchParticle,branchEFlowTrack,branchEFlowPhoton,branchEFlowNeutralHadron,branchTrack);

      //       cout <<"rec electron E: "<< e.E()<<endl;
      //       cout << "true electron 1: "<< pleptonOut.E() <<endl;
      particle = (GenParticle*) electron->Particle.GetObject();
      //       cout <<"truth from rec electron: "<< particle->E<<endl;

       
      //      Kins kinOrig=getKinsFromScatElectron(pleptonIn.E(),pProton.E(),pleptonOut.Px(),pleptonOut.Py(),pleptonOut.Pz(),pleptonOut.E());
      Kins kinOrig=getKinsFromScatElectron(beamEnergy,hadronBeamEnergy,particle->Px,particle->Py,particle->Pz,particle->E);

      Kins kinSmeared=getKinsFromScatElectron(beamEnergy,hadronBeamEnergy,e.Px(),e.Py(),e.Pz(),e.E());
      HadronicVars hvOrig=getOriginalHadronicVars(branchParticle);
      HadronicVars hvSmeared=getHadronicVars(branchElectron, branchEFlowTrack, branchEFlowPhoton,branchEFlowNeutralHadron,branchTrack,branchParticle);
      //
      //	    HadronicVars hvSmeared=getHadronicVars(branchElectron, branchTrack, branchEFlowPhoton,branchEFlowNeutralHadron,branchTrack,branchParticle);

      plots->yRes->Fill((hvOrig.sumEMinusPz-hvSmeared.sumEMinusPz)/(2*beamEnergy));
      //      plots->yRec->Fil1l((hvSmeared.sumEMinusPz)/(2*beamEnergy));
      //      plots->yRec->Fill((hvSmeared.sumEMinusPz));
      plots->yReal->Fill((hvOrig.sumEMinusPz)/(2*beamEnergy));
      plots->yCorr->Fill(hvSmeared.sumEMinusPz/(2*beamEnergy),hvOrig.sumEMinusPz/(2*beamEnergy));

      plots->angRec->Fill(hvSmeared.theta);
      plots->angTrue->Fill(hvOrig.theta);
      plots->angRes->Fill(hvSmeared.theta-hvOrig.theta);
      plots->angCorr->Fill(hvSmeared.theta, hvOrig.theta);



      
      if(useTruth)
	{
	  hvSmeared=hvOrig;
	}
      
      Kins kinJBSmeared=getKinsJB(hvSmeared,beamEnergy,hadronBeamEnergy,sqrtS*sqrtS);
      //Kins kinDASmeared=getKinsDA(e.Px(),e.Py(),e.Pz() ,e.E(),beamEnergy,hvSmeared.theta, sqrtS*sqrtS);
      Kins kinDASmeared=getKinsDA(particle->Px,particle->Py,particle->Pz ,particle->E,beamEnergy,hvSmeared.theta, sqrtS*sqrtS);
      //already have the correct hvSmeared, so just need to set the electron kinematics
      if(useTruth)
	{
	  kinDASmeared=getKinsDA(particle->Px,particle->Py,particle->Pz ,particle->E,beamEnergy,hvSmeared.theta, sqrtS*sqrtS);
	}

      //      cout <<"jb y: "<< kinJBSmeared.y << " sqrtS: "<< sqrtS <<endl;
      float mixedXSmeared=kinSmeared.Q2/(sqrtS*sqrtS*kinJBSmeared.y);
      
      double binlogOrigQ2=log10(kinOrig.Q2);
      double binlogOrigX=log10(kinOrig.x);
      double binRecQ2=log10(kinSmeared.Q2);
      double binRecX=log10(kinSmeared.x);


      TLorentzVector qOrig;
      qOrig.SetPxPyPzE(particle->Px,particle->Py,particle->Pz,particle->E);
      /////////
      TLorentzVector qH;

      //diff is sum since they should be in opposite directions
      //      cout <<"hv smeard pz: "<< hvSmeared.sumPz <<" lepton pz: "<< e.Pz() <<" diff: " << hvSmeared.sumPz+e.Pz()<<endl;
      qH.SetPxPyPzE(-hvSmeared.sumPx,-hvSmeared.sumPy,-hvSmeared.sumPz,hvSmeared.sumE);
      //l direction, h size
      TLorentzVector qLH;
      //
      TVector3 eDir(e.Px(),e.Py(),e.Pz());
      //      cout <<" edir: " << printVect(eDir);
      eDir.SetMag(qH.Vect().Mag());
      //      cout <<" after set mag: " << printVect(eDir);
      qLH.SetVectMag(eDir,sqrt(kinSmeared.Q2));
      //      cout <<"resulting qlh: " << printLVect(qLH);
      //only lepton
      TLorentzVector qL;
      qL.SetPxPyPzE(e.Px(),e.Py(),e.Pz(),e.E());
      //      cout <<" l dir: "<< printLVect(qL);
      //now we have the outgoing lepton direction, q is in-out

      //rec outgoing electron is 'e', true is plepton out
      ///kinOrig.x
      //calculate initial pz of parton
      //      cout <<" hv orig sum: "<< hvOrig.sumPz <<" hv smeared: "<< hvSmeared.sumPz <<" diff: "<< hvOrig.sumPz - hvSmeared.sumPz <<endl;
      //      cout <<" h sum pz: "<< hvSmeared.sumPz << " l pz: " << pleptonOut.Pz() <<endl;
      //      cout <<"initial pz sum: " << kinOrig.x*hadronBeamEnergy-beamEnergy <<" and then " << hvSmeared.sumPz+pleptonOut.Pz()  <<" diff: "<< kinOrig.x*hadronBeamEnergy-beamEnergy -( hvSmeared.sumPz+pleptonOut.Pz()) <<endl;
      //            cout <<"same w/o  x:initial pz sum: " << hadronBeamEnergy-beamEnergy <<" and then " << hvSmeared.sumPz+pleptonOut.Pz()  <<" diff: "<< hadronBeamEnergy-beamEnergy -( hvSmeared.sumPz+pleptonOut.Pz()) <<endl;
      
      
      qH=pleptonIn-qH;
      qL=pleptonIn-qL;
      qLH=pleptonIn-qLH;
      qOrig=pleptonIn-qOrig;


      //////////////////////////
      //-->same for all
      double kappa=pProton.Pz()*pProton.Pz()/(pProton.E()*pProton.E())-1;
      double a=pProton.Dot(pleptonIn);
      double Pz=pProton.Pz();
      //      double Q2=kinOrig.Q2;
      double Pe=pProton.E();
      /////////
      
      double hPx=hvSmeared.sumPx;
      double hPy=hvSmeared.sumPy;
      
      //calculate qTx: (lx-l'x), qTy: (y-l'x). l'x=-hx, l'y=hy
      ///
      double qTx=pleptonIn.Px()+hPx;
      double qTy=pleptonIn.Py()+hPy;
      

      double b=pProton.Px()*qTx+pProton.Py()*qTy;
      double qT2=qTx*qTx+qTy*qTy;
      //      double y=

      ///here there would be differences depending on the y, Q2 used

      double usedQ2[8];
      double usedy[8];
      for(int i=0;i<recTypeEnd;i++)
	{
	  switch(i)
	    {
	      //	      enum RecType{elec, hadronic, da, mixed, mostlyLepton, truth,recTypeEnd};
	    case elec:
	      usedQ2[i]=kinSmeared.Q2;
	      usedy[i]=kinSmeared.y;
	      break;

	    case hadronic:
	      usedQ2[i]=kinJBSmeared.Q2;
	      usedy[i]=kinJBSmeared.y;
	      break;

	    case da:
	      usedQ2[i]=kinDASmeared.Q2;
	      usedy[i]=kinDASmeared.y;
	      break;

	    case mixed:
	      usedQ2[i]=kinSmeared.Q2;
	      usedy[i]=kinJBSmeared.y;
	      break;

	    case mostlyLepton:
	      usedQ2[i]=kinSmeared.Q2;
	      usedy[i]=kinSmeared.y;
	      break;

	    case truth:
	      usedQ2[i]=Q2;
	      usedy[i]=y;
	      break;
	    case gen:
	      usedQ2[i]=Q2;
	      usedy[i]=y;
	      break;
	    case genNoAcc:
	      usedQ2[i]=Q2;
	      usedy[i]=y;
	      break;
	      
	    default: 
	      cout <<"something wrong"<<endl;
	    }
	  
	}
      pair<float,float> res[8];
      double qE[8];
      TLorentzVector q[8];
      TLorentzVector q_P[8];
      TVector3 q_bv[8];
      double W[8];

      TLorentzVector qElec=pleptonIn-e;
      
      for(int i=0;i<recTypeEnd;i++)
	{
	  res[i]=getQz(a,b,kappa,usedy[i],qT2,Pz,Pe,usedQ2[i]);
	  double qz=res[i].second;
	  if(fabs(qElec.Pz()-res[i].second)>fabs(qElec.Pz()-res[i].first))
	    {
	      qz=res[i].first;
	    }
	  qE[i]=getQe(a,b,usedy[i],qT2,Pz,qz,Pe);
	  //	         Kins kinSmeared=getKinsFromScatElectron(beamEnergy,hadronBeamEnergy,e.Px(),e.Py(),e.Pz(),e.E());
	  if(i==elec)
	    {
	      q[i]=qElec;
	    }
	  else
	    {
	      
	      if(i==truth || i==gen || i==genNoAcc)
		{
		  q[i]=pPhoton;
		}
	      else
		{
		  q[i]=TLorentzVector(qTx,qTy,qz,qE[i]);
		}
	    }
	  q_P[i]=q[i]+pProton;
	  q_bv[i]=q_P[i].BoostVector();
	  q_bv[i]=(-1)*q_bv[i];
	  W[i]=q_P[i].M();
	}

      //      cout <<"orig q: "<< printLVect(qOrig) <<endl;
        for(int i=0;i<recTypeEnd;i++)
	{
	  //	  cout <<recTypeNames[i]<<": " << printLVect(q[i]) <<endl;
	}

      
      
      //      pair<float,float> resJB=getQz(a,b,kappa,kinJBSmeared.y,qT2,Pz,Pe,kinJBSmeared.Q2);
      //      pair<float,float> resDA=getQz(a,b,kappa,kinDASmeared.y,qT2,Pz,Pe,kinDASmeared.Q2);
      //      pair<float,float> resMixed=getQz(a,b,kappa,kinJBSmeared.y,qT2,Pz,Pe,kinSmeared.Q2);mix
      //      pair<float,float> resMostlyLepton=getQz(a,b,kappa,kinJBSmeared.y,qT2,Pz,Pe,kinJBSmeared.Q2);



      
      //      double qz=res.second;
      //      double qE=getQe(a,b,y,qT2,Pz,qz,Pe);



      //----> use the other qs here



      
      //      cout <<"qh: "<< printLVect(qH) <<endl;
      //      cout <<"ql: "<< printLVect(qL) <<endl;
      //      cout <<"qlh: "<< printLVect(qLH) <<endl;
      //      cout <<"qOrig: "<< printLVect(qOrig) <<endl;
      
      //calculate the boost to the p-q system

      
      //      TLorentzVector qH_P=qH+pProton;
      //      TLorentzVector qL_P=qL+pProton;
      //      TLorentzVector qLH_P=qLH+pProton;
      TLorentzVector qOrig_P=qOrig+pProton;

      //    float WH=qH_P.M();
      //      float WL=qL_P.M();
      //      float WLH=qLH_P.M();
      float WOrig=qOrig_P.M();
      
      //      TVector3 qH_bv=qH_P.BoostVector();
      //      TVector3 qL_bv=qL_P.BoostVector();
      //      TVector3 qLH_bv=qLH_P.BoostVector();
      TVector3 orig_bv=qOrig_P.BoostVector();

      //      qH_bv=(-1)*qH_bv;
      //      qL_bv=(-1)*qL_bv;
      //      qLH_bv=(-1)*qLH_bv;
      orig_bv=(-1)*orig_bv;
      ///      -->try to get qe not from y but from qz directly (via scattered electron)
      ////////////////////////
      vector<HadronPair> pairsRec[8];
      vector<HadronPair> pairsTrue[8];
      boostVars recBoost;
      boostVars realBoost;

      //recBoost.breitBoost=orig_bv;
      //      recBoost.W=WLH;
      //      recBoost.lv_q=qLH;
      //this is the beam
      //      TLorentzVector lv_beam;
      //      lv_beam.SetPxPyPzE(0.0,0.0,beamEnergy,
      //      recBoost.lv_l=pleptonIn;
      realBoost.breitBoost=orig_bv;
      realBoost.W=WOrig;
      realBoost.lv_q=qOrig;
      realBoost.lv_l=pleptonIn;


      ///get real x, q2 bin
      int xBin=0;
      int q2Bin=0;

      //      cout <<"getting axis " <<endl;
      TAxis *xaxis = plots->Q2VsXSmearNorm->GetXaxis();
      TAxis *yaxis = plots->Q2VsXSmearNorm->GetYaxis();
      xBin = xaxis->FindBin(kinOrig.x);
      q2Bin = yaxis->FindBin(kinOrig.Q2);
      //      cout <<"using xbin: "<< xBin <<" q2Bin: "<< q2Bin << " for x: "<< kinOrig.x <<" q2: " << kinOrig.Q2 <<endl;
      //      cout <<"done " <<endl;


      
      for(int i=0;i<recTypeEnd;i++)
	{
	  recBoost.breitBoost=q_bv[i];
	  recBoost.W=W[i];
	  recBoost.lv_q=q[i];
	  recBoost.lv_l=pleptonIn;
	  switch(i)
	    {
	    case truth:
	      //	      cout <<"truth recboost: "<< printLVect(recBoost.lv_q) <<" realBoost: "<< printLVect(realBoost.lv_q) <<endl;
	      //	      cout <<"rec W: "<< recBoost.W <<" real: "<< realBoost.W <<endl;
	      //	      cout <<"rec bosot : "<< printVect(recBoost.breitBoost) <<" real: "<< printVect(realBoost.breitBoost) <<endl;
	      getHadronPairs(pairsRec[i],pairsTrue[i],branchEFlowTrack,recBoost,realBoost,Pz,true,false,false);
	      break;
	    case gen:
	      getHadronPairs(pairsRec[i],pairsTrue[i],branchParticle,recBoost,realBoost,Pz,true,true,false);
	      break;

	    case genNoAcc:
	      getHadronPairs(pairsRec[i],pairsTrue[i],branchParticle,recBoost,realBoost,Pz,true,true,true);
	      break;
	    default:
	    getHadronPairs(pairsRec[i],pairsTrue[i],branchEFlowTrack,recBoost,realBoost,Pz,false,false,false);
	    }
      ////-----
	  ///	  	  cout <<"looking at rec pairs " <<endl;
	  ////	  	  cout <<" running over " << pairsRec[i].size() << " pairs for comb " << recTypeNames[i] <<endl;

	  //	  cout <<" num pairs " << pairsRec[i].size() <<endl;
	  for(int j=0;j<pairsRec[i].size();j++)
	    {
	      
	      if(!checkSanity(pairsTrue[i][j]))
		{
		  continue;
		}
	      if( !checkSanity(pairsRec[i][j]))
		{
		  continue;
		}
	      
	      int zBin=getBin(plots->zBins,pairsTrue[i][j].z);
	      int yBin=getBin(plots->yBins,y);
	      if(zBin < 0 || yBin < 0)
		{
		  	  cout <<"wrong z, y bin " <<endl;
		  continue;
		}
	      //	      cout <<"zBin: "<< zBin << " yBin: "<< yBin <<" z: "<< pairsTrue[i][j].z<<" y: "<< y<<endl;


	      if(pairsRec[i][j].phi_R < 0 ||  pairsRec[i][j].phi_R> 2*TMath::Pi())
		{
		//		cout <<" phi not in range " <<endl;
		}
	      //		cout <<"filling with index " << i <<endl;
		
	      if(i==truth)
		{
		  plots->phiRHistos[i][zBin][yBin]->Fill(pairsTrue[i][j].phi_R);
		}
	      else
		{
		  plots->phiRHistos[i][zBin][yBin]->Fill(pairsRec[i][j].phi_R);
		}
	      
	      ///	      	      cout <<"rec phiR:" << pairsRec[i][j].phi_R <<" real: "<< pairsTrue[i][j].phi_R <<" % diff: "<< (pairsRec[i][j].phi_R-pairsTrue[i][j].phi_R)/pairsTrue[i][j].phi_R<<endl;
	      plots->ptMeans[i]->SetBinContent(xBin,q2Bin,plots->ptMeans[i]->GetBinContent(xBin,q2Bin)+pairsTrue[i][j].pT);
	      //	      cout <<"filling with pt true: "<< pairsTrue[i][j].pT <<" rec: " << pairsRec[i][j].pT <<endl;
	      plots->ptCounts[i]->SetBinContent(xBin,q2Bin,plots->ptCounts[i]->GetBinContent(xBin,q2Bin)+1);
	      plots->ptDiffs[i]->SetBinContent(xBin,q2Bin,plots->ptDiffs[i]->GetBinContent(xBin,q2Bin)+fabs(pairsRec[i][j].pT-pairsTrue[i][j].pT)/fabs(pairsTrue[i][j].pT));
	      //	      	      cout <<"filling with z1: "<< pairsTrue[i][j].z1 <<" z2: "<< pairsTrue[i][j].z2 <<endl;
	      //      	      cout <<"filling with rec z1: "<< pairsRec[i][j].z1 <<" z2: "<< pairsRec[i][j].z2 <<endl;
	      plots->zMeans[i]->SetBinContent(xBin,q2Bin,plots->zMeans[i]->GetBinContent(xBin,q2Bin)+pairsTrue[i][j].z1);
	      plots->zMeans[i]->SetBinContent(xBin,q2Bin,plots->zMeans[i]->GetBinContent(xBin,q2Bin)+pairsTrue[i][j].z2);
	      plots->zDiffs[i]->SetBinContent(xBin,q2Bin,plots->zDiffs[i]->GetBinContent(xBin,q2Bin)+fabs(pairsRec[i][j].z1-pairsTrue[i][j].z1)/fabs(pairsTrue[i][j].z1));
	      plots->zDiffs[i]->SetBinContent(xBin,q2Bin,plots->zDiffs[i]->GetBinContent(xBin,q2Bin)+fabs(pairsRec[i][j].z2-pairsTrue[i][j].z2)/fabs(pairsTrue[i][j].z2));
	      plots->phiRMeans[i]->SetBinContent(xBin,q2Bin,plots->phiRMeans[i]->GetBinContent(xBin,q2Bin)+pairsTrue[i][j].phi_R);
	      //normalizing phiR doesn't make sense, since we get a lot of larg numbers then for small angles...
	      if(i==truth)
		{
		  //		  cout <<" rec phir: "<< pairsRec[i][j].phi_R <<" real: "<< pairsTrue[i][j].phi_R <<endl;
		}
	      plots->phiRDiffs[i]->SetBinContent(xBin,q2Bin,plots->phiRDiffs[i]->GetBinContent(xBin,q2Bin)+fabs(pairsRec[i][j].phi_R-pairsTrue[i][j].phi_R));
	    }
	  //	  cout <<" done looking at rec pairs.. " << endl;
	}
      
      ////----
      //to get q one can get it only from hadron or use the lepton direction and the Q2 for the length
      ////----
      
      //      cout <<" qH: "<< printLVect(qH);
      //      cout <<" qLH: ";//<< printLVect(qLH);
      //      cout <<" qL: ";//<< printLVect(qL);

      ////      gN.add(lv_target);
      // System.out.println("gN x " + gN.px()+ " y: "+gN.py() + " pz: " +gN.pz() + "
      // e: "+ gN.e() );
      ///		Walt = gN.mass();
      //System.out.println("W: " + W  + " Walt: "+ Walt);
		
      //		gNBoost = gN.boostVector();
      //		gNBoost.negative();
		

      //      cout <<" reconstructed orig x " << kinOrig.x <<" Q2: "<< kinOrig.Q2 <<endl;
      //      cout <<" reconstructed smear x " << kinSmeared.x <<" Q2: "<< kinSmeared.Q2 <<endl;
      //      cout <<" reconstructed JB x " << kinJBSmeared.x <<" Q2: "<< kinJBSmeared.Q2 <<endl;
      //      cout <<" reconstructed DA x " << kinDASmeared.x <<" Q2: "<< kinDASmeared.Q2 <<endl;
      //      cout <<" reconstructed mixed x " << mixedXSmeared <<endl;
      plots->Q2VsXSmearNorm->Fill(kinOrig.x,kinOrig.Q2);
      
      if(fabs(binlogOrigQ2-binRecQ2)<logQ2Range/(2*corrBinsQ2) && fabs(binlogOrigX-binRecX)<logXRange/(2*corrBinsX))
	{
	  plots->Q2VsXSmear->Fill(kinOrig.x,kinOrig.Q2);
	}
      binRecQ2=log10(kinDASmeared.Q2);
      binRecX=log10(kinDASmeared.x);
      
      if(fabs(binlogOrigQ2-binRecQ2)<logQ2Range/(2*corrBinsQ2) && fabs(binlogOrigX-binRecX)<logXRange/(2*corrBinsX))
	{
	  plots->Q2VsXSmearDA->Fill(kinOrig.x,kinOrig.Q2);
	}
      
      binRecQ2=log10(kinJBSmeared.Q2);
      binRecX=log10(kinJBSmeared.x);
      //      cout <<"jb q2: "<< kinJBSmeared.Q2 <<" x: "<< kinJBSmeared.x <<" log q2: "<< binRecQ2 <<" log x: "<< binRecX << " binlogOrig: " << binlogOrigQ2 <<" x: "<< binlogOrigX <<endl;
      if(fabs(binlogOrigQ2-binRecQ2)<logQ2Range/(2*corrBinsQ2) && fabs(binlogOrigX-binRecX)<logXRange/(2*corrBinsX))
	{
	  //	  cout <<" in bin " <<endl;
	  plots->Q2VsXSmearJB->Fill(kinOrig.x,kinOrig.Q2);
	}
      else
	{
	  //	  cout <<" not in bin " <<endl;
	}
      binRecQ2=log10(kinSmeared.Q2);
      binRecX=log10(mixedXSmeared);
      if(fabs(binlogOrigQ2-binRecQ2)<logQ2Range/(2*corrBinsQ2) && fabs(binlogOrigX-binRecX)<logXRange/(2*corrBinsX))
	{
	  plots->Q2VsXSmearMixed->Fill(kinOrig.x,kinOrig.Q2);
	}
	  
      plots->Q2CorrElec->Fill(kinOrig.Q2,kinSmeared.Q2);
      plots->Q2CorrJB->Fill(kinOrig.Q2,kinJBSmeared.Q2);
      plots->Q2CorrDA->Fill(kinOrig.Q2,kinDASmeared.Q2);
      plots->XCorrElec->Fill(kinOrig.x,kinSmeared.x);
      plots->XCorrJB->Fill(kinOrig.x,kinJBSmeared.x);
      plots->XCorrDA->Fill(kinOrig.x,kinDASmeared.x);
      plots->XCorrMixed->Fill(kinOrig.x,mixedXSmeared);

      plots->corrConventionalQ2->Fill(kinOrig.Q2,(kinOrig.Q2-kinSmeared.Q2)/kinOrig.Q2);
      plots->corrConventionalX->Fill(kinOrig.x,(kinOrig.x-kinSmeared.x)/kinOrig.x);
      plots->corrJBQ2_2->Fill(kinOrig.Q2,kinJBSmeared.Q2);
      plots->corrJBQ2->Fill(kinOrig.Q2,(kinOrig.Q2-kinJBSmeared.Q2)/kinOrig.Q2);
      plots->corrJBX->Fill(kinOrig.x,(kinOrig.x-kinJBSmeared.x)/kinOrig.x);
	      
      plots->corrDAQ2->Fill(kinOrig.Q2,(kinOrig.Q2-kinDASmeared.Q2)/kinOrig.Q2);
      plots->corrDAX->Fill(kinOrig.x,(kinOrig.x-kinDASmeared.x)/kinOrig.x);
      //	      cout <<" fill dax with " << (kinOrig.x-kinsDASmeared.x)/kinOrig.x<<endl;
      plots->corrMixedX->Fill(kinOrig.x,(kinOrig.x-mixedXSmeared)/kinOrig.x);

	  
      //get real Q2, x

      //get reconstructed JB etc

      //plot

	  
    }

      
}



int main(int argc, char** argv)
{

  
  int colors[]={kRed,kBlue,kGreen,kBlack,kCyan,kMagenta,kOrange,kYellow};
  int markerStyles[]={20,21,22,23,43,33,34,47};

  char buffer[200];
  gStyle->SetOptStat(0);

  
  srand(time(NULL));
  if(argc<4)
    {
      cout <<"filename and electron+hadron beam energy required " <<endl;
      exit(0);
    }
  //  TFile inputFile(argv[1]);
  
  int beamEnergyI = atoi(argv[2]);
  int hadronBeamEnergyI = atoi(argv[3]);
  double beamEnergy=0.0+beamEnergyI;
  double hadronBeamEnergy=0.0+hadronBeamEnergyI;
  
  cout <<"using beam energy of " << beamEnergy <<" and hadron beam energy of " << hadronBeamEnergy<<endl;
  //      cout <<"original e px: " << ePxOrig <<" py: " << ePyOrig<< " pz: "<< ePzOrig<<endl;
  float sqrtS=0;
  if(hadronBeamEnergy>200)
    {
      sqrtS=140.716; //sqrtS for 18x275;
    }
  else
    {
      if(hadronBeamEnergy<200 && hadronBeamEnergy>50)
	{
	  if(beamEnergy < 9)
	    sqrtS=45;
	  else
	    sqrtS=63.2527;
	}
      else
	{
	  sqrtS=28.6;
	}
    }

  TChain *chain = new TChain("Delphes");
  chain->Add(argv[1]);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;
  
  plots->zBins.push_back(0.2);
  plots->zBins.push_back(0.3);
  plots->zBins.push_back(0.5);
  plots->zBins.push_back(2.0);
  
  plots->yBins.push_back(0.1);
  plots->yBins.push_back(0.2);
  plots->yBins.push_back(0.5);
  plots->yBins.push_back(0.8);
  plots->yBins.push_back(2.0);
  
  cout <<" book histogram: " << endl;
  BookHistograms(result, plots,beamEnergyI,hadronBeamEnergyI);
  cout <<" analyze events: " << endl;
  AnalyzeEvents(treeReader, plots,beamEnergy, hadronBeamEnergy,sqrtS);
  cout <<"done analyzing, dividing... "<<endl;  
  plots->Q2VsXSmear->Divide(plots->Q2VsXSmearNorm);
  plots->Q2VsXSmearDA->Divide(plots->Q2VsXSmearNorm);
  plots->Q2VsXSmearMixed->Divide(plots->Q2VsXSmearNorm);
  plots->Q2VsXSmearJB->Divide(plots->Q2VsXSmearNorm);
  cout <<"plotting... "<<endl;  
  TCanvas c1;
  plots->angCorr->Draw("colz");

  
  c1.SetLogy();
  c1.SetLogx();
  //  cout <<"first... "<<endl;  

  plots->Q2VsXSmear->SetMaximum(1.0);
  plots->Q2VsXSmear->SetMinimum(0.0);
  
  plots->Q2VsXSmearDA->SetMaximum(1.0);
  plots->Q2VsXSmearDA->SetMinimum(0.0);

  plots->Q2VsXSmearJB->SetMaximum(1.0);
  plots->Q2VsXSmearJB->SetMinimum(0.0);
  plots->Q2VsXSmearDA->SetMaximum(1.0);
  plots->Q2VsXSmearDA->SetMinimum(0.0);
  plots->Q2VsXSmearMixed->SetMaximum(1.0);
  plots->Q2VsXSmearMixed->SetMinimum(0.0);



  int d=sqrt(plots->zBins.size());
  
  if(d< sqrt(plots->zBins.size()))
    d++;

  //  TCanvas c2;
  //add one for legend
  c1.Clear();
  c1.Divide(d+1,d);


  c1.SetLogy(false);
  c1.SetLogx(false);

  for(int j=0;j<plots->yBins.size();j++)
    {
	for(int k=0;k<plots->zBins.size();k++)
	{
	  c1.cd(k+1);
	  double max=0;
	  double min=100000000;
	  for(int i=0;i<recTypeEnd;i++)
	    {
	      if(max<plots->phiRHistos[i][k][j]->GetMaximum())
		max=plots->phiRHistos[i][k][j]->GetMaximum();
	      if(min>plots->phiRHistos[i][k][j]->GetMinimum())
		min=plots->phiRHistos[i][k][j]->GetMinimum();
	    }
	  for(int i=0;i<recTypeEnd;i++)
	    {
	      plots->phiRHistos[i][k][j]->SetMaximum(1.2*max);
	      plots->phiRHistos[i][k][j]->SetMinimum(1.2*min);
	    }
	  for(int i=0;i<recTypeEnd;i++)
	    {
	      cout <<"trying to plot i: "<< i <<" j : " << j << " k: "<< k <<endl;
	      plots->phiRHistos[i][k][j]->SetMarkerColor(colors[i]);
	      plots->phiRHistos[i][k][j]->SetLineColor(colors[i]);
	      plots->phiRHistos[i][k][j]->SetMarkerStyle(markerStyles[i]);
	      if(i==0)
		plots->phiRHistos[i][k][j]->Draw("P");
	      else
		plots->phiRHistos[i][k][j]->Draw("SAME P");

	      //save them individually
	      /*	      c1.cd(0);
	      plots->phiRHistos[i][k][j]->Draw();
	      sprintf(buffer,"allPhiRHistos/phiR_%s_zBin%d_yBin%d.png",recTypeNames[i].c_str(),k,j);
	      c1.SaveAs(buffer);
	      cout <<"done " <<endl;*/
	    }
	}
	auto legend=new TLegend(0,0,1.0,1,0);

	for(int i=0;i<recTypeEnd;i++)
	  {
	    legend->AddEntry(plots->phiRHistos[i][0][j],recTypeNames[i].c_str(),"p");

	    TH1* histo=(TH1*)gDirectory->Get(plots->phiRHistos[i][0][j]->GetName());
	    cout <<" pull name for legend: " << plots->phiRHistos[i][0][j]->GetName() <<endl;
	    cout <<"color for marker: "<< plots->phiRHistos[i][0][j]->GetMarkerColor()<<endl;;
	    cout <<" or " << histo->GetMarkerColor()<<endl;
	    
	    cout <<"style for marker: "<< plots->phiRHistos[i][0][j]->GetMarkerStyle()<<endl;
	    	    cout <<" or " << histo->GetMarkerStyle()<<endl;
	  }
	c1.cd(plots->zBins.size()+1);
	legend->Draw();
	legend->Draw();
	sprintf(buffer,"phiRDists_yBin_%d_%d_%d.png",j,beamEnergyI,hadronBeamEnergyI);
	c1.SaveAs(buffer);
	sprintf(buffer,"phiRDists_yBin_%d_%d_%d.root",j,beamEnergyI,hadronBeamEnergyI);
	c1.SaveAs(buffer);

    }
  c1.SetLogy();
  c1.SetLogx();

  c1.cd(0);

  //  c1.SetLogz();
  for(int i=0;i<recTypeEnd;i++)
    {

      plots->ptMeans[i]->Divide(plots->ptCounts[i]);
      plots->ptDiffs[i]->Divide(plots->ptCounts[i]);


      //account for the fact that we fill with z1, z2
      plots->zMeans[i]->Divide(plots->zMeans[i],plots->ptCounts[i],1.0,2.0);
      plots->zDiffs[i]->Divide(plots->zDiffs[i],plots->ptCounts[i],1.0,2.0);


      plots->phiRMeans[i]->Divide(plots->ptCounts[i]);
      plots->phiRDiffs[i]->Divide(plots->ptCounts[i]);

      plots->zMeans[i]->SetMaximum(1.0);
      plots->zDiffs[i]->SetMaximum(1.0);
      plots->zMeans[i]->SetMinimum(0.0);
      plots->zDiffs[i]->SetMinimum(0.0);
      plots->phiRDiffs[i]->SetMaximum(1.0);
      plots->phiRDiffs[i]->SetMinimum(0.0);

      
      sprintf(buffer,"ptMeans_%s_%d_%d.png",recTypeNames[i].c_str(),beamEnergyI,hadronBeamEnergyI);      
      plots->ptMeans[i]->Draw("colz");
      drawYContour(0.01,beamEnergy,hadronBeamEnergy);
      drawYContour(0.05,beamEnergy,hadronBeamEnergy);
      c1.SaveAs(buffer);
      
      sprintf(buffer,"ptDiffs_%s_%d_%d.png",recTypeNames[i].c_str(),beamEnergyI,hadronBeamEnergyI);      
      plots->ptDiffs[i]->Draw("colz");
      drawYContour(0.01,beamEnergy,hadronBeamEnergy);
      drawYContour(0.05,beamEnergy,hadronBeamEnergy);
      c1.SaveAs(buffer);

      sprintf(buffer,"phiRMeans_%s_%d_%d.png",recTypeNames[i].c_str(),beamEnergyI,hadronBeamEnergyI);      
      plots->phiRMeans[i]->Draw("colz");
      drawYContour(0.01,beamEnergy,hadronBeamEnergy);
      drawYContour(0.05,beamEnergy,hadronBeamEnergy);
      c1.SaveAs(buffer);
      sprintf(buffer,"phiRDiffs_%s_%d_%d.png",recTypeNames[i].c_str(),beamEnergyI,hadronBeamEnergyI);      
      plots->phiRDiffs[i]->Draw("colz");
      drawYContour(0.01,beamEnergy,hadronBeamEnergy);
      drawYContour(0.05,beamEnergy,hadronBeamEnergy);
      c1.SaveAs(buffer);


      sprintf(buffer,"zMeans_%s_%d_%d.png",recTypeNames[i].c_str(),beamEnergyI,hadronBeamEnergyI);

  
      plots->zMeans[i]->Draw("colz");
      drawYContour(0.01,beamEnergy,hadronBeamEnergy);
      drawYContour(0.05,beamEnergy,hadronBeamEnergy);
      c1.SaveAs(buffer);
      sprintf(buffer,"zDiffs_%s_%d_%d.png",recTypeNames[i].c_str(),beamEnergyI,hadronBeamEnergyI);      
      plots->zDiffs[i]->Draw("colz");
      drawYContour(0.01,beamEnergy,hadronBeamEnergy);
      drawYContour(0.05,beamEnergy,hadronBeamEnergy);
      c1.SaveAs(buffer);
    }
  c1.SetLogz(false);  


  
  plots->Q2VsXSmearDA->Draw("colz");
  drawYContour(0.01,beamEnergy,hadronBeamEnergy);
  drawYContour(0.05,beamEnergy,hadronBeamEnergy);
  sprintf(buffer,"Q2VsXSmearDA_%dx%d.png",beamEnergyI,hadronBeamEnergyI);
  c1.SaveAs(buffer);
  
  plots->Q2VsXSmearJB->Draw("colz");
  drawYContour(0.01,beamEnergy,hadronBeamEnergy);
  drawYContour(0.05,beamEnergy,hadronBeamEnergy);

  sprintf(buffer,"Q2VsXSmearJB_%dx%d.png",beamEnergyI,hadronBeamEnergyI);
  c1.SaveAs(buffer);
  plots->Q2VsXSmearMixed->Draw("colz");
  drawYContour(0.01,beamEnergy,hadronBeamEnergy);
  drawYContour(0.05,beamEnergy,hadronBeamEnergy);

  sprintf(buffer,"Q2VsXSmearMixed_%dx%d.png",beamEnergyI,hadronBeamEnergyI);
  c1.SaveAs(buffer);
  plots->Q2VsXSmear->Draw("colz");
  drawYContour(0.01,beamEnergy,hadronBeamEnergy);
  drawYContour(0.05,beamEnergy,hadronBeamEnergy);

  sprintf(buffer,"Q2VsXSmear_%dx%d.png",beamEnergyI,hadronBeamEnergyI);
  c1.SaveAs(buffer);

  c1.SetLogz();  
  cout <<" print histo: " << endl;
  PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;

}
