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
#include "studyLambda.h"
using namespace std;


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


  TH1* lambdaMass;
  TH1* lambdaIJMass;
    


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

};

void BookHistograms(ExRootResult *result, TestPlots *plots,int beamEnergyI, int  hadronBeamEnergyI)
{
  TLegend *legend;
  TPaveText *comment;
  char buffer[300];
  sprintf(buffer,"Fraction of events staying in bin (%dx%d)",beamEnergyI,hadronBeamEnergyI);
  plots->lambdaMass=result->AddHist1D("lambdaMass", "lambda mass","mass [GeV]","counts",50,0.8,1.3);
  plots->lambdaIJMass=result->AddHist1D("lambdaInJetMass", "lambda in jet mmass","mass [GeV]","counts",50,0.8,1.3);
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
      bool useTruth=false;
  //    bool useTruth=true;
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
	  //	  cout <<"particle at index " <<i   <<" pid: " << g->PID << " energy: "<< g->E <<endl;
	}
    }


      getLambdas(branchParticle,branchEFlowTrack,plots->lambdaMass);
      getLambdasInJets(branchParticle,branchTrack,branchJet, plots->lambdaIJMass);
      
      
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
      //      cout <<"true q2: " << Q2 <<" x: " << x <<" y: " << y <<endl;
      
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
      studyReconstruction(branchParticle,branchEFlowTrack,branchEFlowPhoton,branchEFlowNeutralHadron,branchTrack);	  
      //      cout <<"rec electron E: "<< e.E()<<endl;
      //      cout << "true electron 1: "<< pleptonOut.E() <<endl;
      particle = (GenParticle*) electron->Particle.GetObject();
      //      cout <<"truth from rec electron: "<< particle->E<<endl;
      //      Kins kinOrig=getKinsFromScatElectron(pleptonIn.E(),pProton.E(),pleptonOut.Px(),pleptonOut.Py(),pleptonOut.Pz(),pleptonOut.E());
            Kins kinOrig=getKinsFromScatElectron(beamEnergy,hadronBeamEnergy,particle->Px,particle->Py,particle->Pz,particle->E);

      Kins kinSmeared=getKinsFromScatElectron(beamEnergy,hadronBeamEnergy,e.Px(),e.Py(),e.Pz(),e.E());
      HadronicVars hvOrig=getOriginalHadronicVars(branchParticle);
      HadronicVars hvSmeared=getHadronicVars(branchElectron, branchEFlowTrack, branchEFlowPhoton,branchEFlowNeutralHadron,branchTrack,branchParticle);

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
  
  //  cout <<"using beam energy of " << beamEnergy <<" and hadron beam energy of " << hadronBeamEnergy<<endl;
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
    cout <<"first... "<<endl;  

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



  
    plots->Q2VsXSmearDA->Draw("colz");
  sprintf(buffer,"Q2VsXSmearDA_%dx%d.png",beamEnergyI,hadronBeamEnergyI);
  c1.SaveAs(buffer);
  plots->Q2VsXSmearJB->Draw("colz");
  sprintf(buffer,"Q2VsXSmearJB_%dx%d.png",beamEnergyI,hadronBeamEnergyI);
  c1.SaveAs(buffer);
  plots->Q2VsXSmearMixed->Draw("colz");
  sprintf(buffer,"Q2VsXSmearMixed_%dx%d.png",beamEnergyI,hadronBeamEnergyI);
  c1.SaveAs(buffer);
  plots->Q2VsXSmear->Draw("colz");
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
