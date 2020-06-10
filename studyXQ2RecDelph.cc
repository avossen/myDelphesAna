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


class ExRootTreeReader;
class ExRootResult;

struct TestPlots
{
  TH2* Q2VsXSmear;
  TH2* Q2VsXSmearDA;
  TH2* Q2VsXSmearJB;
  TH2* Q2VsXSmearMixed;
  TH2* Q2VsXSmearNorm;

};

void BookHistograms(ExRootResult *result, TestPlots *plots,int beamEnergyI, int  hadronBeamEnergyI)
{
  TLegend *legend;
  TPaveText *comment;
  char buffer[300];
  sprintf(buffer,"Fraction of events staying in bin (%dx%d)",beamEnergyI,hadronBeamEnergyI);
  plots->Q2VsXSmear=result->AddHist2D("Q2VsXSmear",buffer,"x","Q^{2}",10,0.0001,1,10,0.5,200,1,1);
  plots->Q2VsXSmearDA=result->AddHist2D("Q2VsXSmearDA",buffer,"x","Q^{2}",10,0.0001,1,10,0.5,200,1,1);
  plots->Q2VsXSmearJB=result->AddHist2D("Q2VsXSmearJB",buffer,"x","Q^{2}",10,0.0001,1,10,0.5,200,1,1);
  plots->Q2VsXSmearMixed=result->AddHist2D("Q2VsXSmearMixed",buffer,"x","Q^{2}",10,0.0001,1,10,0.5,200,1,1);
  plots->Q2VsXSmearNorm=result->AddHist2D("Q2VsXSmearNorm",buffer,"x","Q^{2}",10,0.0001,1,10,0.5,200,1,1);
  
  //  BinLog(plots->Q2VsXSmear->GetXaxis());
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
  float logQ2Range=2.6;
  float logXRange=4;

  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
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
	  cout <<"going through entry " << entry <<endl;
	  cout <<" is there a pointer here? " << branchParticle->At(0) <<endl;
	 TLorentzVector  pProton    =  ((GenParticle*)branchParticle->At(0))->P4();// #these numbers 0 , 3, 5 are hardcoded in Pythia8
	 

	 TLorentzVector	  pleptonIn    = ((GenParticle*)branchParticle->At(3))->P4();
	 cout <<"proton... has energy: " << pProton.E() <<" lepton in: " << pleptonIn.E() <<endl;
	 cout <<"lepton..." <<endl;
	 TLorentzVector	  pleptonOut   = ((GenParticle*)branchParticle->At(5))->P4();
	 cout <<"lep out..." <<endl;
	  TLorentzVector	  pPhoton      = pleptonIn - pleptonOut;
	 cout <<"photon..." <<endl;
	  //Q2, W2, Bjorken x, y, nu.                                                                                                                                                  
	double  Q2 = -pPhoton.M2();
	 double  W2 = (pProton + pPhoton).M2();
	 double  x = Q2 / (2. * pProton.Dot(pPhoton));
	double  y = (pProton.Dot(pPhoton)) / (pProton.Dot(pleptonIn));

	  cout <<"reconstructed: q2: " << Q2 <<" x: " << x <<" y: " << y <<endl;
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
	  
	  cout <<"rec electron E: "<< e.E()<<endl;
	  cout << "true electron 1: "<< pleptonOut.E() <<endl;
	  particle = (GenParticle*) electron->Particle.GetObject();
	  cout <<"truth from rec electron: "<< particle->E<<endl;
	  Kins kinOrig=getKinsFromScatElectron(pleptonIn.E(),pProton.E(),pleptonOut.Px(),pleptonOut.Py(),pleptonOut.Pz(),pleptonOut.E());
	  Kins kinSmeared=getKinsFromScatElectron(beamEnergy,hadronBeamEnergy,e.Px(),e.Py(),e.Pz(),e.E());
	  HadronicVars hvSmeared=getHadronicVars(branchElectron, branchEFlowTrack, branchEFlowPhoton,branchEFlowNeutralHadron);
	  Kins kinJBSmeared=getKinsJB(hvSmeared,beamEnergy,hadronBeamEnergy,sqrtS*sqrtS);
	  Kins kinsDASmeared=getKinsDA(e.Px(),e.Py(),e.Pz() ,e.E(),beamEnergy,hvSmeared.theta, sqrtS*sqrtS);
	  
	  float mixedXSmeared=kinSmeared.Q2/(sqrtS*sqrtS*kinJBSmeared.y);
	  double binlogOrigQ2=log(kinOrig.Q2);
	  double binlogOrigX=log(kinOrig.x);
	  double binRecQ2=log(kinSmeared.Q2);
	  double binRecX=log(kinSmeared.x);
	
	  binRecQ2=log(kinsDASmeared.Q2);
	  binRecX=log(kinsDASmeared.x);

	  if(fabs(binlogOrigQ2-binRecQ2)<logQ2Range/20 && fabs(binlogOrigX-binRecX)<logXRange/20)
	    {

	    }
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
  cout <<" book histogram: " << endl;
  BookHistograms(result, plots,beamEnergyI,hadronBeamEnergyI);
  cout <<" analyze events: " << endl;
  AnalyzeEvents(treeReader, plots,beamEnergy, hadronBeamEnergy,sqrtS);
  cout <<" print histo: " << endl;
  PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;

}
