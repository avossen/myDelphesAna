#include "TMath.h"
#include "classes/DelphesClasses.h"

using namespace std;


float printLVect(TLorentzVector& v)
{
  cout <<"( " << v.Px() <<" " << v.Py() <<" " << v.Pz()<<" " << v.E()<< " )"<<endl;
}

float printVect(TVector3& v)
{
  cout <<"( " << v.Px() <<" " << v.Py() <<" " << v.Pz()<<" )"<<endl;
}
#include "HadronPair.h"
struct Kins
{
  float x;
  float Q2;
  float y;
  float nu;
  float W;
};


pair<double,double> getQz(double a, double b, double kappa, double y, double qT2, double Pz, double Pe,double Q2)
{

  double pHalf=1.0/(Pe*Pe*kappa)*(a*y+b)*Pz;
  double q=-qT2/kappa+Q2/kappa+(a*y+b)*(a*y+b)/(Pe*Pe*kappa);

  double c=pHalf*pHalf-q;
  
  if(c>0)
    {
      //      cout <<"c: "<< c <<endl;
      return pair<double,double>(-pHalf+sqrt(c),-pHalf-sqrt(c));
    }
  else
    {
      //cout <<"qpart is 0, c: " << c <<endl;
      return pair<double,double>(0,0);
    }

}

float getQe(double a, double b, double y, double qT2, double Pz, double qz, double Pe)
{
  return (a*y+b+Pz*qz)/Pe;
}

float drawYContour(float yval,float beamEnergy, float hadronBeamEnergy)
{
  // y contours
  Float_t EdotP; 
  // y-cut contour
  //  yval = 0.01; 
  EdotP = 2*beamEnergy*hadronBeamEnergy; // (scalar product of beams' 4 momenta)
  float Q2min=1.0;
  float xmax=1.0;
  TLine * yline = new TLine(
    Q2min / (2*yval*EdotP), Q2min, xmax, 2*xmax*yval*EdotP);
  yline->SetLineWidth(4);
  yline->SetLineStyle(2);
  yline->Draw();
  // clas12 y=1 line

}

static void BinLog(TAxis * axis)
    {
      Float_t lb = axis->GetXmin();
      Float_t ub = axis->GetXmax();
      cout <<"got lb: " << lb <<" ub: "<< ub <<endl;
      if(lb<=0||ub<=0||lb>=ub) {
        fprintf(stderr,"ERROR: bad axis range for Tools::BinLog\n");
        return;
      };
      lb = TMath::Log10(lb);
      ub = TMath::Log10(ub);
      Int_t nbins = axis->GetNbins();
      Float_t width = (ub-lb)/nbins;
      Float_t * newBins = new Float_t[nbins+1];
      for (int b=0; b<=nbins; b++) newBins[b] = TMath::Power(10,lb+b*width);
      axis->Set(nbins,newBins);
      delete[] newBins;
    } 

struct HadronicVars
{
  double sumPx;
  double sumPy;
  double sumPz;
  double sumE;
  double sumEMinusPz;
  double theta;
};


void studyReconstruction(TClonesArray* branchParticles, TClonesArray* branchEFlowTrack, TClonesArray* branchEFlowPhoton, TClonesArray* branchEFlowNeutralHadron, TClonesArray* branchTrack)
{

  //  cout <<"study rec " <<endl;
  for(int i=6;i<branchParticles->GetEntries();i++)
    {
      GenParticle* mParticle=((GenParticle*)branchParticles->At(i));
      if(mParticle->Status!=1)
	continue;


      //      cout <<"looking at generated particle: "<< mParticle->PID <<" energy: "<< mParticle->E<<" eta: "<< mParticle->Eta <<endl;
      //  cout <<"eflow tracks " <<endl;
        for(int j=0;j<branchEFlowTrack->GetEntries();j++)
	  {
	    GenParticle* p=(GenParticle*)((Track*)branchEFlowTrack->At(j))->Particle.GetObject();
	    if(p==mParticle)
	      {
		//found the particle in the track flow object
		Track* t=(Track*)branchEFlowTrack->At(j);
		//		cout <<"found corresponding flow track with E: "<< t->P4().E()<<endl;
	      }
	  }
	//	      cout <<"tracks " <<endl;
	for(int j=0;j<branchTrack->GetEntries();j++)
	  {
	    GenParticle* p=(GenParticle*)((Track*)branchTrack->At(j))->Particle.GetObject();
	    if(p==mParticle)
	      {
		//found the particle in the track flow object
		Track* t=(Track*)branchTrack->At(j);
		//		cout <<"found corresponding  track with E: "<< t->P4().E()<<endl;
	      }
	  }
	//      cout <<"photon " <<endl;
	for(int j=0;j<branchEFlowPhoton->GetEntries();j++)
	  {
	    int pi=((Tower*)branchEFlowPhoton->At(j))->Particles.IndexOf(mParticle);
	    if(pi>=0)
	      {
		Tower* t=(Tower*)branchEFlowPhoton->At(j);
		//		cout <<"found corresponding ecal tower with E: "<< t->E<<endl;
	      }
	  }
	//	cout <<"neutral hadron? " << branchEFlowNeutralHadron <<endl;
	if(branchEFlowNeutralHadron!=0)
	  {
	    for(int j=0;j<branchEFlowNeutralHadron->GetEntries();j++)
	      {
		int hi=((Tower*)branchEFlowNeutralHadron->At(j))->Particles.IndexOf(mParticle);
		if(hi>=0)
		  {
		    Tower* t=(Tower*)branchEFlowNeutralHadron->At(j);
		//		cout <<"found corresponding hcal tower with E: "<< t->E<<endl;
		  }
	      }
	  }
      
    }


}


//electron goes into the negative direction
Kins getKinsFromScatElectron(double eBeamEnergy, double hadronBeamEnergy,double scatElPx,double scatElPy,double scatElPz,double scatElectronEnergy)
{
  Kins ret;
  double m_e = 0.000510998950;
  double m_p = 0.938272088;
  TLorentzVector lv_e;// = new TLorentzVector();
  TLorentzVector lv_beam;// = new TLorentzVector();
  TLorentzVector lv_target;// = new TLorentzVector();

  //e2 =m2+p2 -->p2=e2-m2
  double longBeamMomentum=sqrt(eBeamEnergy*eBeamEnergy-m_e*m_e);
  double longHadBeamMomentum=sqrt(hadronBeamEnergy*hadronBeamEnergy-m_p*m_p);

  //  cout <<" longbeamMomentum: " << longBeamMomentum << " beam Mom: " << longHadBeamMomentum <<endl;
  lv_e.SetPxPyPzE(scatElPx, scatElPy, scatElPz, scatElectronEnergy);
  lv_beam.SetPxPyPzE(0.0, 0.0, -longBeamMomentum, eBeamEnergy);
  lv_target.SetPxPyPzE(0.0, 0.0, longHadBeamMomentum, hadronBeamEnergy);
  
    TLorentzVector q=lv_beam-lv_e;
  // q.sub(lv_beam);
  ret.Q2 = (-1) * q.M2();
  // need to multiply lorentz vectors... doesn't seem to be implemented in the
  // jlab clas
  // since this is lv_target*q /m_p and the momentum of the target is zero, it is
  // just the product of the energy parts, divided by m_p:
  ret.nu = lv_target.E() * q.E() / m_p;
  //only true for fixed target
  ret.y = (lv_target*q)/(lv_target*lv_beam);
  // System.out.println("target x " + lv_target.px()+ " y: "+lv_target.py() + "
  // pz: " +lv_target.pz() + " e: "+ lv_target.e() );
  ret.x = ret.Q2 / (2 * (lv_target*q));
  ret.W = sqrt(m_p+ret.Q2 * (1 - ret.x) / ret.x);
  return ret;
}


HadronicVars getOriginalHadronicVars(TClonesArray* branchParticles)
{
  HadronicVars ret;
  double EMinusPz=0;
  TLorentzVector vSum;
  for(int i=6;i<branchParticles->GetEntries();i++)
    {
      //at 5 we have the outgoing lepton which we shouldn't include in the hadronic final state
      //I guess dire at 4

      //      if(i==5)
      //      if(i==4)
      //	continue;
      int status=((GenParticle*)branchParticles->At(i))->Status;            
      if(status!=1)
	continue;
      TLorentzVector lvParticle=((GenParticle*)branchParticles->At(i))->P4();
      vSum+=lvParticle;
      EMinusPz+=lvParticle.E()-lvParticle.Pz();
    }
  ret.sumPx=vSum.Px();
  ret.sumPy=vSum.Py();
  ret.sumPz=vSum.Pz();
  ret.sumE=vSum.E();
  ret.sumEMinusPz=EMinusPz;
  ret.theta=2.*TMath::ATan(EMinusPz/vSum.Pt());
  return ret;
}

struct boostVars{
  TVector3 breitBoost;
  double W;
  TLorentzVector lv_q;
  TLorentzVector lv_l;
};


void getHadronPairs(vector<HadronPair>& pairsRec, vector<HadronPair>& pairsTrue,TClonesArray* branchTracks, boostVars recBoost, boostVars realBoost, double targetPz, bool useMatchedTruth, bool useGen, bool useNoAcc)
  {


    
    //      for(int i=6;i<branchParticles->GetEntries();i++)
    //    {
    //      GenParticle* mParticle=((GenParticle*)branchParticles->At(i));
    //      if(mParticle->Status!=1)
    //	continue;

  for(int i=0;i<branchTracks->GetEntries();i++)
    {
     TObject* t=branchTracks->At(i);
      TLorentzVector lv1;

      if(useGen)
	lv1=((GenParticle*)t)->P4();
      else
	lv1=((Track*)t)->P4();
      
      //      Track* t=(Track*)branchTracks->At(i);
      if(!useNoAcc)
	{	
	  if(fabs(lv1.Eta())>4.0 || (fabs(lv1.P())<0.01))
	    {
	      continue;
	    }
	}
      //      TLorentzVector track_mom=t->P4();
      GenParticle* tp;

      //just the same
      if(useGen)
	tp=(GenParticle*)t;
      else
	tp= (GenParticle*) ((Track*)t)->Particle.GetObject();
      if(isnan(lv1.E()))
	 continue;

      //      TLorentzVector genP=tp->P4();
      //      cout <<"looking at rec track: "<< printLVect(track_mom) <<endl;
      //      cout <<"corresponding truth: "<< printLVect(genP)<<endl;
      
      //pi+
      if(tp->PID!=211)
	continue;

        for(int j=0;j<branchTracks->GetEntries();j++)
	  {
	    TObject* t2=branchTracks->At(j);
	    TLorentzVector lv2;
	    
	    if(useGen)
	      lv2=((GenParticle*)t2)->P4();
	    else
	      lv2=((Track*)t2)->P4();
	    
		  //		  Track* t2=(Track*)branchTracks->At(j);
		  //	    TLorentzVector track_mom2=t2->P4();
	    GenParticle* tp2 = 0;
	    if(useGen)
	      tp2=(GenParticle*)t2;
	    else
	      tp2=(GenParticle*)((Track*) t2)->Particle.GetObject();

	    if(!useNoAcc)
	      {
		if(fabs(lv2.Eta())>4.0 || (fabs(lv2.P())<0.01))
		  {
		    continue;
		  }
	      }
	    if(isnan(lv2.E()))
	      continue;
	    //pi-
	    if(tp2->PID!=-211)
	      continue;
	    //should put minimum cuts here, maybe z1, z2> 0.1
	    TLorentzVector real1=tp->P4();
	    TLorentzVector real2=tp2->P4();
	    
	    HadronPair pairRec(lv1,lv2,recBoost.breitBoost,recBoost.W,recBoost.lv_q,recBoost.lv_l,targetPz);
	    HadronPair pairTruth(real1,real2,realBoost.breitBoost,realBoost.W,realBoost.lv_q,realBoost.lv_l,targetPz);

	    //	    //	    cout <<" rec z1: "<< pairRec.z1 <<" truth : "<< pairTruth.z1 <<" z2: "<< pairRec.z2 <<" truth: "<< pairTruth.z2 <<endl;
	    if(useMatchedTruth)
	      {
	      	    if(pairTruth.z1< 0.05 || pairTruth.z2<0.05)
		      continue;
	      }
	    else
	      {
		if(pairRec.z1< 0.05 || pairRec.z2<0.05)
		  {
		    //		cout <<"z 1: "<< pairRec.z1 <<" z2: "<< pairRec.z2 <<endl;
		    continue;
		  }
	      }
	    //	    cout <<"combining " << printLVect(track_mom) <<" and " << printLVect(track_mom2) <<endl;
	    //	    cout <<"combining true " << printLVect(real1) <<" and " << printLVect(real2) <<endl;
	    pairsRec.push_back(pairRec);
	    pairsTrue.push_back(pairTruth);
	    //	    cout <<"ret " <<endl; 
	  }
    }
      //      if(useTruth)
      //	{
      //	  GenParticle* tp = (GenParticle*) t->Particle.GetObject();
      //	  track_mom=tp->P4();
      //	}
  }

    
bool checkSanity(HadronPair& pair)
{
  bool sane=true;
  if(isnan(pair.pT))
    {
      cout <<"nan pt" << endl;
      sane=false;
    }
  if(isnan(pair.phi_R))
    {
      cout <<" nan phiR" <<endl;
            sane=false;

    }
  if(isnan(pair.z1))
    {

      cout <<" nan z1" <<endl;
            sane=false;

    }
  if(isnan(pair.z2))
    {
      cout <<" nan z2 "<<endl;
                  sane=false;

    }
  if(sane)
    return true;
  else
    return false;
}

int getBin(vector<float>& b1, float value)
{
  int coo1=-1;

  for(unsigned int i=0;i<b1.size();i++)
    {
      if(value<=b1[i])
	{
	coo1=i;
	break;
	}
    }
  /*  if(coo1<0)
    {
        cout <<"wrong coo: val: " << value <<endl;
	}*/
  //  cout <<"value: " << value <<" coo: " << coo1 <<endl;
  return coo1;
}

HadronicVars getHadronicVars(TClonesArray* branchElectron,TClonesArray* branchEFlowTrack, TClonesArray* branchEFlowPhoton, TClonesArray* branchEFlowNeutralHadron, TClonesArray* branchTrack, TClonesArray* branchParticles)
{
  bool useTruth=false;
  
  HadronicVars ret;

  double delta_track=0;
  double delta_track_noel=0;
  double delta_photon=0;
  double delta_neut=0;
  
  TLorentzVector vSum(0.0,0.0,0.0,0.0);
  TLorentzVector vSumNoEl(0.0,0.0,0.0,0.0);
  TLorentzVector e(0.0,0.0,0.0,0.0);//only first candidate

  ////use truth for eta > 4 
  bool suppWithTruth=false;
  bool noTruthNoTrack=true;
  for(int i=0;i<branchEFlowTrack->GetEntries();i++)
    {
      Track* t=(Track*)branchEFlowTrack->At(i);

      if(fabs(t->Eta)>4.0 || (fabs(t->P)<0.01))
	{
	  continue;
	}
      TLorentzVector track_mom=t->P4();
      if(isnan(track_mom.E()))
	 continue;
      if(useTruth)
	{
	  GenParticle* tp = (GenParticle*) t->Particle.GetObject();
	  track_mom=tp->P4();
	}
      vSum+=track_mom;
      delta_track+=(track_mom.E()-track_mom.Pz());
    }
  //use tracks for eta>4.0
  for(int i=0;i<branchTrack->GetEntries();i++)
    {

      if(suppWithTruth||noTruthNoTrack)
	continue;
      Track* t=(Track*)branchTrack->At(i);

      if(fabs(t->Eta)<4.0 || (fabs(t->P)<0.01))
	{
	  continue;
	}
      TLorentzVector track_mom=t->P4();
      if(isnan(track_mom.E()))
	 continue;
      if(useTruth)
	{
	  GenParticle* tp = (GenParticle*) t->Particle.GetObject();
	  track_mom=tp->P4();
	}
      vSum+=track_mom;
      delta_track+=(track_mom.E()-track_mom.Pz());
    }
  
  for(int i=0;i<branchEFlowPhoton->GetEntries();i++)
    {
      Tower* to=(Tower*)branchEFlowPhoton->At(i);
      if(fabs(to->Eta)>4.0)
	continue;
      TLorentzVector ph_mom=to->P4();
      if(isnan(ph_mom.E()))
	 continue;

      if(useTruth)
	{
	  //	  GenParticle* tp = (GenParticle*) to->Particle.GetObject();
	  //	  ph_mom=tp->P4();
	}

      vSum+=ph_mom;
      delta_track+=(ph_mom.E()-ph_mom.Pz());
    }

  if(branchEFlowNeutralHadron!=0)
    {    
      for(int i=0;i<branchEFlowNeutralHadron->GetEntries();i++)
	{
	  Tower* nh=(Tower*)branchEFlowNeutralHadron->At(i);
	  if(fabs(nh->Eta)>4.0)
	    continue;
 
	  TLorentzVector neutralHadron=nh->P4();
	  if(isnan(neutralHadron.E()))
	    continue;

	  if(useTruth)
	    {
	      //	  GenParticle* tp = (GenParticle*) nh->Particle.GetObject();
	      //	  neutralHadron=tp->P4();
	    }
	  vSum+=neutralHadron;
	  delta_track+=(neutralHadron.E()-neutralHadron.Pz());
      
	}
    }
    if(suppWithTruth)
      {
	for(int i=6;i<branchParticles->GetEntries();i++)
	  {
	    //at 5 we have the outgoing lepton which we shouldn't include in the hadronic final state
	    //I guess dire at 4

	    //      if(i==5)
	    //      if(i==4)
	    //	continue;
	    GenParticle* part= ((GenParticle*)branchParticles->At(i));
	    int status=((GenParticle*)branchParticles->At(i))->Status;
	    if(fabs(part->Eta)<4)
	      continue;
	    if(status!=1)
	      continue;

	    TLorentzVector lvParticle=((GenParticle*)branchParticles->At(i))->P4();

	    vSum+=lvParticle;
	    delta_track+=lvParticle.E()-lvParticle.Pz();
	  }
      }
    //subtract scattered electron
  delta_track_noel=delta_track;
  vSumNoEl=vSum;
  if(branchElectron->GetEntries()>0)
    {
      Electron* el=(Electron*)branchElectron->At(0);
      TLorentzVector e=el->P4();
      if(useTruth)
	{
	  GenParticle* te = (GenParticle*) el->Particle.GetObject();
	  e=te->P4();
	}
      delta_track_noel = delta_track - (e.E() - e.Pz());
      vSumNoEl-=e;
    }

  ret.sumPx=vSumNoEl.Px();
  ret.sumPy=vSumNoEl.Py();
  ret.sumPz=vSumNoEl.Pz();
  ret.sumE=vSumNoEl.E();
  ret.sumEMinusPz=delta_track_noel;
  ret.theta=2.*TMath::ATan(delta_track_noel/vSumNoEl.Pt());

  return ret;
}



Kins getKinsJB(HadronicVars v, double beamEnergy, double hadronBeamEnergy,double s)
{
  Kins ret;
  double pt2=v.sumPx*v.sumPx+v.sumPy*v.sumPy;
  ret.y=v.sumEMinusPz/(2*beamEnergy); ///-->this is textbook
  //hadron beam pz
  ////  double hadPz=sqrt(hadronBeamEnergy*hadronBeamEnergy-ProtonMass*ProtonMass);
  ////  double elPz=sqrt(beamEnergy*beamEnergy-ElectronMass*ElectronMass);
  
  //apparently there is some improvement to be made in teh presence of rad corrections (see eicsmear code)
  /////-->doesn't seem to work that great  ret.y=(hadronBeamEnergy*v.sumE-hadPz*v.sumPz-ProtonMass*ProtonMass)/(beamEnergy*hadronBeamEnergy-elPz*hadPz);
  
  if(ret.y>0)
    {
      ret.Q2=pt2/(1.-ret.y);
      ret.x=ret.Q2/(s*ret.y);
    }
  return ret;
}
Kins getKinsDA(double scatElPx,double scatElPy,double scatElPz, double scatElectronEnergy,double beamEnergy,double tP, double s)
{
  Kins ret;
  TLorentzVector lv_e;
  lv_e.SetPxPyPzE(scatElPx, scatElPy, scatElPz, scatElectronEnergy);
  double t=lv_e.Theta();
  //  double nominator=4*lv_e.E()*lv_e.E()*cos(t/2)*cos(t/2);
  //  double denom=sin(t/2)*sin(t/2)+sin(t/2)*cos(t/2)*tan(tP/2);
  //      ret.Q2=nominator/denom;
  //  ret.y=1.0-(sin(t/2))/(sin(t/2)+cos(t/2)*tan(tP/2));

  //from the eicsmear sourcecode
  double theta=t;
  double gamma=tP;
  double denominator = tan(theta / 2.) + tan(gamma / 2.);
    if (denominator > 0.) {
      ret.y = tan(gamma / 2.) / denominator;
      ret.Q2 = 4. * pow(beamEnergy, 2.) / tan(theta / 2.) / denominator;
    }  // if             

  //this is from the Bluemlein paper, but doesn't seem to work


  ///Zeus paper:
  /// ret.y=tan(tP/2);
  ///ret.y/=(tan(t/2)+tan(tP/2));
  
  if(ret.y>0)
    {
      ret.x=ret.Q2/(s*ret.y);
    }
  
  
  return ret;
  
}
