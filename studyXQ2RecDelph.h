#include "TMath.h"
using namespace std;
struct Kins
{
  float x;
  float Q2;
  float y;
  float nu;
  float W;
};

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


HadronicVars getOriginalHadronicVar(TClonesArray* branchParticles)
{
  for(int i=0;branchParticles->GetEntries();i++)
    {
      //at 5 we have the outgoing lepton which we shouldn't include in the hadronic final state
      if(i==5)
	continue;
            
      if(status!=1)
	continue;

      TLorentzVector lvParticle=((Particle*)branchParticles->At(i))->P4();
      int status=((Particle*)branchParticles->At(i))->Status();
      
    }
}

HadronicVars getHadronicVars(TClonesArray* branchElectron,TClonesArray* branchEFlowTrack, TClonesArray* branchEFlowPhoton, TClonesArray* branchEFlowNeutralHadron)
{
  HadronicVars ret;
  TVector3 temp_p;
  double delta_track=0;
  double delta_track_noel=0;
  double delta_photon=0;
  double delta_neut=0;
  TLorentzVector vSum;
  TLorentzVector vSumNoEl;
  TLorentzVector e;//only first candidate
  for(int i=0;i<branchEFlowTrack->GetEntries();i++)
    {
      TLorentzVector track_mom=((Track*)branchEFlowTrack->At(i))->P4();
      vSum+=track_mom;
      delta_track+=(track_mom.E()-track_mom.Pz());
      temp_p = temp_p+ track_mom.Vect();

    }
  //subtract scattered electron
  if(branchElectron->GetEntries()>0)
    {
       e = ((Electron*)branchElectron->At(0))->P4();
      delta_track_noel = delta_track - (e.E() - e.Pz());
    }
  
  for(int i=0;i<branchEFlowPhoton->GetEntries();i++)
    {
      TLorentzVector photon = ((Tower*)branchEFlowPhoton->At(i))->P4();
      vSum+=photon;	
    }

  for(int i=0;i<branchEFlowNeutralHadron->GetEntries();i++)
    {
      TLorentzVector neutralHadron=((Tower*)branchEFlowNeutralHadron->At(i))->P4();
      vSum+=neutralHadron;
    }
  ret.sumPx=vSum.Px();
  ret.sumPy=vSum.Py();
  ret.sumPz=vSum.Pz();
  ret.sumE=vSum.E();
  ret.sumEMinusPz=delta_track_noel;
  ret.theta=2.*TMath::ATan((delta_track_noel)/vSum.Pt());
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
