#ifndef HadronPair_H
#define HadronPair_H



class HadronPair
{
 public:
  float phi_R;
  float phi_h;
  int polarization;
  float z;
  float M;
  float pT;
  float xF;
  float m_theta;
  float z1;
  float z2;
  float weight;
  float weightUncert;
  
  //for the theta coverage
  float pt1;
  float pt2;
  
  TLorentzVector m_h1;
  TLorentzVector m_h2;

  TVector3 m_breitBoost;
  bool hasMatchingMC;
  float m_W;
  float m_qE;
  TVector3 vecQ;
  TVector3 vecQLab;
  TVector3 vecL;
  TLorentzVector origQ;
  double pTLab;
  double m_targetPz;
  
  HadronPair(TLorentzVector h1, TLorentzVector h2, TVector3 breitBoost, float W, TLorentzVector& lv_q, TLorentzVector& lv_l, double targetPz)
    {
      //always true in fast simulation
      hasMatchingMC=true;
      m_targetPz=targetPz;
      m_h1 = h1;
      m_h2 = h2;
      pt1=h1.Perp();
      pt2=h2.Perp();
      m_breitBoost = breitBoost;
      m_W = W;
      m_qE = lv_q.E();
      origQ=lv_q;
      TLorentzVector lvQ = lv_q;
      TLorentzVector lvQLab =lv_q;
      TLorentzVector lvL = lv_l;
      lvQ.Boost(breitBoost);
      lvL.Boost(breitBoost);
      vecQ = lvQ.Vect();
      vecL = lvL.Vect();
      vecQLab = lv_q.Vect();
      vecQLab=vecQLab.Unit();
      compute();
    }
  
  void compute()
  {

    TLorentzVector pionPair = (m_h1)+(m_h2);
    M = pionPair.M();
    //    cout <<"pion pair mass: "<< M <<endl;
    TLorentzVector boostedPair(pionPair);
    //    cout <<"boostedPair: "<< printLVect(pionPair)<<endl;
    //    cout <<"boostVect: "<< printVect(m_breitBoost);
    boostedPair.Boost(m_breitBoost);
    //    cout<<"boosted pair:"<< printLVect(boostedPair)<<endl;
    TVector3 vecQUnit;
    vecQUnit.SetMagThetaPhi(1.0, vecQ.Theta(), vecQ.Phi());
    TVector3 vecQLabUnit;
    vecQLabUnit.SetMagThetaPhi(1.0, vecQLab.Theta(), vecQLab.Phi());
    //double otherPt = vecQ.cross(pionPair.vect()).mag();
    double otherPt = vecQUnit.Cross(boostedPair.Vect()).Mag();
		// System.out.println("pt " + pT + " or "+otherPT);
    pT = otherPt;
    pTLab = vecQLabUnit.Cross(pionPair.Vect()).Mag();
    //System.out.println("pt " + pT + " or "+pTLab);
    xF = boostedPair.Pz() / m_W;
    TLorentzVector lTarget;
    lTarget.SetPxPyPzE(0,0,m_targetPz,sqrt(m_targetPz*m_targetPz+0.9315*0.9315));

    z=lTarget*pionPair/(lTarget*origQ);

    //should be in lab system, so that the proton momentum
    //doesn't have to be accounted for
    //    cout <<"mqe: " << m_qE <<" pion1 e: "<< m_h1.E() <<" pion2 e: "<< m_h2.E() <<" pair: "<< pionPair.E()<<endl;
    //only valid for system in which the proton is at rest...
    z1=lTarget*m_h1/(lTarget*origQ);
    z2=lTarget*m_h2/(lTarget*origQ);

    
    TLorentzVector vh1(m_h1.Px(), m_h1.Py(), m_h1.Pz(), m_h1.E());
    TLorentzVector vh1T =vh1;
    TVector3 pairBoostVect = boostedPair.BoostVector();
    pairBoostVect=(-1)*pairBoostVect;
    vh1T.Boost(pairBoostVect);
    m_theta = acos(vh1T.Vect().Dot(boostedPair.Vect()) / (vh1T.Vect().Mag() * boostedPair.Vect().Mag()));
    vh1.Boost(m_breitBoost);
    TLorentzVector vh2(m_h2.Px(), m_h2.Py(), m_h2.Pz(), m_h2.E());
    vh2.Boost(m_breitBoost);
    //so vecR is now in the Breit frame
    TVector3 vecR(vh1.Vect());
  //scale by 1/z1
    vecR.SetMagThetaPhi(vh1.Vect().Mag()/z1, vh1.Vect().Theta(),vh1.Vect().Phi());
    TVector3 scaledVh2;
    scaledVh2.SetMagThetaPhi(vh2.Vect().Mag()/z2, vh2.Vect().Theta(), vh2.Vect().Phi());
    vecR=vecR-scaledVh2;
    //get vecRT:
    TVector3 vecRt;
    TVector3 RAlongQ;
    
    RAlongQ.SetMagThetaPhi(vecR.Dot(vecQUnit), vecQUnit.Theta(), vecQUnit.Phi());
    vecRt=vecR;
    //The transverse part is the original vector minus the one that
    //is along q
    //System.out.println("lenght of r along q : " + RAlongQ.mag());
    vecRt=vecRt-RAlongQ;
    //System.out.println("length of R " + vecRt.mag() + " phi " + vecRt.phi() + " theta: " + vecRt.theta());
    TVector3 vecPh(boostedPair.Vect());
    
    TVector3 PtAlongQ;
    PtAlongQ.SetMagThetaPhi(vecR.Dot(vecQUnit), vecQUnit.Theta(), vecQUnit.Phi());
    TVector3 vecPhT(vecPh);
    vecPhT=vecPhT=PtAlongQ;
    
    TVector3 vT(vecQUnit.Cross(vecL));
    //vT.setMagThetaPhi(1.0, vT.theta(), vT.phi());
    vT=vT.Unit();
    TVector3 vTR(vecQUnit.Cross(vecRt));
    TVector3 vTH(vecQUnit.Cross(vecPhT));
    vTR=vTR.Unit();
    vTH=vTH.Unit();
    double cosPhiR = vT.Dot(vTR);
    double cosPhiH = vT.Dot(vTH);
    
    double sinPhiR = vecL.Cross(vecRt).Dot(vecQUnit);
    ///
  //the scaling for cosPhiR is not necessary anymore since
  //for that quantity we already operate
  //with unit vectors above, scaling again would lead to the wrong value
  //
    double rScale = vecQUnit.Cross(vecL).Mag() * vecQUnit.Cross(vecRt).Mag();
    sinPhiR = sinPhiR / rScale;
    double sinPhiH = vecL.Cross(vecPhT).Dot(vecQUnit);
    double hScale = vecQUnit.Cross(vecL).Mag() * vecQUnit.Cross(vecPh).Mag();
    sinPhiH = sinPhiH / hScale;
    
    phi_h = acos(cosPhiH);
    if (sinPhiH < 0.0) {
      phi_h = 2 * TMath::Pi() - phi_h;
    }
    phi_R = acos(cosPhiR);
    if (sinPhiR < 0.0) {
      phi_R = 2 * TMath::Pi() - phi_R;
      
    }
  //new way:
  // turns out it is the same as the old way with the exception that
  //the angle is in the range -pi - pi instead of 0 to 2pi
  //

    double phi_R2Sign=(vecQUnit.Cross(vecL)).Dot(vecRt)/fabs((vecQUnit.Cross(vecL)).Dot(vecRt));
    double phi_R2=(vecQUnit.Cross(vecL)).Dot(vecQUnit.Cross(vecRt));
    phi_R2=acos(phi_R2/((vecQUnit.Cross(vecL)).Mag()*(vecQUnit.Cross(vecRt)).Mag()));
    phi_R2=phi_R2*phi_R2Sign;
    
  }

};

#endif
