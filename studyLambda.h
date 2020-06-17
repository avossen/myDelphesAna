void getLambdas(TClonesArray* branchParticles, TClonesArray* branchEFlowTrack, TH1* histo,TVector3& boostVect,double W)
{
    for(int i=0;i<branchEFlowTrack->GetEntries();i++)
      {
	//also look at permutations...
	for(int j=0;j<branchEFlowTrack->GetEntries();j++)
	  {
	    Track* t1=(Track*)branchEFlowTrack->At(i);
	    Track* t2=(Track*)branchEFlowTrack->At(j);

	    GenParticle* part1=(GenParticle*)t1->Particle.GetObject();
	    GenParticle* part2=(GenParticle*)t2->Particle.GetObject();

	    int part1Mother=part1->M1;
	    int part2Mother=part2->M1;
	    int  numP=branchParticles->GetEntries();
	    int m1PID=0;
	    int m2PID=0;
	    GenParticle* gen1;
	    GenParticle* gen2;
	    //	    cout <<"PID: "<< t1->PID <<endl;
	    if(t1->PID==2212 && t2->PID==-211)
	      {

		if(numP>part1Mother && part1Mother >=0)
		  {
		     gen1=((GenParticle*)branchParticles->At(part1Mother));
		    //  cout <<"first mother: "<< mParticle->PID<<endl;
		    m1PID=gen1->PID;
		  }
		if(numP>part2Mother && part2Mother>=0)
		  {
		     gen2=((GenParticle*)branchParticles->At(part2Mother));
		    m2PID=gen2->PID;
				    //		    cout <<"second mother: "<< mParticle->PID<<endl;
		  }
		//		if(m1PID==3122 && m2PID==3122 && (part1Mother==part2Mother))
		  {
		    if(m1PID==3122)
		      cout << "found proton " << endl;
		    else
		      cout <<"found pion " <<endl;
		    cout <<"first px : "<< t1->P4().Px() <<" py: " << t1->P4().Py() <<" pz: "<< t1->P4().Pz() <<" E: "<< t1->P4().E()<<endl;
		    cout <<"second px : "<< t2->P4().Px() <<" py: " << t2->P4().Py() <<" pz: "<< t2->P4().Pz() <<" E: "<< t2->P4().E()<<endl;
		    cout <<" t1 mass: "<< t1->P4().M() << " t2 mass: "<< t2->P4().M()<<endl;
		    cout <<"vertex t1: "<< t1->X <<", "<<t1->Y <<", " << t1->Z<<endl;
		    cout <<"vertex t2: "<< t2->X <<", "<<t2->Y <<", " << t2->Z<<endl;
		    cout <<"vertex p1: "<< part1->X <<", "<<part1->Y <<", " << part1->Z<<endl;
		    cout <<"vertex p2: "<< part2->X <<", "<<part2->Y <<", " << part2->Z<<endl;

		    TLorentzVector tl1=t1->P4();
		    TLorentzVector tl2=t2->P4();
		    double mP=0.93828;
		    double mPi=0.139568;
		    tl1.SetE(sqrt(tl1.P()*tl1.P()+mP*mP));
		    tl2.SetE(sqrt(tl2.P()*tl2.P()+mPi*mPi));
		    
		    TLorentzVector lambdaP4=(tl1+tl2);
		    cout <<endl<<" found lambda with mass: "<< lambdaP4.M()<<endl;
		    //get xF of lambda:
		    lambdaP4.Boost(boostVect);
		    if(lambdaP4.Pz()/W>0)
		      histo->Fill(lambdaP4.M());
		  }
		
	      }
	  }
      }
  
}


void getLambdasInJets(TClonesArray* branchParticles, TClonesArray* branchTracks, TClonesArray* jets, TH1* histo,TVector3 boostVect, double W)
{
  //  cout <<" we have " << jets->GetEntries() <<" jets " <<endl;

  for(int i=0;i<jets->GetEntries();i++)
    {

      Jet* jet=(Jet*)jets->At(i);
      {
	for(int j=0;j<jet->Constituents.GetEntries();j++)
	  {
	    TObject* object = jet->Constituents.At(j);
	    // Check if the constituent is accessible
	    if(object == 0)
	      {
		continue;
	      }

	    //could also be a tower
	    if(object->IsA() == Track::Class())
	      {
		Track* track1 = (Track*) object;
		if(track1->PT<0.100) continue;  // if track has less than 200 MeV, then do not consider

		///-->use also pions outside of tracks..
				for(int k=0;k<branchTracks->GetEntries();k++)
		//		  {
		//		for(int k=0;k<jet->Constituents.GetEntries();k++)
		  {
		    /////		 TObject* object2 = jet->Constituents.At(k);
		 ////// Check if the constituent is accessible
		 ////if(object2 == 0)
		 ////  {
		 ////	continue;
		 ////  }
		    //could also be a tower
		    ////		    		    if(object2->IsA() == Track::Class())
		      {
			////						Track* track2 = (Track*) object2;
		   
		    ////		    {
		      //
			Track* track2=(Track*)branchTracks->At(k);

			
		      //	if(track2->PT<0.100) continue;  // if track hasless than 200 MeV, then do not consider
		      GenParticle* part1=(GenParticle*)track1->Particle.GetObject();
		      GenParticle* part2=(GenParticle*)track2->Particle.GetObject();

		      int part1Mother=part1->M1;
		      int part2Mother=part2->M1;
		      int  numP=branchParticles->GetEntries();
		      int m1PID=0;
		      int m2PID=0;
		      GenParticle* gen1;
		      GenParticle* gen2;

		      if(track1->PID==2212 && track2->PID==-211)
			{
			  if(numP>part1Mother && part1Mother >=0)
			    {
			      gen1=((GenParticle*)branchParticles->At(part1Mother));
			      //  cout <<"first mother: "<< mParticle->PID<<endl;
			      m1PID=gen1->PID;
			    }
			  if(numP>part2Mother && part2Mother>=0)
			    {
			      gen2=((GenParticle*)branchParticles->At(part2Mother));
			      m2PID=gen2->PID;
			      //		    cout <<"second mother: "<< mParticle->PID<<endl;
			    }

			  //			    if(m1PID==3122 && m2PID==3122 && (part1Mother==part2Mother))
			  {
			    if(m1PID==3122)
			      cout << "found proton " << endl;
			    else
			      cout <<"found pion " <<endl;
			    //				cout <<"first px : "<< t1->P4().Px() <<" py: " << t1->P4().Py() <<" pz: "<< t1->P4().Pz() <<" E: "<< t1->P4().E()<<endl;
			    //				cout <<"second px : "<< t2->P4().//Px() <<" py: " << t2->P4().Py() <<" pz: "<< t2->P4().Pz() <<" E: "<< t2->P4().E()<<endl;
			    //				cout <<" t1 mass: "<< t1->P4().M() << " t2 mass: "<< t2->P4().M()<<endl;
			    //				cout <<"vertex t1: "<< t1->X <<", "<<t1->Y <<", " << t1->Z<<endl;
			    cout <<"vertex t2: "<< track2->X <<", "<<track2->Y <<", " << track2->Z<<endl;
			    cout <<"vertex p1: "<< part1->X <<", "<<part1->Y <<", " << part1->Z<<endl;
			    cout <<"vertex p2: "<< part2->X <<", "<<part2->Y <<", " << part2->Z<<endl;

			    TLorentzVector tl1=track1->P4();
			    TLorentzVector tl2=track2->P4();
			    double mP=0.93828;
			    double mPi=0.139568;
			    tl1.SetE(sqrt(tl1.P()*tl1.P()+mP*mP));
			    tl2.SetE(sqrt(tl2.P()*tl2.P()+mPi*mPi));
		    
			    TLorentzVector lambdaP4=(tl1+tl2);
			    cout <<endl<<" found lambda with mass: "<< lambdaP4.M()<<endl;
			    histo->Fill(lambdaP4.M());
			  }

			  //			    cout <<"found pair " <<endl;
			  TLorentzVector lambdaP4=track1->P4()+track2->P4();
			  lambdaP4.Boost(boostVect);
			  if(lambdaP4.Pz()/W>0)
			    histo->Fill(lambdaP4.M());
			  //			    cout <<endl<<" foundlambda in jet with mass: "<< lambdaP4.M()<<endl;
			}
		    }

		  }
	      }
  
	  }
      }
    }
}
