#include <math.h>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <iostream>
#include <fstream>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TDatabasePDG.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"

#include "myexampleAnalysis.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/PseudoJet.hh"  
#include "fastjet/tools/Filter.hh"
#include "fastjet/Selector.hh"

#include "fastjet/tools/Recluster.hh"

#include "fastjet/contrib/SoftDrop.hh"

#include "Pythia8/Pythia.h"

using namespace std;

// Constructor 
myexampleAnalysis::myexampleAnalysis(){

	//0.8AK jets
	jetRad = 0.8;
	m_jet_def	= new fastjet::JetDefinition(fastjet::antikt_algorithm, jetRad);
	
	//Constants for use in our algorithm
	//Following paper http://xxx.lanl.gov/pdf/0802.2470v2
	mu = 0.67;
	ycut = 0.15;
}

// Destructor 
myexampleAnalysis::~myexampleAnalysis(){
	delete m_jet_def;
}

// End
void myexampleAnalysis::End(){

	cout << "The number that were cut by PTMISS is " << numPTMISS << endl;
	cout << "The number that were cut by PTCUT is " << numPTCUT << endl;
	cout << "The number that were cut by BTAG is " << numBTAG << endl;
	cout << "The number that were cut by ISOLATED is " << numISOLATED << endl;
	cout << "The number that were cut by RHO is " << numRHO << endl;
	cout << "The number that were cut by N2 is " << numN2 << endl;

	return;
}

fastjet::ClusterSequence *largeRClusterSequence; 
std::vector<fastjet::PseudoJet> *particlesForJets; 
fastjet::PseudoJet *p; 
std::vector <fastjet::PseudoJet> *largeRJets; 
std::vector <fastjet::PseudoJet> *bs;
std::vector <fastjet::PseudoJet> *pieces;
int idmod; //For tagging b-hadrons
float pxMiss = 0; //For missing momenta calculations
float pyMiss = 0;
float success = 0;
fastjet::PseudoJet trimmedJet;
vector <fastjet::PseudoJet> constit;
std::list<int> jetIndices;
std::list<int>::iterator it;
int cSize;
int totalB;
int thisId;

int failed = 0;

//Reculster with AK0.2 jets
fastjet::Recluster recluster_ca_inf = fastjet::Recluster(fastjet::antikt_algorithm, 0.2, fastjet::Recluster::keep_all);

void destroy(){
	delete largeRClusterSequence;
	delete largeRJets;
	delete particlesForJets;
	delete bs;
	delete pieces;
}

// Perform clustering and pre-selection. The first argument is the showered event. 
// The second argument is a target array to write higgs jet information to.
// This target must be a float[4]. It will contain [pt, eta, phi, m] after the method is executed.
void myexampleAnalysis::AnalyzeEvent(Pythia8::Pythia* pythia8, float *target){
	
	//cout << endl;
	//cout << endl;
	//cout << endl;
	//cout << "Event analysis started. Event ID " << ievt << endl;

	if (!pythia8->next()){
		cout << "Failed in Analysis" << endl;
		failed++;
		target[0] = 0; //pT
		target[1] = 0; //eta
		target[2] = 0; //phi
		target[3] = -10; //mass
		return;
	}
	else{
		//cout << "hm" << endl;
	}

	
	pxMiss = 0;
	pyMiss = 0;

	particlesForJets = new std::vector<fastjet::PseudoJet>;
	bs = new std::vector<fastjet::PseudoJet>;

	int ipHiggs=-1;
	// Particle loop
	for (unsigned int ip=0; ip<pythia8->event.size(); ip++){
		p = new fastjet::PseudoJet(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e());
		(*p).set_user_info(new MyUserInfo(pythia8->event[ip].id() , ip , pythia8->event[ip].charge()));

		idmod = abs(pythia8->event[ip].id())%10000;
		if (idmod == 5){
			(*bs).push_back(*p);
		}
		// Checking decay
		if (pythia8->event[ip].idAbs() == 23 || pythia8->event[ip].idAbs() == 24 || pythia8->event[ip].idAbs() == 25){
 		//print PDG ID of two dauguers
			//cout <<ip <<" "<< pythia8->event[ip].id() <<" status " << pythia8->event[ip].status()  <<" decays into "<< pythia8->event[ pythia8->event[ip].daughter1() ].id() <<" and "<< pythia8->event[ pythia8->event[ip].daughter2() ].id()  <<endl;  
			if (pythia8->event[ip].idAbs() == 25 && pythia8->event[ pythia8->event[ip].daughter1() ].id() != 25) {
				ipHiggs=ip;
				//cout << "Final state higgs found, ip = " << ipHiggs <<endl;

			}
		}
		if (!pythia8->event[ip].isFinal())	 continue;
		if (fabs(pythia8->event[ip].id())==12)	 continue; //Neutrinos
		if (fabs(pythia8->event[ip].id())==14)	 continue;
		if (fabs(pythia8->event[ip].id())==16)	 continue;

		pxMiss = pxMiss + pythia8->event[ip].px();
		pyMiss = pyMiss + pythia8->event[ip].py();

		(*particlesForJets).push_back((*p));		
		delete p;	
	} 
	//cout << "Event has " << (*bs).size()<< " b hadrons in total." << endl;


	largeRClusterSequence =  new fastjet::ClusterSequence (*particlesForJets, *m_jet_def);

	// Collection of R=0.8 jets
	largeRJets = new std::vector<fastjet::PseudoJet>;
	*largeRJets =  fastjet::sorted_by_pt((*largeRClusterSequence).inclusive_jets(0)); 
	
	//if too large pTMiss, then reject the event
	//if (pow(pxMiss,2)+pow(pyMiss,2) > 0){
	//	cout << "Missing pT detected! " << pow(pow(pxMiss,2)+pow(pyMiss,2),0.5) << endl;
	//	//cout << "Big pT Miss" << endl;
	//	numPTMISS++;
	//	return (0);
	//}

	
	pieces = new std::vector<fastjet::PseudoJet>;

	fastjet::contrib::SoftDrop sd(0.0,0.1);

	int leadingCentralJetId = 0;
	
	//cout << "# large R jet = " << (*largeRJets).size() << endl;
	for (int iJet = 0; iJet < (*largeRJets).size(); iJet++){
		// Count b in every 0.8 jets
		totalB = 0;
		trimmedJet = sd((*largeRJets)[iJet]);
		fastjet::PseudoJet rec_jet_ca_inf = recluster_ca_inf(trimmedJet);
		*pieces = rec_jet_ca_inf.pieces();
		for (int r = 0; r < (*pieces).size(); r++){
			for (int rr = 0; rr < (*bs).size(); rr++){
				if ((*bs)[rr].delta_R((*pieces)[r]) < 0.2){
					totalB++;
					break;
				}
			}
		}
		
		//cout << "iJet" << iJet << " #b= " << totalB 
                //     << ", (" << pow(pow(trimmedJet.px(), 2) + pow(trimmedJet.py(), 2), 0.5) 
                //     << ", " << trimmedJet.eta() 
                //     << ", " << trimmedJet.phi()
		//     << ", " << trimmedJet.m() 
		//     << ") d eta=" << pythia8->event[ipHiggs].eta() - trimmedJet.eta() 
		//     << ", d phi=" << pythia8->event[ipHiggs].phi() - trimmedJet.phi() << endl;

		// Check if jet is >2 b-tagged and central jet		
		if (totalB >= 2 && (*largeRJets)[iJet].eta() < 2 && (*largeRJets)[iJet].eta() > -2){
			leadingCentralJetId = iJet;
			cout << "Leading central jet chosen with eta=" << (*largeRJets)[iJet].eta() << endl;
			break;
		}else if(iJet == (*largeRJets).size()-1){
			cout << "No leading central b-tagged jet!" << endl;
			numBTAG++;
			destroy();
			target[0] = 0; //pT
			target[1] = 0; //eta
			target[2] = 0; //phi
			target[3] = 0; //mass
			return;
		}			
	}

	trimmedJet = sd((*largeRJets)[leadingCentralJetId]);
	fastjet::PseudoJet rec_jet_ca_inf = recluster_ca_inf(trimmedJet);
	*pieces = trimmedJet.pieces();
	constit = trimmedJet.constituents();
	cSize = constit.size();
	
	// See discarded-b-tagging1.txt
	// pt cut - For now on the base jet (why am I not using trimmed?)
	// if ((*largeRJets)[0].pt() < 400){
	if ((*largeRJets)[leadingCentralJetId].pt() < 450){
		numPTCUT++;
		destroy();
		target[0] = 0; //pT
		target[1] = 0; //eta
		target[2] = 0; //phi
		target[3] = 0; //mass
		return;
	}	


	
	// B-tagging is put before leading jet selection
	//if (totalB != 2){
	//	numBTAG++;
	//	destroy();
	//	return(0);
	//}




	thisId = fabs(constit[0].user_info<MyUserInfo>().pdg_id());

	//Here, we want to check over all the jets
	//And see which are just a single particle, and check if they are isolated or not

	if (cSize == 1){
	if (thisId == 11){
		//we have an ELECTRON 
		if (constit[0].pt() > 10 && fabs(constit[0].eta()) < 2.5){ 
			//we have an isolated electron/muon 
			//cout << "Particle Level : Isolated Electron" << endl;
			numISOLATED++;
			destroy();
			target[0] = 0; //pT
			target[1] = 0; //eta
			target[2] = 0; //phi
			target[3] = 0; //mass
			return;
		}
	}
	if (thisId == 13){
		//we have an muon
		if (constit[0].pt() > 10 && fabs(constit[0].eta()) < 2.4){
			//we have an isolated electron/muon
			//cout << "Particle Level : Isolated Muon" << endl;
			numISOLATED++;
			destroy();
			target[0] = 0; //pT
			target[1] = 0; //eta
			target[2] = 0; //phi
			target[3] = 0; //mass
			return;
		}
	}
	if (thisId ==15){
		//we have a tau
		if (constit[0].pt() > 18 && fabs(constit[0].eta()) < 2.3){
			//we have an isolated tau
			//cout << "Particle Level : Isolated Tau" << endl;
			numISOLATED++;
			destroy();
			target[0] = 0; //pT
			target[1] = 0; //eta
			target[2] = 0; //phi
			target[3] = 0; //mass
			return;
		}
	}
	}

	//Now we do a cut on rho, but we need to get the softdrop working first

	rho = log(pow(trimmedJet.m(),2)/pow((*largeRJets)[leadingCentralJetId].pt(),2)); //Use the ungroomed pt
	if (rho <= -6){
		numRHO++;
		destroy();
		target[0] = 0; //pT
		target[1] = 0; //eta
		target[2] = 0; //phi
		target[3] = 0; //mass
		return;
	}
	if (rho >= -2.1){
		numRHO++;
		destroy();
		target[0] = 0; //pT
		target[1] = 0; //eta
		target[2] = 0; //phi
		target[3] = 0; //mass
		return;
	}

	//Now, lets actually just implement an N12 and see what kind of shenanigans we get
	//Here, beta is set to zero
	//https://arxiv.org/pdf/1609.07483.pdf
	
	float e2 = 0;
	float e3 = 0;
	float sumpt = 0;

	float rs = 0;
	float rt = 0;
	float st = 0;

	for (int r = 0; r < cSize; r++){
		sumpt = sumpt + constit[r].pt();
	for (int s = 0; s < cSize; s++){
		if (r == s) continue;	
		e2 = e2 + constit[r].pt()*constit[s].pt()*constit[r].delta_R(constit[s]);
	for (int t = 0; t < cSize; t++){
		if (s == t) continue;
		rs = constit[r].delta_R(constit[s]);
		rt = constit[r].delta_R(constit[t]);
		st = constit[s].delta_R(constit[t]);
		e3 = e3 + constit[r].pt()*constit[s].pt()*constit[t].pt()*std::min({rs*rt,rs*st,rt*st});
	}
	}
	}

	e2 = e2/(pow(sumpt,2));
	e3 = e3/(pow(sumpt,3));

	float N2 = e3/pow(e2,2);
	
	//N2 cut
	if (N2 > 0.45){
		numN2++;
		destroy();
		target[0] = 0; //pT
		target[1] = 0; //eta
		target[2] = 0; //phi
		target[3] = 0; //mass
		return;
	}
	
	/*
	delete largeRClusterSequence;
	delete largeRJets;
	delete particlesForJets;
	*/
	destroy();	

	//storing higgs jet information to target array
	target[0] = trimmedJet.pow(pow(trimmedJet.px(), 2) + pow(trimmedJet.py(), 2), 0.5); //pT
	target[1] = trimmedJet.eta(); //eta
	target[2] = trimmedJet.phi(); //phi
	target[3] = trimmedJet.m(); //mass
	cout << "higgs output test - analysis ";
	cout << target << endl;
}


