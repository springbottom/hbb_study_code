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

#include "fastjet/contrib/Recluster.hh"

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

fastjet::ClusterSequence *cs; 
std::vector<fastjet::PseudoJet> *particlesForJets; 
fastjet::PseudoJet *p; 
std::vector <fastjet::PseudoJet> *myJets; 
std::vector <fastjet::PseudoJet> *bhadrons;
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
fastjet::contrib::Recluster recluster_ca_inf = fastjet::contrib::Recluster(fastjet::antikt_algorithm, 0.2, false);

void destroy(){
	delete cs;
	delete myJets;
	delete particlesForJets;
	delete bhadrons;
	delete bs;
	delete pieces;
}

// Analyze
float myexampleAnalysis::AnalyzeEvent(int ievt, Pythia8::Pythia* pythia8){
	

	//cout << "call reached" << endl;

	if (!pythia8->next()){
		cout << "Failed in Analysis" << endl;
		failed++;
		return (-10);
	}
	else{
		//cout << "hm" << endl;
	}

	
	pxMiss = 0;
	pyMiss = 0;

	particlesForJets = new std::vector<fastjet::PseudoJet>;
	bhadrons = new std::vector<fastjet::PseudoJet>;
	bs = new std::vector<fastjet::PseudoJet>;


	// Particle loop
	for (unsigned int ip=0; ip<pythia8->event.size(); ip++){
		p = new fastjet::PseudoJet(pythia8->event[ip].px(), pythia8->event[ip].py(), pythia8->event[ip].pz(),pythia8->event[ip].e());
		(*p).set_user_info(new MyUserInfo(pythia8->event[ip].id() , ip , pythia8->event[ip].charge()));

		idmod = abs(pythia8->event[ip].id())%10000;
		if (idmod == 5){
			(*bs).push_back(*p);
		}
		if ((idmod >= 500 && idmod < 600) || (idmod >= 5000 && idmod < 6000)){
			(*bhadrons).push_back(*p);
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

	//cout << "Done Particle Loop" << endl;


	cs =  new fastjet::ClusterSequence (*particlesForJets, *m_jet_def);
	myJets = new std::vector<fastjet::PseudoJet>;
	*myJets =  fastjet::sorted_by_pt((*cs).inclusive_jets(0)); 
	
	//if too large pTMiss, then reject the event
	/*
	if (pow(pxMiss,2)+pow(pyMiss,2) > pow(140,2)){
		cout << pxMiss << " " << pyMiss << endl;
		//cout << "Big pT Miss" << endl;
		numPTMISS++;
		return (0);
	}
	*/
	
	pieces = new std::vector<fastjet::PseudoJet>;

	fastjet::contrib::SoftDrop sd(0.0,0.1);
	trimmedJet = sd((*myJets)[0]);
	*pieces = trimmedJet.pieces();
	constit = trimmedJet.constituents();
	cSize = constit.size();

		
	//Code for looking back in the tree
	/*
	int r = 0;
	totalB = 0;
	while (r < cSize){
		jetIndices.push_back(constit[r].user_info<MyUserInfo>().pythia_id());
		r++;
	}
	for (it = jetIndices.begin(); it != jetIndices.end(); ++it){
		//cout << jetIndices.size() << " " << *it << endl;
		//cout << *it <<" "<<pythia8->event[*it].mother1() <<" "<< pythia8->event[*it].mother2() <<" "<< pythia8->event[*it].id() << endl;
		if (*it == 0){
			continue;
		}
		if (jetIndices.end() == std::find(jetIndices.begin(), jetIndices.end(), pythia8->event[*it].mother1())){
			jetIndices.push_back(pythia8->event[*it].mother1());
		}
		if (jetIndices.end() == std::find(jetIndices.begin(), jetIndices.end(), pythia8->event[*it].mother2())){
			jetIndices.push_back(pythia8->event[*it].mother2());
		}
		if ((pythia8->event[*it].id())==5){
			if (totalB == 0){
				totalB = 1;
			}
			if (totalB == -1){
				totalB = 2;
			}
		}
		else if ((pythia8->event[*it].id()) == -5){
			if (totalB == 0){
				totalB = -1;
			}
			if (totalB == 1){
				totalB = 2;
			}
		}
	}
	jetIndices.clear();
	if (totalB == 2){
		//pythia8->event.list();
	}
	*/
	
	/*
	r = 0;
	while (r < jetIndices.size()){
		if (jetIndices[r] == 0){
			r++;
			continue;
		}
		jetIndices.push_back(pythia8->event[jetIndices[r]].mother1());
		jetIndices.push_back(pythia8->event[jetIndices[r]].mother2());
		if (fabs(pythia8->event[jetIndices[r]].id() == 5){
			totalB++;
		}
	}
	*/
	
	
	/*
	totalB = 0;
	for (int r = 0; r < (*bhadrons).size(); r++){
		if ((*bhadrons)[r].pt() < 5) continue;
		if ((*bhadrons)[r].delta_R(trimmedJet) < jetRad){
			totalB++;
		}
	}
	*/

	//Code for looking at the pieces
	
	/*
	totalB = 0;
	for (int r = 0; r < (*pieces).size(); r++){
	for (int rr = 0; rr < (*bhadrons).size(); rr++){
		if ((*bhadrons)[rr].pt() < 5) continue;
		//if ((*bhadrons)[rr].delta_R((*pieces)[r]) < 0.001){
		//	totalB++;
		//}
		if ((*bhadrons)[rr].delta_R((*pieces)[r]) < 0.2*(*bhadrons)[rr].delta_R((*pieces)[(r+1)%2])){
			totalB++;
		}
	}
	}
	*/

	//cout << "about to B" << endl;
	
	
	
	
	totalB = 0;
	fastjet::PseudoJet rec_jet_ca_inf = recluster_ca_inf(trimmedJet);
	*pieces = rec_jet_ca_inf.pieces();
	for (int r = 0; r < (*pieces).size(); r++){
		/*
		if ((*pieces)[r].pt() < 100){
			continue;
		}
		*/
		
		/*
		for (int rr = 0; rr < (*bhadrons).size(); rr++){
			if ((*bhadrons)[rr].delta_R((*pieces)[r]) < 0.2){
				for (int rrr = 0; rrr < (*pieces).size(); rrr++){
					if (rrr == r) continue;
					if ((*bhadrons)[rr].delta_R((*pieces)[rrr]) < 0.2){
						totalB--;
						break;
					}
				}
				totalB++;
				break;
			}
		}
		*/
		
		
		for (int rr = 0; rr < (*bs).size(); rr++){
			if ((*bs)[rr].delta_R((*pieces)[r]) < 0.2){
				totalB++;
				break;
			}
		}
		
	}

	/*
	for (int r = 0; r < (*pieces).size(); r++){
	for (int rr = 0; rr < (*bhadrons).size(); rr++){
		if ((*bhadrons)[rr].delta_R((*pieces)[r]) < 0.1){
			continue;
		}
	for (int rrr = 0; rrr < (*pieces).size(); rrr++){
		if ((*bhadrons)[rr].delta_R((*pieces)[r]) > (*bhadrons)[rr].delta_R((*pieces)[rrr])){
			break;
		}
		if (rrr == (*pieces).size() - 1){
			totalB++;
		}
		if (rrr == r){continue;}
		
	}}}*/
	
		
	//pt cut - For now on the base jet (why am I not using trimmed?)
	//if ((*myJets)[0].pt() < 400){
	if ((*myJets)[0].pt() < 450){
		numPTCUT++;
		destroy();
		return (0);
	}	


	
	
	if (totalB != 2){
		numBTAG++;
		destroy();
		return(0);
	}




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
			return (0);
		}
	}
	if (thisId == 13){
		//we have an muon
		if (constit[0].pt() > 10 && fabs(constit[0].eta()) < 2.4){
			//we have an isolated electron/muon
			//cout << "Particle Level : Isolated Muon" << endl;
			numISOLATED++;
			destroy();
			return (0);
		}
	}
	if (thisId ==15){
		//we have a tau
		if (constit[0].pt() > 18 && fabs(constit[0].eta()) < 2.3){
			//we have an isolated tau
			//cout << "Particle Level : Isolated Tau" << endl;
			numISOLATED++;
			destroy();
			return (0);
		}
	}
	}

	//Now we do a cut on rho, but we need to get the softdrop working first

	rho = log(pow(trimmedJet.m(),2)/pow((*myJets)[0].pt(),2)); //Use the ungroomed pt
	if (rho <= -6){
		numRHO++;
		destroy();
		return (0);
	}
	if (rho >= -2.1){
		numRHO++;
		destroy();
		return (0);
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
	
	
	if (N2 > 0.45){
		//Here we cut
		numN2++;
		destroy();
		return(0);
	}
	
	

	success = trimmedJet.m();


	/*
	delete cs;
	delete myJets;
	delete particlesForJets;
	delete bhadrons;
	*/
	destroy();	
	return (success);
}


