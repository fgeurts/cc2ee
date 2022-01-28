#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TChain.h"
#include "TSystem.h"
#include "TMath.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TBenchmark.h"
#include "TRandom.h"
#include "TPythia6.h"
#include "TMCParticle.h"
#include <iostream>
#include <vector>
#include <string>
#include "histogram.h"
//#include "TPythia6/TPythia6.h"
//#include "TPythia6/TMCParticle6.h"
using std::cout;
using std::endl;

#define eMass 0.000511

void genPythia(Int_t irun=1, Int_t Nevt=1000, Int_t iseed=789456, Bool_t debug=1)
{
	TStopwatch*   stopWatch = new TStopwatch();
	stopWatch->Start();
	gSystem->Load("libPhysics");
	//gSystem->Load("libEG.so");
	gSystem->Load("libEGPythia6.so");
	gSystem->Load("libPythia6.so");
	//gSystem->Load("TPythia6/TPythia6.so");

	//========================= initial trees ==========================
	TTree *eTree = new TTree("meTree","eTree");
	eTree->SetAutoSave(1000000);

	cout << "Initialize the eTree ... " << endl;

	eTree->Branch("ncString",     &meTree.ncString,    "ncString/I");
	eTree->Branch("ncbarString",  &meTree.ncbarString, "ncbarString/I");

	eTree->Branch("nc",           &meTree.nc,          "nc/I");
	eTree->Branch("isStringC",    meTree.isStringC,    "isStringC[nc]/O");
	eTree->Branch("cPt",          meTree.cPt,          "cPt[nc]/D");
	eTree->Branch("cEta",         meTree.cEta,         "cEta[nc]/D");
	eTree->Branch("cPhi",         meTree.cPhi,         "cPhi[nc]/D");
	eTree->Branch("cY",           meTree.cY,           "cY[nc]/D");

	eTree->Branch("ncbar",        &meTree.ncbar,       "ncbar/I");
	eTree->Branch("isStringCbar", meTree.isStringCbar, "isStringCbar[ncbar]/O");
	eTree->Branch("cbarPt",       meTree.cbarPt,       "cbarPt[ncbar]/D");
	eTree->Branch("cbarEta",      meTree.cbarEta,      "cbarEta[ncbar]/D");
	eTree->Branch("cbarPhi",      meTree.cbarPhi,      "cbarPhi[ncbar]/D");
	eTree->Branch("cbarY",        meTree.cbarY,        "cbarY[ncbar]/D");

	eTree->Branch("nePos",        &meTree.nePos,       "nePos/I");
	eTree->Branch("ePosParentGID",meTree.ePosParentGID,"ePosParentGID[nePos]/I");
	eTree->Branch("ePosPt",       meTree.ePosPt,       "ePosPt[nePos]/D");
	eTree->Branch("ePosEta",      meTree.ePosEta,      "ePosEta[nePos]/D");
	eTree->Branch("ePosPhi",      meTree.ePosPhi,      "ePosPhi[nePos]/D");

	eTree->Branch("neNeg",        &meTree.neNeg,       "neNeg/I");
	eTree->Branch("eNegParentGID",meTree.eNegParentGID,"eNegParentGID[neNeg]/I");
	eTree->Branch("eNegPt",       meTree.eNegPt,       "eNegPt[neNeg]/D");
	eTree->Branch("eNegEta",      meTree.eNegEta,      "eNegEta[neNeg]/D");
	eTree->Branch("eNegPhi",      meTree.eNegPhi,      "eNegPhi[neNeg]/D");

	TH1D *hnStats = new TH1D("hnStats", "hnStats", 10, 0, 10);
	hnStats->GetXaxis()->SetBinLabel(1, "nMBEvts");
	hnStats->GetXaxis()->SetBinLabel(2, "nCEvts");
	hnStats->GetXaxis()->SetBinLabel(3, "nStringCEvts");
	hnStats->GetXaxis()->SetBinLabel(4, "nCbarEvts");
	hnStats->GetXaxis()->SetBinLabel(5, "nStringCbarEvts");
	hnStats->GetXaxis()->SetBinLabel(6, "nTwoStringEvts");
	hnStats->SetLabelSize(0.025,"X");
	hnStats->SetLabelFont(62,   "X");
	hnStats->LabelsOption("d",  "X");

	//======================== initial parameters ==========================
	//Initialize Pythia
	TPythia6 *myPythia6 = new TPythia6();

	//myPythia6->SetMSEL(1); //call PYGIVE('msel = 1')   ! pp min. bias. 
	//myPythia6->SetMSEL(4);//c trigger
	//myPythia6->SetMSEL(5);//b trigger

	//PRC 92, 024912 (2015) tune - use minimum bias trigger and this setting can match published charmed-meson spectrum in pp collision
	myPythia6->SetMSEL(1); //call PYGIVE('msel = 1')   ! pp min. bias. 
	myPythia6->SetPARP(91,1.0);//<kt>
	myPythia6->SetPARP(67,1.0);//mstp32*4 high pT tuned parameter

	//CDFA tuned parameters: STAR default tune
	//myPythia6->SetMSEL(4);//c trigger
	//myPythia6->SetMSTP(51,7);
	//myPythia6->SetMSTP(82,4);
	//myPythia6->SetPARP(82,2.0);
	//myPythia6->SetPARP(83,0.5);
	//myPythia6->SetPARP(84,0.4);
	//myPythia6->SetPARP(85,0.9);
	//myPythia6->SetPARP(86,0.95);
	//myPythia6->SetPARP(89,1800);
	//myPythia6->SetPARP(90,0.25);
	//myPythia6->SetPARP(91,1.0);
	//myPythia6->SetPARP(67,4.0);

	////Xin's tune: for charm trigger
	//myPythia6->SetMSEL(4);//c trigger
	//myPythia6->SetMSTP(51,7);
	//myPythia6->SetMSTP(82,4);
	//myPythia6->SetPARP(82,2.0);
	//myPythia6->SetPARP(83,0.5);
	//myPythia6->SetPARP(84,0.4);
	//myPythia6->SetPARP(85,0.9);
	//myPythia6->SetPARP(86,0.95);
	//myPythia6->SetPARP(89,1800);
	//myPythia6->SetPARP(90,0.25);  
	//myPythia6->SetPARP(91,1.);//<kt>
	//myPythia6->SetPARP(67,1);//mstp32*4 high pT tuned parameter

	//====================== particle decay mode ==========================
	//switch off non-electron decays
	//---------------------------------------------------------------------
	myPythia6->SetMDME(818,1,0); //Ds
	for(int i=684;  i<=735;  i++) myPythia6->SetMDME(i,1,0); //D+
	for(int i=755;  i<=807;  i++) myPythia6->SetMDME(i,1,0); //D0
	for(int i=824;  i<=850;  i++) myPythia6->SetMDME(i,1,0); //Ds
	for(int i=857;  i<=862;  i++) myPythia6->SetMDME(i,1,0); //eta_c,J/psi,chi_2c
	for(int i=1097; i<=1165; i++) myPythia6->SetMDME(i,1,0); //Lc
	//====================== particle decay mode ==========================

	//============================ initial run =============================
	myPythia6->SetMRPY(1,iseed);
	myPythia6->Initialize("CMS","p","p",62.4);
	//myPythia6->Initialize("CMS","p","p",200);
	myPythia6->Pystat(2);
	//============================ initial run =============================
	
	
	//============================ run events ==============================
	Int_t nCharmEvents = 0;
	for(Int_t i = 1; i<=Nevt; i++)
	{
		myPythia6->GenerateEvent();

		hnStats->Fill(0.5);

		TObjArray *particles = myPythia6->GetListOfParticles();
		Int_t     nParticles = particles->GetEntries();

		if(i%10000==0)
		{
			cout<<"Woring on "<<i<<"-th event ..."<<endl;
		}

		memset(&meTree, 0, sizeof(meTree));

		Int_t nc          = 0;
		Int_t ncbar       = 0;
		Int_t ncString    = 0;
		Int_t ncbarString = 0;
		Int_t nePos       = 0;
		Int_t neNeg       = 0;

		for(Int_t l=0; l<nParticles; l++)
		{
			TMCParticle *mParticle = (TMCParticle*)particles->At(l);
			if(!mParticle) continue;

			Int_t pid = mParticle->GetKF(); // particle id
			Int_t pks = mParticle->GetKS(); // particle status, stable?, 1 means final state

			if(pid==4 && pks==12) //c quark and status is incoming beam
			{
				Double_t cpx = mParticle->GetPx();
				Double_t cpy = mParticle->GetPy();
				Double_t cpz = mParticle->GetPz();
				Double_t cmass = mParticle->GetMass();
			
				TLorentzVector cFourMom(0,0,0,0);
				cFourMom.SetXYZM(cpx,cpy,cpz,cmass);
				
				meTree.cPt[nc]       = cFourMom.Pt();
				meTree.cEta[nc]      = cFourMom.PseudoRapidity();
				meTree.cPhi[nc]      = cFourMom.Phi();
				meTree.cY[nc]        = cFourMom.Rapidity();
				meTree.isStringC[nc] = kFALSE;;

				Int_t firstChildIdx = mParticle->GetFirstChild() - 1;
				Int_t lastChildIdx  = mParticle->GetLastChild()  - 1;

				for(Int_t j=firstChildIdx; j<=lastChildIdx; j++)
				{
					TMCParticle *mString = (TMCParticle*)particles->At(j);
					if(mString)
					{
						Int_t stringId = mString->GetKF();
						
						if(stringId == 92)
						{
							ncString++;
							meTree.isStringC[nc] = kTRUE;
						}
					}
				}

				nc++;
			}

			if(pid==-4 && pks==12)//cbar quark and status is incoming beam
			{
				Double_t cpx   = mParticle->GetPx();
				Double_t cpy   = mParticle->GetPy();
				Double_t cpz   = mParticle->GetPz();
				Double_t cmass = mParticle->GetMass();
				
				TLorentzVector cbarFourMom(0,0,0,0);
				cbarFourMom.SetXYZM(cpx,cpy,cpz,cmass);
				
				meTree.cbarPt[ncbar]       = cbarFourMom.Pt();
				meTree.cbarEta[ncbar]      = cbarFourMom.PseudoRapidity();
				meTree.cbarPhi[ncbar]      = cbarFourMom.Phi();
				meTree.cbarY[ncbar]        = cbarFourMom.Rapidity();
				meTree.isStringCbar[ncbar] = kFALSE;

				Int_t firstChildIdx = mParticle->GetFirstChild() - 1;
				Int_t lastChildIdx  = mParticle->GetLastChild()  - 1;

				for(Int_t j=firstChildIdx; j<=lastChildIdx; j++)
				{
					TMCParticle *mString = (TMCParticle*)particles->At(j);
					if(mString)
					{
						Int_t stringId = mString->GetKF();
						
						if(stringId == 92)
						{
							ncbarString++;
							meTree.isStringCbar[ncbar] = kTRUE;
						}
					}
				}

				ncbar++;
			}

			if(TMath::Abs(pid) != 11 || pks != 1) continue; //the final particle must be e+: -11, or  e-: 11

			Int_t eParentIdx = mParticle->GetParent()-1;
			TMCParticle *meParent = (TMCParticle*)particles->At(eParentIdx);
			if(!meParent) continue;

			Int_t eParentId = TMath::Abs(meParent->GetKF());
			
			//only let D, D0. Ds, Lambdc pass here
			if(
					eParentId != 411 &&
					eParentId != 421 &&
					eParentId != 431 &&
					eParentId != 4122
			  )  continue;

			Double_t px = mParticle->GetPx();
			Double_t py = mParticle->GetPy();
			Double_t pz = mParticle->GetPz();

			TVector3 Mom(px,py,pz);
			Double_t pt  = Mom.Perp();
			Double_t eta = Mom.PseudoRapidity();
			Double_t phi = Mom.Phi();

			if(pid == -11) 
			{
				meTree.ePosParentGID[nePos] = eParentId;
				meTree.ePosPt[nePos]        = pt;
				meTree.ePosEta[nePos]       = eta;
				meTree.ePosPhi[nePos]       = phi;
				nePos++;
			}

			if(pid == 11)
			{
				meTree.eNegParentGID[neNeg] = eParentId;
				meTree.eNegPt[neNeg]        = pt;
				meTree.eNegEta[neNeg]       = eta;
				meTree.eNegPhi[neNeg]       = phi;
				neNeg++;
			}
		}//loop all particles

		if(nc>0)                          hnStats->Fill(1.5);
		if(ncString>0)                    hnStats->Fill(2.5);
		if(ncbar>0)                       hnStats->Fill(3.5);
		if(ncbarString>0)                 hnStats->Fill(4.5);
		if(ncString==1 && ncbarString==1) hnStats->Fill(5.5);

		meTree.ncString    = ncString;
		meTree.ncbarString = ncbarString;
		meTree.nc          = nc;
		meTree.ncbar       = ncbar;
		meTree.nePos       = nePos;
		meTree.neNeg       = neNeg;

		if ( (nc+ncbar) == 0) continue;

		if(debug)
		{
			cout<<endl;
			cout<<"Woring on "<<i<<"-th event ..."<<endl;
			myPythia6->Pylist(1);
		}
		else if(nCharmEvents<10) myPythia6->Pylist(1);

		nCharmEvents++;

		eTree->Fill();
	}//event

	//============================= output ==============================
	char rootfilename[100];
	if(debug)
	{
		sprintf(rootfilename, "outtest/pythiaevent_test.root");
	}
	else
	{
		sprintf(rootfilename,"output/pythiaevent%d.root",irun);
	}

	TFile* file = new TFile(rootfilename,"RECREATE");
	file->cd();
	hnStats->Write();
	eTree->Write();
	file->Close();

	myPythia6->Pystat(1);

	stopWatch->Stop();
	stopWatch->Print();
}
