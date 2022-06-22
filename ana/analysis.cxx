#define Mkaon 0.493677
#define Mproton 0.93827231
#define M_el 0.00051099907

#include <iostream>
#include <fstream>
#include <cstdlib>
#include "sys/types.h"
#include "dirent.h"

#include "math.h"
#include "string.h"
//Add the data structure
#include <iomanip>
#include "meTree.h"
//#include "StRefMultCorr.h"
//#include "Histograms.h"
//#include <stdio.h>

#ifndef __CINT__
#include "TROOT.h"
#include "TFile.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TH3.h"
#include "TF1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TVector2.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TUnixSystem.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using std::cout;
using std::endl;
using std::setw;
#endif


//fg Mass,pT bin ranges copied from RiceU-HeavyIons/eeAuAu27/charming.cpp
const Int_t MassBins = 47;
const Double_t LowEdges[48]={0, 0.004, 0.008, 0.012, 0.016, 0.02, 0.024, 0.028, 0.032, 0.036, 0.04, 0.044, 0.048, 0.052, 0.056, 0.06, 0.064, 0.068, 0.072, 0.076, 0.08, 0.084, 0.088, 0.092, 0.096, 0.1, 0.25, 0.4, 0.5164, 0.6332, 0.75, 0.7716, 0.7832, 0.8072, 0.9958, 1.0108, 1.0258, 1.0408, 1.2008, 1.5508, 3.011, 3.081, 3.097, 3.113, 3.129, 3.145, 3.225, 3.305};
const Int_t HpTBins = 12;
const Double_t HlowEdgespT[13]={0.0,0.2,0.4,0.6,0.8,1.0,1.4,1.8,2.2,2.6,3.0,3.5,4.0};
//-fg


double PI;
Double_t GetBR(Int_t id);
Double_t pTSmear(double *x, double *p);
Double_t DoubleCrystalBall(double *x, double *p);

int main(int argc, char** argv){
    PI = TMath::Pi();

  // process command line arguments
  if(argc!=1 && argc!=4) {
    cout << "Usage: " << argv[0] << " <energy (27,39,64)> <input.list> <output.root>" << endl;
    return -1; 
  }
  // defaults
    int ENERGY = 27; //  27, 39, 62
    char *inFile = (char*)"test.list";
    char outFile[100];
    sprintf(outFile,"test.root");

    if(argc==4){
      ENERGY = atoi(argv[1]);
	inFile = argv[2];
	sprintf(outFile,"%s",argv[3]);
    }

    if (ENERGY!=27 && ENERGY!=39 && ENERGY!=62){
      cout << "Energy " << ENERGY << " is not supported." << endl;
      return -1;
    }
    cout << "Applying run specific parameters for beam energy " << ENERGY << " GeV" << endl;

    TFile *fout = new TFile(outFile,"recreate");
    Int_t NPX = 10000;
    TF1* fDCB = new TF1("fDCB",DoubleCrystalBall,-1.,1.,7);

    // Set run-specific DCB and smearing parameters
    if(ENERGY==39)
      fDCB->SetParameters(1.746,1.673,1.531,8.267,5.3e-5,0.0092,1);  // source:  STAR Analysis Note PSN0656 (fig. 26, page 21)
    else if(ENERGY==62)
      fDCB->SetParameters(1.727,1.665,1.571,7.839,-2.5e-5,0.0092,1); // source:  STAR Analysis Note PSN0656 (fig. 27, page 22)
    else // ENERGY==27
      fDCB->SetParameters(1.812, 2.145, 1.224, 4.325, -3.3E-4, 0.0093, 1.); // source: STAR Analysis Note PSN0656 (fig. 40, page 31)
    fDCB->SetNpx(NPX);
    TH1D *hDoubleCrystalBall = (TH1D*) fDCB->GetHistogram();

    TF1 *fpTSmear = new TF1("fpTSmear",pTSmear,0.,10.,2);//15.,0); 
    if(ENERGY==39)
      fpTSmear->SetParameters(0.001973, 0.008535); // source:  STAR Analysis Note PSN0656 (fig. 26, page 21)
    else if(ENERGY==62)
      fpTSmear->SetParameters(0.002296, 0.008505); // source:  STAR Analysis Note PSN0656 (fig. 27, page 22)
    else // ENERGY==27
      fpTSmear->SetParameters(0.003011, 0.007934); // source: STAR Analysis Note PSN0656 (fig. 39, page 31)
    fpTSmear->SetNpx(NPX);//Added 1/26/15

    // Histograms
    TH1D *hCounter = new TH1D("hCounter","hCounter",10,0,10);
    hCounter->GetXaxis()->SetBinLabel(1, "nMBEvts");
    hCounter->GetXaxis()->SetBinLabel(2, "nCEvts");
    hCounter->GetXaxis()->SetBinLabel(3, "nStringCEvts");
    hCounter->GetXaxis()->SetBinLabel(4, "nCbarEvts");
    hCounter->GetXaxis()->SetBinLabel(5, "nStringCbarEvts");
    hCounter->GetXaxis()->SetBinLabel(6, "nTwoStringEvts");

    TH2D *h0pTInvMassAccccbar      = new TH2D("h0pTInvMassAccccbar"      , "No smearing ccbar acc"   , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *h0pTInvMassccbar         = new TH2D("h0pTInvMassccbar"         , "No smearing ccbar "      , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *hpTInvMassAccccbar       = new TH2D("hpTInvMassAccccbar"       , "Smearing ccbar acc"      , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *hpTInvMassccbar          = new TH2D("hpTInvMassccbar"          , "Smearing ccbar "         , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *hpTInvMassRmDPhiAccccbar = new TH2D("hpTInvMassRmDPhiAccccbar" , "RmDPhi ccbar acc"        , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *hpTInvMassRmDPhiccbar    = new TH2D("hpTInvMassRmDPhiccbar"    , "RmDPhi ccbar "           , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *hpTInvMassRmPhiAccccbar  = new TH2D("hpTInvMassRmPhiAccccbar"  , "RmPhi ccbar acc"         , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *hpTInvMassRmPhiccbar     = new TH2D("hpTInvMassRmPhiccbar"     , "RmPhi ccbar "            , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *hpTInvMassRmEtaAccccbar  = new TH2D("hpTInvMassRmEtaAccccbar"  , "RmEta and Phi ccbar acc" , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *hpTInvMassRmEtaccbar     = new TH2D("hpTInvMassRmEtaccbar"     , "RmEta and Phi ccbar "    , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *hpTInvMassRmAllAccccbar  = new TH2D("hpTInvMassRmAllAccccbar"  , "RmAll ccbar acc"         , 1800 , 0. , 3.6 , 200 , 0. , 5.);
    TH2D *hpTInvMassRmAllccbar     = new TH2D("hpTInvMassRmAllccbar"     , "RmAll ccbar "            , 1800 , 0. , 3.6 , 200 , 0. , 5.);
 //fg - Define LE-based histograms
    TH2D *h0pTInvMassAccccbarLE      = new TH2D("h0pTInvMassAccccbarLE"      , "No smearing ccbar acc"   , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *h0pTInvMassccbarLE         = new TH2D("h0pTInvMassccbarLE"         , "No smearing ccbar "      , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *hpTInvMassAccccbarLE       = new TH2D("hpTInvMassAccccbarLE"       , "Smearing ccbar acc"      , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *hpTInvMassccbarLE          = new TH2D("hpTInvMassccbarLE"          , "Smearing ccbar "         , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *hpTInvMassRmDPhiAccccbarLE = new TH2D("hpTInvMassRmDPhiAccccbarLE" , "RmDPhi ccbar acc"        , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *hpTInvMassRmDPhiccbarLE    = new TH2D("hpTInvMassRmDPhiccbarLE"    , "RmDPhi ccbar "           , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *hpTInvMassRmPhiAccccbarLE  = new TH2D("hpTInvMassRmPhiAccccbarLE"  , "RmPhi ccbar acc"         , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *hpTInvMassRmPhiccbarLE     = new TH2D("hpTInvMassRmPhiccbarLE"     , "RmPhi ccbar "            , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *hpTInvMassRmEtaAccccbarLE  = new TH2D("hpTInvMassRmEtaAccccbarLE"  , "RmEta and Phi ccbar acc" , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *hpTInvMassRmEtaccbarLE     = new TH2D("hpTInvMassRmEtaccbarLE"     , "RmEta and Phi ccbar "    , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *hpTInvMassRmAllAccccbarLE  = new TH2D("hpTInvMassRmAllAccccbarLE"  , "RmAll ccbar acc"         , MassBins, LowEdges, HpTBins, HlowEdgespT);
    TH2D *hpTInvMassRmAllccbarLE     = new TH2D("hpTInvMassRmAllccbarLE"     , "RmAll ccbar "            , MassBins, LowEdges, HpTBins, HlowEdgespT);
 //-fg
    TH1D *hPosPt  = new TH1D("hPosPt", "hPosPt",500,0,5);
    TH1D *hPosEta = new TH1D("hPosEta","hPosEta",400,-4,4);
    TH1D *hPosPhi = new TH1D("hPosPhi","hPosPhi",360,-PI,PI);

    TH1D *hNegPt  = new TH1D("hNegPt", "hNegPt",500,0,5);
    TH1D *hNegEta = new TH1D("hNegEta","hNegEta",400,-4,4);
    TH1D *hNegPhi = new TH1D("hNegPhi","hNegPhi",360,-PI,PI);

    //---------------------------------------------------
    // open files and add to the chain
    //---------------------------------------------------
    TChain *chain = new TChain("meTree");

    int ifile = 0;
    char filename[512];
    ifstream *inputStream = new ifstream;
    inputStream->open(inFile);
    if(!(inputStream)){
	cout<<"Can't open list file!!"<<endl;
	return 0;
    }
    while(inputStream->good()){
	inputStream->getline(filename,512);
	if(inputStream->good()){
	    TFile *ftmp = new TFile(filename);
	    if(!ftmp||!(ftmp->IsOpen())||!(ftmp->GetNkeys())){
		cout<<"Can't open file or file has no key!!"<<endl;
	    }
	    else {
		cout<<"read in "<<ifile<<"th file: "<<filename<<endl;
		TH1F *htmp = (TH1F *)ftmp->Get("hnStats");
		hCounter->Add(htmp);
		chain->Add(filename);
		ifile++;
		delete htmp;
	    }
	    delete ftmp;
	}
    }
    delete inputStream;
    meTree *tree = new meTree(chain);

    //---------------------------------------------------------
    // loop events
    //---------------------------------------------------------
    int n = chain->GetEntries();
    cout<<n<<" events"<<endl;

    TRandom3 *rnu = new TRandom3();
    rnu->SetSeed(12345);

    cout<<"### get pt eta phi distribution"<<endl;
    for(int i=0;i<n;i++){
      if(i%100000==0) cout<<i<<" events"<<endl;
      chain->GetEntry(i);
      Int_t nePos = tree->nePos;
      Int_t neNeg = tree->neNeg;

      for (Int_t ii=0; ii<nePos; ii++){
	Float_t ePosPt = tree->ePosPt[ii];
	Float_t ePosPhi = tree->ePosPhi[ii];
	Float_t ePosEta = tree->ePosEta[ii];
	hPosPt->Fill(ePosPt);
	hPosEta->Fill(ePosEta);
	hPosPhi->Fill(ePosPhi);
      }

      for (Int_t ii=0; ii<neNeg; ii++){
	Float_t eNegPt = tree->eNegPt[ii];
	Float_t eNegPhi = tree->eNegPhi[ii];
	Float_t eNegEta = tree->eNegEta[ii];
	hNegPt->Fill(eNegPt);
	hNegEta->Fill(eNegEta);
	hNegPhi->Fill(eNegPhi);
      }
    }
    cout<<"### get pt eta phi distribution--- done!!!"<<endl;


    cout<<"### begin to get pairs"<<endl;
    for(int i=0;i<n;i++){
	if(i%100000==0) cout<<i<<" events"<<endl;
	chain->GetEntry(i);

	// get the number of electrons and positrons for this event
	Int_t nePos = tree->nePos;
	Int_t neNeg = tree->neNeg;

	for (Int_t ii=0; ii<nePos; ii++){
	  Float_t ePosPt = tree->ePosPt[ii];
	  Float_t ePosPhi = tree->ePosPhi[ii];
	  Float_t ePosEta = tree->ePosEta[ii];
	  Int_t   ePosParentGID = tree->ePosParentGID[ii];
	  double BRPos = GetBR(ePosParentGID);

	  for (Int_t jj=0; jj<neNeg; jj++){
	    Float_t eNegPt = tree->eNegPt[jj];
	    Float_t eNegPhi = tree->eNegPhi[jj];
	    Float_t eNegEta = tree->eNegEta[jj];
	    Int_t   eNegParentGID = tree->eNegParentGID[jj];
	    double BRNeg = GetBR(eNegParentGID);

	    TLorentzVector d01;
	    d01.SetPtEtaPhiM(ePosPt,ePosEta,ePosPhi,M_el);
	    TLorentzVector d02;
	    d02.SetPtEtaPhiM(eNegPt,eNegEta,eNegPhi,M_el);
	    TLorentzVector pair0;
	    pair0 = d01+d02;
	    bool acc0 = fabs(pair0.Rapidity())<1
	                 && d01.Pt()>0.2 && fabs(d01.Eta()<1)
	                 && d02.Pt()>0.2 && fabs(d02.Eta()<1);

	    h0pTInvMassccbar->Fill(pair0.M(),pair0.Pt(),BRNeg*BRPos);
	    h0pTInvMassccbarLE->Fill(pair0.M(),pair0.Pt(),BRNeg*BRPos);
	    if(acc0) {
	      h0pTInvMassAccccbar->Fill(pair0.M(),pair0.Pt(),BRNeg*BRPos);
	      h0pTInvMassAccccbarLE->Fill(pair0.M(),pair0.Pt(),BRNeg*BRPos);
	    }

	    Double_t ptRes1 = fpTSmear->Eval(d01.Pt())*hDoubleCrystalBall->GetRandom()/0.01;
	    Double_t ptRes2 = fpTSmear->Eval(d02.Pt())*hDoubleCrystalBall->GetRandom()/0.01;

	    Double_t newpt1 =  d01.Pt() + d01.Pt()*ptRes1;
	    Double_t newpt2 =  d02.Pt() + d02.Pt()*ptRes2;
	    Double_t neweta1=-999.,newphi1=-999.;
	    Double_t neweta2=-999.,newphi2=-999.;
	    neweta1 = d01.Eta();
	    neweta2 = d02.Eta();
	    newphi1 = d01.Phi();
	    newphi2 = d02.Phi();

	    TLorentzVector d1;
	    TLorentzVector d2;
	    TLorentzVector pair;

	    d1.SetPtEtaPhiM(newpt1,neweta1,newphi1,M_el);
	    d2.SetPtEtaPhiM(newpt2,neweta2,newphi2,M_el);
	    pair = d1+d2;
	    bool acc;

	    acc = fabs(pair.Rapidity())<1
	          && d1.Pt()>0.2 && fabs(d1.Eta())<1
	          && d2.Pt()>0.2 && fabs(d2.Eta())<1;

	    hpTInvMassccbar->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    hpTInvMassccbarLE->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    if(acc) {
	      hpTInvMassAccccbar->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	      hpTInvMassAccccbarLE->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    }

	    // Random deltPhi
	    double dphi = rnu->Rndm()*2.*PI;
	    newphi2 = newphi1 + dphi;
	    newphi2 = newphi2>2.*PI ? newphi2-2.*PI :newphi2;
	    d1.SetPtEtaPhiM(newpt1,neweta1,newphi1,M_el);
	    d2.SetPtEtaPhiM(newpt2,neweta2,newphi2,M_el);
	    pair = d1+d2;

	    acc = fabs(pair.Rapidity())<1
	          && d1.Pt()>0.2 && fabs(d1.Eta())<1
	          && d2.Pt()>0.2 && fabs(d2.Eta())<1;
 	    hpTInvMassRmDPhiccbar->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
 	    hpTInvMassRmDPhiccbarLE->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    if(acc) {
	      hpTInvMassRmDPhiAccccbar->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	      hpTInvMassRmDPhiAccccbarLE->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    }

	    // Random Phi
	    newphi1 = hPosPhi->GetRandom();
	    newphi2 = hNegPhi->GetRandom();
	    d1.SetPtEtaPhiM(newpt1,neweta1,newphi1,M_el);
	    d2.SetPtEtaPhiM(newpt2,neweta2,newphi2,M_el);
	    pair = d1+d2;

	    acc = fabs(pair.Rapidity())<1
	          && d1.Pt()>0.2 && fabs(d1.Eta())<1
	          && d2.Pt()>0.2 && fabs(d2.Eta())<1;
	    hpTInvMassRmPhiccbar->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    hpTInvMassRmPhiccbarLE->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    if(acc) {
	      hpTInvMassRmPhiAccccbar->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	      hpTInvMassRmPhiAccccbarLE->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    }

	    // Random Eta and Phi
	    newphi1 = hPosEta->GetRandom();
	    newphi2 = hNegEta->GetRandom();
	    d1.SetPtEtaPhiM(newpt1,neweta1,newphi1,M_el);
	    d2.SetPtEtaPhiM(newpt2,neweta2,newphi2,M_el);
	    pair = d1+d2;

	    acc = fabs(pair.Rapidity())<1
	          && d1.Pt()>0.2 && fabs(d1.Eta())<1
	          && d2.Pt()>0.2 && fabs(d2.Eta())<1;
	    hpTInvMassRmEtaccbar->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    hpTInvMassRmEtaccbarLE->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    if(acc) {
	      hpTInvMassRmEtaAccccbar->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	      hpTInvMassRmEtaAccccbarLE->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    }

	    // Random All
	    newpt1 = hPosPt->GetRandom()*(1+ptRes1);
	    newpt2 = hNegPt->GetRandom()*(1+ptRes2);
	    d1.SetPtEtaPhiM(newpt1,neweta1,newphi1,M_el);
	    d2.SetPtEtaPhiM(newpt2,neweta2,newphi2,M_el);
	    pair = d1+d2;

	    acc = fabs(pair.Rapidity())<1
                  && d1.Pt()>0.2 && fabs(d1.Eta())<1
	          && d2.Pt()>0.2 && fabs(d2.Eta())<1;
 	    hpTInvMassRmAllccbar->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
 	    hpTInvMassRmAllccbarLE->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    if(acc) {
	      hpTInvMassRmAllAccccbar->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	      hpTInvMassRmAllAccccbarLE->Fill(pair.M(),pair.Pt(),BRNeg*BRPos);
	    }
	  } //eNeg loop
	} // ePos loop
    }
    fout->Write();
    fout->Close();
}

Double_t GetBR(Int_t id){
    Double_t BR =0.;
    if(abs(id)==411) BR = 1.607E-1;//0.172;//1.607E-1;
    if(abs(id)==421) BR = 0.0649;//0.0671;////0.0649;
    if(abs(id)==431) BR = 0.065;//0.08;//0.065;
    if(abs(id)==4122) BR = 0.045;

    return BR;
}
Double_t pTSmear(double *x, double *p){
    Double_t pT = x[0];
    Double_t a = p[0];
    Double_t b = p[1];

    return TMath::Sqrt(a*a*pT*pT+b*b);
}//pTSmear
Double_t DoubleCrystalBall(double *x, double *p){//you can see twice as far into the future

    Double_t xx = x[0];
    Double_t alpha1 = p[0];
    Double_t alpha2 = p[1];
    Double_t n1 = p[2];
    Double_t n2 = p[3];
    Double_t xxmean = p[4];
    Double_t sigma = p[5];
    Double_t N = p[6];


    Double_t jmb = (xx - xxmean)/sigma;

    Double_t A = TMath::Power(n1/TMath::Abs(alpha1),n1)*TMath::Exp(-alpha1*alpha1/2.);
    Double_t B = n1/TMath::Abs(alpha1)-TMath::Abs(alpha1);
    Double_t C = TMath::Power(n2/TMath::Abs(alpha2),n2)*TMath::Exp(-alpha2*alpha2/2.);
    Double_t D = n2/TMath::Abs(alpha2)-TMath::Abs(alpha2);

    if(jmb < -alpha1){
	return N*A*TMath::Power(B-jmb,-n1);
    }//alpha1
    else if(jmb < alpha2){
	return N*TMath::Exp(-jmb*jmb/2.);
    }//alpha2
    else{
	return N*C*TMath::Power((D+jmb),-n2);
    }//else
}//doublecrystalball
