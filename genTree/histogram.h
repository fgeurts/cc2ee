const Int_t mMax = 100;

struct trkTree
{
	Int_t    ncString;
	Int_t    ncbarString;

	Int_t    nc;
	Bool_t   isStringC[mMax];
	Double_t cPt[mMax];
	Double_t cEta[mMax];
	Double_t cPhi[mMax];
	Double_t cY[mMax];

	Int_t    ncbar;
	Bool_t   isStringCbar[mMax];
	Double_t cbarPt[mMax];
	Double_t cbarEta[mMax];
	Double_t cbarPhi[mMax];
	Double_t cbarY[mMax];

	Int_t    nePos;
	Int_t    ePosParentGID[mMax];
	Double_t ePosPt[mMax];
	Double_t ePosEta[mMax];
	Double_t ePosPhi[mMax];

	Int_t    neNeg;
	Int_t    eNegParentGID[mMax];
	Double_t eNegPt[mMax];
	Double_t eNegEta[mMax];
	Double_t eNegPhi[mMax];

}; trkTree meTree;
