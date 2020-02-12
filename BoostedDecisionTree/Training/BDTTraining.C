// --------------------------------- Both settings
#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <stdio.h>
#include "TChain.h"
#include "TFile.h"
#include "TGraph.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TGClient.h"

// --------------------------------- ROOT6 settings
#if defined(__CLING__)

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"
#include "TMVA/MethodBase.h"
#include "TMVA/Types.h"
#include "TMVA/IMethod.h"
#include "TMVA/mvas.h"

// --------------------------------- ROOT5 settings
#elif defined(__CINT__)

#endif
// --------------------------------- XXXXXXXXXXXXXX

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */

enum treeType {kSgn=0,kBkg,kUndef};
enum variableCategories {kVariable=0,kSpectator,kNotUsed};

// PtBins - settings
const Int_t nPtBins = 6;
Float_t 	ptBins[nPtBins+1]     = {1.,    2.,     4.,     6.,     8.,    12.,           24.};
Int_t 		ptBinsInt[nPtBins]    =    {     1,      2,      3,      4,      5,             6};
Bool_t 		boolTrainBin[nPtBins] =    { kFALSE,  kFALSE,  kFALSE,  kTRUE,   kFALSE,   kFALSE};

// Tree settings
const Int_t  nTrees = 3;
Long64_t     nSgnEntries[nTrees][nPtBins];
TTree*       inputSgnTrees[nTrees][nPtBins];
Bool_t		 flagUseSgnTrees[nTrees][nPtBins] = {kFALSE};
Long64_t     nSgnEntriesTot[nPtBins] = {0};

Long64_t     nBkgEntries[nTrees][nPtBins];
TTree*       inputBkgTrees[nTrees][nPtBins];
Bool_t		 flagUseBkgTrees[nTrees][nPtBins] = {kFALSE};
Long64_t     nBkgEntriesTot[nPtBins] = {0};

// Initialize Variables settings
const Int_t nVar = 35;
TString  variableTreeNames[nVar];
TString  variableGraphNames[nVar];
TString  variableUnits[nVar];
Bool_t	 isTrainVariable[nVar];
Int_t	 variableCategories[nVar];

// TMVA Objects
TString 	dataSetNames[nPtBins];		// Name of Dataloader
TString     nVarTitle;

// Other global objects
Int_t 	multiplicityLow = 0;	    // Multiplicity/Centrality percentile up
Int_t 	multiplicityUp  = 0;	    // Multiplicity/Centrality percentile low
Double_t 	ptLow       = 0;		// pT lower range
Double_t 	ptUp        = 0;		// pT upper range

// BDT Parameters
int NTrees = 850; 				    // default: 800
int MaxDepth = 3; 					// default: 3
const char* MinNodeSize = "5.0%"; 	// default: '5.0%'
int nCuts = 20; 					// default: 20 
const char* BoostType = "AdaBoost";
float AdaBoostBeta = 0.5; 			// default: 0.5
float BaggedSampleFraction = 0.6;
Bool_t UseBaggedBoost = kFALSE;
const char* SeparationType = "GiniIndex";

// Training and Test Parameters
Int_t maxTrainEntries = 500000;
Double_t ratioTrainTest = 0.7; 				// (0.6 equals 60% training, 40% testing)
const char* NormMode = "EqualNumEvents";
const char* SplitMode = "Random";			// default: "Random"
const char* Mixmode  = "SameAsSplitMode";	// default: "SameAsSplitMode"
const int SplitSeed = 100;					// default: 100

void SetVariables(){

	// Variable name in the tree
    variableTreeNames[0] = "massLc2K0Sp";
    variableTreeNames[1] = "alphaArm";
    variableTreeNames[2] = "massK0S";
    variableTreeNames[3] = "massLambda";
    variableTreeNames[4] = "massLambdaBar";
    variableTreeNames[5] = "cosPAK0S";
    variableTreeNames[6] = "dcaV0";
    variableTreeNames[7] = "tImpParBach";
    variableTreeNames[8] = "tImpParV0";
    variableTreeNames[9] = "nSigmaTPCpr";
    variableTreeNames[10] = "nSigmaTOFpr";
    variableTreeNames[11] = "bachelorPt";
    variableTreeNames[12] = "V0positivePt";
    variableTreeNames[13] = "V0negativePt";
    variableTreeNames[14] = "dcaV0pos";
    variableTreeNames[15] = "dcaV0neg";
    variableTreeNames[16] = "v0Pt";
    variableTreeNames[17] = "bachTPCmom";
    variableTreeNames[18] = "LcPt";
    variableTreeNames[19] = "combinedProtonProb";
    variableTreeNames[20] = "V0positiveEta";
    variableTreeNames[21] = "bachelorP";
    variableTreeNames[22] = "bachelorEta";
    variableTreeNames[23] = "v0P";
    variableTreeNames[24] = "DecayLengthK0S";
    variableTreeNames[25] = "nSigmaTPCpi";
    variableTreeNames[26] = "nSigmaTPCka";
    variableTreeNames[27] = "NtrkRaw";
    variableTreeNames[28] = "NtrkCorr";
    variableTreeNames[29] = "CosThetaStar";
    variableTreeNames[30] = "signd0";        
    variableTreeNames[31] = "centrality"; 
    variableTreeNames[32] = "NtrkAll";
    variableTreeNames[33] = "origin";
    variableTreeNames[34] = "ptArm";

    // // Variable train of spectector, (only if kIsUsed)
    // isTrainVariable[0]=kFALSE;
    // isTrainVariable[1]=kTRUE;
    // isTrainVariable[2]=kTRUE;
    // isTrainVariable[3]=kFALSE;
    // isTrainVariable[4]=kFALSE;
    // isTrainVariable[5]=kTRUE;
    // isTrainVariable[6]=kTRUE;
    // isTrainVariable[7]=kTRUE;
    // isTrainVariable[8]=kTRUE;
    // isTrainVariable[9]=kTRUE;
    // isTrainVariable[10]=kTRUE;
    // isTrainVariable[11]=kTRUE;
    // isTrainVariable[12]=kFALSE;
    // isTrainVariable[13]=kFALSE;
    // isTrainVariable[14]=kFALSE;
    // isTrainVariable[15]=kFALSE;
    // isTrainVariable[16]=kFALSE;
    // isTrainVariable[17]=kFALSE;
    // isTrainVariable[18]=kFALSE;
    // isTrainVariable[19]=kTRUE;
    // isTrainVariable[20]=kFALSE;
    // isTrainVariable[21]=kFALSE;
    // isTrainVariable[22]=kFALSE;
    // isTrainVariable[23]=kFALSE;
    // isTrainVariable[24]=kFALSE;
    // isTrainVariable[25]=kFALSE;
    // isTrainVariable[26]=kFALSE;
    // isTrainVariable[27]=kFALSE;
    // isTrainVariable[28]=kFALSE;
    // isTrainVariable[29]=kTRUE;
    // isTrainVariable[30]=kTRUE;
    // isTrainVariable[31]=kFALSE;
    // isTrainVariable[32]=kFALSE;
    // isTrainVariable[33]=kFALSE;
    // isTrainVariable[34]=kTRUE;

    // // Is the variable used (reduces memory)
    // variableCategories[0]  = kNotUsed;
    // variableCategories[1]  = kNotUsed;
    // variableCategories[2]  = kVariable;
    // variableCategories[3]  = kNotUsed;
    // variableCategories[4]  = kNotUsed;
    // variableCategories[5]  = kVariable;
    // variableCategories[6]  = kNotUsed;
    // variableCategories[7]  = kVariable;
    // variableCategories[8]  = kVariable;
    // variableCategories[9]  = kVariable;
    // variableCategories[10] = kVariable;
    // variableCategories[11] = kVariable;
    // variableCategories[12] = kNotUsed;
    // variableCategories[13] = kNotUsed;
    // variableCategories[14] = kNotUsed;
    // variableCategories[15] = kNotUsed;
    // variableCategories[16] = kNotUsed;
    // variableCategories[17] = kVariable;
    // variableCategories[18] = kNotUsed;
    // variableCategories[19] = kNotUsed;
    // variableCategories[20] = kNotUsed;
    // variableCategories[21] = kVariable;
    // variableCategories[22] = kNotUsed;
    // variableCategories[23] = kNotUsed;
    // variableCategories[24] = kNotUsed;
    // variableCategories[25] = kVariable;
    // variableCategories[26] = kVariable;
    // variableCategories[27] = kNotUsed;
    // variableCategories[28] = kNotUsed;
    // variableCategories[29] = kVariable;
    // variableCategories[30] = kVariable;
    // variableCategories[31] = kNotUsed;
    // variableCategories[32] = kNotUsed;
    // variableCategories[33] = kNotUsed;
    // variableCategories[34] = kNotUsed;

    // Variable names for graph title and x-axjs
    variableGraphNames[0] = "M(K^{0}_{s}p)";
    variableGraphNames[1] = "#alpha_{Arm}";
    variableGraphNames[2] = "M(K^{0}_{s})";
    variableGraphNames[3] = "M(#bar{#Lambda})";
    variableGraphNames[4] = "M(#Lambda)";
    variableGraphNames[5] = "cos(PA_{V0})";
    variableGraphNames[6] = "dca_{V0}";
    variableGraphNames[7] = "d_{0,bachelor}";
    variableGraphNames[8] = "d_{0,V0}";
    variableGraphNames[9] = "n#sigma_{TPC,proton}";
    variableGraphNames[10] = "n#sigma_{TOF,proton}";
    variableGraphNames[11] = "p_{T,bachelor}";
    variableGraphNames[12] = "p_{T,daught+}";
    variableGraphNames[13] = "p_{T,daught-}";
    variableGraphNames[14] = "d_{0,daught+}";
    variableGraphNames[15] = "d_{0,daught-}";
    variableGraphNames[16] = "p_{T,V0}";
    variableGraphNames[17] = "p_{TPC,bachelor}";
    variableGraphNames[18] = "p_{T,#Lambda_{c}}";
    variableGraphNames[19] = "combProtonProb";
    variableGraphNames[20] = "#eta_{daught+}";
    variableGraphNames[21] = "p_{bachelor}";
    variableGraphNames[22] = "#eta_{bachelor}";
    variableGraphNames[23] = "p_{V0}";
    variableGraphNames[24] = "L_{V0}";
    variableGraphNames[25] = "n#sigma_{TPC,#pi}";
    variableGraphNames[26] = "n#sigma_{TPC,K}";
    variableGraphNames[27] = "Ntrk_{Raw}";
    variableGraphNames[28] = "Ntrk_{Corr}";
    variableGraphNames[29] = "Cos(#theta*)";
    variableGraphNames[30] = "signed d_{0}";        
    variableGraphNames[31] = "Centrality"; 
    variableGraphNames[32] = "Ntrk_{All}";
    variableGraphNames[33] = "Origin";
    variableGraphNames[34] = "p_{T}^{Arm}";

    // Variable unit
    variableUnits[0] = "GeV/c^{2}";
    variableUnits[1] = " ";
    variableUnits[2] = "GeV/c^{2}";
    variableUnits[3] = "GeV/c^{2}";
    variableUnits[4] = "GeV/c^{2}";
    variableUnits[5] = " ";
    variableUnits[6] = "cm";
    variableUnits[7] = "cm";
    variableUnits[8] = "cm";
    variableUnits[9] = " ";
    variableUnits[10] = " ";
    variableUnits[11] = "GeV/c";
    variableUnits[12] = "GeV/c";
    variableUnits[13] = "GeV/c";
    variableUnits[14] = "cm";
    variableUnits[15] = "cm";
    variableUnits[16] = "GeV/c";
    variableUnits[17] = " ";
    variableUnits[18] = "GeV/c";
    variableUnits[19] = " ";
    variableUnits[20] = " ";
    variableUnits[21] = "GeV/c";
    variableUnits[22] = " ";
    variableUnits[23] = "GeV/c";
    variableUnits[24] = "cm";
    variableUnits[25] = " ";
    variableUnits[26] = " ";
    variableUnits[27] = "#";
    variableUnits[28] = "#";
    variableUnits[29] = " ";
    variableUnits[30] = "cm";        
    variableUnits[31] = "%"; 
    variableUnits[32] = "#";
    variableUnits[33] = " ";
    variableUnits[34] = " ";

    Printf("-> Variables setted");

}

void BDTTrainingFactory(TFile* outputFile,Int_t nvars){

	Printf("\n========= BDTTrainingFactory started =========");

	// Timer
	TStopwatch sw;
	sw.Start();

	// // Training and Test parameters
	// char ttParameters[300];
	// sprintf(ttParameters,"!V:nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d:SplitMode=%s:"\
	// 		"Mixmode=%s:SplitSeed=%d:NormMode=%s",
	// 		nTrain_Signal,nTrain_Background,nTest_Signal,nTest_Background,SplitMode,Mixmode,SplitSeed,NormMode);
	// Printf("* Training and test parameters");
	// Printf("%s",ttParameters);

	// BDT Parameters
	char bdtParameters [300];
	sprintf(bdtParameters,"!H:!V:NTrees=%d:MinNodeSize=%s:MaxDepth=%d:BoostType=%s:AdaBoostBeta=%.2f:SeparationType=%s:nCuts=%d%s",
			NTrees,MinNodeSize,MaxDepth,BoostType,AdaBoostBeta,SeparationType,nCuts,UseBaggedBoost ? Form(":UseBaggedBoost:BaggedSampleFraction=%0.2f",BaggedSampleFraction) : "");
	Printf("* BDT parameters");
	Printf("%s\n",bdtParameters);
	// Create the factory object. Later you can choose the methods
	// whose performance you'd like to investigate. The factory is
	// the only TMVA object you have to interact with
	// - https://root.cern.ch/download/doc/tmva/TMVAUsersGuide.pdf#page=20
	TMVA::Factory *factory = new TMVA::Factory( "LHC17j4b_TMVAClassification_BDT", outputFile,
		"!V:!Silent:Color:DrawProgressBar:AnalysisType=Classification" ); //Transformations=I;D;P;U;G,D

	// Create dataloaders
	TMVA::DataLoader *dataloader;

	// Later used objects
	Int_t nTrain_Signal[nPtBins];
	Int_t nTest_Signal[nPtBins];
	Int_t nTrain_Background[nPtBins];
	Int_t nTest_Background[nPtBins];

	for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){

		// pT Ranges
		ptLow = ptBins[iPtBin];
		ptUp = ptBins[iPtBin+1];

		Printf("\n%s---------------------------------\n",boolTrainBin[iPtBin] ? GREEN : RED);
		Printf("Multiplicity = %d - %d, pT = %.2f-%.2f",multiplicityLow,multiplicityUp,ptLow,ptUp);
		Printf("     %s",boolTrainBin[iPtBin] ? "> Trained" : "> NOT Trained");
		Printf("---------------------------------%s\n",RESET);

		if(!boolTrainBin[iPtBin]) continue;

		char* dataSetName = Form("LcTraining_%s_Mult%d_%d_pT%.0f_%.0f",nVarTitle.Data(),multiplicityLow,multiplicityUp,ptLow,ptUp);
		dataSetNames[iPtBin]=dataSetName;

		TMVA::DataLoader *dataloader=new TMVA::DataLoader(dataSetName); // TMVA::DataLoader("<name>");

		char ttParameters[300];
		nTrain_Signal[iPtBin]     = ratioTrainTest*Int_t(nSgnEntriesTot[iPtBin]);
		if(nTrain_Signal[iPtBin] > maxTrainEntries) 
			  nTrain_Signal[iPtBin] 	= maxTrainEntries;
		nTest_Signal[iPtBin] 	    = (nTrain_Signal[iPtBin]/ratioTrainTest)*(1. - ratioTrainTest);

		nTrain_Background[iPtBin] = ratioTrainTest*Int_t(nBkgEntriesTot[iPtBin]);
		if(nTrain_Background[iPtBin] > maxTrainEntries) 
			  nTrain_Background[iPtBin] = maxTrainEntries;
		nTest_Background[iPtBin] 	= (nTrain_Background[iPtBin]/ratioTrainTest)*(1. - ratioTrainTest);

		sprintf(ttParameters,"!V:nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d:nTest_Background=%d:"\
				"SplitMode=%s:Mixmode=%s:SplitSeed=%d:NormMode=%s",
				nTrain_Signal[iPtBin],nTrain_Background[iPtBin],nTest_Signal[iPtBin],nTest_Background[iPtBin],SplitMode,Mixmode,SplitSeed,NormMode);
		// Printf("* TT parameters");
		// Printf("%s\n",ttParameters);

          if (nvars == 14){
            dataloader->AddVariable("massK0S","massK0S","GeV/c^{2}",'F');
            dataloader->AddVariable("tImpParBach","tImpParBach","cm",'F');
            dataloader->AddVariable("tImpParV0","tImpParV0","cm",'F');
            dataloader->AddVariable("bachelorPt","bachelorPt","GeV/c",'F');
            dataloader->AddVariable("CtK0S := DecayLengthK0S*0.497/v0P","CtK0S","cm",'F');
            dataloader->AddVariable("cosPAK0S","cosPAK0S","",'F');
            dataloader->AddVariable("CosThetaStar","CosThetaStar","",'F');
            dataloader->AddVariable("signd0","signd0","cm",'F');
            dataloader->AddVariable("bachelorP","bachelorP","GeV/c",'F');
            dataloader->AddVariable("nSigmaTOFpr","nSigmaTOFpr","",'F');
            dataloader->AddVariable("nSigmaTPCpr","nSigmaTPCpr","",'F');
            dataloader->AddVariable("nSigmaTPCpi","nSigmaTPCpi","",'F');
            dataloader->AddVariable("nSigmaTPCka","nSigmaTPCka","",'F');
            dataloader->AddVariable("bachTPCmom","bachTPCmom","",'F');
          }
          else if (nvars == 11) {
            dataloader->AddVariable("massK0S","massK0S","GeV/c^{2}",'F');
            dataloader->AddVariable("tImpParBach","tImpParBach","cm",'F');
            dataloader->AddVariable("tImpParV0","tImpParV0","cm",'F');
            dataloader->AddVariable("CtK0S := DecayLengthK0S*0.497/v0P","CtK0S","cm",'F');
            dataloader->AddVariable("cosPAK0S","cosPAK0S","",'F');
            dataloader->AddVariable("CosThetaStar","CosThetaStar","",'F');
            dataloader->AddVariable("signd0","signd0","cm",'F');
            dataloader->AddVariable("nSigmaTOFpr","nSigmaTOFpr","",'F');
            dataloader->AddVariable("nSigmaTPCpr","nSigmaTPCpr","",'F');
            dataloader->AddVariable("nSigmaTPCpi","nSigmaTPCpi","",'F');
            dataloader->AddVariable("nSigmaTPCka","nSigmaTPCka","",'F');

            dataloader->AddSpectator("massLc2K0Sp","massLc2K0Sp","GeV/c^{2}",'F');
            dataloader->AddSpectator("LcPt","LcPt","GeV/c",'F');
            // dataloader->AddSpectator("massLc2Lambdapi","massLc2Lambdapi","GeV/c^{2}",'F');
            dataloader->AddSpectator("centrality","massLc2Lambdapi","units",'F');
            dataloader->AddSpectator("massLambda","massLambda","GeV/c^{2}",'F');
            dataloader->AddSpectator("massLambdaBar","massLambdaBar","GeV/c^{2}",'F');
            dataloader->AddSpectator("cosPAK0S","cosPAK0S","",'F');
            dataloader->AddSpectator("V0positivePt","V0positivePt","GeV/c",'F');
            dataloader->AddSpectator("V0negativePt","V0negativePt","GeV/c",'F');
            dataloader->AddSpectator("dcaV0pos","dcaV0pos","cm",'F');
            dataloader->AddSpectator("dcaV0neg","dcaV0neg","cm",'F');
            dataloader->AddSpectator("v0Pt","v0Pt","GeV/c",'F');
            dataloader->AddSpectator("dcaV0","dcaV0","cm",'F');
            dataloader->AddSpectator("V0positiveEta","V0positiveEta","",'F');
            dataloader->AddSpectator("bachelorEta","bachelorEta","",'F');
            dataloader->AddSpectator("centrality","centrality","",'F');
          }
          else if (nvars == 10){
            dataloader->AddVariable("massK0S","massK0S","GeV/c^{2}",'F');
            dataloader->AddVariable("tImpParBach","tImpParBach","cm",'F');
            dataloader->AddVariable("tImpParV0","tImpParV0","cm",'F');
            dataloader->AddVariable("CtK0S := DecayLengthK0S*0.497/v0P","CtK0S","cm",'F');
            dataloader->AddVariable("cosPAK0S","cosPAK0S","",'F');
            dataloader->AddVariable("signd0","signd0","cm",'F');
            dataloader->AddVariable("nSigmaTOFpr","nSigmaTOFpr","",'F');
            dataloader->AddVariable("nSigmaTPCpr","nSigmaTPCpr","",'F');
            dataloader->AddVariable("nSigmaTPCpi","nSigmaTPCpi","",'F');
            dataloader->AddVariable("nSigmaTPCka","nSigmaTPCka","",'F');
          }
          else if (nvars == 7) {
            dataloader->AddVariable("massK0S","massK0S","GeV/c^{2}",'F');
            dataloader->AddVariable("tImpParBach","tImpParBach","cm",'F');
            dataloader->AddVariable("tImpParV0","tImpParV0","cm",'F');
            dataloader->AddVariable("CtK0S := DecayLengthK0S*0.497/v0P","CtK0S","cm",'F');
            dataloader->AddVariable("cosPAK0S","cosPAK0S","",'F');
            dataloader->AddVariable("CosThetaStar","CosThetaStar","",'F');
            dataloader->AddVariable("signd0","signd0","cm",'F');

            dataloader->AddSpectator("massLc2K0Sp","massLc2K0Sp","GeV/c^{2}",'F');
            dataloader->AddSpectator("LcPt","LcPt","GeV/c",'F');
            // dataloader->AddSpectator("massLc2Lambdapi","massLc2Lambdapi","GeV/c^{2}",'F');
            dataloader->AddSpectator("centrality","massLc2Lambdapi","units",'F');
            dataloader->AddSpectator("massLambda","massLambda","GeV/c^{2}",'F');
            dataloader->AddSpectator("massLambdaBar","massLambdaBar","GeV/c^{2}",'F');
            dataloader->AddSpectator("cosPAK0S","cosPAK0S","",'F');
            dataloader->AddSpectator("V0positivePt","V0positivePt","GeV/c",'F');
            dataloader->AddSpectator("V0negativePt","V0negativePt","GeV/c",'F');
            dataloader->AddSpectator("dcaV0pos","dcaV0pos","cm",'F');
            dataloader->AddSpectator("dcaV0neg","dcaV0neg","cm",'F');
            dataloader->AddSpectator("v0Pt","v0Pt","GeV/c",'F');
            dataloader->AddSpectator("dcaV0","dcaV0","cm",'F');
            dataloader->AddSpectator("V0positiveEta","V0positiveEta","",'F');
            dataloader->AddSpectator("bachelorEta","bachelorEta","",'F');
            dataloader->AddSpectator("centrality","centrality","",'F');
          }


  //       //_________________________________________________________________________________

  //       TString fullVariableString = "";
  //       TString fullSpectatorString = "";

		// Printf("\nType     :%20.20s \t %20.20s \t %10.10s","TreeName","Name","Unit");
		// Printf("---------------------------------------------------------------------------");
		// for(Int_t iVar=0;iVar<nVar;iVar++){

		// 	if     (variableCategories[iVar]==kNotUsed){ continue;}
  //           else if(variableCategories[iVar]==kVariable){
  //               fullVariableString.Append(Form("%s,",variableTreeNames[iVar].Data()));
  //               Printf("Variable : %20.20s \t %20.20s \t %10.10s",variableTreeNames[iVar].Data(),variableGraphNames[iVar].Data(),variableUnits[iVar].Data());
  //               dataloader->AddVariable(  variableTreeNames[iVar].Data(),variableGraphNames[iVar].Data(),variableUnits[iVar].Data(),'F');
  //           }
  //           else if(variableCategories[iVar]==kSpectator){
  //               fullSpectatorString.Append(Form("%s,",variableTreeNames[iVar].Data()));
  //               Printf("Spectator: %20.20s \t %20.20s \t %10.10s",variableTreeNames[iVar].Data(),variableGraphNames[iVar].Data(),variableUnits[iVar].Data());
  //               dataloader->AddSpectator( variableTreeNames[iVar].Data(),variableGraphNames[iVar].Data(),variableUnits[iVar].Data(),'F');
  //           }
		// }

  //       fullVariableString.Append("ctau");
  //       dataloader->AddVariable(  "ctau := DecayLengthK0S*0.497/v0P","ctau","cm",'F');
  //       Printf("Variable : %20.20s \t %20.20s \t %10.10s","ctau","ctau","cm");
		// Printf("\n");
  //       Printf("Variables : %s",fullVariableString.Data());
  //       Printf("Spectators: %s",fullSpectatorString.Data());
  //       Printf("\n");

  //       //_________________________________________________________________________________

		// You can add an arbitrary number of signal or background trees
		for(Int_t iT = 0; iT < nTrees;iT++){
			if(flagUseSgnTrees[iT][iPtBin]) dataloader->AddSignalTree ( inputSgnTrees[iT][iPtBin], Double_t(nSgnEntries[iT][iPtBin])/Double_t(nSgnEntriesTot[iPtBin]) );
		}
		for(Int_t iT = 0; iT < nTrees;iT++){
			if(flagUseBkgTrees[iT][iPtBin]) dataloader->AddBackgroundTree ( inputBkgTrees[iT][iPtBin], Double_t(nBkgEntries[iT][iPtBin])/Double_t(nBkgEntriesTot[iPtBin]) );
		}

		// Apply additional cuts on the signal and background samples (can be different)
		TCut mycuts = "";

		// Training and Test parameters in Dataloader
		dataloader->PrepareTrainingAndTestTree( mycuts, ttParameters);

		// // BDT Parameters in Bookmethod
		TMVA::MethodBase* method = factory->BookMethod(dataloader, TMVA::Types::kBDT, Form("%s_Mult%d_%d_pT%.0f_%.0f",nVarTitle.Data(),multiplicityLow,multiplicityUp,ptLow,ptUp),bdtParameters);

	}

	// Loop over pT bins
	Printf("\n\n\nptBin \t ptLow \t  ptUp \t nTrain_Signal \t nTrain_Backgr \t  nTest_Signal \t  nTest_Backgr \t NormMode");
	for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
		if(!boolTrainBin[iPtBin]) continue;
		Printf("%5.2d \t %5.2f \t %5.2f \t %13d \t %13d \t %13d \t %13d \t %8.8s",ptBinsInt[iPtBin],ptBins[iPtBin],ptBins[iPtBin+1],nTrain_Signal[iPtBin],nTrain_Background[iPtBin],nTest_Signal[iPtBin],nTest_Background[iPtBin],NormMode);
	}
	Printf("\n\n");

	Printf("\n%s---------------------------------\n",MAGENTA);
	Printf("  === Prepration done -> Start Training");
	Printf("\n---------------------------------%s\n",RESET);

    // Train MVAs using the set of training events
	factory->TrainAllMethods();

	Printf("\n%s---------------------------------\n",MAGENTA);
	Printf("  === Training done -> Start Testing");
	Printf("\n---------------------------------%s\n",RESET);

	// Evaluate all MVAs using the set of test events
	factory->TestAllMethods();

	Printf("\n%s---------------------------------\n",MAGENTA);
	Printf("  === Testing done -> Start Evaluating");
	Printf("\n---------------------------------%s\n",RESET);

	// Evaluate and compare performance of all configured MVAs
	factory->EvaluateAllMethods();

	Printf("\n%s---------------------------------\n",MAGENTA);
	Printf("  === Evaluating done...");
	Printf("\n---------------------------------%s\n",RESET);

	// Cd to main directory of outputfile
	outputFile->cd();

	// Get elapsed time
	sw.Stop();
	std::cout<< "BDTFactory: "; sw.Print();
	std::cout<<std::endl;

	return;
}

void BDTTraining(Int_t multLow,
		         Int_t multUp,
				 TString  inputSgnBkgTreesFile,
                 Int_t nvars = 14,
				 TString  fileNameAddition = ""
				 ){


	std::cout << std::endl << std::flush;
	Printf("\n========= BDTTraining.C Started =========");

	// This loads the library
	TMVA::Tools::Instance();

	multiplicityLow = multLow;
	multiplicityUp  = multUp;

	//================================================
  	//         Input
  	//================================================

	TFile* inputFile(0);
	if (!gSystem->AccessPathName( inputSgnBkgTreesFile )) {
	  inputFile = TFile::Open( inputSgnBkgTreesFile ); // check if file in local directory exists
	}
	else {
	  Printf("ERROR: could not open: %s",inputSgnBkgTreesFile.Data());
	  exit(1);
	}

	for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){

		for(Int_t iT = 0; iT < nTrees; iT++){
			
			TTree* tTree = 0x0;
			// pT Bin Ranges
			ptLow = ptBins[iPtBin];
			ptUp  = ptBins[iPtBin+1];

			tTree = (TTree*)inputFile->Get(Form("Sgn[%d]Pt[%.0f-%.0f]Mult[%d-%d]",iT,ptLow,ptUp,multiplicityLow,multiplicityUp));
			if(!tTree){
				flagUseSgnTrees[iT][iPtBin]=kFALSE;
			}
			else{
				inputSgnTrees[iT][iPtBin]   = tTree;
				flagUseSgnTrees[iT][iPtBin] = kTRUE;
				nSgnEntries[iT][iPtBin]     = tTree->GetEntries();
				nSgnEntriesTot[iPtBin]	   += nSgnEntries[iT][iPtBin];
			}

			tTree = (TTree*)inputFile->Get(Form("Bkg[%d]Pt[%.0f-%.0f]Mult[%d-%d]",iT,ptLow,ptUp,multiplicityLow,multiplicityUp));
			if(!tTree){
				flagUseBkgTrees[iT][iPtBin]=kFALSE;
			}
			else{
				inputBkgTrees[iT][iPtBin]   = tTree;
				flagUseBkgTrees[iT][iPtBin] = kTRUE;
				nBkgEntries[iT][iPtBin]     = tTree->GetEntries();
				nBkgEntriesTot[iPtBin]	   += nBkgEntries[iT][iPtBin];
			}
		}
	}

    if(nvars==14)       nVarTitle = "all";
    else if(nvars==11)  nVarTitle = "noP";
    else if(nvars==10)  nVarTitle = "noPCts";
    else if(nvars==7)   nVarTitle = "noNsigma";

	for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
		// pT Bin Ranges
		ptLow = ptBins[iPtBin];
		ptUp  = ptBins[iPtBin+1];

		Printf("\nptBin = %2d \t ptLow = %4.2f \t ptUp = % 4.2f \t Train = %s",ptBinsInt[iPtBin],ptLow,ptUp,boolTrainBin[iPtBin] ? Form("%sYes%s",GREEN,RESET) : Form("%sNo%s",RED,RESET));
		Printf("  iT \t Type \t flagUse \t %10.10s \t Weight","Entries");
		for(Int_t iT = 0; iT < nTrees; iT++){
			Printf("  %2d \t Sign \t %7d \t %10lld \t %6.2f", iT, Int_t(flagUseSgnTrees[iT][iPtBin]),nSgnEntries[iT][iPtBin],Double_t(nSgnEntries[iT][iPtBin])/Double_t(nSgnEntriesTot[iPtBin]));
		}
		for(Int_t iT = 0; iT < nTrees; iT++){
			Printf("  %2d \t Bkg \t %7d \t %10lld \t %6.2f", iT, Int_t(flagUseBkgTrees[iT][iPtBin]),nBkgEntries[iT][iPtBin],Double_t(nBkgEntries[iT][iPtBin])/Double_t(nBkgEntriesTot[iPtBin]));
		}

	}

	// ================================================
 	//  	        Output
 	// ================================================

  	// Save output
  	TString outputFileName;
	TFile *outputFile(0);
	if(fileNameAddition==""){
		outputFileName = Form("TrainingResults_%s_Mult%d_%d.root",nVarTitle.Data(),multiplicityLow,multiplicityUp);
		outputFile = new TFile(outputFileName.Data(),"RECREATE");
	}
	else{
		outputFileName = Form("TrainingResults_%s_Mult%d_%d_%s.root",nVarTitle.Data(),multiplicityLow,multiplicityUp,fileNameAddition.Data());
		outputFile = new TFile(outputFileName.Data(),"RECREATE");
	}
  	if(!outputFile){
		Printf("ERROR: outputFile not created");
		exit(1);
	}
	outputFile->cd();

	//================================================
  	//         Start training
  	//================================================

	SetVariables();

	BDTTrainingFactory(outputFile, nvars);

	//================================================
  	//         Evaluating output
  	//================================================

  	// Save output
	outputFile->Close();

	for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
		if(boolTrainBin[iPtBin]){

			TMVA::variables(dataSetNames[iPtBin],outputFileName,"InputVariables_Id","Input Variables (training sample)");
			TMVA::correlations(dataSetNames[iPtBin],outputFileName);
			TMVA::mvas(dataSetNames[iPtBin],outputFileName,TMVA::kMVAType);
			TMVA::mvas(dataSetNames[iPtBin],outputFileName,TMVA::kCompareType);
			TMVA::efficiencies(dataSetNames[iPtBin],outputFileName);
			TMVA::efficiencies(dataSetNames[iPtBin],outputFileName,3);
			TMVA::mvaeffs(dataSetNames[iPtBin],outputFileName); // Off on submit mode
			TMVA::BDTControlPlots(dataSetNames[iPtBin],outputFileName);
			// TMVA::TMVAGui( outputFileName ); // Off on submit mode
			// if(!gROOT->IsBatch()) gPad->Close();
		}

	}

	return;

}




