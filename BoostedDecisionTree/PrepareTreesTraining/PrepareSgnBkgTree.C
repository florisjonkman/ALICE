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

// --------------------------------- ROOT5 settings
#elif defined(__CINT__)

#endif
// --------------------------------- XXXXXXXXXXXXXX

enum treeType {kSgn=0,kBkg,kUndef};
enum variableCategories {kIsUsed=0,kNotUsed};

// PtBins - settings
const Int_t nPtBins = 6;
const Double_t ptBins[nPtBins+1] = {1.,2.,4.,6.,8.,12.,24.}; // pT bin ranges

// MultiplicityBins - settings
const Int_t nMultBins = 1;
Int_t multBinsLow[nMultBins] = {1   }; // Low included
Int_t multBinsUp[nMultBins] =  {9999}; // Up excluded

// Initialize Variables settings
const Int_t nVar = 35;
TString  variableTreeNames[nVar];
TString  variableGraphNames[nVar];
Int_t	 variableCategories[nVar];
TString	 variableUnits[nVar];

// Tree settings
Int_t nTrees = 0;
Long64_t nTreeEntries[6];

// Other global objects
TTree* outputTrees[6][nPtBins][nMultBins]; // iTree = 0 -> Signal, iTree = 1 ->Background
TH2F* hSgnEntriesPtVsMult;
TH2F* hBkgEntriesPtVsMult;

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
    // isTrainVariable[11]=kFALSE;
    // isTrainVariable[12]=kFALSE;
    // isTrainVariable[13]=kFALSE;
    // isTrainVariable[14]=kTRUE;
    // isTrainVariable[15]=kTRUE;
    // isTrainVariable[16]=kFALSE;
    // isTrainVariable[17]=kFALSE;
    // isTrainVariable[18]=kFALSE;
    // isTrainVariable[19]=kTRUE;
    // isTrainVariable[20]=kFALSE;
    // isTrainVariable[21]=kFALSE;
    // isTrainVariable[22]=kFALSE;
    // isTrainVariable[23]=kFALSE;
    // isTrainVariable[24]=kTRUE;
    // isTrainVariable[25]=kFALSE;
    // isTrainVariable[26]=kFALSE;
    // isTrainVariable[27]=kFALSE;
    // isTrainVariable[28]=kFALSE;
    // isTrainVariable[29]=kTRUE;
    // isTrainVariable[30]=kTRUE;
    // isTrainVariable[31]=kFALSE;
    // isTrainVariable[32]=kFALSE;
    // isTrainVariable[33]=kFALSE;
    // isTrainVariable[34]=kFALSE;

    // Is the variable used (reduces memory)
    variableCategories[0]  = kIsUsed;
    variableCategories[1]  = kIsUsed;
    variableCategories[2]  = kIsUsed;
    variableCategories[3]  = kIsUsed;
    variableCategories[4]  = kIsUsed;
    variableCategories[5]  = kIsUsed;
    variableCategories[6]  = kIsUsed;
    variableCategories[7]  = kIsUsed;
    variableCategories[8]  = kIsUsed;
    variableCategories[9]  = kIsUsed;
    variableCategories[10] = kIsUsed;
    variableCategories[11] = kIsUsed;
    variableCategories[12] = kIsUsed;
    variableCategories[13] = kIsUsed;
    variableCategories[14] = kIsUsed;
    variableCategories[15] = kIsUsed;
    variableCategories[16] = kIsUsed;
    variableCategories[17] = kIsUsed;
    variableCategories[18] = kIsUsed;
    variableCategories[19] = kIsUsed;
    variableCategories[20] = kIsUsed;
    variableCategories[21] = kIsUsed;
    variableCategories[22] = kIsUsed;
    variableCategories[23] = kIsUsed;
    variableCategories[24] = kIsUsed;
    variableCategories[25] = kIsUsed;
    variableCategories[26] = kIsUsed;
    variableCategories[27] = kIsUsed;
    variableCategories[28] = kIsUsed;
    variableCategories[29] = kIsUsed;
    variableCategories[30] = kIsUsed;
    variableCategories[31] = kIsUsed;
    variableCategories[32] = kIsUsed;
    variableCategories[33] = kIsUsed;
    variableCategories[34] = kIsUsed;

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

void GenerateSignalBackgroundTrees(TTree **inputTrees, TString* treeTitles, Int_t* treeTypes){

	// Variables needed
	Float_t pTBinLow;
    Float_t pTBinUp;
    Int_t multLow;
    Int_t multUp;

    Float_t variables[nVar];

    // Input Tree variables initializing
    for (Int_t iT = 0; iT < nTrees; iT++)
	{
		for(Int_t iVar=0;iVar<nVar;iVar++){

			if(variableCategories[iVar]==kNotUsed) continue;
		
			inputTrees[iT]->SetBranchAddress(variableTreeNames[iVar].Data(),&variables[iVar]);

		}	// End iVar - loop

	}	// End iT-loop

	// Output Tree variables initializing
	for (Int_t iT = 0; iT < nTrees; iT++)
	{
		for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
			// pT Bin Ranges
	    	pTBinLow = ptBins[iPtBin];
	    	pTBinUp  = ptBins[iPtBin+1];
		
			for(Int_t iMult = 0; iMult < nMultBins; iMult++){
				// Mult Bin Ranges
	    		multLow = multBinsLow[iMult];
	    		multUp  = multBinsUp [iMult];

				outputTrees[iT][iPtBin][iMult] = new TTree(Form("%sPt[%.0f-%.0f]Mult[%d-%d]",treeTitles[iT].Data(),pTBinLow,pTBinUp,multLow,multUp), Form("%sPt[%.0f-%.0f]Mult[%d-%d]",treeTitles[iT].Data(),pTBinLow,pTBinUp,multLow,multUp));
				
				for(Int_t iVar=0;iVar<nVar;iVar++){

					if(variableCategories[iVar]==kNotUsed) continue;

					 outputTrees[iT][iPtBin][iMult]->Branch(variableTreeNames[iVar].Data(),&variables[iVar],Form("%s/F",variableTreeNames[iVar].Data()));

				} // End iVar - loop

			} // End iMult-loop

		}	// End iPt-loop		

	} // End iT-loop

	// Fill the trees
	hSgnEntriesPtVsMult    = new TH2F("hSgnEntriesPtVsMult", "Number of Signal Entries / p_{T}-bin / mult-bin; p_{T}; Multiplicity",nPtBins,ptBins,nMultBins,-0.5,nMultBins-0.5);
	hBkgEntriesPtVsMult    = new TH2F("hBkgEntriesPtVsMult", "Number of Backg. Entries / p_{T}-bin / mult-bin; p_{T}; Multiplicity",nPtBins,ptBins,nMultBins,-0.5,nMultBins-0.5);
    for (Int_t iMult = 0; iMult < nMultBins; iMult++){
      hSgnEntriesPtVsMult->GetYaxis()->SetBinLabel(iMult+1, Form("[%d - %d]",multBinsLow[iMult],multBinsUp[iMult]));
      hBkgEntriesPtVsMult->GetYaxis()->SetBinLabel(iMult+1, Form("[%d - %d]",multBinsLow[iMult],multBinsUp[iMult]));
    }

	for (Int_t iT = 0; iT < nTrees; iT++){

		// Fill Histogram
		Long64_t nEntries = inputTrees[iT]->GetEntries(); //100000;//
		// if( maxEntries!=0 && nEntries > maxEntries) nEntries = maxEntries;

		Printf("\n--- %10.10s: %10lld candidates ",treeTitles[iT].Data(),nEntries);

		for(Long64_t iLc = 0; iLc <nEntries; iLc++){
			
			inputTrees[iT]->GetEntry(iLc);

			Float_t pTLc = variables[18];
			Float_t NtrkCorr = variables[28];

            // if(variables[10] <= -998) continue;
			
			for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){

				// pT Bin Ranges
				pTBinLow = ptBins[iPtBin];
				pTBinUp  = ptBins[iPtBin+1];

				if(pTLc > pTBinLow && pTLc <= pTBinUp){

					for(Int_t iMult = 0; iMult < nMultBins; iMult++){
						// Mult Bin Ranges
		    			multLow = multBinsLow[iMult];
		    			multUp  = multBinsUp [iMult];


		    			if(NtrkCorr >= multLow && NtrkCorr < multUp){
                            // if(iMult==2) Printf("NtrkCorr = %.0f",NtrkCorr);

		    				if (iLc % (nEntries/10) == 0) Printf("... Processing: %10lld \t %.0f%% \t pTLc = %6.2f (%6.2f-%6.2f), mult = %4.0f (%3.1d-%3.1d)",iLc,Double_t(iLc)/Double_t(nEntries)*100.,pTLc,pTBinLow,pTBinUp,NtrkCorr,multLow,multUp);

		    				outputTrees[iT][iPtBin][iMult]->Fill();
		    				if(treeTypes[iT]==kSgn) hSgnEntriesPtVsMult->Fill(pTLc,iMult);
		    				if(treeTypes[iT]==kBkg) hBkgEntriesPtVsMult->Fill(pTLc,iMult);

		    				// break;
		    				
		    			} // End: If multLow < pTLc <= multUp

	    			} // End iMult - loop

	    			break;

				} // End :If pTlow < pTLc <= pTUp

			} // End iPt-loop

		} // End iLc-loop

	} // End iT-loop


}


void PrepareSgnBkgTree( TString  inputSgnFile1,
				  		TString  inputBkgFile1,
				  		TString  inputSgnFile2 = "",
				  		TString  inputBkgFile2 = "",
				  		TString  inputSgnFile3 = "",
				  		TString  inputBkgFile3 = "",
				  		TString  fileNameAddition = ""
				  	   ){

	std::cout << std::endl << std::flush;
	Printf("========= PrepareSgnBkgTree.C Started =========");

	//================================================
  	//         Input
  	//================================================

	TTree*  inputTrees[6];
	TString treeTitles[6];
	Int_t  treeTypes[6] = {kUndef};

	Int_t nSgnTrees = 1;
	TString inputSgnDirs[3];
	
	Int_t nBkgTrees = 1;
	TString inputBkgDirs[3];
	
	inputSgnDirs[0]=inputSgnFile1;
	if( !(inputSgnFile2.EqualTo("")) ){inputSgnDirs[1]=inputSgnFile2; nSgnTrees++;}
	if( !(inputSgnFile3.EqualTo("")) ){inputSgnDirs[2]=inputSgnFile3; nSgnTrees++;}

	inputBkgDirs[0]=inputBkgFile1;
	if( !(inputBkgFile2.EqualTo("")) ){inputBkgDirs[1]=inputBkgFile2; nBkgTrees++;}
	if( !(inputBkgFile3.EqualTo("")) ){inputBkgDirs[2]=inputBkgFile3; nBkgTrees++;}

	for(Int_t i = 0; i < nSgnTrees;i++){
		TFile* inputFile(0);
		if (!gSystem->AccessPathName( inputSgnDirs[i].Data() )) {
		  inputFile = TFile::Open( inputSgnDirs[i].Data() ); // check if file in local directory exists
		}
		else {
		  Printf("ERROR: could not open input signal file %d: %s",i,inputSgnDirs[i].Data());
		  exit(1);
		}
		TTree* tTree = (TTree*)inputFile->Get("treeList_PfCuts_Sgn");
		if(!tTree){
			Printf("ERROR: could not open signal tree %d",i);
			exit(1);
		}
		inputTrees[i]   = tTree;
		treeTypes[i]    = kSgn;
		treeTitles[i]   = Form("Sgn[%d]",i);

		Printf("* Input signal file%d: %s", i, inputFile->GetName());
		Printf("  Tree: %s (%10lld entries)", tTree->GetName(), tTree->GetEntries());

	}

	for(Int_t i = 0; i < nBkgTrees;i++){
		TFile* inputFile(0);
		if (!gSystem->AccessPathName( inputBkgDirs[i].Data() )) {
		  inputFile = TFile::Open( inputBkgDirs[i].Data() ); // check if file in local directory exists
		}
		else {
		  Printf("ERROR: could not open input background file %d: %s",i,inputBkgDirs[i].Data());
		  exit(1);
		}
		TTree* tTree = (TTree*)inputFile->Get("treeList_Bkg");
		if(!tTree){
			Printf("ERROR: could not open background tree %d",i);
			exit(1);
		}
		inputTrees[nSgnTrees+i]   = tTree;
		treeTypes[nSgnTrees+i]    = kBkg;
		treeTitles[nSgnTrees+i]   = Form("Bkg[%d]",i);

		Printf("* Input background file%d: %s",i,inputFile->GetName());
		Printf("  Tree: %s (%10lld entries)",tTree->GetName(),tTree->GetEntries());
	}

	nTrees=nSgnTrees+nBkgTrees;

	//================================================
  	//         Start anaylsis
  	//================================================

	// Save output
	TFile *outputFile(0);
	if(fileNameAddition==""){
		outputFile = new TFile("SgnBkgTrees.root","RECREATE");
	}
	else{
		outputFile = new TFile(Form("SgnBkgTrees_%s.root",fileNameAddition.Data()),"RECREATE");
	}
  	if(!outputFile){
		Printf("ERROR: outputFile not created");
		exit(1);
	}
	outputFile->cd();

	// Set variables of trees
	SetVariables();

  	// Generate Trees
	GenerateSignalBackgroundTrees(inputTrees,treeTitles,treeTypes);

	hSgnEntriesPtVsMult->Write();
	hBkgEntriesPtVsMult->Write();
	for (Int_t iT = 0; iT < nTrees; iT++){
		for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
			for(Int_t iMult = 0; iMult < nMultBins; iMult++){
				outputTrees[iT][iPtBin][iMult]->Write();
			}
		}
	}
	
	outputFile->Close();

	return;
}