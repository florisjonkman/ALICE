#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TMath.h"

#include <fstream>
#include <iostream>
#include <cstdio>

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */

enum variableCategories {kIsUsed=0,kNotUsed};

// Initialize Variables settings
const Int_t nVar = 36;
TString  variableTreeNames[nVar];
TString  variableGraphNames[nVar];
Int_t    variableCategories[nVar];
TString  variableUnits[nVar];
Double_t variableRanges[nVar][2];
Int_t    variableNBins[nVar];
Bool_t	 variableFlagLogScale[nVar];

// Bins settings
const Int_t nPtBins = 6;
Double_t ptBins[nPtBins+1]= {1.,2.,4.,6.,8.,12.,24.};; // pT bin ranges

// Other global objects
// TH1F** hVariablePtTree = new TH1F*[nVar][nPtBins][3];
TH1F* hVariablePtTree[nVar][nPtBins][2];
TH1F* hVariablePtTreeNorm[nVar][nPtBins][2];
Int_t treeColors[2] = {kRed,kBlue};

TH1F* hVariablePtTreeSideBands[nVar][nPtBins][2];
TH1F* hVariablePtTreeNormSideBands[nVar][nPtBins][2];
Int_t sideBandsColors[2] = {kGreen+2,kOrange+2};
Int_t nSideBands = 2;

void SetVariables(){

	Printf("Variable settings");

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
    variableTreeNames[35] = "ctau";

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
    variableCategories[31] = kNotUsed; // Not used
    variableCategories[32] = kIsUsed;
    variableCategories[33] = kIsUsed; // Not used
    variableCategories[34] = kIsUsed;
    variableCategories[35] = kIsUsed;

    // Variable names for graph title and x-axjs
    variableGraphNames[0] = "M(K^{0}_{S}p)";
    variableGraphNames[1] = "#alpha_{Arm}";
    variableGraphNames[2] = "M(K^{0}_{S})";
    variableGraphNames[3] = "M(#bar{#Lambda})";
    variableGraphNames[4] = "M(#Lambda)";
    variableGraphNames[5] = "cos(PA)";
    variableGraphNames[6] = "DCA(V0)";
    variableGraphNames[7] = "d_{0}(bach.)";
    variableGraphNames[8] = "d_{0}(V0)";
    variableGraphNames[9] = "n#sigma_{TPC}(p)";
    variableGraphNames[10] = "n#sigma_{TOF}(proton)";
    variableGraphNames[11] = "p_{T}(bach.)";
    variableGraphNames[12] = "p_{T}(#pi+)";
    variableGraphNames[13] = "p_{T}(#pi-)";
    variableGraphNames[14] = "d_{0}(#pi+)";
    variableGraphNames[15] = "d_{0}(#pi-)";
    variableGraphNames[16] = "p_{T}(V0)";
    variableGraphNames[17] = "p_{TPC}(bach.)";
    variableGraphNames[18] = "p_{T}(#Lambda_{c})";
    variableGraphNames[19] = "combProtonProb";
    variableGraphNames[20] = "#eta(#pi+)";
    variableGraphNames[21] = "p(bach.)";
    variableGraphNames[22] = "#eta(bach.)";
    variableGraphNames[23] = "p(V0)";
    variableGraphNames[24] = "L(V0)";
    variableGraphNames[25] = "n#sigma_{TPC}(#pi)";
    variableGraphNames[26] = "n#sigma_{TPC}(K)";
    variableGraphNames[27] = "Ntrk_{Raw}";
    variableGraphNames[28] = "Ntrk_{Corr}";
    variableGraphNames[29] = "Cos(#theta*)";
    variableGraphNames[30] = "signed d_{0}";        
    variableGraphNames[31] = "Centrality"; 
    variableGraphNames[32] = "Ntrk_{All}";
    variableGraphNames[33] = "Origin";
    variableGraphNames[34] = "p_{T}^{Arm}";
    variableGraphNames[35] = "c#tau";

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
    variableUnits[35] = "cm";

    // Variable ranges for plot
    // Lowerlimit						// Upperlimit
	variableRanges[0][0]  =  2.050;		variableRanges[0][1]  =  2.500;
    variableRanges[1][0]  = -1.500;		variableRanges[1][1]  =  1.500;	
    variableRanges[2][0]  =  0.450;		variableRanges[2][1]  =  0.550;
    variableRanges[3][0]  =  1.000;		variableRanges[3][1]  =  3.000;
    variableRanges[4][0]  =  1.000;		variableRanges[4][1]  =  3.000;
    variableRanges[5][0]  =  0.997;		variableRanges[5][1]  =  1.0005;
    variableRanges[6][0]  = -0.005;		variableRanges[6][1]  =  0.805;
    variableRanges[7][0]  = -0.200;		variableRanges[7][1]  =  0.200;
    variableRanges[8][0]  = -1.200;		variableRanges[8][1]  =  1.200;
    variableRanges[9][0]  =-10.000;		variableRanges[9][1]  = 10.000;
    variableRanges[10][0] =-41.000;		variableRanges[10][1] = 10.000;
    variableRanges[11][0] =  0.000;		variableRanges[11][1] = 36.000;
    variableRanges[12][0] =  0.000;		variableRanges[12][1] = 16.000;
    variableRanges[13][0] =  0.000;		variableRanges[13][1] = 16.000;
    variableRanges[14][0] =  0.000;		variableRanges[14][1] = 25.000;
    variableRanges[15][0] =  0.000;		variableRanges[15][1] = 25.000;
    variableRanges[16][0] =  0.000;		variableRanges[16][1] = 24.000;
    variableRanges[17][0] =  0.000;		variableRanges[17][1] = 30.000;
    variableRanges[18][0] =  0.000;		variableRanges[18][1] = 36.000;
    variableRanges[19][0] = -0.020;		variableRanges[19][1] =  1.020;
    variableRanges[20][0] = -0.840;		variableRanges[20][1] =  0.840;
    variableRanges[21][0] =  0.000;		variableRanges[21][1] = 25.000;
    variableRanges[22][0] = -0.840;		variableRanges[22][1] =  0.840;
    variableRanges[23][0] =  0.000;		variableRanges[23][1] =  30.00;
    variableRanges[24][0] =-20.000;		variableRanges[24][1] = 200.00;
    variableRanges[25][0] =-20.000;		variableRanges[25][1] = 100.00;
    variableRanges[26][0] = -5.000;		variableRanges[26][1] = 65.000;
    variableRanges[27][0] =  0.000;		variableRanges[27][1] = 250.00;
    variableRanges[28][0] =  0.000;		variableRanges[28][1] = 250.00;
    variableRanges[29][0] = -1.000;		variableRanges[29][1] =  1.000;
    variableRanges[30][0] = -0.200;		variableRanges[30][1] =  0.200;
    variableRanges[31][0] =  0.000;		variableRanges[31][1] = 100.00;
    variableRanges[32][0] =  0.000;		variableRanges[32][1] = 250.00;
    variableRanges[33][0] = -2.000;		variableRanges[33][1] =  6.000;
    variableRanges[34][0] =  0.000;		variableRanges[34][1] =  0.240;
    variableRanges[35][0] =  0.000;     variableRanges[35][1] =  100. ;

    // Variable number of bins
	variableNBins[0]  = 100;
    variableNBins[1]  = 40;
    variableNBins[2]  = 40;
    variableNBins[3]  = 40;
    variableNBins[4]  = 40;
    variableNBins[5]  = 40;
    variableNBins[6]  = 40;
    variableNBins[7]  = 50;
    variableNBins[8]  = 60;
    variableNBins[9]  = 60;
    variableNBins[10] = 60;
    variableNBins[11] = 36;
    variableNBins[12] = 16;
    variableNBins[13] = 16;
    variableNBins[14] = 40;
    variableNBins[15] = 40;
    variableNBins[16] = 24;
    variableNBins[17] = 40;
    variableNBins[18] = 36;
    variableNBins[19] = 64;
    variableNBins[20] = 50;
    variableNBins[21] = 25;
    variableNBins[22] = 50;
    variableNBins[23] = 30;
    variableNBins[24] = 40;
    variableNBins[25] = 40;
    variableNBins[26] = 40;
    variableNBins[27] = 40;
    variableNBins[28] = 40;
    variableNBins[29] = 40;
    variableNBins[30] = 60;
    variableNBins[31] = 40;
    variableNBins[32] = 40;
    variableNBins[33] = 8;
    variableNBins[34] = 40;
    variableNBins[35] = 50;

    // flag if y-axis is LogScale
	variableFlagLogScale[0]  = kFALSE;
    variableFlagLogScale[1]  = kFALSE;
    variableFlagLogScale[2]  = kFALSE;
    variableFlagLogScale[3]  = kFALSE;
    variableFlagLogScale[4]  = kFALSE;
    variableFlagLogScale[5]  = kTRUE;
    variableFlagLogScale[6]  = kFALSE;
    variableFlagLogScale[7]  = kTRUE;
    variableFlagLogScale[8]  = kTRUE;
    variableFlagLogScale[9]  = kFALSE;
    variableFlagLogScale[10] = kFALSE;
    variableFlagLogScale[11] = kFALSE;
    variableFlagLogScale[12] = kFALSE;
    variableFlagLogScale[13] = kFALSE;
    variableFlagLogScale[14] = kFALSE;
    variableFlagLogScale[15] = kFALSE;
    variableFlagLogScale[16] = kFALSE;
    variableFlagLogScale[17] = kFALSE;
    variableFlagLogScale[18] = kFALSE;
    variableFlagLogScale[19] = kTRUE;
    variableFlagLogScale[20] = kFALSE;
    variableFlagLogScale[21] = kFALSE;
    variableFlagLogScale[22] = kFALSE;
    variableFlagLogScale[23] = kFALSE;
    variableFlagLogScale[24] = kTRUE;
    variableFlagLogScale[25] = kFALSE;
    variableFlagLogScale[26] = kFALSE;
    variableFlagLogScale[27] = kFALSE;
    variableFlagLogScale[28] = kFALSE;
    variableFlagLogScale[29] = kFALSE;
    variableFlagLogScale[30] = kTRUE;
    variableFlagLogScale[31] = kFALSE;
    variableFlagLogScale[32] = kFALSE;
    variableFlagLogScale[33] = kFALSE;
    variableFlagLogScale[34] = kFALSE;
    variableFlagLogScale[35] = kTRUE;


    Printf("%10.10s \t %10.10s \t %5.5s \t %5.5s \t %5.5s \t %10.10s","treeName","varGraphName","Lower","Upper","nBins","flagLogScale");
    Printf("------------------------------------------------------------------------------------");
    for(Int_t iVar=0;iVar<nVar;iVar++){
    	Printf("%10.10s \t %10.10s \t %5.3f \t %5.3f \t %5d \t %10d",variableTreeNames[iVar].Data(),variableGraphNames[iVar].Data(),variableRanges[iVar][0],variableRanges[iVar][1],variableNBins[iVar],variableFlagLogScale[iVar]);
    }
    Printf("");

}


void PlotVariables(TTree **trees,TString* treeNames,TString* treeTitles,const Int_t nTrees, Bool_t flagDrawCanvas, Bool_t flagKeepCanvasOpen, Bool_t flagSideBandsSeperate){

	Printf("\nFill histograms (and draw canvasses)");
	Printf("-----------------");

	// Variables needed
	Double_t pTBinLow;
    Double_t pTBinUp;

	// Create Canvasses
    TCanvas* canvasVar[nVar];

  	// Initializing all histograms, 1 hist / pTbin / var / tree
  	for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
		// pT Bin Ranges
    	pTBinLow = ptBins[iPtBin];
    	pTBinUp  = ptBins[iPtBin+1];

    	for(Int_t iVar=0;iVar<nVar;iVar++){

    		if(variableCategories[iVar]==kNotUsed) continue;

    		for(Int_t iT=0; iT<nTrees; iT++){

    			hVariablePtTree[iVar][iPtBin][iT] =new TH1F(Form("h%sPt[%.2f-%.2f]T[%s]",variableTreeNames[iVar].Data(),pTBinLow,pTBinUp,treeTitles[iT].Data()),
    											   Form("h%sPt[%.2f-%.2f]T[%s]",variableTreeNames[iVar].Data(),pTBinLow,pTBinUp,treeTitles[iT].Data()),
    											   variableNBins[iVar],
    											   variableRanges[iVar][0],
    											   variableRanges[iVar][1]);

    		} // End iT-loop

            if(flagSideBandsSeperate){
                for(Int_t iSB=0; iSB<nSideBands; iSB++){
                    hVariablePtTreeSideBands[iVar][iPtBin][iSB] =new TH1F(Form("h%sPt[%.2f-%.2f]T[BkgSB%d]",variableTreeNames[iVar].Data(),pTBinLow,pTBinUp,iSB),
                                                   Form("h%sPt[%.2f-%.2f]T[BkgSB%d]",variableTreeNames[iVar].Data(),pTBinLow,pTBinUp,iSB),
                                                   variableNBins[iVar],
                                                   variableRanges[iVar][0],
                                                   variableRanges[iVar][1]);

                } // End Side Band loop

            } // End if SideBands

    	} // End iVar-loop

	} // End iPt-loop


	// Settings branches of Trees
	Float16_t pTLc;
    Float16_t massLc2K0Sp;
	Float16_t variables[nVar];

	for (Int_t iT = 0; iT < nTrees; iT++)
	{
		for(Int_t iVar=0;iVar<nVar;iVar++){
            if(iVar!=35) trees[iT]->SetBranchAddress(variableTreeNames[iVar].Data(),&variables[iVar]);
		} // End iVar-loop

	} // End iT-loop

	// Filling the histograms
	for (Int_t iT = 0; iT < nTrees; iT++)
	{
		// Fill Histogram
		Long64_t nEntries = trees[iT]->GetEntries(); //100000;// 

		for(Int_t iLc = 0; iLc <nEntries; iLc++){
			
			trees[iT]->GetEntry(iLc);
            pTLc=variables[18];
            massLc2K0Sp=variables[0];

			for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){

				// pT Bin Ranges
				pTBinLow = (Double_t) ptBins[iPtBin];
				pTBinUp  = (Double_t) ptBins[iPtBin+1];

				if(pTLc > pTBinLow && pTLc <= pTBinUp){

					if(iLc % Int_t(nEntries/5) == 0){
    						Printf("pTLc = %.2f (%.2f-%.2f)   %s",pTLc,pTBinLow,pTBinUp,treeTitles[iT].Data());
    				}

					for(Int_t iVar=0;iVar<nVar;iVar++){

						if(variableCategories[iVar]==kNotUsed) continue;
					
                        if(iVar!=35){
                            hVariablePtTree[iVar][iPtBin][iT]->Fill( (Float_t) variables[iVar] );
                            if(iLc % Int_t(nEntries/5) == 0){
                                Printf("   %10.10s = %.3f",variableTreeNames[iVar].Data(),variables[iVar]);
                            }

                            if(iT==0 && flagSideBandsSeperate){
                                if(massLc2K0Sp < 2.289){ hVariablePtTreeSideBands[iVar][iPtBin][0]->Fill( (Float_t) variables[iVar] );}
                                else{                    hVariablePtTreeSideBands[iVar][iPtBin][1]->Fill( (Float_t) variables[iVar] );}
                            }
                        }
						else{
                            Float_t ctau = ( Float_t(variables[24])*0.497 / Float_t(variables[23])  );
                            
                            hVariablePtTree[iVar][iPtBin][iT]->Fill( ctau );
                            if(iLc % Int_t(nEntries/5) == 0){
                                // Printf("DecayLength = %f, v0P = %f, ctau = %f",variables[24],variables[23],ctau);
                                Printf("   %10.10s = %.3f",variableTreeNames[iVar].Data(),ctau);

                            }

                            if(iT==0 && flagSideBandsSeperate){
                                if(massLc2K0Sp < 2.289){ hVariablePtTreeSideBands[iVar][iPtBin][0]->Fill( ctau );}
                                else{                    hVariablePtTreeSideBands[iVar][iPtBin][1]->Fill( ctau );}
                            }

                        }

					} // End iVar loop

					break;

    			} // End if statement pTLc range

			} // End iPt-loop

		} // End iLc-loop

	} // End iT-loop

	Printf("All histograms filled!");

	// Drawing the histograms
	for(Int_t iVar=0;iVar<nVar;iVar++){

		if(variableCategories[iVar]==kNotUsed) continue;

		if(flagDrawCanvas){
			canvasVar[iVar] = new TCanvas(Form("c%s",variableTreeNames[iVar].Data()), Form("Canvas %s",variableTreeNames[iVar].Data()), 1200, 550);
			if(nPtBins<=6){
  				canvasVar[iVar]->Divide(3,2);
  			}
  			else{
  				canvasVar[iVar]->Divide(4,2);
  			}
	  	}

        Bool_t isLegDrawn = kFALSE;

		for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
			// pT Bin Ranges
    		pTBinLow = ptBins[iPtBin];
    		pTBinUp  = ptBins[iPtBin+1];

			TLegend *leg;
            if(flagSideBandsSeperate){ leg = new TLegend(0.0982102, 0.779231, 0.623276, 0.900355);}
            else{leg = new TLegend(0.0982102, 0.811531, 0.620994, 0.900355);}
			leg->SetFillColor(gPad->GetFillColor());
			leg->SetNColumns(2);

			Double_t lmax = 0;
			Double_t gmax = 0;
			Bool_t flagDrawHist[2]; 

			for (Int_t iT = 0; iT < nTrees; iT++){

				hVariablePtTree[iVar][iPtBin][iT]->SetMarkerColor(treeColors[iT]);
				hVariablePtTree[iVar][iPtBin][iT]->SetFillColorAlpha(treeColors[iT], 0.1);
				hVariablePtTree[iVar][iPtBin][iT]->SetLineColor(treeColors[iT]);
				hVariablePtTree[iVar][iPtBin][iT]->SetLineWidth(1);
				hVariablePtTree[iVar][iPtBin][iT]->SetTitle(Form("%s | %.2f < p_{T} < %.2f",variableGraphNames[iVar].Data(),pTBinLow,pTBinUp));
				hVariablePtTree[iVar][iPtBin][iT]->GetXaxis()->SetTitle(Form("%s [%s]",variableGraphNames[iVar].Data(),variableUnits[iVar].Data()));
  				hVariablePtTree[iVar][iPtBin][iT]->GetYaxis()->SetTitle("Entries");

				Double_t hEntries = hVariablePtTree[iVar][iPtBin][iT]->GetEntries();
				
				if(hEntries > 0){
					Double_t integral = hVariablePtTree[iVar][iPtBin][iT]->Integral();
					hVariablePtTree[iVar][iPtBin][iT]->Sumw2();
					if(integral>0){
						hVariablePtTreeNorm[iVar][iPtBin][iT] = (TH1F*) hVariablePtTree[iVar][iPtBin][iT]->Clone();
						hVariablePtTreeNorm[iVar][iPtBin][iT]->SetName(Form("h%sPt[%.2f-%.2f]T[%s]Norm",variableTreeNames[iVar].Data(),pTBinLow,pTBinUp,treeTitles[iT].Data()));
						hVariablePtTreeNorm[iVar][iPtBin][iT]->Scale(1/integral);
						hVariablePtTreeNorm[iVar][iPtBin][iT]->GetYaxis()->SetTitle("Normalized counts");                        
						lmax = hVariablePtTreeNorm[iVar][iPtBin][iT]->GetMaximum();
						if(lmax > gmax ) gmax = lmax;
						flagDrawHist[iT]=kTRUE;

                        if(iT==0 && flagSideBandsSeperate){
                            for(Int_t iSB=0; iSB<nSideBands; iSB++){
                                Double_t hEntriesSB = hVariablePtTreeSideBands[iVar][iPtBin][iSB]->GetEntries();
                                if(hEntries > 0){
                                    Double_t integralSB = hVariablePtTreeSideBands[iVar][iPtBin][iSB]->Integral();
                                    hVariablePtTreeSideBands[iVar][iPtBin][iSB]->Sumw2();
                                if(integralSB > 0){
                                        hVariablePtTreeNormSideBands[iVar][iPtBin][iSB] = (TH1F*) hVariablePtTreeSideBands[iVar][iPtBin][iSB]->Clone();
                                        hVariablePtTreeNormSideBands[iVar][iPtBin][iSB]->SetName(Form("h%sPt[%.2f-%.2f]T[BkgSB%d]Norm",variableTreeNames[iVar].Data(),pTBinLow,pTBinUp,iSB));
                                        if(iVar > 0){
                                            hVariablePtTreeNormSideBands[iVar][iPtBin][iSB]->Scale(1/integralSB);
                                        }
                                        else{
                                            hVariablePtTreeNormSideBands[iVar][iPtBin][iSB]->Scale(1/integral);
                                        }
                                    }
                                }
                            }
                        }// End Sidebands flag

					}
					else{
						flagDrawHist[iT]=kFALSE;
					}

				} // End if entries

			} // End iT-loop	

			Int_t nDrawn = 0;
			if(flagDrawCanvas){

				for (Int_t iT = 0; iT < nTrees; iT++){

					canvasVar[iVar]->cd(iPtBin+1);

                    if(!flagDrawHist[iT]) continue;

					if(variableFlagLogScale[iVar]){
						gPad-> SetLogy(variableFlagLogScale[iVar]);
						hVariablePtTreeNorm[iVar][iPtBin][iT]->GetYaxis()->SetRangeUser(0.00001,gmax*15.00);
					}
					else{
						hVariablePtTreeNorm[iVar][iPtBin][iT]->GetYaxis()->SetRangeUser(0.0,gmax*1.35);
					}
					 

					if(nDrawn==0){
                        gStyle->SetOptStat("ne");
						hVariablePtTreeNorm[iVar][iPtBin][iT]->Draw("Hist");
						hVariablePtTreeNorm[iVar][iPtBin][iT]->Draw("E0 same");
                        // gPad->Update();
                        // TPaveStats *st = (TPaveStats*)hVariablePtTreeNorm[iVar][iPtBin][iT]->FindObject("stats");
                        // st->SetX1NDC(0.568244);
                        // st->SetX2NDC(0.705709);
                        // gPad->Update();

                        if(iT==0 && flagSideBandsSeperate){
                            for(Int_t iSB=0; iSB<nSideBands; iSB++){
                                hVariablePtTreeNormSideBands[iVar][iPtBin][iSB]->SetMarkerColor(sideBandsColors[iSB]);
                                // hVariablePtTreeNormSideBands[iVar][iPtBin][iT]->SetFillColorAlpha(treeColors[iT], 0.1);
                                hVariablePtTreeNormSideBands[iVar][iPtBin][iSB]->SetLineColor(sideBandsColors[iSB]);
                                hVariablePtTreeNormSideBands[iVar][iPtBin][iSB]->SetMarkerStyle(8);
                                hVariablePtTreeNormSideBands[iVar][iPtBin][iSB]->SetMarkerSize(0.5);
                                hVariablePtTreeNormSideBands[iVar][iPtBin][iSB]->Draw("E0SAME");
                            }   
                        }

					}
					else{
                        gStyle->SetOptStat("ne");
						hVariablePtTreeNorm[iVar][iPtBin][iT]->Draw("SAMES hist");
						hVariablePtTreeNorm[iVar][iPtBin][iT]->Draw("E0 same");
						gPad->Update();
						TPaveStats *st = (TPaveStats*)hVariablePtTreeNorm[iVar][iPtBin][iT]->FindObject("stats");
						// st->SetX1NDC(0.768244);
						// st->SetX2NDC(0.955709);
						st->SetY1NDC(0.755725);
						st->SetY2NDC(0.853355);
						gPad->Update();
						
					}
					leg->AddEntry(hVariablePtTreeNorm[iVar][iPtBin][iT], treeTitles[iT].Data(), "f");
					nDrawn++;
					
				} // End iT-loop

                if(flagSideBandsSeperate){
                    leg->AddEntry(hVariablePtTreeNormSideBands[iVar][iPtBin][0], "Bkg, Data, Left SB", "P");
                    leg->AddEntry(hVariablePtTreeNormSideBands[iVar][iPtBin][1], "Bkg, Data, Right SB", "P");
                }
                
				if(nDrawn > 0 && !isLegDrawn){
					leg->Draw();
					gPad->Update();
                    isLegDrawn=kTRUE;
				}

			} // End if flagDrawCanvas
	
		} // End iPt-loop

		if(flagDrawCanvas){
			canvasVar[iVar]->SaveAs(Form("Output/%s.pdf",variableTreeNames[iVar].Data()));
			if(!flagKeepCanvasOpen) canvasVar[iVar]->Close();
		}

	} // End iVar loop

}


void DrawVariablesFromTree(
	const char*  dirTreeFile1,
	const char*  treeName1,
	const char*  treeTitle1,
	const char*  dirTreeFile2 = "",
	const char*  treeName2 = "",
	const char*  treeTitle2 = "",
	const char*  outputFileName = "output",
	const Bool_t flagSaveHist = kTRUE,
	const Bool_t flagDrawCanvas = kTRUE,
	const Bool_t flagKeepCanvasOpen =kFALSE,
    const Bool_t flagSideBandsSeperate = kFALSE
	)
{

	Int_t   treeCounter = 0;
	TTree*  inputTrees[2];
	TString treeNames [2];
	TString treeTitles[2];

	//===========================================================
	//     Reading all information and creating output
	//===========================================================

	// Read input file 1, must always be given
	TFile *inputFile1(0);

	if (!gSystem->AccessPathName(dirTreeFile1 )){
		inputFile1  	= new TFile(dirTreeFile1);
	}
	else{
		Printf("ERROR: could not open input root file 1");
		exit(1);
	}

	TTree *inputTree1 = (TTree*)inputFile1->Get(treeName1);
	if(!inputTree1){
		Printf("ERROR: %s not in inputFile1",treeName1);
		exit(1);
	}
	inputTrees[0] = inputTree1;
	treeNames[0]  = treeName1;
	treeTitles[0] = treeTitle1;
	treeCounter++;

	// Read input file 2, only if given
	TTree *inputTree2;
	TFile *inputFile2(0);

	if(!(dirTreeFile2=="")){
		if (!gSystem->AccessPathName(dirTreeFile2)){
			inputFile2 = new TFile(dirTreeFile2);
		}
		else{
			Printf("ERROR: could not open input root file 2");
			exit(1);
		}
		inputTree2 = (TTree*)inputFile2->Get(treeName2);
		if(!inputTree2){
			Printf("ERROR: %s not in inputFile2",treeName2);
			exit(1);
		}
		else{
			inputTrees[1] = inputTree2;
			treeNames[1]  = treeName2;
			treeTitles[1] = treeTitle2;
			treeCounter++;
		}
	}
	
	// Open output file 
	gSystem->Exec("rm -r Output/");
	gSystem->Exec("mkdir Output");
 	TFile *outputFile(0);
  	outputFile = new TFile(Form("Output/%s.root",outputFileName),"RECREATE");
  	if(!outputFile){
		Printf("ERROR: outputFile not created");
		exit(1);
	}

	
	Printf("Information");
	Printf("\n-----------------");
	Printf("* nTrees = %d",treeCounter);
	for (Int_t iT = 0; iT < treeCounter; iT++){
		Printf("* iT = %d, %s, %s",iT, treeTitles[iT].Data(), treeNames[iT].Data());
		if (iT==0) Printf("  Path: %s",dirTreeFile1);
		if (iT==1) Printf("  Path: %s",dirTreeFile2);
	}
	Printf("* Outputfile: %s.root",outputFileName);
	Printf("\n");

	// Set all settings for the histograms
	SetVariables();

	// Draw all variables
	PlotVariables(inputTrees,treeNames,treeTitles,treeCounter,flagDrawCanvas,flagKeepCanvasOpen,flagSideBandsSeperate);	
	
	inputFile1->Close();
	inputFile2->Close();
	
	// Saveing output
	TDirectory* dirVar[nVar];
	TDirectory* dirVarTrees[nVar][2];

	for(Int_t iVar=0;iVar<nVar;iVar++){

		if(variableCategories[iVar]==kNotUsed) continue;
		
		outputFile->cd();
		for (Int_t iT = 0; iT < treeCounter; iT++){
			dirVarTrees[iVar][iT] = outputFile->mkdir(Form("%s/%s",variableTreeNames[iVar].Data(),treeTitles[iT].Data()));
			dirVarTrees[iVar][iT]->cd(Form("%s",treeTitles[iT].Data()));
			for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
                if(hVariablePtTree[iVar][iPtBin][iT]) hVariablePtTree[iVar][iPtBin][iT]->Write();
				if(hVariablePtTreeNorm[iVar][iPtBin][iT]) hVariablePtTreeNorm[iVar][iPtBin][iT]->Write();
			}
        }
	}

    for(Int_t iVar=0;iVar<nVar;iVar++){

        if(variableCategories[iVar]==kNotUsed) continue;
        outputFile->cd();
        dirVarTrees[iVar][0]->cd(Form("%s",treeTitles[0].Data()));

        for(Int_t iSB=0; iSB<nSideBands; iSB++){
            for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
                if(hVariablePtTreeNormSideBands[iVar][iPtBin][iSB]) hVariablePtTreeNormSideBands[iVar][iPtBin][iSB]->Write();
            }
        }
    }

	outputFile->cd();

	outputFile->Close();

}