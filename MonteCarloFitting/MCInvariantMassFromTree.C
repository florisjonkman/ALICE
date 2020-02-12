// --------------------------------- Both settings
#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TMath.h"

#include <TDatabasePDG.h>

#include <fstream>
#include <iostream>
#include <cstdio>

// --------------------------------- ROOT6 settings
#if defined(__CLING__)

#include "AliVertexingHFUtils.h"
#include "AliHFMassFitter.h"

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

enum {kExpo=0, kLinear, kPol2,kPol3,kPow,kPowEx};
enum {kGaus=0, kDoubleGaus};

//===========================================================
//     Fit Parameters 
//===========================================================
	
// Bins: pT
const Int_t nPtBins = 6;
Double_t ptBins[nPtBins+1]= {1.,2.,4.,6.,8.,12.,24.};; // pT bin ranges

const Int_t nFiles = 1;
TTree* inputTrees[nFiles];
TH1D*  hMass[nPtBins];
TString inputFiles[nFiles];
Double_t ratios[nFiles];


void SetInputFiles(){

  inputFiles[0]="../../Results/GRID/Trees/20191025_GenerateTrees_LHC19h4c2_2016_kINT7_AllCuts_kNtrk10_refM12_25.root";
  // inputFiles[1]="/home/alidock/Lc/Results/GRID/Efficiencies/CorrectTPCCuts/Lc2pK0S_CutEffAllCuts_TMVAEffAllCuts_LHC18f4a_16k.root";
  // inputFiles[2]="/home/alidock/Lc/Results/GRID/Efficiencies/CorrectTPCCuts/Lc2pK0S_CutEffAllCuts_TMVAEffAllCuts_LHC18f4a_16l.root";
  // inputFiles[3]="/home/alidock/Lc/Results/GRID/Efficiencies/CorrectTPCCuts/Lc2pK0S_CutEffAllCuts_TMVAEffAllCuts_LHC18l4a_2017.root";
  // inputFiles[4]="/home/alidock/Lc/Results/GRID/Efficiencies/CorrectTPCCuts/Lc2pK0S_CutEffAllCuts_TMVAEffAllCuts_LHC18l4b_2018.root";
  // // inputFiles[0]="/home/alidock/Lc/Results/GRID/Efficiencies/TMVAEfficiency/Lc2pK0S_TMVAEff_allCuts_LHC17h8a.root";
  // inputFiles[1]="/home/alidock/Lc/Results/GRID/Efficiencies/TMVAEfficiency/Lc2pK0S_TMVAEff_allCuts_LHC18f4a_16k.root";
  // inputFiles[2]="/home/alidock/Lc/Results/GRID/Efficiencies/TMVAEfficiency/Lc2pK0S_TMVAEff_allCuts_LHC18f4a_16l.root";
  // inputFiles[3]="/home/alidock/Lc/Results/GRID/Efficiencies/TMVAEfficiency/Lc2pK0S_TMVAEff_allCuts_CutEff_PfCuts_LHC18l4a.root";
  // inputFiles[4]="/home/alidock/Lc/Results/GRID/Efficiencies/TMVAEfficiency/Lc2pK0S_TMVAEff_allCuts_CutEff_PfCuts_LHC18l4b.root";
  ratios[0] = 1.00;
  // ratios[1] = 1.00;
  // ratios[2] = 1.00;
  // ratios[3] = 1.00;
  // ratios[4] = 1.00;

 //  inputFiles[0]="/home/alidock/Lc/Results/GRID/GenerateTrees/20191025_GenerateTrees_LHC19h4c2_2016_kINT7_AllCuts_kNtrk10_refM12_25.root";
	// inputFiles[1]="/home/alidock/Lc/Results/GRID/GenerateTrees/20191025_GenerateTrees_LHC19h4b2_2017_kINT7_AllCuts_kNtrk10_refM12_25.root";
	// inputFiles[2]="/home/alidock/Lc/Results/GRID/GenerateTrees/20191025_GenerateTrees_LHC19h4a2_2018_kINT7_AllCuts_kNtrk10_refM12_25.root";
 //  ratios[0] = 0.56;
 //  ratios[1] = 0.26;
 //  ratios[2] = 0.17;
}

void SetInvariantMassHistograms(){

  for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
    // pT Bin Ranges
    Double_t pTBinLow = ptBins[iPtBin];
    Double_t pTBinUp  = ptBins[iPtBin+1];

    hMass[iPtBin] = new TH1D(Form("hMass_pT_%0.f_%0.f",pTBinLow,pTBinUp),Form("hMass_pT_%0.f_%0.f",pTBinLow,pTBinUp),1000, 2.05, 2.55);
    hMass[iPtBin]->SetTitle(Form("%.0f < #it{p}_{T} < %0.f",pTBinLow,pTBinUp));
    hMass[iPtBin]->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    hMass[iPtBin]->GetYaxis()->SetTitle(Form("Entries/(%.1f MeV/c^{2})",(hMass[iPtBin]->GetBinWidth(1)*1000)));
    hMass[iPtBin]->SetLineColor(kBlue+1);
    hMass[iPtBin]->SetMarkerColor(kBlack);
    hMass[iPtBin]->SetMarkerStyle(8);
    hMass[iPtBin]->SetMarkerSize(0.5);
    hMass[iPtBin]->GetYaxis()->SetTitleOffset(1.6);
  }

}

void FillInvariantMassPlotsFromTree(TTree* inputTree,TString treeName,Double_t ratio,Int_t multLow, Int_t multUp){

	// Initialize tree variables
  	const Int_t nVar = 4;
  	TString fCandidateVariableNames[nVar];
  	Float16_t fCandidateVariables[nVar];

  	fCandidateVariableNames[0] = "massLc2K0Sp";
    fCandidateVariableNames[1] = "LcPt";
    fCandidateVariableNames[2] = "NtrkRaw";
    fCandidateVariableNames[3] = "origin";

    for(Int_t ivar=0; ivar<nVar; ivar++){
    	inputTree->SetBranchAddress(Form("%s",fCandidateVariableNames[ivar].Data()), &fCandidateVariables[ivar]);
    }

    // Fill histograms
    Long64_t nEntries = inputTree->GetEntries();// 100000; //inputTree->GetEntries();
    Printf("\n--- Total number of candidates: %10lld candidates",nEntries);

    nEntries = nEntries*ratio;
    Printf("\n--- Ratio number of candidates: %10lld candidates",nEntries);

	for(Long64_t iLC = 0; iLC <nEntries; iLC++){
		// Get Entry
		inputTree->GetEntry(iLC);
		if (iLC % (nEntries/10) == 0) Printf("... Processing: %10lld \t %.0f%%",iLC,Double_t(iLC)/Double_t(nEntries)*100.);
		
    Float_t massLc  = (Float_t) fCandidateVariables[0];
    Float_t pTLc    = (Float_t) fCandidateVariables[1];
    Float_t nTrkRaw = (Float_t) fCandidateVariables[2];
    Float_t origin  = (Float_t) fCandidateVariables[3];

    if(nTrkRaw >= multLow && nTrkRaw < multUp){
      for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
        // pT Bin Ranges
        Double_t pTBinLow = ptBins[iPtBin];
        Double_t pTBinUp  = ptBins[iPtBin+1];

        if(pTLc>=pTBinLow && pTLc <=pTBinUp){

          if(origin == 4 || origin == 5) 	
            hMass[iPtBin]->Fill(massLc);
          break;
        }
      }
    }

		
	}

	return;

}

void MCInvariantMassFromTree(
  	Int_t    multLow = 0,
  	Int_t    multUp = 999,
	TString outputDirName = "Output"
	)
{

	//===========================================================
	//     Reading all information and creating output
	//===========================================================

	SetInputFiles();

	// Read input file
	for (Int_t i = 0; i < nFiles; ++i)
	{
		TFile *inputFile(0);
		if ( !gSystem->AccessPathName(inputFiles[i].Data()) ){
			inputFile  	= new TFile( inputFiles[i].Data() );
		}
		else{
			Printf("ERROR: could not open input root file %d",i);
			exit(1);
		}

		// Get Tree from inputfile
		inputTrees[i] = (TTree*)inputFile->Get("treeList_PfCuts_Sgn");
		if(!inputTrees[i]){
		Printf("ERROR: treeList_Sgn not in inputFile");
		exit(1);
		}

    // inputFile->Close();
	}
	
	// // Get Tree from inputfile
	// TTree *inputTreeBkg = (TTree*)inputFile->Get("treeList_stdCuts_Bkg");
	// if(!inputTreeBkg){
	// Printf("ERROR: treeList_Bkg not in inputFile");
	// exit(1);
	// }

	// Open output for Invariantmass file 
	// gSystem->Exec(Form("rm -r %s/",outputDirName.Data()));
	// gSystem->Exec(Form("mkdir %s",outputDirName.Data()));
	TFile *outputFile(0);
	outputFile = new TFile(Form("MCInvariantMass_mult_%d_%d.root",multLow,multUp),"RECREATE");
	if(!outputFile){
		Printf("ERROR: outputFile not created");
	exit(1);
	}
	outputFile->cd();

  SetInvariantMassHistograms();

	FillInvariantMassPlotsFromTree(inputTrees[0],"treeList_PfCuts_Sgn",ratios[0],multLow,multUp);
  // FillInvariantMassPlotsFromTree(inputTrees[1],"treeListTMVAEff_stdCuts_Sgn",ratios[1],multLow,multUp);
  // FillInvariantMassPlotsFromTree(inputTrees[2],"treeListTMVAEff_stdCuts_Sgn",ratios[2],multLow,multUp);
  // FillInvariantMassPlotsFromTree(inputTrees[3],"treeListTMVAEff_stdCuts_Sgn",ratios[3],multLow,multUp);
  // FillInvariantMassPlotsFromTree(inputTrees[4],"treeListTMVAEff_stdCuts_Sgn",ratios[4],multLow,multUp);


	outputFile->cd();

  for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
    hMass[iPtBin]->Write();
  }
	outputFile->Close();


}