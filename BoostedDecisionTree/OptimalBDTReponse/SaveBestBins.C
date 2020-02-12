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

enum {kExpo=0, kLinear, kPol2, kNoBk,kPowEx,kThrExpo};
enum {kGaus=0, kDoubleGaus};

// Bins: pT
const Int_t nPtBins = 6;
Double_t ptBins[nPtBins+1]= {1.,2.,4.,6.,8.,12.,24.};; // pT bin ranges
Int_t bestBDTbin[nPtBins];

// Fit Paramters
Int_t factor4refl = 0;                // 0 = gaus; 1 = gaus+gaus broadened
Int_t nSigmaCount = 3;

Double_t minMassForPlot = 2.050;
Double_t maxMassForPlot = 2.550;

// Intialize other objects
Double_t massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
TString gPlotDir;
Int_t gPtLow;
Int_t gPtUp;
Int_t gMultLow;
Int_t gMultUp;
Int_t gBinMin;
Int_t gBinMax;
Double_t gSigmaMC;
Double_t gMeanMC;

TH2D* hListBDTHistoTMVA;

TH1F* GetInvariantMassHistogram(TH2D* hBDTHistoTMVA, Int_t binBdtMin,Double_t &bdtValue, Int_t rebin = 1){

  TH1F* hInvariantMass = (TH1F*) hBDTHistoTMVA->ProjectionY(Form("hInvariantMass_pT_%d_%d_mult_%d_%d",gPtLow,gPtUp,gMultLow,gMultUp),binBdtMin,-1);
  TH1F* hRebinned   = (TH1F*) AliVertexingHFUtils::RebinHisto(hInvariantMass,rebin,-1);
  // hRebinned->SetAxisRange(minMassForPlot,maxMassForPlot);
  hRebinned->SetName(Form("hRebinned_pT_%d_%d_mult_%d_%d_binbdt_%d",gPtLow,gPtUp,gMultLow,gMultUp,binBdtMin));
  bdtValue = (hBDTHistoTMVA->GetXaxis())->GetBinLowEdge(binBdtMin);
  hRebinned->SetTitle(Form("bdt (bin) #geq %.4f (%d)",bdtValue,binBdtMin));
  hRebinned->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
  hRebinned->GetYaxis()->SetTitle(Form("Entries/(%.1f MeV/c^{2})",(hRebinned->GetBinWidth(1)*1000)));
  hRebinned->SetLineColor(kBlue+1);
  hRebinned->SetMarkerColor(kBlack);
  hRebinned->SetMarkerStyle(8);
  hRebinned->SetMarkerSize(0.5);
  hRebinned->GetYaxis()->SetTitleOffset(1.6);

  return hRebinned;

}

void SaveBestBins(Int_t multLow,
                  Int_t multUp,
                  TString fileName
                  ){

  std::cout << std::endl << std::flush;
  Printf("========= SaveBestBins.C Started =========");

  gMultLow = multLow;
  gMultUp = multUp;
  gMeanMC = 2.289;
  gSigmaMC = 0.010;

  bestBDTbin[0] = 5002;
  bestBDTbin[1] = 4800;
  bestBDTbin[2] = 4900;
  bestBDTbin[3] = 5100;
  bestBDTbin[4] = 5150;
  bestBDTbin[5] = 4570;
  // bestBDTbin[6] = 5000;
  // bestBDTbin[7] = 5082;

  //================================================
  //         Input
  //================================================

  // Read input file
  TFile *inputFile(0);
  if (!gSystem->AccessPathName(fileName.Data() )){
    inputFile   = new TFile(fileName.Data());
  }
  else{
    Printf("ERROR: could not open input root file");
    return NULL;
  }

  TFile *outputFile(0);
  outputFile = new TFile(Form("InvariantMassOptimalBDT_mult_%d_%d.root",multLow,multUp),"RECREATE");
  if(!outputFile){
    Printf("ERROR: outputFile not created");
    exit(1);
  }
  outputFile->cd();

  TList* list;
  TH1D* hBDTBinCut = new TH1D("hBDTBinCut","hBDTBinCut; #it{p}_{T}; BDTBinCut;",nPtBins,ptBins);
  TH1D* hBDTValueCut = new TH1D("hBDTValueCut","hBDTValueCut; #it{p}_{T}; BDTValueCut;",nPtBins,ptBins);

  for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
    // pT Bin Ranges
    Double_t pTBinLow = ptBins[iPtBin];
    Double_t pTBinUp  = ptBins[iPtBin+1];

    TString listName = Form("treeListxml_PfCuts_MB_pT_%.0f_%.0f",pTBinLow,pTBinUp);
    TList* list = (TList*)inputFile->Get(listName.Data());
    if(!list){
      Printf("ERROR: List %s not in inputFile",listName.Data());
      Printf("> Check object listName and compare it with inputfile");
      return NULL;
    }

    TH2D* hListBDTHistoTMVA = (TH2D*)list->FindObject(Form("fBDTHistoTMVA"));
    hListBDTHistoTMVA->SetName(Form("fBDTHistoTMVA_pT_%.0f_%.0f",pTBinLow,pTBinUp));
    if(!hListBDTHistoTMVA){
      Printf("ERROR: hBDTHistoTMVA[%d] not in inputFile",iPtBin);
      return NULL;
    }
    Printf("fBDTHistoTMVA_pT_%.0f_%.0f",pTBinLow,pTBinUp);
  
    Printf("* Multiplicity = [%d - %d]", multLow,multUp);
    Printf("* Input file: %s", inputFile->GetName());
    Printf("* Output file: %s", outputFile->GetName());

    outputFile->cd();

    Double_t bdtValue = 0;

    TH1F* hMass;
    hMass = GetInvariantMassHistogram(hListBDTHistoTMVA,bestBDTbin[iPtBin],bdtValue,1);
    // hMass->SetAxisRange(minMassForPlot,maxMassForPlot);
    hMass->SetName(Form("hMass_pT_%.0f_%.0f",pTBinLow,pTBinUp));
    hMass->SetTitle(Form("%.1f < #it{p}_{T} < %.1f",pTBinLow,pTBinUp));
    hMass->Write();

    TH1F* hMass_1216_1624_merge;
    if(iPtBin == 5){
      hMass_1216_1624_merge = (TH1F*) hMass->Clone("hMass_1216_1624_merge");
    }
    if(iPtBin == 6){
      hMass_1216_1624_merge->Add(hMass);
      hMass_1216_1624_merge->SetName("hMass_pT_12_24");
      hMass_1216_1624_merge->SetTitle("12.0 < #it{p}_{T} < 24.0");
      hMass_1216_1624_merge->Write();
    }

    hBDTBinCut->SetBinContent(iPtBin+1,bestBDTbin[iPtBin]);
    hBDTValueCut->SetBinContent(iPtBin+1,bdtValue);

    delete hListBDTHistoTMVA;
    delete list;
  }
  
  hBDTBinCut->Write();
  hBDTValueCut->Write();

  outputFile->Close();


  return;
}