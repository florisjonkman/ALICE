
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

// Parameters
// --------------------------------------
Int_t    typeB    = kPol2;
Int_t    typeS    = kGaus;
Double_t minFit   = 2.150; // Save = 2.150
Double_t maxFit   = 2.440; // Save = 2.440
TString fileName     = "../../../Results/GRID/Application/AnalysisResults_Lc2pK0S_TMVA_PfCuts_pp2016_2017_2018.root";
TString  tmvaEffFile = "../../Efficiencies/TMVAEfficiencies/TMVAEffOutput_PfCuts_OptimalBDTResponse/TMVAEff_PfCuts_OptimalBDTResponse.root";
// --------------------------------------

// Bins: pT
const Int_t nPtBins = 6;
Double_t ptBins[nPtBins+1]= {1.,2.,4.,6.,8.,12.,24.};; // pT bin ranges
Int_t bestBDTbin[nPtBins];

// NoPCTs - 6 Bins!
Double_t sigmaMC[nPtBins]    = {0.007520,0.007553,0.008798,0.010374,0.012303,0.015439}; // Sigma
Double_t sigmaMCErr[nPtBins] = {0.000048,0.000034,0.000051,0.000093,0.000121,0.000270}; // SigmaErr
Double_t meanMC[nPtBins]     = {2.285820,2.286810,2.287827,2.287972,2.288575,2.288733}; // Mean 
Double_t meanMCErr[nPtBins]  = {0.000066,0.000047,0.000069,0.000124,0.000165,0.000369}; // MeanErr

// CorrectNTPCCuts
// Double_t sigmaMC[nPtBins]    = {0.007520,0.007553,0.008798,0.010374,0.012303,0.016364}; // Sigma
// Double_t sigmaMCErr[nPtBins] = {0.000048,0.000034,0.000052,0.000092,0.000121,0.000244}; // SigmaErr
// Double_t meanMC[nPtBins]     = {2.285820,2.286810,2.287827,2.287972,2.288575,2.288647}; // Mean 
// Double_t meanMCErr[nPtBins]  = {0.000066,0.000047,0.000069,0.000124,0.000165,0.000328}; // MeanErr

// // IncorrectNTPCCuts
// Double_t sigmaMC[nPtBins]    = {0.007501,0.007546,0.008982,0.010531,0.012178,0.016148}; // Sigma
// Double_t sigmaMCErr[nPtBins] = {0.000049,0.000035,0.000052,0.000094,0.000122,0.000246}; // SigmaErr
// Double_t meanMC[nPtBins]     = {2.285841,2.286819,2.287857,2.287810,2.288578,2.288788}; // Mean 
// Double_t meanMCErr[nPtBins]  = {0.000067,0.000048,0.000072,0.000128,0.000167,0.000332}; // MeanErr

// Fit Paramters
Int_t factor4refl = 0;                // 0 = gaus; 1 = gaus+gaus broadened
Int_t nSigmaCount = 3;

Double_t minMassForPlot;
Double_t maxMassForPlot;

// Intialize other objects
Double_t massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
TString gPlotDir;
Int_t gPtLow;
Int_t gPtUp;
Int_t gMultLow;
Int_t gMultUp;
Int_t gBinMin;
Int_t gBinMax;
Double_t gFixSigma;
Double_t gSigmaMC;
Double_t gMeanMC;
Double_t gSigmaMCErr;
Double_t gMeanMCErr;

// Global variables
Double_t gFitMassMean;
Double_t gFitMassMeanErr;
Double_t gFitMassSigma;
Double_t gFitMassSigmaErr;
Double_t gChiSquare;
Double_t gSignificance;
Double_t gSignificanceErr;
Double_t gSignal;
Double_t gSignalErr;
Double_t gBackground;
Double_t gBackgroundErr;
Double_t gBdtValue;
Double_t gRawyield;
Double_t gRawyieldErr;
Double_t gPromptEff;
Double_t gPromptEffErr;

// list of Global variables
Bool_t   lIsFit[10000] = {kFALSE};
Int_t    lBDTBin[10000] = {1};
Double_t lBDTValue[10000] = {0.0};
Double_t lFitMassMean[10000] = {0.0};
Double_t lFitMassMeanErr[10000] = {0.0};
Double_t lFitMassSigma[10000] = {0.0};
Double_t lFitMassSigmaErr[10000] = {0.0};
Double_t lChiSquare[10000] = {0.0};
Double_t lSignificance[10000] = {0.0};
Double_t lSignificanceErr[10000] = {0.0};
Double_t lSignal[10000] = {0.0};
Double_t lSignalErr[10000] = {0.0};
Double_t lBackground[10000] = {0.0};
Double_t lBackgroundErr[10000] = {0.0};
Double_t lRawyield[10000] = {0.0};
Double_t lRawyieldErr[10000] = {0.0};
Double_t lPromptEff[10000] = {0.0};
Double_t lPromptEffErr[10000] = {0.0};

// Global objects
TH2D* gBDTHistoTMVA;
TH2D* hListBDTHistoTMVA[nPtBins];
TH2D* hPtBDTLcPrompt;
TH2D* hPtBDTLcBfd;

// Function to compute the significance
void ComputeSignificance(Double_t signal, Double_t  errsignal, Double_t  background, Double_t  errbackground, Double_t &significance,Double_t &errsignificance){
  Double_t errSigSq=errsignal*errsignal;
  Double_t errBkgSq=errbackground*errbackground;
  Double_t sigPlusBkg=signal+background;
  if(sigPlusBkg>0. && signal>0.){
    significance =  signal/TMath::Sqrt(signal+background);
    errsignificance = significance*TMath::Sqrt((errSigSq+errBkgSq)/(4.*sigPlusBkg*sigPlusBkg)+(background/sigPlusBkg)*errSigSq/signal/signal);
  }else
  {
    significance=0.;
    errsignificance=0.;
  }
  return;
}

// Modified version of AliHFMassFitter
Bool_t MassFitter(TH1F* hist,Bool_t flagDraw = kTRUE, Bool_t debug = kTRUE){

  TH1F* hMass = (TH1F*) hist->Clone();
	Int_t nBins = hMass->GetNbinsX();

  TF1* funBckStore1=0x0;
  TF1* funBckStore2=0x0;
  TF1* funBckStore3=0x0;

	// Define minimum fit range
    Double_t hMin=TMath::Max(minFit,hMass->GetBinLowEdge(1));
  	Double_t hMax=TMath::Min(maxFit,hMass->GetBinLowEdge(hMass->GetNbinsX()+1));

  // Get minimum and maximum bin content within fit range
  Double_t contentMax = 0;
  Double_t contentMin = TMath::Power(10,99);
  for(Int_t iBin = 1; iBin <= nBins; iBin++){
      if(hMass->GetBinCenter(iBin) >= hMin && hMass->GetBinCenter(iBin) <= hMax){
          Double_t content = hMass->GetBinContent(iBin);
          if(content > contentMax) contentMax = content;
          if(content < contentMin) contentMin = content;
      }
  }

  // Draw the histogram the invariant mass
  Double_t plotYmin = contentMin*0.95;
  Double_t plotYmax = contentMax*1.05;

  if(flagDraw){
  	gStyle->SetOptStat(0);
    hMass->GetXaxis()->SetRangeUser(minMassForPlot,maxMassForPlot);
    hMass->GetYaxis()->SetRangeUser(plotYmin,plotYmax);
    hMass->Draw();
  }

  // Creating AliHFMassFitter and setting values
  AliHFMassFitter* fitter=new AliHFMassFitter( hMass,hMin, hMax,1/* 1=no rebin */,typeB,typeS);
  Double_t rebin=nBins/fitter->GetBinN();
  fitter->SetReflectionSigmaFactor(factor4refl);
  fitter->SetInitialGaussianMean(gMeanMC);
  fitter->SetInitialGaussianSigma(gSigmaMC);
  if(gFixSigma > 0) fitter->SetFixGaussianSigma(gFixSigma);
  fitter->SetMinRangeFit(hMin);
  fitter->SetMaxRangeFit(hMax);
  // fitter->SetFitOption("0,R,M,Q");

  // Print set values
  if(debug){
    Printf("");
    Printf("=== Values set for fitting");
    Printf("* Initial mean: %.5f GeV/c^2", gMeanMC);
    Printf("* Min fit range: %.5f GeV/c^2", hMin);
    Printf("* Max fit range: %.5f GeV/c^2", hMax);
    Printf("* Sigma, %s: %.4f", gFixSigma > 0 ? "free" : "fixed", gFixSigma);
    Printf("");
  }

  // Fit the function
  Bool_t out;
  out=fitter->MassFitter(0);
  if(!out){
    std::cout << RED << "FAILED: Full fitfuction did not converge" << RESET << std::endl;
    Printf("  > Trying to fit background function");
    fitter->RefitWithBkgOnly(kTRUE);
    return kFALSE;
  }

  // Extracting variables
  Double_t fitMassMean    = fitter->GetMean();
  Double_t fitMassMeanErr = fitter->GetMeanUncertainty();
  Double_t fitMassSigma   = fitter->GetSigma(); 
  Double_t fitMassSigmaErr= fitter->GetSigmaUncertainty();                   
  Double_t dChiSquare     = fitter->GetReducedChiSquare();

  TF1* fB1               = fitter->GetBackgroundFullRangeFunc();
  TF1* fB2               = fitter->GetBackgroundRecalcFunc();
  TF1* fM                = fitter->GetMassFunc();
  
  if(fB1)  funBckStore1=new TF1(*fB1);
  if(fB2)  funBckStore2=new TF1(*fB2);
  if(fM)   funBckStore3=new TF1(*fM);

  Double_t dSignal          = 0;
  Double_t dSignalErr       = 0;
  Double_t dBackground      = 0;
  Double_t dBackgroundErr   = 0;
  Double_t dSignificance    = 0;
  Double_t dSignificanceErr = 0;
  if(typeB!=kNoBk) fitter->Signal(nSigmaCount,dSignal,dSignalErr);
  if(typeB!=kNoBk) fitter->Background(nSigmaCount,dBackground,dBackgroundErr);
  if(typeB!=kNoBk) fitter->Significance(nSigmaCount,dSignificance,dSignificanceErr);
  Double_t dRawyield     = fitter->GetRawYield();
  Double_t dRawyieldErr  = fitter->GetRawYieldError();

  if(debug){
    Printf("\n> Extract parameters");
    Printf("  * Rawyield signal:\t %.0f ± %.0f Entries", dRawyield,dRawyieldErr);
    Printf("  * Mean mass:\t\t %.3f ± %.3f GeV/c^2", fitMassMean,fitMassMeanErr);
    Printf("  * Mean sigma:\t\t %.3f ± %.3f GeV/c^2", fitMassSigma,fitMassSigmaErr);
    Printf("  * Reduced chisq:\t %.2f ", dChiSquare);
    Printf("  * Signal:\t\t %.0f ± %.0f", dSignal,dSignalErr);
    Printf("  * Background:\t\t %.0f ± %f", dBackground,dBackgroundErr);
    Printf("  * Significance:\t %.2f ± %.2f", dSignificance,dSignificanceErr);
    // Printf("  * integralS = %.1f (sqrt %.1f)",integralS,dSignalErr);
    // Printf("  * errIntS (corr) = %.1f",TMath::Sqrt((errIntS*errIntS)/(hRebinned[iPtBin]->GetBinWidth(4))));
    // Printf("  * integralB = %.1f (sqrt %.1f)",integralB,dBackgroundErr);
    // Printf("  * errIntB (corr) = %.1f",TMath::Sqrt((errIntB*errIntB)/(hRebinned[iPtBin]->GetBinWidth(4))));
    Printf("");
  }
  if(flagDraw){
    	gStyle->SetOptStat(0);
	    // if(minForPlot > fitMassMean-9*fitMassSigma) minForPlot = fitMassMean-9*fitMassSigma;
	    // if(maxForPlot < fitMassMean+9*fitMassSigma) maxForPlot = fitMassMean+9*fitMassSigma;

    	// Draw 3 sigma line
	    TLine* line3sMin = new TLine(fitMassMean-nSigmaCount*fitMassSigma,plotYmin,fitMassMean-nSigmaCount*fitMassSigma,plotYmax);
	    TLine* line3sMax = new TLine(fitMassMean+nSigmaCount*fitMassSigma,plotYmin,fitMassMean+nSigmaCount*fitMassSigma,plotYmax);
	    line3sMin->SetLineColorAlpha(kGray,1.0); line3sMin->SetLineWidth(1); line3sMin->SetLineStyle(2);
	    line3sMax->SetLineColorAlpha(kGray,1.0); line3sMax->SetLineWidth(1); line3sMax->SetLineStyle(2);
	    line3sMin->Draw("Same");
	    line3sMax->Draw("Same");

	    // Plot the Background
      if(typeB!=kNoBk ){
        if(fB1){
          fB1->SetLineColor(kGray+1);
          fB1->SetLineWidth(1);
          fB1->Draw("Same");
        }
        if(fB2){
          fB2->SetLineColor(kRed);
          fB2->SetLineWidth(1);
          fB2->Draw("Same");
        }
      }

      if(fM){
        fM->SetLineColor(kBlue);
        fM->SetLineWidth(1);
        fM->Draw("Same");
      }

	    // Draw Info arround plot
	    if(dChiSquare>0 && dChiSquare<1000){
	      TPaveText *paveChiSquared=new TPaveText(0.721805,0.66899,0.889724,0.834495,"NDC");
	      paveChiSquared->SetBorderSize(0);
	      paveChiSquared->SetFillStyle(0);
	      paveChiSquared->AddText(Form("#chi_{red}^{2} = %.2f",dChiSquare));
	      paveChiSquared->SetTextColor(kGreen+3);
	      paveChiSquared->SetTextAlign(32);
	      paveChiSquared->Draw();

	      TPaveText *paveMeanSigma=new TPaveText(0.594612,0.785714,0.906642,0.900697,"NDC");
	      paveMeanSigma->SetBorderSize(0);
	      paveMeanSigma->SetFillStyle(0);
	      paveMeanSigma->AddText(Form("#mu= %.3f #pm %.3f",fitMassMean,fitMassMeanErr));
	      paveMeanSigma->AddText(Form("#sigma = %.3f #pm %.3f",fitMassSigma,fitMassSigmaErr));
	      paveMeanSigma->SetTextAlign(32);
	      paveMeanSigma->SetTextColor(kBlue);
	      paveMeanSigma->Draw();

	      TPaveText *paveSigBG=new TPaveText(0.110276,0.118467,0.446115,0.283972,"NDC");
	      paveSigBG->SetBorderSize(0);
	      paveSigBG->SetFillStyle(0);
	      paveSigBG->AddText(Form("S (%d#sigma) %.0f #pm %.0f",nSigmaCount,dSignal,dSignalErr));
	      paveSigBG->AddText(Form("B (%d#sigma) %.0f #pm %.0f",nSigmaCount,dBackground,dBackgroundErr));
	      paveSigBG->AddText(Form("S/B (%d#sigma) %.4f",nSigmaCount,dSignal/dBackground));
	      paveSigBG->SetTextAlign(12);
	      paveSigBG->Draw();

	      TPaveText *paveSign=new TPaveText(0.117794,0.792683,0.389724,0.921603,"NDC");
	      paveSign->SetBorderSize(0);
	      paveSign->SetFillStyle(0);
	      paveSign->AddText(Form("Sign (%d#sigma) %.2f #pm %.2f",nSigmaCount,dSignificance,dSignificanceErr));
	      paveSign->SetTextAlign(12);
	      paveSign->SetTextColor(kRed);
	      paveSign->Draw();
	    }
    }

    gFitMassMean     = fitMassMean;
    gFitMassMeanErr  = fitMassMeanErr;
    gFitMassSigma    = fitMassSigma;
    gFitMassSigmaErr = fitMassSigmaErr;
    gChiSquare       = dChiSquare;
    gSignificance    = dSignificance;
    gSignificanceErr = dSignificanceErr;
    gSignal          = dSignal;
    gSignalErr       = dSignalErr;
    gBackground      = dBackground;
    gBackgroundErr   = dBackgroundErr;
    gRawyield        = dRawyield;
    gRawyieldErr     = dRawyieldErr;

	  return kTRUE;

}


TH1F* GetInvariantMassHistogram(TH2D* hBDTHistoTMVA, Int_t binBdtMin,Double_t &bdtValue, Int_t rebin = 1){

	TH1F* hInvariantMass = (TH1F*) hBDTHistoTMVA->ProjectionY(Form("hInvariantMass_pT_%d_%d_mult_%d_%d",gPtLow,gPtUp,gMultLow,gMultUp),binBdtMin,-1);
	TH1F* hRebinned		= (TH1F*) AliVertexingHFUtils::RebinHisto(hInvariantMass,rebin,-1);
	hRebinned->SetAxisRange(minMassForPlot,maxMassForPlot);
	hRebinned->SetName(Form("hRebinned_pT_%d_%d_mult_%d_%d_binbdt_%d",gPtLow,gPtUp,gMultLow,gMultUp,binBdtMin));
	bdtValue = (hBDTHistoTMVA->GetXaxis())->GetBinLowEdge(binBdtMin);
  hRebinned->SetTitle(Form("BDT Response (BDT Bin) #geq %.4f (%d)",bdtValue,binBdtMin));
  hRebinned->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
	hRebinned->GetYaxis()->SetTitle(Form("Entries/(%.1f MeV/c^{2})",(hRebinned->GetBinWidth(1)*1000)));
	hRebinned->SetLineColor(kBlue+1);
  hRebinned->SetMarkerColor(kBlack);
  hRebinned->SetMarkerStyle(8);
  hRebinned->SetMarkerSize(0.5);
  hRebinned->GetYaxis()->SetTitleOffset(1.6);

  return hRebinned;

}

void MakeQAPlots(Int_t nSteps, Int_t groupSize){

  Printf("\n--- MakeQAPlots ---\n");

  // lIsFit[i] = isFit;
  // lBDTBin[i] = binLow;
  // lBDTValue[i] = gBdtValue;
  // lFitMassMean[i] = gFitMassMean;
  // lFitMassMeanErr[i] = gFitMassMeanErr;
  // lFitMassSigma[i] = gFitMassSigma;
  // lFitMassSigmaErr[i] = gFitMassSigmaErr;
  // lChiSquare[i] = gChiSquare;
  // lSignificance[i] = gSignificance; 
  // lSignificanceErr[i] = gSignificanceErr; 
  // lSignal[i] = gSignal;
  // lSignalErr[i] = gSignalErr;
  // lBackground[i] = gBackground;
  // lBackgroundErr[i] = gBackgroundErr;

  Int_t leftOver = nSteps % groupSize;
  if(leftOver > 0) Printf("NOTE: Please choose nSteps MOD groupSize = 0, currently leftOver = %d",leftOver);

  Double_t grFitMassMean      [nSteps/groupSize];
  Double_t grFitMassMeanErr   [nSteps/groupSize];
  Double_t grFitMassSigma     [nSteps/groupSize];
  Double_t grFitMassSigmaErr  [nSteps/groupSize];
  Double_t grChiSquare        [nSteps/groupSize];
  Double_t grSignificance     [nSteps/groupSize];
  Double_t grSignificanceErr  [nSteps/groupSize];
  Double_t grSoverB           [nSteps/groupSize];
  Double_t grBDTValue         [nSteps/groupSize];
  Double_t grRawyield         [nSteps/groupSize];
  Double_t grRawyieldErr      [nSteps/groupSize];
  Double_t grPromptEff        [nSteps/groupSize];
  Double_t grPromptEffErr     [nSteps/groupSize];

  Double_t binCenter          [nSteps/groupSize];
  Double_t binErr             [nSteps/groupSize];

  Int_t nGroup = 0;
  Double_t massMean         = 0.;
  Double_t massMeanErr      = 0.;
  Double_t massSigma        = 0.;
  Double_t massSigmaErr     = 0.;
  Double_t chiSquare        = 0.;
  Double_t significance     = 0.;
  Double_t significanceErr  = 0.;
  Double_t signalOverBackgr = 0.;
  Double_t bdtValue         = 0.;
  Double_t rawyield         = 0.;
  Double_t rawyieldErr      = 0.;
  Double_t promptEff        = 0.;
  Double_t promptEffErr     = 0.;

  Int_t    nContribute      = 0;
  Int_t    binLow           = 0;
  Int_t    binUp            = 0;

  FILE *fp = fopen(Form("%s/FitGroupResults_pT_%d_%d_mult_%d_%d.csv",gPlotDir.Data(),gPtLow,gPtLow,gMultLow,gMultUp), "w");
  fprintf(fp,"nGroup,nContr,binLow,binUp,binCentr,binErr,Mean,MeanErr,Sigma,SigmaErr,ChiSq,Sign,SignErr,Yield,YieldErr,PrEff,PrEffErr\n");
  Printf("%7.7s  %6.6s  %7.7s  %7.7s  %7.7s  %7.7s %6.6s %6.6s  %8.8s  %6.6s  %10.10s  %6.6s  %6.6s  %8.8s  %8.8s  %8.8s  %8.8s  %8.8s","nGroup","nContr","binLow","binUp","binCentr","binErr","BDT","Mean","MeanErr","Sigma","SigmaErr","ChiSq","Sign","SignErr","Yield","YieldErr","PrEff","PrEffErr");

  for (Int_t i = 0; i < nSteps; ++i)
  {
      if(i % groupSize == 0){
        binLow = lBDTBin[i];
      }
      if(i % groupSize == groupSize-1){
        binUp = lBDTBin[i];
      }

      if(lIsFit[i]){
        nContribute++;
        massMean        += lFitMassMean[i];
        massMeanErr     += TMath::Power(lFitMassMeanErr[i],2);
        massSigma       += lFitMassSigma[i];
        massSigmaErr    += TMath::Power(lFitMassSigmaErr[i],2);
        chiSquare       += lChiSquare[i];
        significance    += lSignificance[i];
        significanceErr += TMath::Power(lSignificanceErr[i],2);
        signalOverBackgr+= lBackground[i] > 0 ? lSignal[i]/lBackground[i] : 0;
        bdtValue        += lBDTValue[i];
        rawyield        += lRawyield[i];
        rawyieldErr     += lRawyieldErr[i];
        promptEff       += lPromptEff[i];
        promptEffErr    += lPromptEffErr[i];
      }

      if(i % groupSize == groupSize-1 && nContribute > 0.0 ){

          grFitMassMean      [nGroup] = massMean/Double_t(nContribute);
          grFitMassMeanErr   [nGroup] = TMath::Sqrt(massMeanErr)/Double_t(nContribute);
          grFitMassSigma     [nGroup] = massSigma/Double_t(nContribute);
          grFitMassSigmaErr  [nGroup] = TMath::Sqrt(massSigmaErr)/Double_t(nContribute);
          grChiSquare        [nGroup] = chiSquare/Double_t(nContribute);
          grSignificance     [nGroup] = significance/Double_t(nContribute);
          grSignificanceErr  [nGroup] = TMath::Sqrt(significanceErr)/Double_t(nContribute);
          grSoverB           [nGroup] = signalOverBackgr/Double_t(nContribute);
          grBDTValue         [nGroup] = bdtValue/Double_t(nContribute);
          grRawyield         [nGroup] = rawyield/Double_t(nContribute);
          grRawyieldErr      [nGroup] = rawyieldErr/Double_t(nContribute);
          grPromptEff        [nGroup] = promptEff/Double_t(nContribute);
          grPromptEffErr     [nGroup] = promptEffErr/Double_t(nContribute);

          binCenter          [nGroup] = ( Double_t(binUp)+Double_t(binLow) ) / 2.;
          binErr             [nGroup] = ( Double_t(binUp)-Double_t(binLow) ) / 2.;
      
          Printf("%7.1d  %6.1d  %7.1d  %7.1d  %7.0f  %7.0f  %6.4f  %6.3f  %8.3f  %6.3f  %10.3f  %6.2f  %6.3f  %8.3f  %8.0f  %8.0f  %8.3f  %8.3f",nGroup,nContribute,binLow,binUp,binCenter[nGroup],binErr[nGroup],grBDTValue[nGroup],grFitMassMean[nGroup],grFitMassMeanErr[nGroup],grFitMassSigma[nGroup],grFitMassSigmaErr[nGroup],grChiSquare[nGroup],grSignificance[nGroup],grSignificanceErr[nGroup],grRawyield[nGroup],grRawyieldErr[nGroup],grPromptEff[nGroup],grPromptEffErr[nGroup]);
          fprintf(fp,"%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",nGroup,nContribute,binLow,binUp,binCenter[nGroup],binErr[nGroup],grBDTValue[nGroup],grFitMassMean[nGroup],grFitMassMeanErr[nGroup],grFitMassSigma[nGroup],grFitMassSigmaErr[nGroup],grChiSquare[nGroup],grSignificance[nGroup],grSignificanceErr[nGroup],grRawyield[nGroup],grRawyieldErr[nGroup],grPromptEff[nGroup],grPromptEffErr[nGroup]);

          nGroup++;
          massMean         = 0.;
          massMeanErr      = 0.;
          massSigma        = 0.;
          massSigmaErr     = 0.;
          chiSquare        = 0.;
          significance     = 0.;
          significanceErr  = 0.;
          signalOverBackgr = 0.;
          bdtValue         = 0.;
          rawyield         = 0.;
          rawyieldErr      = 0.;
          promptEff        = 0.;
          promptEffErr     = 0.;
          nContribute      = 0;
      }
  }
  fclose(fp);

  // Canvas of Mean Sigma and Chi-Squared
  TCanvas *cPerformance=new TCanvas("cPerformance","Performance",1800,1800);
  cPerformance->Divide(3,3,0.0000001,0.0000001);
  cPerformance->cd(1);
  gPad->SetGridx(1);
  TGraphErrors* graphMean=new TGraphErrors(nGroup,binCenter,grFitMassMean,binErr,grFitMassMeanErr);
  graphMean->SetName("graphMean");
  graphMean->SetTitle("Mean from Fit;BDT Bin;Mean");
  graphMean->SetMarkerStyle(21);
  graphMean->SetMarkerSize(1);
  graphMean->SetLineColor(kGray+1);
  graphMean->Draw("AP");
  TLine* lineMeanMC  = new TLine(binCenter[0],gMeanMC,binCenter[nGroup-1],gMeanMC);
  TLine* lineMeanMCErr1  = new TLine(binCenter[0],gMeanMC-gMeanMCErr,binCenter[nGroup-1],gMeanMC-gMeanMCErr);
  TLine* lineMeanMCErr2  = new TLine(binCenter[0],gMeanMC+gMeanMCErr,binCenter[nGroup-1],gMeanMC+gMeanMCErr);
  lineMeanMC->SetLineColorAlpha(kRed,1.0); lineMeanMC->SetLineWidth(1); lineMeanMC->SetLineStyle(1);
  lineMeanMCErr1->SetLineColorAlpha(kRed,1.0); lineMeanMCErr1->SetLineWidth(1); lineMeanMCErr1->SetLineStyle(2);
  lineMeanMCErr2->SetLineColorAlpha(kRed,1.0); lineMeanMCErr2->SetLineWidth(1); lineMeanMCErr2->SetLineStyle(2);
  lineMeanMC->Draw(); lineMeanMCErr1->Draw(); lineMeanMCErr2->Draw();

  cPerformance->cd(2);
  gPad->SetGridx(1);
  TGraphErrors* graphSigma=new TGraphErrors(nGroup,binCenter,grFitMassSigma,binErr,grFitMassSigmaErr);
  graphSigma->SetName("graphSigma");
  graphSigma->SetTitle("Sigma from Fit;BDT Bin;Sigma");
  graphSigma->SetMarkerStyle(21);
  graphSigma->SetLineColor(kGray+1);
  graphSigma->SetMarkerSize(1);
  graphSigma->GetYaxis()->SetRangeUser(gSigmaMC-0.006,gSigmaMC+0.006);
  graphSigma->Draw("AP");
  TLine* lineSigmaMC  = new TLine(binCenter[0],gSigmaMC,binCenter[nGroup-1],gSigmaMC);
  TLine* lineSigmaMCErr1  = new TLine(binCenter[0],gSigmaMC-gSigmaMCErr,binCenter[nGroup-1],gSigmaMC-gSigmaMCErr);
  TLine* lineSigmaMCErr2  = new TLine(binCenter[0],gSigmaMC+gSigmaMCErr,binCenter[nGroup-1],gSigmaMC+gSigmaMCErr);
  lineSigmaMC->SetLineColorAlpha(kRed,1.0); lineSigmaMC->SetLineWidth(1); lineSigmaMC->SetLineStyle(1);
  lineSigmaMCErr1->SetLineColorAlpha(kRed,1.0); lineSigmaMCErr1->SetLineWidth(1); lineSigmaMCErr1->SetLineStyle(2);
  lineSigmaMCErr2->SetLineColorAlpha(kRed,1.0); lineSigmaMCErr2->SetLineWidth(1); lineSigmaMCErr2->SetLineStyle(2);
  lineSigmaMC->Draw(); lineSigmaMCErr1->Draw(); lineSigmaMCErr2->Draw();

  cPerformance->cd(3);
  gPad->SetGridx(1);
  TGraphErrors* graphChiSq=new TGraphErrors(nGroup,binCenter,grChiSquare,binErr,0);
  graphChiSq->SetName("graphChiSq");
  graphChiSq->SetTitle("Reduced Chi-squared;BDT Bin;Reduced Chi-squared");
  graphChiSq->SetMarkerStyle(21);
  graphChiSq->SetLineColor(kGray+1);
  graphChiSq->SetMarkerSize(1);
  graphChiSq->GetYaxis()->SetRangeUser(0.2,1.8);
  graphChiSq->Draw("AP");

  // cMeanSigmaChi->SaveAs(Form("%s/MeanSigmaChi.pdf",gPlotDir.Data()));

  // Canvas of Significance and S/B
  // TCanvas *cSignificanceSB=new TCanvas("cSignificanceSB","Significance, S/B",1200,600);
  // cSignificanceSB->Divide(2,1);

  TGraphErrors* graphSignificance=new TGraphErrors(nGroup,binCenter,grSignificance,binErr,grSignificanceErr);
  cPerformance->cd(4);
  gPad->SetGridx(1);
  graphSignificance->SetName("graphSignificance");
  graphSignificance->SetTitle("Significance;BDT Bin;Significance");
  graphSignificance->SetMarkerStyle(21);
  graphSignificance->SetMarkerSize(1);
  graphSignificance->SetLineColor(kGray+1);
  graphSignificance->Draw("AP");

  TGraphErrors* graphSB=new TGraphErrors(nGroup,binCenter,grSoverB,binErr,0);
  cPerformance->cd(5);
  gPad->SetGridx(1);
  graphSB->SetName("graphSoverB");
  graphSB->SetTitle("Signal/Background ;BDT Bin; ");
  graphSB->SetMarkerStyle(21);
  graphSB->SetMarkerSize(1);
  graphSB->SetLineColor(kGray+1);
  graphSB->Draw("AP");

  TGraphErrors* graphYield=new TGraphErrors(nGroup,binCenter,grRawyield,binErr,grRawyieldErr);
  cPerformance->cd(6);
  gPad->SetGridx(1);
  graphYield->SetName("graphYield");
  graphYield->SetTitle("Raw yield ;BDT Bin;");
  graphYield->SetMarkerStyle(21);
  graphYield->SetMarkerSize(1);
  graphYield->SetLineColor(kGray+1);
  graphYield->Draw("AP");

  TGraphErrors* graphEff=new TGraphErrors(nGroup,binCenter,grPromptEff,binErr,grPromptEffErr);
  cPerformance->cd(7);
  gPad->SetGridx(1);
  graphEff->SetName("graphEff");
  graphEff->SetTitle("TMVA Efficiency;BDT Bin;");
  graphEff->SetMarkerStyle(21);
  graphEff->SetMarkerSize(1);
  graphEff->SetLineColor(kGray+1);
  graphEff->Draw("AP");

  TGraphErrors* graphBDT=new TGraphErrors(nGroup,binCenter,grBDTValue,binErr,0);
  cPerformance->cd(8);
  gPad->SetGridx(1);
  gPad->SetGridy(1);
  graphBDT->SetName("graphSoverB");
  graphBDT->SetTitle("BDT Response vs BDT Bin ;BDT Bin;BDT Response");
  graphBDT->GetYaxis()->SetTitleOffset(1.1);
  graphBDT->SetMarkerStyle(21);
  graphBDT->SetMarkerSize(1);
  graphBDT->SetLineColor(kGray+1);
  graphBDT->Draw("AP");

  cPerformance->cd(9);
  gStyle->SetOptStat(0);
  TH1F* hBDT = (TH1F*) gBDTHistoTMVA->ProjectionX(Form("hBDT%d_%d_mult_%d_%d",gPtLow,gPtUp,gMultLow,gMultUp),0,-1);
  hBDT->SetTitle("BDT Response + Ranges");
  hBDT->GetXaxis()->SetTitle("BDT Response");
  hBDT->GetYaxis()->SetTitle("Entries");
  hBDT->GetYaxis()->SetTitleOffset(1.1);
  hBDT->Draw("hist");

  TLine* lineLeft  = new TLine(lBDTValue[0],0,lBDTValue[0],hBDT->GetMaximum());
  TLine* lineRight = new TLine(lBDTValue[nSteps-1],0,lBDTValue[nSteps-1],hBDT->GetMaximum());
  lineLeft->SetLineColorAlpha(kRed,1.0); lineLeft->SetLineWidth(1); lineLeft->SetLineStyle(3);
  lineRight->SetLineColorAlpha(kRed,1.0); lineRight->SetLineWidth(1); lineRight->SetLineStyle(3);
  lineLeft->Draw("Same");
  lineRight->Draw("Same");

  cPerformance->SaveAs(Form("%s/Performance.pdf",gPlotDir.Data()));

}

void HFMassFitterBDT(
        Int_t pTBin,
        Int_t multLow,            // Lower multiplicity value, only for output
        Int_t multUp,             // Upper multiplicity value, only for output
		    Int_t binMin,             // Minimum BDT bin, minimum = 1, default: 4000
		    Int_t binMax,             // Maximum BDT bin, maximum = 10000, default: 5500, if crashes, this value could be to high
		    Int_t nSteps,             // Minimum (binMax-binMin)/nSteps > 1, maximum (binMax-binMin)
	      Int_t reb,                // Rebin factor
        Bool_t fixSigma,          // kTRUE, take sigma from sigmaMC, kFALSE, do not
        Bool_t drawTMVAEff = kFALSE, // drawTMVAEff=kTRUE, tmvaEffFile must be set to the correct path
	      Int_t  groupSize = 10,     // Group size to prevent statistical fluctuations
	      Bool_t flagDrawMassHist = kTRUE, // Draw histograms, if kFALSE fit will not be stored
	      Bool_t debug = kTRUE     // Print all output
				  ){

	std::cout << std::endl << std::flush;
	Printf("========= HFMassFitterBDT.C Started =========");

	gPtLow = ptBins[pTBin];
	gPtUp = ptBins[pTBin+1];
	gMultLow = multLow;
	gMultUp = multUp;
  gFixSigma = fixSigma ? sigmaMC[pTBin] : 0.0 ;
  gMeanMC = meanMC[pTBin];
  gMeanMCErr = meanMCErr[pTBin];
  gSigmaMC = sigmaMC[pTBin];
  gSigmaMCErr = sigmaMCErr[pTBin];

	//================================================
  	//         Input
  	//================================================

  	// Read input file
	TFile *inputFile(0);
	if (!gSystem->AccessPathName(fileName.Data() )){
		inputFile  	= new TFile(fileName.Data());
	}
	else{
		Printf("ERROR: could not open input root file");
		return NULL;
	}

  TString listName = Form("treeListxml_PfCuts_MB_pT_%d_%d",gPtLow,gPtUp);
	TList *list = (TList*)inputFile->Get(listName.Data()); // treeListxml_stdCuts_MB_pT_%d_%d
	if(!list){
		Printf("ERROR: List %s not in inputFile",listName.Data());
		Printf("> Check object listName and compare it with inputfile");
		return NULL;
	}

	TH2D* hBDTHistoTMVA = (TH2D*)list->FindObject(Form("fBDTHistoTMVA"));
	if(!hBDTHistoTMVA){
		Printf("ERROR: hBDTHistoTMVA not in inputFile");
		return NULL;
	}
  gBDTHistoTMVA=hBDTHistoTMVA;

	Int_t binStep = ((binMax-binMin)/nSteps);
	if(binStep < 1){
		Printf("ERROR: binStep = %d, this is not allowed",binStep);
		return NULL;
	}
  gBinMin = binMin;
  gBinMax = binMin+nSteps*binStep;

  // TMVA Eff
  // Read the input bdt cuts file
  TFile *input2(0);
  if(drawTMVAEff){
    if (!gSystem->AccessPathName( tmvaEffFile.Data() )) {
    input2 = TFile::Open( tmvaEffFile.Data() ); // check if file in local directory exists
    }
    else {
      Printf("ABORT: could not open %s\nFill in correct path, or run with drawTMVAEff=kFALSE",tmvaEffFile.Data());
      return NULL;
    }
    Printf("* Input TMVA eff file: %s",input2->GetName());

    hPtBDTLcPrompt = (TH2D*) input2->Get("hPtBDTLcPrompt");
    if(!hPtBDTLcPrompt){
      Printf("ABORT: hPtBDTLcPrompt not found in tmvaEffFile"); 
      return NULL;
    }
  }

  // Output directory
  gPlotDir = Form("Output_pT_%d_%d_reb_%d_minMax_%d_%d_nSteps_%d_grSize_%d",gPtLow,gPtUp,reb,binMin,binMax,nSteps,groupSize);
  gSystem->Exec(Form("rm -r %s",gPlotDir.Data()));
  gSystem->Exec(Form("mkdir %s",gPlotDir.Data()));
  gSystem->Exec(Form("mkdir %s/MassFits",gPlotDir.Data()));
  Printf("* Files saved in: %s/",gPlotDir.Data());

	Printf("* pT = [%d - %d]", gPtLow,gPtUp);
	Printf("* Multiplicity = [%d - %d]", multLow,multUp);
	Printf("* BinRange = [%d - %d], nSteps = %d,binStep = %d", binMin,binMax,nSteps,binStep);
	Printf("* Input file: %s", inputFile->GetName());
	Printf("* TH2D: %s", hBDTHistoTMVA->GetName());
	Printf("  - Entries: %10.0f",hBDTHistoTMVA->GetEntries());
	Printf("  - X (%12s): Nbins = %5d, Low = %7.3f, Up = %7.3f, Width = %6.4f","Bdt",hBDTHistoTMVA->GetNbinsX(),(hBDTHistoTMVA->GetXaxis())->GetBinLowEdge(1),(hBDTHistoTMVA->GetXaxis())->GetBinLowEdge(hBDTHistoTMVA->GetNbinsX()+1),(hBDTHistoTMVA->GetXaxis())->GetBinWidth(1));
	Printf("  - Y (%12s): Nbins = %5d, Low = %7.3f, Up = %7.3f, Width = %6.4f","Inv. Mass Lc",hBDTHistoTMVA->GetNbinsY(),(hBDTHistoTMVA->GetYaxis())->GetBinLowEdge(1),(hBDTHistoTMVA->GetYaxis())->GetBinLowEdge(hBDTHistoTMVA->GetNbinsY()+1),(hBDTHistoTMVA->GetYaxis())->GetBinWidth(1));

	// if(minMassForPlot > minFit) minMassForPlot = minFit;
 //  if(maxMassForPlot < maxFit) maxMassForPlot = maxFit;
  minMassForPlot = minFit;
  maxMassForPlot = maxFit;

 //    Int_t binMin = 3200;
	// Int_t binMax = 5500;
	// Int_t nSteps = 100;

	// TCanvas* cInvariantMass = new TCanvas("cInvariantMass", "InvariantMass", 1200, 1000);
	// cInvariantMass->Divide(10,10);
	const Int_t nCanvases = ceil(Double_t(nSteps) / 25.);
	TCanvas* canvases[nCanvases];

  if(flagDrawMassHist){
    for (Int_t n = 0; n < nCanvases; ++n){
      canvases[n] = new TCanvas(Form("cMassFit%d",n), "InvariantMass + Fit", 1200, 1000);
      canvases[n]->Divide(5,5);
    }
  }

  Int_t binCanvasRanges[nCanvases+1];

	for (Int_t i = 0; i < nSteps; ++i){

		if(debug) Printf("nStep = %d",i);

		Int_t binLow = binMin+binStep*i;

    gFitMassMean     = 0.;
    gFitMassMeanErr  = 0.;
    gFitMassSigma    = 0.;
    gFitMassSigmaErr = 0.;
    gChiSquare       = 0.;
    gSignificance    = 0.;
    gSignificanceErr = 0.;
    gSignal          = 0.;
    gSignalErr       = 0.;
    gBackground      = 0.;
    gBackgroundErr   = 0.;
    gBdtValue        = 0.;
    gRawyield        = 0.;
    gRawyieldErr     = 0.;
    gPromptEff       = 0.;
    gPromptEffErr    = 0.;

		TH1F* hRebinned = GetInvariantMassHistogram(hBDTHistoTMVA,binLow,gBdtValue,reb);

		if(flagDrawMassHist) canvases[i/ 25]->cd((i % 25)+1);

		Bool_t isFit = MassFitter(hRebinned,flagDrawMassHist,debug);

    if((i % 25) == 0)   binCanvasRanges[i/25] = binLow;
    if(i == nSteps-1  ) binCanvasRanges[nCanvases] = binLow;

    if(drawTMVAEff){
      Double_t nAccErr = 0;
      Double_t nAcc = hPtBDTLcPrompt->IntegralAndError(pTBin+1,pTBin+1,binLow,10000,nAccErr);
      Double_t nRejErr = 0;
      Double_t nRej = hPtBDTLcPrompt->IntegralAndError(pTBin+1,pTBin+1,1,binLow-1,nRejErr);
      gPromptEff = nAcc/(nAcc+nRej);
      gPromptEffErr = gPromptEff*TMath::Sqrt((nAcc > 0. ? TMath::Power(nAccErr/nAcc,2) : 0. ) + (nRej > 0 ? TMath::Power(nRejErr/nRej,2) : 0.) );
      Printf("nAcc = %.0f , nRej = %.0f",nAcc,nRej);
      Printf("Eff = %.3f ± %.3f",gPromptEff,gPromptEffErr);
    }
    
    lIsFit[i]           = isFit;
    lBDTBin[i]          = binLow;
    lBDTValue[i]        = gBdtValue;
    lFitMassMean[i]     = gFitMassMean;
    lFitMassMeanErr[i]  = gFitMassMeanErr;
    lFitMassSigma[i]    = gFitMassSigma;
    lFitMassSigmaErr[i] = gFitMassSigmaErr;
    lChiSquare[i]       = gChiSquare;
    lSignificance[i]    = gSignificance; 
    lSignificanceErr[i] = gSignificanceErr; 
    lSignal[i]          = gSignal;
    lSignalErr[i]       = gSignalErr;
    lBackground[i]      = gBackground;
    lBackgroundErr[i]   = gBackgroundErr;
    lRawyield[i]        = gRawyield;
    lRawyieldErr[i]     = gRawyieldErr;
    lPromptEff[i]       = gPromptEff;
    lPromptEffErr[i]    = gPromptEffErr;
	}

  
  MakeQAPlots(nSteps,groupSize);

	// Saving fitparameters in file
	FILE *fp = fopen(Form("%s/FitResults_pT_%d_%d_mult_%d_%d_minMax_%d_%d.csv",gPlotDir.Data(),gPtLow,gPtUp,multLow,multUp,binMin,binMax), "w");
	fprintf(fp,"IsFit,Bin,BDT,Mean,MeanErr,Sigma,SigmaErr,ChiSq,Sign,SignErr,Signal,SignalErr,Backg,BackgErr,Yield,YieldErr,PromptEff,PromptEffErr\n");
	for(Int_t i = 0; i < nSteps; i++){
	fprintf(fp,"%d,%d,%.4f,%.3f,%.3f,%.3f,%.3f,%.2f,%.2f,%.2f,%.0f,%.0f,%0.f,%.0f,%.0f,%.0f,%.3f,%.3f\n",
		lIsFit[i],lBDTBin[i],lBDTValue[i],lFitMassMean[i],lFitMassMeanErr[i],lFitMassSigma[i],lFitMassSigmaErr[i],lChiSquare[i],lSignificance[i],lSignificanceErr[i],lSignal[i],lSignalErr[i],lBackground[i],lBackgroundErr[i],lRawyield[i],lRawyieldErr[i],lPromptEff[i],lPromptEffErr[i]);
	}
	fclose(fp);

  if(flagDrawMassHist){
    for (Int_t n = 0; n < nCanvases; ++n){
      canvases[n]->SaveAs(Form("%s/MassFits/MassFit%d_pT_%d_%d_mult_%d_%d_minMax_%d_%d.pdf",gPlotDir.Data(),n,gPtLow,gPtUp,multLow,multUp,binCanvasRanges[n],binCanvasRanges[n+1]));
      canvases[n]->Close();
    }
  }

	Printf("\nDone! Parameters:");
	Printf("* pT = [%d - %d]", gPtLow,gPtUp);
	Printf("* Multiplicity = [%d - %d]", multLow,multUp);
	Printf("* Bin = [%d - %d], step = %d",binMin,binMax,binStep);
	Printf("* Input file: %s", inputFile->GetName());
	Printf("* list: %s",list->GetName());
	Printf("* TH2D: %s", hBDTHistoTMVA->GetName());
	Printf("* Rebin = %d",reb);
	Printf("* GroupSize = %d",groupSize);
	Printf("* MeanMC = %.3f",gMeanMC);
	Printf("* SigmaMC = %.3f",gSigmaMC);
	Printf("* SigmaFix = %s",gFixSigma > 0 ? "True" : "False");
	Printf("* MinFit = %.3f",minFit);
	Printf("* MaxFit = %.3f",maxFit);
	Printf("* TypeB = %d",typeB);
	Printf("* TypeS = %d",typeS);
	Printf("\n");
	
	return;
}

void FitSingleBin(Int_t pTBin,
                  Int_t multLow,            // Lower multiplicity value, only for output
                  Int_t multUp,             // Upper multiplicity value, only for output
                  Int_t bin,                // BDT bin, min 1, maximum 10000
                  Int_t reb,                // Rebin factor
                  Bool_t fixSigma          // kTRUE, take sigma from sigmaMC, kFALSE, do not
                  ){

  std::cout << std::endl << std::flush;
  Printf("========= FitSingleBin.C Started =========");

  gPtLow = ptBins[pTBin];
  gPtUp = ptBins[pTBin+1];
  gMultLow = multLow;
  gMultUp = multUp;
  gFixSigma = fixSigma ? sigmaMC[pTBin] : 0.0 ;
  gMeanMC = meanMC[pTBin];
  gSigmaMC = sigmaMC[pTBin];

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

  TString listName = Form("treeListxml_PfCuts_MB_pT_%d_%d",gPtLow,gPtUp);
  TList *list = (TList*)inputFile->Get(listName.Data());
  if(!list){
    Printf("ERROR: List %s not in inputFile",listName.Data());
    Printf("> Check object listName and compare it with inputfile");
    return NULL;
  }

  TH2D* hBDTHistoTMVA = (TH2D*)list->FindObject(Form("fBDTHistoTMVA"));
  if(!hBDTHistoTMVA){
    Printf("ERROR: hBDTHistoTMVA not in inputFile");
    return NULL;
  }

  Printf("* pT = [%d - %d]", gPtLow,gPtUp);
  Printf("* Multiplicity = [%d - %d]", multLow,multUp);
  Printf("* Bin = %d",bin);
  Printf("* Input file: %s", inputFile->GetName());
  Printf("* TH2D: %s", hBDTHistoTMVA->GetName());
  Printf("  - Entries: %10.0f",hBDTHistoTMVA->GetEntries());
  Printf("  - X (%12s): Nbins = %5d, Low = %7.3f, Up = %7.3f, Width = %6.4f","Bdt",hBDTHistoTMVA->GetNbinsX(),(hBDTHistoTMVA->GetXaxis())->GetBinLowEdge(1),(hBDTHistoTMVA->GetXaxis())->GetBinLowEdge(hBDTHistoTMVA->GetNbinsX()+1),(hBDTHistoTMVA->GetXaxis())->GetBinWidth(1));
  Printf("  - Y (%12s): Nbins = %5d, Low = %7.3f, Up = %7.3f, Width = %6.4f","Inv. Mass Lc",hBDTHistoTMVA->GetNbinsY(),(hBDTHistoTMVA->GetYaxis())->GetBinLowEdge(1),(hBDTHistoTMVA->GetYaxis())->GetBinLowEdge(hBDTHistoTMVA->GetNbinsY()+1),(hBDTHistoTMVA->GetYaxis())->GetBinWidth(1));

  if(minMassForPlot > minFit) minMassForPlot = minFit;
  if(maxMassForPlot < maxFit) maxMassForPlot = maxFit;
  minMassForPlot = minFit;
  maxMassForPlot = maxFit;

  TCanvas* cInvariantMass = new TCanvas("cMassFit", "InvariantMass", 800 , 600);

  gFitMassMean = 0.;
  gFitMassMeanErr = 0.;
  gFitMassSigma = 0.;
  gFitMassSigmaErr = 0.;
  gChiSquare = 0.;
  gSignificance = 0.;
  gSignificanceErr = 0.;
  gSignal = 0.;
  gSignalErr = 0.;
  gBackground = 0.;
  gBackgroundErr = 0.;
  gBdtValue = 0.;

  TH1F* hRebinned = GetInvariantMassHistogram(hBDTHistoTMVA,bin,gBdtValue,reb);

  Bool_t isFit = MassFitter(hRebinned,kTRUE,kTRUE);

  cInvariantMass->SaveAs(Form("MassFit_pT_%d_%d_reb_%d_mult_%d_%d_bin_%d.pdf",gPtLow,gPtUp,reb,multLow,multUp,bin));
  // cInvariantMass->Close();

  FILE *fp = fopen(Form("Parameters_pT_%d_%d_reb_%d_mult_%d_%d_bin_%d.txt",gPtLow,gPtUp,reb,multLow,multUp,bin), "w");

  fprintf(fp,"* listName = %s",listName.Data());
	fprintf(fp,"* pT = [%d - %d]\n", gPtLow,gPtUp);
	fprintf(fp,"* Multiplicity = [%d - %d]\n", multLow,multUp);
	fprintf(fp,"* Bin = %d\n",bin);
	fprintf(fp,"* Input file: %s\n", inputFile->GetName());
	fprintf(fp,"* list: %s\n",list->GetName());
	fprintf(fp,"* TH2D: %s\n", hBDTHistoTMVA->GetName());
	fprintf(fp,"* Rebin: %d\n",reb);
	fprintf(fp,"* MeanMC: %.3f\n",gMeanMC);
	fprintf(fp,"* SigmaMC: %.3f\n",gSigmaMC);
	fprintf(fp,"* SigmaFix: %.3f\n",gFixSigma);
	fprintf(fp,"* MinFit: %.3f\n",minFit);
	fprintf(fp,"* MaxFit: %.3f\n",maxFit);
	fprintf(fp,"* TypeB = %d\n",typeB);
	fprintf(fp,"* TypeS = %d\n",typeS);

	fclose(fp);

	Printf("\nDone! Parameters:");
  Printf("* listName = %s",listName.Data());
	Printf("* pT = [%d - %d]", gPtLow,gPtUp);
	Printf("* Multiplicity = [%d - %d]", multLow,multUp);
	Printf("* Bin = %d",bin);
	Printf("* Input file: %s", inputFile->GetName());
	Printf("* list: %s",list->GetName());
	Printf("* TH2D: %s", hBDTHistoTMVA->GetName());
	Printf("* Rebin = %d",reb);
	Printf("* MeanMC = %.3f",gMeanMC);
	Printf("* SigmaMC = %.3f",gSigmaMC);
	Printf("* SigmaFix = %s",gFixSigma > 0 ? "True" : "False");
	Printf("* MinFit = %.3f",minFit);
	Printf("* MaxFit = %.3f",maxFit);
	Printf("* TypeB = %d",typeB);
	Printf("* TypeS = %d",typeS);
	Printf("\n");

  return;
}