#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TChain.h"
#include "TList.h"
#include "TFile.h"
#include "TH1F.h"
#include "TROOT.h"
#include "TBrowser.h"
#include "TMath.h"
#include <TDatabasePDG.h>

#include "AliVertexingHFUtils.h"
#include "AliHFMassFitter.h"

#include <fstream>
#include <iostream>
#include <cstdio>
#endif

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */

// Debugging
Bool_t debug = kTRUE;
Int_t debugNentries = 0;

// Keep canvasses open
Bool_t KeepOpen = kFALSE;

enum {K0Short, Lambda, AntiLambda,LambdaC};
enum {kBoth, kParticleOnly, kAntiParticleOnly};
enum {kExpo=0, kLinear, kPol2, kNoBk,kPowEx,kThrExpo};
enum {kGaus=0, kDoubleGaus};

// Bins: pT
const Int_t nPtBins = 6;
Double_t ptBins[nPtBins+1] =  {  1,  2,  4,  6,  8, 12, 24}; // pT bin ranges

TString gOutputDirname;
TDirectory* directories[nPtBins];

Int_t    isFit[nPtBins];
Int_t    rebin[nPtBins];
Int_t    typeB[nPtBins];
Int_t    typeS[nPtBins];
Double_t minFit[nPtBins];
Double_t maxFit[nPtBins];
Double_t minPlot[nPtBins];
Double_t maxPlot[nPtBins];
Double_t fixSigma[nPtBins];

TH1D* hMass[nPtBins];
TH1F* hRebinned[nPtBins];

void SetFitParameters(){

  // Fit parameters per Pt Bin
  // Fit?             Rebin             BackgroudFunc        SignalFunct           MinimumFitRange       MaxFitRange         MinimumPlotRange       MaxPlotRange
  isFit[0] = kTRUE;   rebin[0]  = 3;    typeB[0]  = kNoBk;   typeS[0]  = kGaus;    minFit[0]  = 2.260;   maxFit[0]  = 2.310; minPlot[0]  = 2.230;   maxPlot[0]  = 2.345;  fixSigma[0] = 0.000;
  isFit[1] = kTRUE;   rebin[1]  = 2;    typeB[1]  = kNoBk;   typeS[1]  = kGaus;    minFit[1]  = 2.260;   maxFit[1]  = 2.310; minPlot[1]  = 2.230;   maxPlot[1]  = 2.345;  fixSigma[1] = 0.000;
  isFit[2] = kTRUE;   rebin[2]  = 3;    typeB[2]  = kNoBk;   typeS[2]  = kGaus;    minFit[2]  = 2.262;   maxFit[2]  = 2.315; minPlot[2]  = 2.230;   maxPlot[2]  = 2.345;  fixSigma[2] = 0.000;
  isFit[3] = kTRUE;   rebin[3]  = 4;    typeB[3]  = kNoBk;   typeS[3]  = kGaus;    minFit[3]  = 2.258;   maxFit[3]  = 2.320; minPlot[3]  = 2.230;   maxPlot[3]  = 2.345;  fixSigma[3] = 0.000;
  isFit[4] = kTRUE;   rebin[4]  = 4;    typeB[4]  = kNoBk;   typeS[4]  = kGaus;    minFit[4]  = 2.250;   maxFit[4]  = 2.330; minPlot[4]  = 2.190;   maxPlot[4]  = 2.382;  fixSigma[4] = 0.000;
  isFit[5] = kTRUE;   rebin[5]  = 4;    typeB[5]  = kNoBk;   typeS[5]  = kGaus;    minFit[5]  = 2.240;   maxFit[5]  = 2.340; minPlot[5]  = 2.190;   maxPlot[5]  = 2.382;  fixSigma[5] = 0.000;
  // isFit[6] = kTRUE;   rebin[6]  = 6;    typeB[6]  = kNoBk;   typeS[6]  = kGaus;    minFit[6]  = 2.240;   maxFit[6]  = 2.340; minPlot[6]  = 2.190;   maxPlot[6]  = 2.382;  fixSigma[6] = 0.000;
  // isFit[6] = kTRUE;   rebin[6]  = 2;   typeB[6]  = kNoBk;   typeS[6]  = kGaus;    minFit[6]  = 2.150;   maxFit[6]  = 2.440;
  // isFit[7] = kTRUE;    rebin[7]  = 2;   typeB[7]  = kNoBk;   typeS[7]  = kGaus;    minFit[7]  = 2.150;   maxFit[7]  = 2.440;
  
  // Saving fitparameters in file
  FILE *fp = fopen(Form("%s/SetParameters.csv",gOutputDirname.Data()), "w");
  fprintf(fp, "pTBin, pTLow, pTUp, isFit, Rebin, BgFunc, SigFunc, MinRange, MaxRange\n");
  for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
    // pT Bin Ranges
    Double_t pTBinLow = ptBins[iPtBin];
    Double_t pTBinUp  = ptBins[iPtBin+1];
    fprintf(fp,"%d, %.1f, %.1f, %d, %d, %d, %d, %.3f, %.3f \n",iPtBin+1,pTBinLow,pTBinUp,isFit[iPtBin],rebin[iPtBin],typeB[iPtBin],typeS[iPtBin],minFit[iPtBin],maxFit[iPtBin]);
  }
  fclose(fp);

}

void MassFitter(TFile* outputFile, Bool_t flagDraw = kTRUE
				         ){

	//===========================================================
	//     Fit Parameters per PT
	//===========================================================

	Int_t factor4refl = 0;				        // 0 = gaus; 1 = gaus+gaus broadened
 	Int_t nSigmaCount = 3;

  //===========================================================
  //     Creating histograms
  //===========================================================  


  // // Original Invariant Mass Histograms
  // TH1D** hMass = new TH1D*[nPtBins];
  // for(Int_t i=0;i<nPtBins;i++) hMass[i]=0x0;

  TH1D* hCntSig1            =new TH1D("hCntSig1","hCntSig1",nPtBins,ptBins);
  TH1D* hCntSig2            =new TH1D("hCntSig2","hCntSig2",nPtBins,ptBins);
  TH1D* hNDiffCntSig1       =new TH1D("hNDiffCntSig1","hNDiffCntSig1",nPtBins,ptBins);
  TH1D* hNDiffCntSig2       =new TH1D("hNDiffCntSig2","hNDiffCntSig2",nPtBins,ptBins);
  TH1D* hSignal             =new TH1D("hSignal","hSignal",nPtBins,ptBins);
  TH1D* hRawyield           =new TH1D("hRawyield","hRawyield",nPtBins,ptBins);
  TH1D* hRelErrSig          =new TH1D("hRelErrSig","hRelErrSig",nPtBins,ptBins);
  TH1D* hInvSignif          =new TH1D("hInvSignif","hInvSignif",nPtBins,ptBins);
  TH1D* hBackground         =new TH1D("hBackground","hBackground",nPtBins,ptBins);
  TH1D* hSignalBackground   =new TH1D("hSignalBackground","hSignalBackground",nPtBins,ptBins);  
  TH1D* hBackgroundNormSigma=new TH1D("hBackgroundNormSigma","hBackgroundNormSigma",nPtBins,ptBins);
  TH1D* hSignificance       =new TH1D("hSignificance","Signifi",nPtBins,ptBins);
  TH1D* hMassMean           =new TH1D("hMassMean","hMassMean",nPtBins,ptBins);
  TH1D* hSigma              =new TH1D("hSigma","hSigma",nPtBins,ptBins);

  //===========================================================
  //     Prepare before looping over pT
  //=========================================================== 

  Double_t gaussSet      = 2.289;
  Double_t sigmaSet      = 0.010;
  Double_t minMassForFit = 2.150;

  // Initialize variables for fitting
  AliHFMassFitter** fitter=new AliHFMassFitter*[nPtBins];

  //===========================================================
  //     Create canvasses
  //=========================================================== 

  TCanvas* cInvariantMass; 
  if(flagDraw){
    cInvariantMass = new TCanvas("cInvariantMass", "InvariantMass", 1200, 550);
    cInvariantMass->Divide(3,2);
  }
 
	Double_t dSignificance,dSignificanceErr,dSignal,dSignalErr,dBackground,dBackgroundErr;
  TF1* funBckStore1=0x0;
  TF1* funBckStore2=0x0;
  TF1* funBckStore3=0x0;

   Double_t dChiSquare[nPtBins] = {0.0};
  Double_t ptBinsCenter[nPtBins] = {0.0};
  Double_t ptBinsError[nPtBins] = {0.0};

  Double_t fitMassMeanTotal =0;
  Double_t fitMassMeanErrTotal =0;
  Double_t fitMassSigmaTotal =0;
  Double_t fitMassSigmaErrTotal =0;
  Double_t nFit=0;

	for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
	
    if(flagDraw) cInvariantMass->cd(iPtBin+1);

    // pT Bin Ranges
    Double_t pTBinLow = ptBins[iPtBin];
    Double_t pTBinUp  = ptBins[iPtBin+1];
    Printf("\n%s---------------------------------\n",MAGENTA);
    Printf("=== %.1f < pT < %.1f ===",pTBinLow,pTBinUp);
    Printf("\n---------------------------------%s\n",RESET);

    // pT Center and Error
    ptBinsCenter[iPtBin] = (pTBinLow+pTBinUp)/2.;
    ptBinsError[iPtBin]  = TMath::Abs((pTBinUp-pTBinLow)/2.);

    Double_t minMassForPlot = minPlot[iPtBin];
    Double_t maxMassForPlot = maxPlot[iPtBin];

    if(minMassForPlot > minFit[iPtBin]){
      minMassForPlot = minFit[iPtBin];
      Printf("NOTE: Changed minMassForPlot from %.3f to %.3f",minMassForPlot,minPlot[iPtBin]);
    }
    if(maxMassForPlot < maxFit[iPtBin]){
      minMassForPlot = minFit[iPtBin];
      Printf("NOTE: Changed minMassForPlot from %.3f to %.3f",minMassForPlot,maxPlot[iPtBin]);
    }

    // Rebin histogram
		Int_t orignalNBins	= hMass[iPtBin]->GetNbinsX();
		hRebinned[iPtBin]		= (TH1F*) AliVertexingHFUtils::RebinHisto(hMass[iPtBin],rebin[iPtBin],-1);
		hRebinned[iPtBin] ->SetAxisRange(minMassForPlot,maxMassForPlot);
    hRebinned[iPtBin] ->SetName(Form("hRebinned_pT_%.0f_%.0f",pTBinLow,pTBinLow));
    hRebinned[iPtBin] ->SetTitle(Form("%.1f < #it{p}_{T} < %.1f",pTBinLow,pTBinUp));
    hRebinned[iPtBin] ->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    hRebinned[iPtBin] ->GetYaxis()->SetTitle(Form("Entries/(%.1f MeV/c^{2})",(hRebinned[iPtBin] ->GetBinWidth(1)*1000)));
    hRebinned[iPtBin] ->SetLineColor(kBlue+1);
    hRebinned[iPtBin] ->SetMarkerColor(kBlack);
    hRebinned[iPtBin] ->SetMarkerStyle(8);
    hRebinned[iPtBin] ->SetMarkerSize(0.5);
    hRebinned[iPtBin] ->GetYaxis()->SetTitleOffset(1.6);
    Int_t nBinsRebinned = hRebinned[iPtBin]->GetNbinsX();

    // Define minimum fit range
  	Double_t hMin=TMath::Max(minFit[iPtBin],hRebinned[iPtBin]->GetBinLowEdge(2));
  	Double_t hMax=TMath::Min(maxFit[iPtBin],hRebinned[iPtBin]->GetBinLowEdge(hRebinned[iPtBin]->GetNbinsX()));

    // Get minimum and maximum bin content within fit range
    Double_t contentMax = 0;
    Double_t contentMin = TMath::Power(10,99);
    for(Int_t iBin = 1; iBin <= nBinsRebinned; iBin++){
        if( hRebinned[iPtBin]->GetBinCenter(iBin) >= hMin &&  hRebinned[iPtBin]->GetBinCenter(iBin) <= hMax){
            Double_t content = hRebinned[iPtBin]->GetBinContent(iBin);
            if(content > contentMax) contentMax = content;
            if(content < contentMin) contentMin = content;
        }
    }

    // Draw the histogram the invariant mass
    Double_t plotYmin = contentMin*0.95;
    Double_t plotYmax = contentMax*1.05;
    
    if(flagDraw){
      gStyle->SetOptStat(0);
      hRebinned[iPtBin]->GetXaxis()->SetRangeUser(minMassForPlot,maxMassForPlot);
      hRebinned[iPtBin]->GetYaxis()->SetRangeUser(0,plotYmax);
      hRebinned[iPtBin]->Draw();
    }  

    // If not fitted, only plot the histogram
    if(!isFit[iPtBin]){
      continue;
    }      

    // Creating AliHFMassFitter and setting values
  	fitter[iPtBin]=new AliHFMassFitter( hRebinned[iPtBin],hMin, hMax,1/* 1=no rebin */,typeB[iPtBin],typeS[iPtBin]);
  	rebin[iPtBin]=orignalNBins/fitter[iPtBin]->GetBinN();
    fitter[iPtBin]->SetReflectionSigmaFactor(factor4refl);
  	fitter[iPtBin]->SetInitialGaussianMean(gaussSet);
    fitter[iPtBin]->SetInitialGaussianSigma(sigmaSet);
    if(fixSigma[iPtBin] > 0) fitter[iPtBin]->SetFixGaussianSigma(fixSigma[iPtBin]);
    fitter[iPtBin]->SetMinRangeFit(hMin);
    fitter[iPtBin]->SetMaxRangeFit(hMax);

    // Print set values
    Printf("=== Values set for fitting, %.0f < pT < %.0f",pTBinLow,pTBinUp);
    Printf("* Initial mean: %.5f GeV/c^2", gaussSet);
    Printf("* Min fit range: %.5f GeV/c^2", hMin);
    Printf("* Max fit range: %.5f GeV/c^2", hMax);

    // Fit the function
    Bool_t out;
    out=fitter[iPtBin]->MassFitter(0);
    if(!out){
      std::cout << RED << "FAILED: Full fitfuction did not converge" << RESET << std::endl;
  		Printf("  > Trying to fit background function");
      // fitter[iPtBin]->RefitWithBkgOnly(kTRUE);
      continue;
  	}
    nFit++;

    // Extracting variables
  	Double_t fitMassMean   = fitter[iPtBin]->GetMean();
    Double_t fitMassMeanErr= fitter[iPtBin]->GetMeanUncertainty();
  	Double_t fitMassSigma  = fitter[iPtBin]->GetSigma(); 
    Double_t fitMassSigmaErr= fitter[iPtBin]->GetSigmaUncertainty();                   
  	dChiSquare[iPtBin]     = fitter[iPtBin]->GetReducedChiSquare();

    TF1* fB1               = fitter[iPtBin]->GetBackgroundFullRangeFunc();
  	TF1* fB2               = fitter[iPtBin]->GetBackgroundRecalcFunc();
  	TF1* fM                = fitter[iPtBin]->GetMassFunc();
    
    if(iPtBin==0 && fB1)  funBckStore1=new TF1(*fB1);
    if(iPtBin==0 && fB2)  funBckStore2=new TF1(*fB2);
    if(iPtBin==0 && fM)   funBckStore3=new TF1(*fM);

    fitMassMeanTotal     += fitMassMean;
    fitMassMeanErrTotal  += fitMassMeanErr*fitMassMeanErr;
    fitMassSigmaTotal    += fitMassSigma;
    fitMassSigmaErrTotal += fitMassSigmaErr*fitMassSigmaErr;

    if(typeB[iPtBin]!=kNoBk) fitter[iPtBin]->Signal(nSigmaCount,dSignal,dSignalErr);
    if(typeB[iPtBin]!=kNoBk) fitter[iPtBin]->Background(nSigmaCount,dBackground,dBackgroundErr);
    if(typeB[iPtBin]!=kNoBk) fitter[iPtBin]->Significance(nSigmaCount,dSignificance,dSignificanceErr);
    Double_t dRawyield     = fitter[iPtBin]->GetRawYield();
    Double_t dRawyieldErr  = fitter[iPtBin]->GetRawYieldError();

    Float_t minBinSum = hMass[iPtBin]->FindBin(fitMassMean-nSigmaCount*fitMassSigma);
    Float_t maxBinSum = hMass[iPtBin]->FindBin(fitMassMean+nSigmaCount*fitMassSigma);

    Float_t cntSig1=0.;
    Float_t cntSig2=0.;
    Float_t cntErr=0.;
    for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
      Float_t bkg1=fB1 ? fB1->Eval(hMass[iPtBin]->GetBinCenter(iMB))/rebin[iPtBin] : 0;
      Float_t bkg2=fB2 ? fB2->Eval(hMass[iPtBin]->GetBinCenter(iMB))/rebin[iPtBin] : 0;
      cntSig1+=(hMass[iPtBin]->GetBinContent(iMB)-bkg1);
      cntSig2+=(hMass[iPtBin]->GetBinContent(iMB)-bkg2);
      cntErr+=(hMass[iPtBin]->GetBinContent(iMB));
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
        if(typeB[iPtBin]!=kNoBk){
          fB1->SetLineColor(kGray+1);
          fB1->SetLineWidth(1);
          fB1->Draw("Same");

          fB2->SetLineColor(kRed);
          fB2->SetLineWidth(1);
          fB2->Draw("Same");
        }

        fM->SetLineColor(kBlue);
        fM->SetLineWidth(1);
        fM->Draw("Same");

        // Draw Info arround plot
        if(dChiSquare[iPtBin]>0 && dChiSquare[iPtBin]<1000){
          TPaveText *paveChiSquared=new TPaveText(0.721805,0.66899,0.889724,0.834495,"NDC");
          paveChiSquared->SetBorderSize(0);
          paveChiSquared->SetFillStyle(0);
          paveChiSquared->AddText(Form("#chi_{red}^{2} = %.2f",dChiSquare[iPtBin]));
          paveChiSquared->SetTextColor(kGreen+3);
          paveChiSquared->SetTextAlign(32);
          paveChiSquared->Draw();

          TPaveText *paveMeanSigma=new TPaveText(0.594612,0.785714,0.906642,0.900697,"NDC");
          paveMeanSigma->SetBorderSize(0);
          paveMeanSigma->SetFillStyle(0);
          paveMeanSigma->AddText(Form("#mu = %.3f #pm %.3f",fitMassMean,fitMassMeanErr));
          paveMeanSigma->AddText(Form("#sigma = %.3f #pm %.3f",fitMassSigma,fitMassSigmaErr));
          paveMeanSigma->SetTextAlign(32);
          paveMeanSigma->SetTextColor(kBlue);
          paveMeanSigma->Draw();

          TPaveText *paveSigBG=new TPaveText(0.110276,0.118467,0.446115,0.283972,"NDC");
          paveSigBG->SetBorderSize(0);
          paveSigBG->SetFillStyle(0);
          paveSigBG->AddText(Form("S (%d#sigma) %.0f #pm %.0f",nSigmaCount,dSignal,dSignalErr));
          paveSigBG->AddText(Form("B (%d#sigma) %.0f #pm %.0f",nSigmaCount,dBackground,dBackgroundErr));
          if(typeB[iPtBin]!=kNoBk) paveSigBG->AddText(Form("S/B (%d#sigma) %.4f",nSigmaCount,dSignal/dBackground));
          else paveSigBG->AddText(Form("S/B (%d#sigma) N/A",nSigmaCount));
          paveSigBG->SetTextAlign(12);
          paveSigBG->Draw();

          TPaveText *paveSign=new TPaveText(0.109278,0.774332,0.418304,0.933365,"NDC");
          paveSign->SetBorderSize(0);
          paveSign->SetFillStyle(0);
          paveSign->AddText(Form("Sign (%d#sigma) %.2f #pm %.2f",nSigmaCount,dSignificance,dSignificanceErr));
          paveSign->SetTextAlign(12);
          paveSign->SetTextColor(kRed);
          paveSign->Draw();

        }
      }

    // Fill histograms
    hCntSig1->SetBinContent(iPtBin+1,cntSig1);
    hCntSig1->SetBinError(iPtBin+1,TMath::Sqrt(cntErr));
    hNDiffCntSig1->SetBinContent(iPtBin+1,dSignal > 0 ? (dSignal-cntSig1)/dSignal : 0);
    hNDiffCntSig1->SetBinError(iPtBin+1,dSignal > 0 ? TMath::Sqrt(cntErr)/dSignal : 0);
    hCntSig2->SetBinContent(iPtBin+1,cntSig2);
    hNDiffCntSig2->SetBinContent(iPtBin+1,dSignal > 0 ? (dSignal-cntSig2)/dSignal : 0);
    hNDiffCntSig2->SetBinError(iPtBin+1,dSignal > 0 ? TMath::Sqrt(cntErr)/dSignal : 0);
    hCntSig2->SetBinError(iPtBin+1,TMath::Sqrt(cntErr));
    hRawyield->SetBinContent(iPtBin+1,dRawyield);
    hRawyield->SetBinError(iPtBin+1,dRawyieldErr);
    hRelErrSig->SetBinContent(iPtBin+1,dSignal > 0 ? dSignalErr/dSignal : 0);
    hInvSignif->SetBinContent(iPtBin+1,dSignificance > 0 ? 1/dSignificance : 0);
    hInvSignif->SetBinError(iPtBin+1,dSignificance > 0 ? dSignificanceErr/(dSignificance*dSignificance) : 0);
    hSignal->SetBinContent(iPtBin+1,dSignal); //consider sigma
    hSignal->SetBinError(iPtBin+1,dSignalErr);
    hBackground->SetBinContent(iPtBin+1,dBackground); //consider sigma
    hBackground->SetBinError(iPtBin+1,dBackgroundErr);
    // hSignalBackground->SetBinContent(iPtBin+1,dSignal/dBackground); //consider sigma
    // hSignalBackground->SetBinError(iPtBin+1,1/(dBackground*dBackground)*((dSignalErr*dSignalErr+dSignal*dSignal*dBackgroundErr*dBackgroundErr/(dBackground*dBackground))));
    hBackgroundNormSigma->SetBinContent(iPtBin+1,dBackground/(3*fitter[iPtBin]->GetSigma())*(3*0.012)); //consider sigma
    hBackgroundNormSigma->SetBinError(iPtBin+1,dBackgroundErr);
    hSignificance->SetBinContent(iPtBin+1,dSignificance);
    hSignificance->SetBinError(iPtBin+1,dSignificanceErr);
    hMassMean->SetBinContent(iPtBin+1,fitMassMean);
    hMassMean->SetBinError(iPtBin+1,fitter[iPtBin]->GetMeanUncertainty());
    hSigma->SetBinContent(iPtBin+1,fitMassSigma);
    hSigma->SetBinError(iPtBin+1,fitter[iPtBin]->GetSigmaUncertainty());

	}

  hSignalBackground = (TH1D*) hSignal->Clone("hSignalBackground");
  hSignalBackground->Divide(hBackground);

   // Plotting the MassMean and Sigma
  TCanvas *cMassSigma=new TCanvas("cMassSigma","Fit params",1400,600);
  cMassSigma->Divide(2,1);
  cMassSigma->cd(1);
  hMassMean->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hMassMean->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
  hMassMean->GetYaxis()->SetTitleOffset(1.6);
  hMassMean->GetYaxis()->SetRangeUser(2.280,2.300);
  hMassMean->SetTitle("Mean from gaussian");
  hMassMean->SetMarkerStyle(20);
  hMassMean->Draw("PE");
  fitMassMeanTotal = fitMassMeanTotal/Double_t(nFit);
  fitMassMeanErrTotal = TMath::Sqrt(fitMassMeanErrTotal)/Double_t(nFit);
  TLine* lineMassMean = new TLine(ptBins[0],fitMassMeanTotal,ptBins[nPtBins],fitMassMeanTotal);
  TLine* lineMassMeanUp = new TLine(ptBins[0],fitMassMeanTotal+fitMassMeanErrTotal,ptBins[nPtBins],fitMassMeanTotal+fitMassMeanErrTotal);
  TLine* lineMassMeanLow = new TLine(ptBins[0],fitMassMeanTotal-fitMassMeanErrTotal,ptBins[nPtBins],fitMassMeanTotal-fitMassMeanErrTotal);
  lineMassMean->SetLineColorAlpha(kGray+1,1.0); lineMassMean->SetLineWidth(1); lineMassMean->SetLineStyle(2); lineMassMean->Draw("Same");
  lineMassMeanUp->SetLineColorAlpha(kGray+1,1.0); lineMassMeanUp->SetLineWidth(1); lineMassMeanUp->SetLineStyle(3); lineMassMeanUp->Draw("Same");
  lineMassMeanLow->SetLineColorAlpha(kGray+1,1.0); lineMassMeanLow->SetLineWidth(1); lineMassMeanLow->SetLineStyle(3);lineMassMeanLow->Draw("Same");
  TPaveText *paveMean=new TPaveText(0.470941,0.814665,0.907576,0.917763,"NDC");
  paveMean->SetBorderSize(0);
  paveMean->SetFillStyle(0);
  paveMean->AddText(Form("Mean (uncorr.) = %.4f #pm %.4f",fitMassMeanTotal,fitMassMeanErrTotal));
  paveMean->SetTextAlign(32);
  paveMean->Draw();

  cMassSigma->cd(2);
  hSigma->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hSigma->GetYaxis()->SetTitle("Sigma (GeV/c^{2})");
  hSigma->GetYaxis()->SetTitleOffset(1.6);
  hSigma->SetTitle("Sigma from gaussian");
  hSigma->SetMarkerStyle(20);
  hSigma->Draw("PE");
  fitMassSigmaTotal = fitMassSigmaTotal/Double_t(nFit);
  fitMassSigmaErrTotal = TMath::Sqrt(fitMassSigmaErrTotal)/Double_t(nFit);
  TLine* lineMassSigma = new TLine(ptBins[0],fitMassSigmaTotal,ptBins[nPtBins],fitMassSigmaTotal);
  TLine* lineMassSigmaUp = new TLine(ptBins[0],fitMassSigmaTotal+fitMassSigmaErrTotal,ptBins[nPtBins],fitMassSigmaTotal+fitMassSigmaErrTotal);
  TLine* lineMassSigmaLow = new TLine(ptBins[0],fitMassSigmaTotal-fitMassSigmaErrTotal,ptBins[nPtBins],fitMassSigmaTotal-fitMassSigmaErrTotal);
  lineMassSigma->SetLineColorAlpha(kGray+1,1.0); lineMassSigma->SetLineWidth(1); lineMassSigma->SetLineStyle(2); lineMassSigma->Draw("Same");
  lineMassSigmaUp->SetLineColorAlpha(kGray+1,1.0); lineMassSigmaUp->SetLineWidth(1); lineMassSigmaUp->SetLineStyle(3); lineMassSigmaUp->Draw("Same");
  lineMassSigmaLow->SetLineColorAlpha(kGray+1,1.0); lineMassSigmaLow->SetLineWidth(1); lineMassSigmaLow->SetLineStyle(3);lineMassSigmaLow->Draw("Same");
  TPaveText *paveSigma=new TPaveText(0.470941,0.814665,0.907576,0.917763,"NDC");
  paveSigma->SetBorderSize(0);
  paveSigma->SetFillStyle(0);
  paveSigma->AddText(Form("Mean (uncorr.) = %.4f #pm %.4f",fitMassSigmaTotal,fitMassSigmaErrTotal));
  paveSigma->SetTextAlign(32);
  paveSigma->Draw();


// Plot of Yield
  TCanvas* cYield=new TCanvas("cYield","Results",1400,600);
  cYield->Divide(2,1);
  cYield->cd(1);
  gPad->SetLogy(1);
  hRawyield->SetMarkerStyle(20);
  hRawyield->SetMarkerColor(4);
  hRawyield->SetLineColor(4);
  hRawyield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hRawyield->GetYaxis()->SetTitle("Yield");
  hRawyield->GetYaxis()->SetTitleOffset(1.4);
  hRawyield->SetTitle(Form("Yield (%d#sigma)",nSigmaCount));
  hRawyield->Draw("P");
  hCntSig1->SetMarkerStyle(26);
  hCntSig1->SetMarkerColor(2);
  hCntSig1->SetLineColor(2);
  hCntSig1->Draw("PSAME");
  hCntSig2->SetMarkerStyle(29);
  hCntSig2->SetMarkerColor(kGray+1);
  hCntSig2->SetLineColor(kGray+1);
  hCntSig2->Draw("PSAME");
  TLegend* leg=new TLegend(0.5576,0.7422,0.8997,0.9007);
  leg->SetFillColor(0);
  TLegendEntry* ent=leg->AddEntry(hRawyield,"From Fit","PL");
  ent->SetTextColor(hRawyield->GetMarkerColor());
  ent=leg->AddEntry(hCntSig1,"From Counting1","PL");
  ent->SetTextColor(hCntSig1->GetMarkerColor());
  ent=leg->AddEntry(hCntSig2,"From Counting2","PL");
  ent->SetTextColor(hCntSig2->GetMarkerColor());
  leg->Draw();
  cYield->cd(2);
  hNDiffCntSig1->SetMarkerStyle(26);
  hNDiffCntSig1->SetMarkerColor(2);
  hNDiffCntSig1->SetLineColor(2);
  hNDiffCntSig1->SetTitle("Cmp Fit-Count;#it{p}_{T} (GeV/c);(S_{fit}-S_{count})/S_{fit}");
  hNDiffCntSig1->GetYaxis()->SetTitleOffset(1.2);  
  hNDiffCntSig1->Draw("P");
  hNDiffCntSig2->SetMarkerStyle(29);
  hNDiffCntSig2->SetMarkerColor(kGray+1);
  hNDiffCntSig2->SetLineColor(kGray+1);
  hNDiffCntSig2->Draw("PSAME");
  TLegend* legYield=new TLegend(0.5878,0.7951,0.8991,0.9000);
  legYield->SetFillColor(0);
  TLegendEntry* entYield=legYield->AddEntry(hNDiffCntSig1,"From Counting1","PL");
  entYield->SetTextColor(hNDiffCntSig1->GetMarkerColor());
  entYield=legYield->AddEntry(hNDiffCntSig2,"From Counting2","PL");
  entYield->SetTextColor(hNDiffCntSig2->GetMarkerColor());
  legYield->Draw();
  TLine* line0 = new TLine(ptBins[0],0.,ptBins[nPtBins],0.);
  line0->SetLineColorAlpha(kGray+1,1.0); line0->SetLineWidth(1); line0->SetLineStyle(2); line0->Draw("Same");


  /// Plot of Signal,Background, Signal/Background, Significance
  TCanvas* cResults=new TCanvas("cResults","Results",1200,1200);
  cResults->Divide(2,2);
  cResults->cd(1);
  gPad->SetLogy(1);
  hBackground->SetMarkerStyle(20);
  hBackground->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hBackground->GetYaxis()->SetTitle("Yield");
  hBackground->SetTitle(Form("Yield BG from fit (%d#sigma)",nSigmaCount));
  hBackground->Draw("P");
  cResults->cd(2);
  gPad->SetLogy(1);
  hSignal->SetMarkerStyle(20);
  hSignal->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hSignal->GetYaxis()->SetTitle("Yield");
  hSignal->SetTitle(Form("Yield signal from fit (%d#sigma)",nSigmaCount));
  hSignal->Draw("P");
  cResults->cd(3);
  hSignalBackground->SetMarkerStyle(20);
  hSignalBackground->Draw("P");
  //hSignalBackground->GetYaxis()->SetRangeUser(0.0,20); 
  hSignalBackground->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hSignalBackground->GetYaxis()->SetTitle("Signal/Background");
  hSignalBackground->SetTitle(Form("Ratio signal/BG from fit (%d#sigma)",nSigmaCount));
  cResults->cd(4);
  hSignificance->SetMarkerStyle(20);
  hSignificance->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hSignificance->GetYaxis()->SetTitle("Significance");
  hSignificance->SetTitle(Form("Significance from fit (%d#sigma)",nSigmaCount));
  hSignificance->Draw("P");

    TCanvas* cDiffS=new TCanvas("cDiffS","Difference",600,600);
  hRelErrSig->SetMarkerStyle(20); //fullcircle
  hRelErrSig->SetTitleOffset(1.4);  
  hRelErrSig->SetTitle("Rel Error from Fit;#it{p}_{T} (GeV/c);Signal Relative Error");
  hRelErrSig->GetYaxis()->SetTitleOffset(1.3);
  hRelErrSig->Draw("P");
  hInvSignif->SetMarkerStyle(21); //fullsquare
  hInvSignif->SetMarkerColor(kMagenta+1);
  hInvSignif->SetLineColor(kMagenta+1);
  hInvSignif->Draw("PSAMES");
  TLegend* leg2=new TLegend(0.5878,0.7951,0.8991,0.9000);
  leg2->SetFillColor(0);
  TLegendEntry* ent2=leg2->AddEntry(hRelErrSig,"From Fit","P");
  ent2->SetTextColor(hRelErrSig->GetMarkerColor());
  ent2=leg2->AddEntry(hInvSignif,"1/Significance","PL");
  ent2->SetTextColor(hInvSignif->GetMarkerColor());
  leg2->Draw();

  TCanvas *cChi2=new TCanvas("cChi2","Reduced chi squared",800,600);
  cChi2->cd();
  TGraph* grReducedChiSquare=new TGraphErrors(nPtBins,ptBinsCenter,dChiSquare,ptBinsError,0);
  // gPad->SetLogy(1);
  grReducedChiSquare->SetName("grReducedChiSquare");
  grReducedChiSquare->SetTitle("Reduced chi squared;#it{p}_{T} (GeV/c);#chi^{2}_{red}");
  grReducedChiSquare->SetMarkerStyle(21);
  grReducedChiSquare->Draw("AP");

  TCanvas* cbkgNormSigma=new TCanvas("cbkgNormSigma","Background normalized to sigma",400,500);
  cbkgNormSigma->cd();
  gPad->SetLogy(1); 
  hBackgroundNormSigma->SetMarkerStyle(20);
  hBackgroundNormSigma->SetTitle("Background normalized to sigma");
  hBackgroundNormSigma->GetXaxis()->SetTitle("#it{p}_{T} (GeV/c)");
  hBackgroundNormSigma->GetYaxis()->SetTitle("Background #times 3 #times 0.012/ (3 #times #sigma)");
  hBackgroundNormSigma->GetYaxis()->SetTitleOffset(1.4); 
  hBackgroundNormSigma->Draw("P");

  outputFile->cd();

  // Save Canvases
  cInvariantMass   ->SaveAs(Form("%s/MassFit_Pt.pdf"             ,gOutputDirname.Data()));
  cMassSigma    ->SaveAs(Form("%s/MassSigma.pdf"              ,gOutputDirname.Data()));
  cResults      ->SaveAs(Form("%s/Results.pdf"                ,gOutputDirname.Data()));
  cYield        ->SaveAs(Form("%s/Yield.pdf"                  ,gOutputDirname.Data()));
  cDiffS        ->SaveAs(Form("%s/Difference.pdf"             ,gOutputDirname.Data()));
  cChi2         ->SaveAs(Form("%s/ReducedChi.pdf"             ,gOutputDirname.Data()));
  cbkgNormSigma ->SaveAs(Form("%s/BackgroundNormSigma.pdf"    ,gOutputDirname.Data()));
  
  hMassMean->Write();
  hSigma->Write();
  hCntSig1->Write();
  hCntSig2->Write();
  hNDiffCntSig1->Write();
  hNDiffCntSig2->Write();
  hRawyield->Write();
  hSignalBackground->Write();
  hRelErrSig->Write();
  hInvSignif->Write();
  hSignal->Write();
  hBackground->Write();
  hBackgroundNormSigma->Write();
  hSignificance->Write();
  grReducedChiSquare->Write();

  Printf("\nSigma + SigmaErr");
  printf("{");
  for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
    if(iPtBin != (nPtBins-1)) printf("%.6f,",hSigma->GetBinContent(iPtBin+1));
    else                      printf("%.6f", hSigma->GetBinContent(iPtBin+1));
  }
  printf("}\n{");
  for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
    if(iPtBin != (nPtBins-1)) printf("%.6f,",hSigma->GetBinError(iPtBin+1));
    else                      printf("%.6f", hSigma->GetBinError(iPtBin+1));
  }
  printf("}");

  Printf("\n\nMean + MeanErr");
  printf("{");
  for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
    if(iPtBin != (nPtBins-1)) printf("%.6f,",hMassMean->GetBinContent(iPtBin+1));
    else                      printf("%.6f", hMassMean->GetBinContent(iPtBin+1));
  }
  printf("}\n{");
  for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
    if(iPtBin != (nPtBins-1)) printf("%.6f,",hMassMean->GetBinError(iPtBin+1));
    else                      printf("%.6f", hMassMean->GetBinError(iPtBin+1));
  }
  printf("}");
  std::cout << std::endl << std::flush;

}

void HFMassFitter(TString fileName = "MCInvariantMass_D2H_stdCuts_mult_0_999_cbOnly.root",
                  TString outputDirName = "OutputHF",
                  Bool_t keepCanvasOpen = kTRUE
                  ){

  std::cout << std::endl << std::flush;
  Printf("========= HFMassFitter.C Started =========");

  gOutputDirname = outputDirName;

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

  TList* list;
  for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
    // pT Bin Ranges
    Double_t pTBinLow = ptBins[iPtBin];
    Double_t pTBinUp  = ptBins[iPtBin+1];
    hMass[iPtBin] = (TH1D*)inputFile->Get(Form("hMass_pT_%.0f_%.0f",pTBinLow,pTBinUp));
    if(!hMass[iPtBin]){
      Printf("ERROR: hMass_pT_%.0f_%.0f not in inputFile",pTBinLow,pTBinUp);
      return NULL;
    }
  }

  // Open output file 
  gSystem->Exec(Form("rm -r %s/",outputDirName.Data()));
  gSystem->Exec(Form("mkdir %s",outputDirName.Data()));
  TFile *outputFile(0);
  outputFile = new TFile(Form("%s/FitResults.root",outputDirName.Data()),"RECREATE");
  if(!outputFile){
    Printf("ERROR: outputFile not created");
    exit(1);
  }

  for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
    // pT Bin Ranges
    Double_t pTBinLow = ptBins[iPtBin];
    Double_t pTBinUp  = ptBins[iPtBin+1];

    outputFile->cd();
    directories[iPtBin] = outputFile->mkdir(Form("pT_%.1f_%.1f",pTBinLow,pTBinUp));
  }

  Printf("* Input file: %s", inputFile->GetName());
  Printf("* Output file: %s", outputFile->GetName());

  SetFitParameters();

  MassFitter(outputFile);

  if(!keepCanvasOpen) outputFile->Close();

  return;
}


