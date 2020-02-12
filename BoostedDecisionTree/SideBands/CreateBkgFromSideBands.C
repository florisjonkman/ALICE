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

// Int_t typeB = kPowEx;              // Type of background fit function
// Int_t typeS = kGaus;                 // Type of signal fit function
Int_t factor4refl = 0;                // 0 = gaus; 1 = gaus+gaus broadened
Int_t nSigmaCount = 3;

Int_t 	 rebin[nPtBins];
Bool_t 	 isFit[nPtBins];
Int_t 	 typeB[nPtBins];
Int_t 	 typeS[nPtBins];
// Int_t f4ref[nPtBins];
Double_t minFit[nPtBins];
Double_t maxFit[nPtBins];

// Standard cuts
Double_t meanStdCut[nPtBins];
Double_t sigmaStdCut[nPtBins];

// Plot range
Double_t minMassForPlot = 2.149;
Double_t maxMassForPlot = 2.441;

// Fit parameters
// Mean of gaussian
Double_t gaussSet = 2.289;
Double_t gaussLow = 2.282;   // Lower bound for fit [GeV]
Double_t gaussUp  = 2.295;   // Upper bound fot fit [GeV]

// Sigma of gaussian
Double_t sigmaSet = 0.009;
Double_t sigmaLow = 0.003;
Double_t sigmaUp =  0.018;

Double_t massPi = TDatabasePDG::Instance()->GetParticle(211)->Mass();
TString gOutputDirname;

// Fit settings
// Int_t typeB = kPowEx;					    // Type of background fit function
// Int_t typeS = kGaus;					    // Type of signal fit function
// Int_t factor4refl = 0;				        // 0 = gaus; 1 = gaus+gaus broadened
// Double_t minMassForFit = 2.2;
// Double_t maxMassForFit = 2.38;
// Double_t initialGaussMean = 2.289;
// Double_t initialGaussSigma = 0.010;

Double_t dSignificance,dSignificanceErr,dSignal,dSignalErr,dSignalErr2,dBackground,dBackgroundErr;
Double_t ptBinsCenter[nPtBins] = {0.0};
Double_t ptBinsError[nPtBins] = {0.0};
Double_t dChiSquare[nPtBins] = {0.0};

TH1F** hMass = new TH1F*[nPtBins];
TH1F** hRebinned = new TH1F*[nPtBins];
TDirectory* directories[nPtBins];

Bool_t flagIsFitted[nPtBins];
Double_t cutMassMean[nPtBins];
Double_t cutMassSigma[nPtBins];

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
TH1D* hSignificance       =new TH1D("hSignificance","hSignificance",nPtBins,ptBins);
TH1D* hMassMean           =new TH1D("hMassMean","hMassMean",nPtBins,ptBins);
TH1D* hSigma              =new TH1D("hSigma","hSigma",nPtBins,ptBins);

void SetFitParameters(){

	// Fit parameters per Pt Bin
	// Fit?             Rebin            BackgroudFunc        SignalFunct           MinimumFitRange       MaxFitRange
	isFit[0] = kFALSE;  rebin[0]  = 8;   typeB[0]  = kPol2;   typeS[0]  = kGaus;    minFit[0]  = 2.150;   maxFit[0]  = 2.440;  meanStdCut[0] = 2.288;  sigmaStdCut[0] = 0.008;
	isFit[1] = kFALSE;   rebin[1]  = 9;   typeB[1]  = kPol2;   typeS[1]  = kGaus;    minFit[1]  = 2.150;   maxFit[1]  = 2.440; meanStdCut[1] = 2.288;  sigmaStdCut[1] = 0.008;
	isFit[2] = kFALSE;   rebin[2]  = 16;  typeB[2]  = kPol2;   typeS[2]  = kGaus;    minFit[2]  = 2.150;   maxFit[2]  = 2.440; meanStdCut[2] = 2.290;  sigmaStdCut[2] = 0.010;
	isFit[3] = kFALSE;   rebin[3]  = 18;  typeB[3]  = kPol2;   typeS[3]  = kGaus;    minFit[3]  = 2.150;   maxFit[3]  = 2.440; meanStdCut[3] = 2.290;  sigmaStdCut[3] = 0.011;
	isFit[4] = kFALSE;   rebin[4]  = 18;  typeB[4]  = kPol2;   typeS[4]  = kGaus;    minFit[4]  = 2.150;   maxFit[4]  = 2.440; meanStdCut[4] = 2.292;  sigmaStdCut[4] = 0.014;
	isFit[5] = kFALSE;   rebin[5]  = 32;  typeB[5]  = kPol2;   typeS[5]  = kGaus;    minFit[5]  = 2.150;   maxFit[5]  = 2.440; meanStdCut[5] = 2.295;  sigmaStdCut[5] = 0.016;
	
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

void SetStandardValues(){

	for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
		flagIsFitted[iPtBin]=kFALSE;
		cutMassMean[iPtBin] = meanStdCut[iPtBin];
		cutMassSigma[iPtBin] = sigmaStdCut[iPtBin];
	}
}

void GenerateInvariantMassHistograms(TTree *tree, Bool_t flagDrawCanvas, Bool_t flagKeepCanvasOpen){

  	// Create Canvas
  	TCanvas* cInvariantMass_Pt;
  	if(flagDrawCanvas){
  		cInvariantMass_Pt = new TCanvas("cInvariantMass_Pt", "InvariantMass", 1200, 550);
  		if(nPtBins<=6){
  			cInvariantMass_Pt->Divide(3,2);
  		}
  		else{
  			cInvariantMass_Pt->Divide(4,2);
  		}
  	}

  	// Bins
	const Int_t nMassBins = 1000;
	Double_t dMassLow = 2.050;		// [GeV]
	Double_t dMassUp = 2.550;		// [GeV]
	Double_t dMassBins[nMassBins+1] = {0.0};
	for(Int_t i = 0; i <= nMassBins; i++) dMassBins[i]=((dMassUp-dMassLow)/nMassBins)*i+dMassLow;

	// Array of Histograms, x = mass

	for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
		// pT Bin Ranges
    	Double_t pTBinLow = ptBins[iPtBin];
    	Double_t pTBinUp  = ptBins[iPtBin+1];

    	hMass[iPtBin] = new TH1F(Form("hMass_pT_%.0f_%.0f",pTBinLow,pTBinUp),Form("hMass_pT_%.0f_%.0f",pTBinLow,pTBinUp),nMassBins,dMassBins);
	}

  	// TH2F* hPtMass = new TH2F("hPtMassK0Short","Mass per pT;pT [GeV/c];M [GeV/c^{2}]",nPtBins,ptBins,nMassBins,dMassBins);

  	// Initialize tree variables
  	Float16_t massLc;
  	Float16_t pTLc;

	tree->SetBranchAddress("massLc2K0Sp", &massLc);
	tree->SetBranchAddress("LcPt", &pTLc);

    // Fill histograms
    Long64_t nEntries = tree->GetEntries(); // 100000; //

    Printf("\n--- Processing: %10lld candidates",nEntries);

	for(Long64_t iLC = 0; iLC <nEntries; iLC++){
		// Get Entry
		tree->GetEntry(iLC);

		if (iLC % (nEntries/10) == 0) Printf("... Processing: %10lld \t %.0f%%",iLC,Double_t(iLC)/Double_t(nEntries)*100.);
		

		for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
			// pT Bin Ranges
	    	Double_t pTBinLow = ptBins[iPtBin];
	    	Double_t pTBinUp  = ptBins[iPtBin+1];

	    	if(pTLc>pTBinLow && pTLc <=pTBinUp){
	    		hMass[iPtBin]->Fill( (Float_t) massLc );
	    		break;
	    	}
		}
	}

	// Draw the histograms
	for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
		// pT Bin Ranges
    	Double_t pTBinLow = ptBins[iPtBin];
    	Double_t pTBinUp  = ptBins[iPtBin+1];

		hMass[iPtBin]->SetTitle(Form("%.0f < #it{p}_{T} < %.0f",pTBinLow,pTBinUp));
		hMass[iPtBin]->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
  		hMass[iPtBin]->GetYaxis()->SetTitle(Form("Entries/(%.1f MeV/c^{2})",(hMass[iPtBin]->GetBinWidth(1)*1000)));
  		hMass[iPtBin]->SetLineColor(kBlue+1);
    	hMass[iPtBin]->SetMarkerColor(kBlack);
    	hMass[iPtBin]->SetMarkerStyle(8);
    	hMass[iPtBin]->SetMarkerStyle(8);
    	hMass[iPtBin]->SetMarkerSize(0.5);
    	hMass[iPtBin]->GetYaxis()->SetTitleOffset(1.6);
    	if(flagDrawCanvas){
    		cInvariantMass_Pt->cd(iPtBin+1); // Open canvas
			hMass[iPtBin]->Draw("E0");

    	}
	}

	if(flagDrawCanvas){
		cInvariantMass_Pt->SaveAs(Form("%s/hMass.pdf",gOutputDirname.Data()));
		if(!flagKeepCanvasOpen) cInvariantMass_Pt->Close();
	}
}

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

void MassFitter(TFile* outputFile, Bool_t flagDraw = kTRUE){

  Printf("\n================================================");
  Printf("         Massfitter started");
  Printf("================================================\n");

  // Create Canvas
  TCanvas* cInvariantMass;
  if(flagDraw){
    cInvariantMass = new TCanvas("cInvariantMass", "InvariantMass", 1200, 550);
    cInvariantMass->Divide(3,2);
  }
  
  TF1* funBckStore1=0x0;
  TF1* funBckStore2=0x0;
  TF1* funBckStore3=0x0;

  Double_t fitMassMeanTotal =0;
  Double_t fitMassMeanErrTotal =0;
  Double_t fitMassSigmaTotal =0;
  Double_t fitMassSigmaErrTotal =0;

  outputFile->cd();

  for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){

    if(flagDraw) cInvariantMass->cd(iPtBin+1); // Open canvas 

    // pT Bin Ranges
    Double_t pTBinLow = ptBins[iPtBin];
    Double_t pTBinUp  = ptBins[iPtBin+1];
    Printf("\n%s---------------------------------\n",MAGENTA);
    Printf("=== %.1f < pT < %.1f ===",pTBinLow,pTBinUp);
    Printf("\n---------------------------------%s\n",RESET);

    // pT Center and Error
    ptBinsCenter[iPtBin] = (pTBinLow+pTBinUp)/2.;
    ptBinsError[iPtBin]  = TMath::Abs((pTBinUp-pTBinLow)/2.);

    if(minMassForPlot > minFit[iPtBin]) minMassForPlot = minFit[iPtBin];
    if(maxMassForPlot < maxFit[iPtBin]) maxMassForPlot = maxFit[iPtBin];

    hRebinned[iPtBin]   = (TH1F*) AliVertexingHFUtils::RebinHisto(hMass[iPtBin],rebin[iPtBin],-1);
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
    Double_t hMin=TMath::Max(minFit[iPtBin],hRebinned[iPtBin]->GetBinLowEdge(1));
    Double_t hMax=TMath::Min(maxFit[iPtBin],hRebinned[iPtBin]->GetBinLowEdge(hRebinned[iPtBin]->GetNbinsX()+1));

    // Get minimum and maximum bin content within fit range
    Double_t contentMax = 0;
    Double_t contentMin = TMath::Power(10,99);
    for(Int_t iBin = 1; iBin <= nBinsRebinned; iBin++){
        if(hRebinned[iPtBin]->GetBinCenter(iBin) >= hMin && hRebinned[iPtBin]->GetBinCenter(iBin) <= hMax){
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
      hRebinned[iPtBin]->GetYaxis()->SetRangeUser(plotYmin,plotYmax);
      hRebinned[iPtBin]->Draw();
    }    
    
    // If not fitted, only plot the histogram
    if(!isFit[iPtBin]) continue;

    TString signalFunc;
    TString bgFuncOnly;
    TString bgFunc;
    Int_t nParsB;
    Int_t nParsS;
    Int_t nPars;
    Double_t* x;
    Double_t* par;

    // Signal fit fuction
    TF1 *fitFunc;
    if(typeS[iPtBin]==kGaus){
      //signalFunc = "gaus(0)";
      signalFunc = "[0]*TMath::Gaus(x,[1],[2])";
      nParsS=3;
    }
    else if(typeS[iPtBin]==kDoubleGaus){
      //signalFunc = "gaus(0)+gaus(3)";
      signalFunc = "[0]*TMath::Gaus(x,[1],[2])+[3]*TMath::Gaus(x,[4],[5])";
      nParsS=6;
    }
    else{
      Printf("ABORT FIT: Unknown signal fit function");
      continue;
    }

    // BG fit function 
    if(typeB[iPtBin]==kExpo){
      bgFunc = Form("expo(%d)",nParsS);
      nParsB=2;
    }
    else if(typeB[iPtBin]==kLinear){
      bgFuncOnly = "pol1";
      bgFunc = Form("pol1(%d)",nParsS);
      nParsB=2;
    }
    else if(typeB[iPtBin]==kPol2){
      bgFuncOnly = "pol2";
      bgFunc = Form("pol2(%d)",nParsS);
      nParsB=3;
    }
    else if(typeB[iPtBin]==kPol3){
      bgFuncOnly = "pol3";
      bgFunc = Form("pol3(%d)",nParsS);
      nParsB=4; 
    }
    else if(typeB[iPtBin]==kPow){
      bgFuncOnly = Form("[0]*TMath::Power((x-%f),[1])",massPi);
      bgFunc = Form("[%d]*TMath::Power((x-%f),[%d])",nParsS,massPi,nParsS+1);
      nParsB=2;
    }
    else if(typeB[iPtBin]==kPowEx){
      bgFuncOnly = Form("[0]*TMath::Sqrt(x-%f)*TMath::Exp(-[1]*(x-%f))",massPi,massPi);
      bgFunc = Form("[%d]*TMath::Sqrt(x-%f)*TMath::Exp(-[%d]*(x-%f))",nParsS,massPi,nParsS+1,massPi);
      nParsB=2;
    }
    else{
      Printf("ABORT FIT: Unknown background fit function");
      continue;
    }

    // Total function
    char* totalFunc = Form("%s+%s",signalFunc.Data(),bgFunc.Data());
    fitFunc = new TF1("fitFunc",totalFunc,hMin,hMax);
    TF1* bgFitFunc = new TF1("bgFitFunc",bgFuncOnly,hMin,hMax);
    nPars = nParsS + nParsB;
    Printf("\n> Fit function");
    Printf("  %s",totalFunc);

    // Set parameters
    // Height of gaussian
    Double_t setPeak = (contentMax-contentMin);
    Double_t peakUp  = contentMax+((contentMax-contentMin)/2.); //contentMax+((contentMax-contentMin)/2.)
    Double_t peakLow = 0;
    // fitFunc->SetParameter(0,setPeak);         // Set peak: Max content
    fitFunc->SetParLimits(0,0,peakUp);           // Set Limits: [contentMin,contentMax+difference]
    if(typeS[iPtBin]==kDoubleGaus){
      // fitFunc->SetParameter(3,setPeak);       // Set mean: mass PDG particle
      fitFunc->SetParLimits(3,0,peakUp);       // Set Limits
    }

    // Mean of gaussian
    fitFunc->SetParameter(1,gaussSet);         // Set mean: mass PDG particle
    fitFunc->SetParLimits(1,gaussLow,gaussUp);   // Set Limits
    if(typeS[iPtBin]==kDoubleGaus){
      fitFunc->SetParameter(4,gaussSet);         // Set mean: mass PDG particle
      fitFunc->SetParLimits(4,gaussLow,gaussUp); // Set Limits
    }

    // Sigma of gaussian
    fitFunc->SetParameter(2,sigmaSet);           
    fitFunc->SetParLimits(2, sigmaLow, sigmaUp);
    if(typeS[iPtBin]==kDoubleGaus){
      fitFunc->SetParameter(5,sigmaSet);       // Set mean: mass PDG particle
      fitFunc->SetParLimits(5,sigmaLow,sigmaUp);   // Set Limits
    }

    // Print set values
    Printf("");
    Printf("> Values set for fitting");
    Printf("  * Peak guassian: %.0f [%.0f,%.0f] Entries", setPeak,peakUp,peakLow);
    Printf("  * Mean gaussian: %.5f [%.5f,%.5f] GeV/c^2", gaussSet,gaussLow,gaussUp);
    Printf("  * Sigma gaussian: %.5f [%.5f,%.5f] GeV/c^2", sigmaSet,sigmaLow,sigmaUp);
    Printf("");

    // Fit function
    hRebinned[iPtBin]->Fit("fitFunc","0RM");
    TString statusFit = (gMinuit->fCstatu);
    if(statusFit.Contains("FAILED") || statusFit.Contains("CALL LIMIT")){
      std::cout << RED << "FAILED: Full fitfuction did not converge, status:" << statusFit.Data() << RESET << std::endl;
      Printf("  > Trying to fit background function");
      hRebinned[iPtBin]->Fit(bgFitFunc,"R");
      TString statusFit = (gMinuit->fCstatu);
      if(statusFit.Contains("FAILED") || statusFit.Contains("CALL LIMIT")){
        std::cout << RED << "  > FAILED: Background function did not converge, status:" << statusFit.Data() << RESET << std::endl;
      }
      continue;
    }
    // else{
    //  Printf("%sSUCCES: Full fitfuction fitted, status: %s%s",GREEN,statusFit.Data(),RESET);
    // }    

    // Getting the fit results
    TF1 *fitResult = hRebinned[iPtBin]->GetFunction("fitFunc");
    TVirtualFitter * fitter = TVirtualFitter::GetFitter();
    assert(fitter != 0);
    Double_t* covMatrix = fitter->GetCovarianceMatrix();
    Printf("\n> Total Covariance matrix");
    for (Int_t i = 0; i < nPars; ++i) printf("\t%d",i);
    Printf("");
    for (Int_t i = 0; i < nPars; ++i)
    {
       printf("%d: ",i);
       for (Int_t j = 0; j < nPars; ++j)
       {
          printf(" %4f ",TMath::Sqrt(TMath::Abs(covMatrix[j+i*nPars])));
          if(j == nPars-1) printf("\n");
       }
    }

    // Parameters from Fit
    Double_t dRawyieldS         = fitResult->GetParameter(0);
    Double_t dRawyieldErrS      = fitResult->GetParError(0);
    Double_t fitMassMean        = fitResult->GetParameter(1);
    Double_t fitMassMeanErr     = fitResult->GetParError(1);
    Double_t fitMassSigma       = fitResult->GetParameter(2);
    Double_t fitMassSigmaErr    = fitResult->GetParError(2);
    Double_t dRawyieldS2, dRawyieldErrS2, fitMassMean2, fitMassMeanErr2, fitMassSigma2,fitMassSigmaErr2;
    if(typeS[iPtBin]==kDoubleGaus){
      dRawyieldS2         = fitResult->GetParameter(3);
      dRawyieldErrS2      = fitResult->GetParError(3);
      fitMassMean2        = fitResult->GetParameter(4);
      fitMassMeanErr2     = fitResult->GetParError(4);
      fitMassSigma2       = fitResult->GetParameter(5);
      fitMassSigmaErr2    = fitResult->GetParError(5);
      fitMassMean = (fitMassMean + fitMassMean2)/2.;
      fitMassMeanErr = (1./2.)*TMath::Sqrt(fitMassMeanErr*fitMassMeanErr+fitMassMeanErr2+fitMassMeanErr2);
      fitMassSigma = (fitMassSigma + fitMassSigma2)/2.;
      fitMassSigmaErr = (1./2.)*TMath::Sqrt(fitMassSigmaErr*fitMassSigmaErr+fitMassSigmaErr2+fitMassSigmaErr2);
    }
    dChiSquare[iPtBin]          = fitResult->GetChisquare() / Double_t(fitResult->GetNDF());

    flagIsFitted[iPtBin]=kTRUE;
    cutMassMean[iPtBin]=fitMassMean;
    cutMassSigma[iPtBin]=fitMassSigma;
   	Printf("%sSUCCES: background and signal fitted%s",GREEN,RESET);
   	Printf("%sMean = %.3f, sigma = %.3f%s",GREEN,cutMassMean[iPtBin],cutMassSigma[iPtBin],RESET);


    // Extract the Background
    TF1* bgFit;
    if(typeB[iPtBin]==kExpo)  bgFit = new TF1("bgFit", "expo(0)", hMin, hMax);
    if(typeB[iPtBin]==kLinear)bgFit = new TF1("bgFit", "pol1(0)", hMin, hMax);
    if(typeB[iPtBin]==kPol2)  bgFit = new TF1("bgFit", "pol2(0)", hMin, hMax);
    if(typeB[iPtBin]==kPol3)  bgFit = new TF1("bgFit", "pol3(0)", hMin, hMax);
    if(typeB[iPtBin]==kPow)   bgFit = new TF1("bgFit", Form("[0]*TMath::Power((x-%f),[1])",massPi), hMin, hMax);
    if(typeB[iPtBin]==kPowEx) bgFit = new TF1("bgFit", Form("[0]**TMath::Sqrt(x-%f)*TMath::Exp(-[1]*(x-%f))",massPi,massPi), hMin, hMax);

    // Parameters of BG
    Double_t pariBG;
    Double_t erriBG;
    Double_t covMatrixB[nParsB*nParsB];

    for(Int_t iPar=0; iPar < nParsB; iPar++){
      pariBG = fitResult->GetParameter(nParsS+iPar);
      erriBG = fitResult->GetParError(nParsS+iPar);
      bgFit->SetParameter(iPar,pariBG);
      bgFit->SetParError(iPar,erriBG);
      for(Int_t jPar=0; jPar < nParsB; jPar++){
        covMatrixB[jPar+iPar*nParsB]=covMatrix[(nPars*nParsS)+nParsS + jPar + iPar*nPars];
      }
    }
    Printf("\n> Bkg Covariance matrix");
    for (Int_t i = 0; i < nParsB; ++i) printf("\t%d",i);
    Printf("");
    for (Int_t i = 0; i < nParsB; ++i)
    {
       printf("%d: ",i);
       for (Int_t j = 0; j < nParsB; ++j)
       {
          printf(" %4f ",TMath::Sqrt(TMath::Abs(covMatrixB[j+i*nParsB])));
          if(j == nParsB-1) printf("\n");
       }
    }

    // Signal and background
    Double_t minInt    = fitMassMean - (nSigmaCount * fitMassSigma);
    Double_t maxInt    = fitMassMean + (nSigmaCount * fitMassSigma);
    
    Double_t integralS = fitResult->Integral(minInt, maxInt)/(hRebinned[iPtBin]->GetBinWidth(4));
    // Double_t errIntS   = fitResult->IntegralError(minInt,maxInt);
    Double_t integralB = bgFit->Integral(minInt, maxInt)/(hRebinned[iPtBin]->GetBinWidth(4));
    // Double_t errIntB   = bgFit->IntegralError(minInt,maxInt,bgFit->GetParameters(),covMatrixB);
    dBackground        = integralB;
    dBackgroundErr     = TMath::Sqrt(integralB);
    dSignal            = integralS - integralB;
    dSignalErr         = TMath::Sqrt(integralS+integralB);

    ComputeSignificance(dSignal,dBackgroundErr,dBackground,dBackgroundErr,dSignificance,dSignificanceErr);

    Printf("\n> Extract parameters");
    Printf("  * Rawyield signal:\t %.0f ± %.0f Entries", dRawyieldS,dRawyieldErrS);
    Printf("  * Mean mass:\t\t %.3f ± %.3f GeV/c^2", fitMassMean,fitMassMeanErr);
    Printf("  * Mean sigma:\t\t %.3f ± %.3f GeV/c^2", fitMassSigma,fitMassSigmaErr);
    if(typeS[iPtBin]==kDoubleGaus){
      Printf("  * Rawyield2 signal:\t %.0f ± %.0f Entries", dRawyieldS2,dRawyieldErrS2);
    Printf("  * Mean mass2:\t\t %.3f ± %.3f GeV/c^2", fitMassMean2,fitMassMeanErr2);
    Printf("  * Mean sigma2:\t %.3f ± %.3f GeV/c^2", fitMassSigma2,fitMassSigmaErr2);
    }
    Printf("  * Reduced chisq:\t %.2f ", dChiSquare[iPtBin]);
    Printf("  * Signal:\t\t %.0f ± %.0f", dSignal,dSignalErr);
    Printf("  * Background:\t\t %.0f ± %f", dBackground,dBackgroundErr);
    Printf("  * Significance:\t %.0f ± %.1f", dSignificance,dSignificanceErr);
    // Printf("  * integralS = %.1f (sqrt %.1f)",integralS,dSignalErr);
    // Printf("  * errIntS (corr) = %.1f",TMath::Sqrt((errIntS*errIntS)/(hRebinned[iPtBin]->GetBinWidth(4))));
    // Printf("  * integralB = %.1f (sqrt %.1f)",integralB,dBackgroundErr);
    // Printf("  * errIntB (corr) = %.1f",TMath::Sqrt((errIntB*errIntB)/(hRebinned[iPtBin]->GetBinWidth(4))));
    Printf("");

    fitMassMeanTotal     += fitMassMean;
    fitMassMeanErrTotal  += fitMassMeanErr*fitMassMeanErr;
    fitMassSigmaTotal    += fitMassSigma;
    fitMassSigmaErrTotal += fitMassSigmaErr*fitMassSigmaErr;

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

        // // Draw fit sigma line
        // TLine* linefitMin = new TLine(hMin,plotYmin,hMin,plotYmin+(plotYmax-plotYmin)*0.25);
        // TLine* linefitMax = new TLine(hMax,plotYmin,hMax,plotYmin+(plotYmax-plotYmin)*0.25);
        // linefitMin->SetLineColorAlpha(kAzure-9,0.1); linefitMin->SetLineStyle(2);
        // linefitMax->SetLineColorAlpha(kAzure-9,0.1); linefitMax->SetLineStyle(2);
        // linefitMin->Draw("Same");
        // linefitMax->Draw("Same");

        // Plot the Background
        bgFit->SetLineColor(2);
        bgFit->SetLineWidth(1);
        bgFit->Draw("Same");

        // Plot the signal
        fitResult->SetRange(hMin,hMax);
        fitResult->SetLineColor(4);
        fitResult->SetLineWidth(1);
        fitResult->Draw("Same");

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

      directories[iPtBin]->cd();
      hMass[iPtBin]->Write();
      hRebinned[iPtBin]->Write();
      fitResult->Write();
      bgFit->Write();

      // Different counting method
      if(iPtBin==0 && bgFit)  funBckStore1=new TF1(*bgFit);
      //  if(iPtBin==0 && fB2)  funBckStore2=new TF1(*fB2);
      //  if(iPtBin==0 && fM)   funBckStore3=new TF1(*fM);

      Float_t cntSig1=0.;
      //  Float_t cntSig2=0.;
      Float_t cntErr=0.;

      Float_t minBinSum = hMass[iPtBin]->FindBin(fitMassMean-nSigmaCount*fitMassSigma);
      Float_t maxBinSum = hMass[iPtBin]->FindBin(fitMassMean+nSigmaCount*fitMassSigma);

      for(Int_t iMB=minBinSum; iMB<=maxBinSum; iMB++){
        Float_t bkg1=bgFit ? bgFit->Eval(hMass[iPtBin]->GetBinCenter(iMB))/rebin[iPtBin] : 0;
        // Float_t bkg2=fB2 ? fB2->Eval(hMass[iPtBin]->GetBinCenter(iMB))/rebin[iPtBin] : 0;
        cntSig1+=(hMass[iPtBin]->GetBinContent(iMB)-bkg1);
        // cntSig2+=(hMass[iPtBin]->GetBinContent(iMB)-bkg2);
        cntErr+=(hMass[iPtBin]->GetBinContent(iMB));
      }

      // Fill histograms
      hCntSig1->SetBinContent(iPtBin+1,cntSig1);
      hCntSig1->SetBinError(iPtBin+1,TMath::Sqrt(cntErr));
      hNDiffCntSig1->SetBinContent(iPtBin+1,(dSignal-cntSig1)/dSignal);
      hNDiffCntSig1->SetBinError(iPtBin+1,TMath::Sqrt(cntErr)/dSignal);
      // hCntSig2->SetBinContent(iPtBin+1,cntSig2);
      // hCntSig2->SetBinContent(iPtBin+1,TMath::Sqrt(cntErr));
      // hNDiffCntSig2->SetBinContent(iPtBin+1,(dSignal-cntSig2)/dSignal);
      // hNDiffCntSig2->SetBinError(iPtBin+1,TMath::Sqrt(cntErr)/dSignal);
      hRawyield->SetBinContent(iPtBin+1,dSignal);
      hRawyield->SetBinError(iPtBin+1,dSignalErr);
      hRelErrSig->SetBinContent(iPtBin+1,dSignalErr/dSignal);
      hInvSignif->SetBinContent(iPtBin+1,1/dSignificance);
      hInvSignif->SetBinError(iPtBin+1,dSignificanceErr/(dSignificance*dSignificance));
      hSignal->SetBinContent(iPtBin+1,dSignal); //consider sigma
      hSignal->SetBinError(iPtBin+1,dSignalErr);
      hBackground->SetBinContent(iPtBin+1,dBackground); //consider sigma
      hBackground->SetBinError(iPtBin+1,dBackgroundErr);
      // hSignalBackground->SetBinContent(iPtBin+1,dSignal/dBackground); //consider sigma
      // hSignalBackground->SetBinError(iPtBin+1,1/(dBackground*dBackground)*((dSignalErr*dSignalErr+dSignal*dSignal*dBackgroundErr*dBackgroundErr/(dBackground*dBackground))));
      hBackgroundNormSigma->SetBinContent(iPtBin+1,dBackground/(3*fitMassSigma)*(3*0.012)); //consider sigma
      hBackgroundNormSigma->SetBinError(iPtBin+1,dBackgroundErr);
      hSignificance->SetBinContent(iPtBin+1,dSignificance);
      hSignificance->SetBinError(iPtBin+1,dSignificanceErr);
      hMassMean->SetBinContent(iPtBin+1,fitMassMean);
      hMassMean->SetBinError(iPtBin+1,fitMassMeanErr);
      hSigma->SetBinContent(iPtBin+1,fitMassSigma);
      hSigma->SetBinError(iPtBin+1,fitMassSigmaErr);
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
  hMassMean->SetTitle("Mean from gaussian");
  hMassMean->SetMarkerStyle(20);
  hMassMean->Draw("PE");
  fitMassMeanTotal = fitMassMeanTotal/Double_t(nPtBins);
  fitMassMeanErrTotal = TMath::Sqrt(fitMassMeanErrTotal)/Double_t(nPtBins);
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
  fitMassSigmaTotal = fitMassSigmaTotal/Double_t(nPtBins);
  fitMassSigmaErrTotal = TMath::Sqrt(fitMassSigmaErrTotal)/Double_t(nPtBins);
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
  // hCntSig2->SetMarkerStyle(29);
  // hCntSig2->SetMarkerColor(kGray+1);
  // hCntSig2->SetLineColor(kGray+1);
  // hCntSig2->Draw("PSAME");
  TLegend* leg=new TLegend(0.5576,0.7422,0.8997,0.9007);
  leg->SetFillColor(0);
  TLegendEntry* ent=leg->AddEntry(hRawyield,"From Fit","PL");
  ent->SetTextColor(hRawyield->GetMarkerColor());
  ent=leg->AddEntry(hCntSig1,"From Counting1","PL");
  ent->SetTextColor(hCntSig1->GetMarkerColor());
  // ent=leg->AddEntry(hCntSig2,"From Counting2","PL");
  // ent->SetTextColor(hCntSig2->GetMarkerColor());
  leg->Draw();
  cYield->cd(2);
  hNDiffCntSig1->SetMarkerStyle(26);
  hNDiffCntSig1->SetMarkerColor(2);
  hNDiffCntSig1->SetLineColor(2);
  hNDiffCntSig1->SetTitle("Cmp Fit-Count;#it{p}_{T} (GeV/c);(S_{fit}-S_{count})/S_{fit}");
  hNDiffCntSig1->GetYaxis()->SetTitleOffset(1.2);  
  hNDiffCntSig1->Draw("P");
  // hNDiffCntSig2->SetMarkerStyle(29);
  // hNDiffCntSig2->SetMarkerColor(kGray+1);
  // hNDiffCntSig2->SetLineColor(kGray+1);
  // hNDiffCntSig2->Draw("PSAME");
  TLegend* leg1=new TLegend(0.5878,0.7951,0.8991,0.9000);
  leg1->SetFillColor(0);
  TLegendEntry* ent1=leg1->AddEntry(hNDiffCntSig1,"From Counting1","PL");
  ent1->SetTextColor(hNDiffCntSig1->GetMarkerColor());
  // ent1=leg1->AddEntry(hNDiffCntSig2,"From Counting2","PL");
  // ent1->SetTextColor(hNDiffCntSig2->GetMarkerColor());
  leg1->Draw();
  TLine* line0 = new TLine(ptBins[0],0.,ptBins[nPtBins],0.);
  line0->SetLineColorAlpha(kGray+1,1.0); line0->SetLineWidth(1); line0->SetLineStyle(2); line0->Draw("Same");


  // Plot of Signal,Background, Signal/Background, Significance
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
  cInvariantMass->SaveAs(Form("%s/MassFit_Pt.pdf"             ,gOutputDirname.Data()));
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

  return;	
}

TTree* CreateTreeFromSideBands(TTree* inputTree,const Int_t  nsigmaLow,const Int_t  nsigmaUp){

	// Intitalize outputTree
	TTree *outputTree = new TTree("treeList_Bkg", "Candidates variables tree, Background");

	// Initialize tree variables
  	const Int_t nVar = 35;
  	TString fCandidateVariableNames[nVar];
  	Float16_t fCandidateVariables[nVar];

  	fCandidateVariableNames[0] = "massLc2K0Sp";
    fCandidateVariableNames[1] = "alphaArm";
    fCandidateVariableNames[2] = "massK0S";
    fCandidateVariableNames[3] = "massLambda";
    fCandidateVariableNames[4] = "massLambdaBar";
    fCandidateVariableNames[5] = "cosPAK0S";
    fCandidateVariableNames[6] = "dcaV0";
    fCandidateVariableNames[7] = "tImpParBach";
    fCandidateVariableNames[8] = "tImpParV0";
    fCandidateVariableNames[9] = "nSigmaTPCpr";
    fCandidateVariableNames[10] = "nSigmaTOFpr";
    fCandidateVariableNames[11] = "bachelorPt";
    fCandidateVariableNames[12] = "V0positivePt";
    fCandidateVariableNames[13] = "V0negativePt";
    fCandidateVariableNames[14] = "dcaV0pos";
    fCandidateVariableNames[15] = "dcaV0neg";
    fCandidateVariableNames[16] = "v0Pt";
    fCandidateVariableNames[17] = "bachTPCmom";
    fCandidateVariableNames[18] = "LcPt";
    fCandidateVariableNames[19] = "combinedProtonProb";
    fCandidateVariableNames[20] = "V0positiveEta";
    fCandidateVariableNames[21] = "bachelorP";
    fCandidateVariableNames[22] = "bachelorEta";
    fCandidateVariableNames[23] = "v0P";
    fCandidateVariableNames[24] = "DecayLengthK0S";
    fCandidateVariableNames[25] = "nSigmaTPCpi";
    fCandidateVariableNames[26] = "nSigmaTPCka";
    fCandidateVariableNames[27] = "NtrkRaw";
    fCandidateVariableNames[28] = "NtrkCorr";
    fCandidateVariableNames[29] = "CosThetaStar";
    fCandidateVariableNames[30] = "signd0";        
    fCandidateVariableNames[31] = "centrality"; 
    fCandidateVariableNames[32] = "NtrkAll";
    fCandidateVariableNames[33] = "origin";
    fCandidateVariableNames[34] = "ptArm";

    for(Int_t ivar=0; ivar<nVar; ivar++){
    	inputTree->SetBranchAddress(Form("%s",fCandidateVariableNames[ivar].Data()), &fCandidateVariables[ivar]);
    	outputTree->Branch(Form("%s",fCandidateVariableNames[ivar].Data()),&fCandidateVariables[ivar],Form("%s/f",fCandidateVariableNames[ivar].Data()));
    }

    // Fill histograms
    Long64_t nEntries = inputTree->GetEntries();// 100000; //inputTree->GetEntries();

    Printf("\n--- Processing: %10lld candidates",nEntries);

	for(Long64_t iLC = 0; iLC <nEntries; iLC++){
		// Get Entry
		inputTree->GetEntry(iLC);
		if (iLC % (nEntries/10) == 0) Printf("... Processing: %10lld \t %.0f%%",iLC,Double_t(iLC)/Double_t(nEntries)*100.);
		
		Float16_t massLc = (Double_t) fCandidateVariables[0];
		Float16_t pTLc = (Double_t) fCandidateVariables[18];

		for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
			// pT Bin Ranges
	    	Double_t pTBinLow = ptBins[iPtBin];
	    	Double_t pTBinUp  = ptBins[iPtBin+1];

	    	if(pTLc>pTBinLow && pTLc <=pTBinUp){

	    		Double_t nSigmaOfMean = TMath::Abs(massLc-cutMassMean[iPtBin])/cutMassSigma[iPtBin];
	    		
	    		if(nSigmaOfMean >= nsigmaLow && nSigmaOfMean <= nsigmaUp){
          // Printf("Filled");
          outputTree->Fill();
					break;
				}
	    	}
		}
	}

	return outputTree;

}

void CreateBkgFromSideBands(
	TString fileName,
	TString dir,
	const Int_t  nsigmaLow = 7,
	const Int_t  nsigmaUp = 9,
	const Bool_t useStdCuts = kTRUE,
	const Bool_t flagKeepCanvasOpen = kFALSE,
	TString outputDirName = "OutputSideBands"
	)
{

	gOutputDirname = outputDirName;

	//===========================================================
	//     Reading all information and creating output
	//===========================================================

	// Read input file
	TFile *inputFile(0);
	if (!gSystem->AccessPathName( Form("%s%s",dir.Data(),fileName.Data()) )){
		inputFile  	= new TFile( Form("%s%s",dir.Data(),fileName.Data()) );
	}
	else{
		Printf("ERROR: could not open input root file");
		exit(1);
	}

	// Get Tree from inputfile
	TTree *inputTree = (TTree*)inputFile->Get("treeList_PfCuts_Sgn");
	if(!inputTree){
		Printf("ERROR: treeList_Sgn not in inputFile");
		exit(1);
	}

	// Open output for Invariantmass file 
	gSystem->Exec(Form("rm -r %s/",outputDirName.Data()));
	gSystem->Exec(Form("mkdir %s",outputDirName.Data()));
	TFile *outputFileFitResults(0);
	outputFileFitResults = new TFile(Form("%s/FitResults.root",outputDirName.Data()),"RECREATE");
	if(!outputFileFitResults){
		Printf("ERROR: outputFileFitResults not created");
	exit(1);
	}
	outputFileFitResults->cd();

	for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++){
		// pT Bin Ranges
		Double_t pTBinLow = ptBins[iPtBin];
		Double_t pTBinUp  = ptBins[iPtBin+1];

		outputFileFitResults->cd();
		directories[iPtBin] = outputFileFitResults->mkdir(Form("pT_%.1f_%.1f",pTBinLow,pTBinUp));
	}

	TTree *outputTree;

	SetFitParameters();

	SetStandardValues();

  outputFileFitResults->cd();

	if(!useStdCuts){
		
		GenerateInvariantMassHistograms(inputTree,kTRUE,flagKeepCanvasOpen);

		MassFitter(outputFileFitResults,kTRUE);

		outputTree = CreateTreeFromSideBands(inputTree,nsigmaLow,nsigmaUp);

		Printf("\nptBin \t ptLow \t ptUp \t isFitted \t Mean \t\t Sigma");
		for(Int_t iPtBin=0;iPtBin<nPtBins;iPtBin++){
			if(flagIsFitted[iPtBin]) Printf("%s%d \t %.1f \t %.1f \t %d \t\t %.3f \t %.3f%s",GREEN,iPtBin+1,ptBins[iPtBin],ptBins[iPtBin+1],flagIsFitted[iPtBin],cutMassMean[iPtBin],cutMassSigma[iPtBin],RESET);
			else Printf("%s%d \t %.1f \t %.1f \t %d \t\t %.3f \t %.3f%s",RED,iPtBin+1,ptBins[iPtBin],ptBins[iPtBin+1],flagIsFitted[iPtBin],cutMassMean[iPtBin],cutMassSigma[iPtBin],RESET);
		}
	}
	else{
		outputTree = CreateTreeFromSideBands(inputTree,nsigmaLow,nsigmaUp);
	}

	inputFile->Close();

	outputFileFitResults->cd();
	outputTree->Write();

	if(!flagKeepCanvasOpen){
		outputFileFitResults->Close();
	}

}