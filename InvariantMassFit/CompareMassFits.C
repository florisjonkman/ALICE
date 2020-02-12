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

enum {k1999=0, k19, k1029, k3059};

Int_t nFiles;
Int_t ratioBy;
TString inputFiles[5];
TString inputTitles[5];
Int_t   colors[5] = {kBlack,kBlue+1,kRed+1,kOrange+1,kRed+1};
Int_t   markers[5]= {20,21,22,23,29};
// Int_t   colors[5] = {kBlack,kBlue+1,kCyan,kRed-9,kRed+1};
// Int_t   markers[5]= {20,22,22,23,23};

// Int_t   colors[4] = {kBlack,kBlue+1,kBlack,kBlue+1};
// Int_t   markers[4]= {20,21,24,25};

void SetInputFiles1_999(){

  nFiles = 3;

  inputFiles[0]  = "HFMassFitter_AllPt_1_999_stdCuts_FreeSigma/FitResults.root";
  inputFiles[1]  = "HFMassFitter_AllPt_1_999_stdCutsTMVA_FreeSigma/FitResults.root";
  inputFiles[2]  = "HFMassFitter_AllPt_1_999_PfCutsTMVA_FreeSigma/FitResults.root";

  inputTitles[0] = "Standard";
  inputTitles[1] = "Standard + MVA";
  inputTitles[2] = "Prefiltering + MVA";
  // inputTitles[3] = "Std. + MVA, 1 - 999, Free";

  // inputFiles[0]  = "HFMassFitter_AllPt_1_999_stdCuts_FreeSigma/FitResults.root";
  // inputFiles[1]  = "HFMassFitter_AllPt_1_999_stdCutsTMVA_FreeSigma/FitResults.root";
  // inputFiles[2]  = "HFMassFitter_AllPt_1_999_PfCutsTMVA_FreeSigma/FitResults.root";
  // inputFiles[3]  = "HFMassFitter_AllPt_1_999_stdCutsTMVA_FreeSigma/FitResults.root";

  // inputFiles[0]  = "HFMassFitter_AllPt_1_999_PfCutsTMVA_FixSigma/FitResults.root";
  // inputFiles[1]  = "HFMassFitter_AllPt_1_999_PfCutsTMVA_FixSigma_loose2/FitResults.root";
  // inputFiles[2]  = "HFMassFitter_AllPt_1_999_PfCutsTMVA_FixSigma_loose1/FitResults.root";
  // inputFiles[3]  = "HFMassFitter_AllPt_1_999_PfCutsTMVA_FixSigma_tight1/FitResults.root";
  // inputFiles[4]  = "HFMassFitter_AllPt_1_999_PfCutsTMVA_FixSigma_tight2/FitResults.root";

  // inputTitles[0] = "PfCuts. + MVA (Central), 1 - 999";
  // inputTitles[1] = "Loose2";
  // inputTitles[2] = "Loose1";
  // inputTitles[3] = "Tight1";
  // inputTitles[4] = "Tight2";

}

void SetInputFiles30_59(){

  nFiles = 3;

  // inputFiles[0]  = "IncorrectTPCCuts/HFMassFitter_AllPt_30_59_stdCuts_SigmaFix/FitResults.root";
  // inputFiles[1]  = "IncorrectTPCCuts/HFMassFitter_AllPt_30_59_stdCutsTMVA_SigmaFix/FitResults.root";
  // inputFiles[2]  = "IncorrectTPCCuts/HFMassFitter_AllPt_30_59_PfCutsTMVA_SigmaFix/FitResults.root";

  // inputTitles[0] = "Standard, 30 - 59";
  // inputTitles[1] = "Standard + MVA, 30 - 59";
  // inputTitles[2] = "Prefiltering + MVA, 30 - 59";

  inputFiles[0]  = "CorrectTPCCuts/HFMassFitter_AllPt_30_59_stdCuts_SigmaFix/FitResults.root";
  inputFiles[1]  = "CorrectTPCCuts/HFMassFitter_AllPt_30_59_stdCutsTMVA_SigmaFix/FitResults.root";
  inputFiles[2]  = "CorrectTPCCuts/HFMassFitter_AllPt_30_59_PfCutsTMVA_SigmaFix/FitResults.root";

  inputTitles[0] = "Standard, 30 - 59";
  inputTitles[1] = "Standard + MVA, 30 - 59";
  inputTitles[2] = "Prefiltering + MVA, 30 - 59";

}

void CompareMassFits(Int_t multRange, TString outputDirName = "CompareMassFits"){

	if     (multRange == k1999) SetInputFiles1_999();
	else if(multRange == k3059) SetInputFiles30_59();
	else   {Printf("Multrange not yet available"); return NULL;}

	// Bins: pT
	const Int_t nPtBins = 6;
	Double_t ptBins[nPtBins+1] =  {  1,  2,  4,  6,  8, 12, 24}; // pT bin ranges

	Double_t sigmaMC[nPtBins]    = {0.007520,0.007553,0.008798,0.010374,0.012303,0.016364};
	Double_t sigmaMCErr[nPtBins] = {0.000048,0.000034,0.000052,0.000092,0.000121,0.000244}; // SigmaErr
	Double_t meanMC[nPtBins]     = {2.285820,2.286810,2.287827,2.287972,2.288575,2.288647};
	Double_t meanMCErr[nPtBins]  = {0.000066,0.000047,0.000069,0.000124,0.000165,0.000328}; // MeanErr
	TH1D* hSigmaMC = new TH1D("hSigmaMC","hSigmaMC",nPtBins,ptBins);
	TH1D* hMeanMC  = new TH1D("hMeanMC","hMeanMC",nPtBins,ptBins); 

	for (Int_t iPt = 0; iPt < nPtBins; ++iPt){
		hSigmaMC->SetBinContent(iPt+1,sigmaMC[iPt] );
		hSigmaMC->SetBinError(iPt+1,sigmaMCErr[iPt] );
		hMeanMC ->SetBinContent(iPt+1,meanMC[iPt]  );
		hMeanMC ->SetBinError(iPt+1,meanMCErr[iPt]  );
	}
	hSigmaMC->SetLineColor(kMagenta+1);
	hSigmaMC->SetMarkerStyle(6);
	hSigmaMC->SetMarkerColor(kMagenta+1);
	hSigmaMC->SetFillColorAlpha(kMagenta+1, 0.3);
	hSigmaMC->SetLineStyle(8);
	hMeanMC->SetLineColor(kMagenta+1);
	hMeanMC->SetMarkerStyle(6);
	hMeanMC->SetMarkerColor(kMagenta+1);
	hMeanMC->SetFillColorAlpha(kMagenta+1, 0.3);
	hMeanMC->SetLineStyle(8);

	TCanvas *cMassSigma=new TCanvas("cMassSigma","cMassSigma",1200,600);
	cMassSigma->Divide(2,1);

	TCanvas *cYield=new TCanvas("cYield","cYield",1200,600);
	cYield->Divide(2,1);

	TCanvas *cSignalBackgroundSoverBSignificance=new TCanvas("cSignalBackgroundSoverBSignificance","cSignalBackgroundSoverBSignificance",1800,1200);
	cSignalBackgroundSoverBSignificance->Divide(4,2);

	TCanvas *cYieldSoverBSignificance =new TCanvas("cYieldSoverBSignificance","cYieldSoverBSignificance",1800,800);
	cYieldSoverBSignificance->Divide(3,1,0.0000001,0.0000001);
	TCanvas *cYieldSoverBSignificanceRatio =new TCanvas("cYieldSoverBSignificanceRatio","cYieldSoverBSignificanceRatio",1800,400);
	cYieldSoverBSignificanceRatio->Divide(3,1,0.0000001,0.0000001);

	TH1D* hYieldDen;
	TH1D* hSignalDen;
	TH1D* hBackgroundDen;
	TH1D* hSoverBDen;
	TH1D* hSignificanceDen;

	TH1D* hMass;
	TH1D* hSigma;
	TH1D* hYield;
	TH1D* hSignal;
	TH1D* hBackground;
	TH1D* hSoverB;
	TH1D* hSignificance;

	TH1D* hYieldRatio;
	TH1D* hSignalRatio;
	TH1D* hBackgroundRatio;
	TH1D* hSoverBRatio;
	TH1D* hSignificanceRatio;

	TLegend* legMeanSigma = new TLegend(0.128721 , 0.132013, 0.653902, 0.26356);
	TLegend* legYield = new TLegend(0.384356 ,0.681327, 0.885191, 0.883986);

	for (Int_t i = 0; i < nFiles; ++i)
	{
		TFile *inputFile(0);
		if (!gSystem->AccessPathName(inputFiles[i].Data() )){
			inputFile   = new TFile(inputFiles[i].Data());
		}
		else{
			Printf("ERROR: could not open input root file %d, %s",i,inputFiles[i].Data());
			return NULL;
		}

		Printf("i = %d",i);

		gStyle->SetOptStat(0);
		// Getting all the histograms
		hMass = (TH1D*) inputFile->Get("hMassMean");
		hMass->SetName(Form("hMass%d",i));
		hMass->SetMarkerStyle(markers[i]);
		hMass->SetMarkerColor(colors[i]);
		hMass->SetLineColor(colors[i]);
		cMassSigma->cd(1);
		if(i==0) hMass->Draw("PL");
		else	 hMass->Draw("Same");
		legMeanSigma->AddEntry(hMass,inputTitles[i],"PL");

		hSigma = (TH1D*) inputFile->Get("hSigma");
		hSigma->SetName(Form("hSigma%d",i));
		hSigma->SetMarkerStyle(markers[i]);
		hSigma->SetMarkerColor(colors[i]);
		hSigma->SetLineColor(colors[i]);
		hSigma->GetYaxis()->SetRangeUser(0.004,0.025);
		cMassSigma->cd(2);
		if(i==0) hSigma->Draw("P");
		else	 hSigma->Draw("PSAME");

		hYield = (TH1D*) inputFile->Get("hRawyield");
		hYield->SetName(Form("hYield%d",i));
		hYield->SetMarkerStyle(markers[i]);
		hYield->SetMarkerColor(colors[i]);
		hYield->SetLineColor(colors[i]);
		cYield->cd(1);
		gPad-> SetLogy(1);
		// hYield->GetYaxis()->SetRangeUser(0,2);
		if(i==0){
			hYieldDen = (TH1D*) hYield->Clone("hYieldDen");
			hYield->Draw("PL");
		}
		else	 hYield->Draw("Same");
		legYield->AddEntry(hMass,inputTitles[i],"PL");
		cYieldSoverBSignificance->cd(1);
		if(i==0){
			hYield->Draw("PL");
		}
		else hYield->Draw("Same");
		cYield->cd(2);
		hYieldRatio = (TH1D*) hYield->Clone(Form("hYieldRatio%d",i));
		// hYieldRatio->Divide(hYieldDen); // Uncorrelated
		hYieldRatio->Divide(hYieldRatio,hYieldDen,1,1,"B");	// Correlated
		hYieldRatio->GetYaxis()->SetRangeUser(0,5);
		hYieldRatio->SetTitle(Form("Yield / ( %s )",inputTitles[0].Data()));
		hYieldRatio->GetYaxis()->SetTitle(Form("Yield / ( %s )",inputTitles[0].Data()));
		if(i==0) hYieldRatio->Draw("PL");
		else	 hYieldRatio->Draw("Same");
		cYieldSoverBSignificanceRatio->cd(1);
		if(i==0){
			hYieldRatio->GetYaxis()->SetTitle("Ratio");
			hYieldRatio->SetTitle("");
			hYieldRatio->Draw("PL");
		}
		else hYieldRatio->Draw("Same");
		
		hSignal = (TH1D*) inputFile->Get("hSignal");
		hSignal->SetName(Form("hSignal%d",i));
		hSignal->SetMarkerStyle(markers[i]);
		hSignal->SetMarkerColor(colors[i]);
		hSignal->SetLineColor(colors[i]);
		cSignalBackgroundSoverBSignificance->cd(1);
		gPad-> SetLogy(1);
		hSignal->GetYaxis()->SetRangeUser(10,1e5);
		hSignal->GetYaxis()->SetTitle("Signal");
		if(i==0){
			hSignalDen = (TH1D*) hYield->Clone("hSignalDen");
			hSignal->Draw("PL");
		}
		else	 hSignal->Draw("Same");
		cSignalBackgroundSoverBSignificance->cd(1+4);
		hSignalRatio = (TH1D*) hSignal->Clone(Form("hSignalRatio%d",i));
		// hSignalRatio->Divide(hSignalDen); // Uncorrelated
		hSignalRatio->Divide(hSignalRatio,hSignalDen,1,1,"B");	// Correlated
		hSignalRatio->GetYaxis()->SetRangeUser(0,5);
		hSignalRatio->SetTitle(Form("Signal / ( %s )",inputTitles[0].Data()));
		hSignalRatio->GetYaxis()->SetTitle(Form("Signal / ( %s )",inputTitles[0].Data()));
		if(i==0) hSignalRatio->Draw("PL");
		else	 hSignalRatio->Draw("Same");

		TH1D* hBackground = (TH1D*) inputFile->Get("hBackground");
		hBackground->SetName(Form("hBackground%d",i));
		hBackground->SetMarkerStyle(markers[i]);
		hBackground->SetMarkerColor(colors[i]);
		hBackground->SetLineColor(colors[i]);
		cSignalBackgroundSoverBSignificance->cd(2);
		gPad-> SetLogy(1);
		hBackground->GetYaxis()->SetRangeUser(10,1e7);
		hBackground->GetYaxis()->SetTitle("Background");
		if(i==0){
			hBackgroundDen = (TH1D*) hBackground->Clone("hSignalDen");
			hBackground->Draw("PL");
		}
		else	 hBackground->Draw("Same");
		cSignalBackgroundSoverBSignificance->cd(2+4);
		hBackgroundRatio = (TH1D*) hBackground->Clone(Form("hBackgroundRatio%d",i));
		// hBackgroundRatio->Divide(hBackgroundDen); // Uncorrelated
		hBackgroundRatio->Divide(hBackgroundRatio,hBackgroundDen,1,1,"B");	// Correlated
		hBackgroundRatio->GetYaxis()->SetRangeUser(0,5);
		hBackgroundRatio->SetTitle(Form("Background / ( %s )",inputTitles[0].Data()));
		hBackgroundRatio->GetYaxis()->SetTitle(Form("Background / ( %s )",inputTitles[0].Data()));
		if(i==0) hBackgroundRatio->Draw("PL");
		else	 hBackgroundRatio->Draw("Same");

		hSoverB = (TH1D*) inputFile->Get("hSignalBackground");
		hSoverB->SetName(Form("hSoverB%d",i));
		hSoverB->SetMarkerStyle(markers[i]);
		hSoverB->SetMarkerColor(colors[i]);
		hSoverB->SetLineColor(colors[i]);
		hSoverB->GetYaxis()->SetRangeUser(0,0.25);
		cSignalBackgroundSoverBSignificance->cd(3);
		if(i==0){
			hSoverBDen = (TH1D*) hSoverB->Clone("hSoverBDen");
			hSoverB->Draw("PL");
		}
		else	 hSoverB->Draw("Same");
		cYieldSoverBSignificance->cd(2);
		if(i==0){
			hSoverB->Draw("PL");
		}
		else hSoverB->Draw("Same");
		cSignalBackgroundSoverBSignificance->cd(3+4);
		hSoverBRatio = (TH1D*) hSoverB->Clone(Form("hSoverBRatio%d",i));
		// hSoverBRatio->Divide(hSoverBDen); // Uncorrelated
		hSoverBRatio->Divide(hSoverBRatio,hSoverBDen,1,1,"B");	// Correlated
		hSoverBRatio->GetYaxis()->SetRangeUser(0,5);
		hSoverBRatio->SetTitle(Form("(S/B) / ( %s )",inputTitles[0].Data()));
		hSoverBRatio->GetYaxis()->SetTitle(Form("(S/B) / ( %s )",inputTitles[0].Data()));
		if(i==0) hSoverBRatio->Draw("PL");
		else	 hSoverBRatio->Draw("Same");
		cYieldSoverBSignificanceRatio->cd(2);
		if(i==0){
			hSoverBRatio->GetYaxis()->SetTitle("Ratio");
			hSoverBRatio->SetTitle("");
			hSoverBRatio->Draw("PL");
		}
		else hSoverBRatio->Draw("Same");

		hSignificance = (TH1D*) inputFile->Get("hSignificance");
		hSignificance->SetName(Form("hSignificance%d",i));
		hSignificance->SetMarkerStyle(markers[i]);
		hSignificance->SetMarkerColor(colors[i]);
		hSignificance->SetLineColor(colors[i]);
		hSignificance->GetYaxis()->SetRangeUser(0,17);
		cSignalBackgroundSoverBSignificance->cd(4);
		if(i==0){
			hSignificanceDen = (TH1D*) hSignificance->Clone("hSignificanceDen");
			hSignificance->Draw("PL");
		}
		else	 hSignificance->Draw("Same");
		cYieldSoverBSignificance->cd(3);
		if(i==0){
			hSignificance->Draw("PL");
		}
		else hSignificance->Draw("Same");
		cSignalBackgroundSoverBSignificance->cd(4+4);
		hSignificanceRatio = (TH1D*) hSignificance->Clone(Form("hSignificanceRatio%d",i));
		// hSignificanceRatio->Divide(hSignificanceDen); // Uncorrelated
		hSignificanceRatio->Divide(hSignificanceRatio,hSignificanceDen,1,1,"B");	// Correlated
		hSignificanceRatio->GetYaxis()->SetRangeUser(0,5);
		hSignificanceRatio->SetTitle(Form("Significance / ( %s )",inputTitles[0].Data()));
		hSignificanceRatio->GetYaxis()->SetTitle(Form("Significance / ( %s )",inputTitles[0].Data()));
		if(i==0) hSignificanceRatio->Draw("PL");
		else	 hSignificanceRatio->Draw("Same");
		cYieldSoverBSignificanceRatio->cd(3);
		if(i==0){
			hSignificanceRatio->GetYaxis()->SetTitle("Ratio");
			hSignificanceRatio->SetTitle("");
			hSignificanceRatio->Draw("PL");
		}
		else hSignificanceRatio->Draw("Same");

	}
	cMassSigma->cd(1);
	hMeanMC->Draw("E4Same");
	cMassSigma->cd(2);
	hSigmaMC->Draw("E4Same");

	// legend
	legMeanSigma->AddEntry(hMeanMC,"Monte Carlo","L");
	cMassSigma->cd(1);
	legMeanSigma->Draw();
	cYield->cd(1);
	legYield->Draw();
	cSignalBackgroundSoverBSignificance->cd(1);
	legYield->Draw();
	cYieldSoverBSignificance->cd(1);
	legYield->Draw();

	// Open output file 
	gSystem->Exec(Form("rm -r %s/",outputDirName.Data()));
	gSystem->Exec(Form("mkdir %s",outputDirName.Data()));

	cMassSigma   				->SaveAs(Form("%s/MassSigma.pdf"   			 	,outputDirName.Data()));
	cYield       			  	->SaveAs(Form("%s/Yield.pdf"    				,outputDirName.Data()));
	cSignalBackgroundSoverBSignificance 	->SaveAs(Form("%s/SignalBackgroundSoverB.pdf"   ,outputDirName.Data()));
	cYieldSoverBSignificance    ->SaveAs(Form("%s/YieldSoverBSignificance.pdf"    				,outputDirName.Data()));
	cYieldSoverBSignificanceRatio->SaveAs(Form("%s/YieldSoverBSignificanceRatio.pdf"    				,outputDirName.Data()));
	// For thesis

	TCanvas *cThesis1   =new TCanvas("cThesis1","cThesis1",1200,400);
	TCanvas *cThesis2   =new TCanvas("cThesis2","cThesis2",800,800);
	
	cThesis1->Divide(3,1,0.0000001,0.0000001);
	cThesis2->Divide(2,2,0.0000001,0.0000001);
	
	TLegend *legThesis1 = new TLegend(0.406511 ,0.679144, 0.899833,  0.90107);
	TLegend *legThesis2 = new TLegend(0.408521 ,0.697674, 0.902256, 0.899225);

	for (Int_t i = 0; i < nFiles; ++i)
	{
		TFile *inputFile1(0);
		if (!gSystem->AccessPathName(inputFiles[i].Data() )){
			inputFile1   = new TFile(inputFiles[i].Data());
		}
		else{
			Printf("ERROR: could not open input root file %d, %s",i,inputFiles[i].Data());
			return NULL;
		}

		gStyle->SetPadBorderMode(0);
		gStyle->SetOptStat(0);
		// Getting all the histograms
		hMass = (TH1D*) inputFile1->Get("hMassMean");
		hMass->SetName(Form("hMass%d",i));
		hMass->SetMarkerStyle(markers[i]);
		hMass->SetMarkerColor(colors[i]);
		hMass->SetLineColor(colors[i]);
		hMass->SetTitle("Mean");
		hMass->GetYaxis()->SetTitle("Mean (GeV/#it{c}^{2})");
		hMass->GetYaxis()->SetLabelSize(0.03);
		hMass->GetYaxis()->SetTitleOffset(1.52);
		hMass->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
		cThesis1->cd(1);
		if(i==0) hMass->Draw("PL");
		else	 hMass->Draw("Same");
		legThesis1->AddEntry(hMass,inputTitles[i],"PL");
		legThesis2->AddEntry(hMass,inputTitles[i],"PL");

		hSigma = (TH1D*) inputFile1->Get("hSigma");
		hSigma->SetName(Form("hSigma%d",i));
		hSigma->SetMarkerStyle(markers[i]);
		hSigma->SetMarkerColor(colors[i]);
		hSigma->SetLineColor(colors[i]);
		hSigma->SetTitle("Sigma (Unconstrained)");
		hSigma->GetYaxis()->SetTitle("Sigma (GeV/#it{c}^{2})");
		hSigma->GetYaxis()->SetLabelSize(0.03);
		hSigma->GetYaxis()->SetTitleOffset(1.52);
		hSigma->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
		hSigma->GetYaxis()->SetRangeUser(0.005,0.018);
		cThesis1->cd(2);
		if(i==0) hSigma->Draw("PL");
		else	 hSigma->Draw("PSame");

		hYield = (TH1D*) inputFile1->Get("hRawyield");
		hYield->SetName(Form("hYield%d",i));
		hYield->SetMarkerStyle(markers[i]);
		hYield->SetMarkerColor(colors[i]);
		hYield->SetLineColor(colors[i]);
		hYield->SetTitle("Raw yield");
		hYield->GetYaxis()->SetTitle("Raw yield");
		hYield->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
		cThesis1->cd(3);
		gPad-> SetLogy(1);
		if(i==0){
			hYield->Draw("PL");
		}
		else	 hYield->Draw("Same");


		cThesis2->cd(1);
		gPad-> SetLogy(0);
		hSoverB = (TH1D*) inputFile1->Get("hSignalBackground");
		hSoverB->SetName(Form("hSoverB%d",i));
		hSoverB->SetMarkerStyle(markers[i]);
		hSoverB->SetMarkerColor(colors[i]);
		hSoverB->SetLineColor(colors[i]);
		hSoverB->GetYaxis()->SetRangeUser(0,0.16);
		hSoverB->GetYaxis()->SetTitle("S/B");
		hSoverB->SetTitle("Signal/Background");
		hSoverB->GetYaxis()->SetTitleOffset(1.45);
		// hSoverB->SetTitle("S/B");
		if(i==0){
			hSoverBDen = (TH1D*) hSoverB->Clone("hSoverBDen");
			hSoverB->Draw("PL");
		}
		else	 hSoverB->Draw("Same");
		cYieldSoverBSignificance->cd(2);
		if(i==0){
			hSoverB->Draw("PL");
		}
		else hSoverB->Draw("Same");
		cThesis2->cd(3);
		hSoverBRatio = (TH1D*) hSoverB->Clone(Form("hSoverBRatio2_%d",i));
		// hSoverBRatio->Divide(hSoverBDen); // Uncorrelated
		hSoverBRatio->Divide(hSoverBRatio,hSoverBDen,1,1,"B");	// Correlated
		hSoverBRatio->GetYaxis()->SetRangeUser(0.5,5.5);
		hSoverBRatio->SetTitle("Ratio Signal/Background");
		hSoverBRatio->GetYaxis()->SetTitle("(S/B) / Standard");
		hSoverBRatio->GetYaxis()->SetTitleOffset(1.45);
		if(i==0) hSoverBRatio->Draw("PL");
		else	 hSoverBRatio->Draw("Same");

		hSignificance = (TH1D*) inputFile1->Get("hSignificance");
		hSignificance->SetName(Form("hSignificance%d",i));
		hSignificance->SetMarkerStyle(markers[i]);
		hSignificance->SetMarkerColor(colors[i]);
		hSignificance->SetLineColor(colors[i]);
		hSignificance->GetYaxis()->SetTitle("Significance");
		hSignificance->SetTitle("Significance");
		hSignificance->GetYaxis()->SetRangeUser(0,16);
		cThesis2->cd(2);
		if(i==0){
			hSignificanceDen = (TH1D*) hSignificance->Clone("hSignificanceDen2");
			hSignificance->Draw("PL");
		}
		else	 hSignificance->Draw("Same");
		cThesis2->cd(4);
		hSignificanceRatio = (TH1D*) hSignificance->Clone(Form("hSignificanceRatio2_%d",i));
		// hSignificanceRatio->Divide(hSignificanceDen); // Uncorrelated
		hSignificanceRatio->Divide(hSignificanceRatio,hSignificanceDen,1,1,"B");	// Correlated
		hSignificanceRatio->GetYaxis()->SetRangeUser(0.8,2.8);
		hSignificanceRatio->SetTitle("Ratio Significance");
		hSignificanceRatio->GetYaxis()->SetTitle("Significance / Standard");
		if(i==0)hSignificanceRatio->Draw("PL");
		else	 hSignificanceRatio->Draw("Same");
		
	}

	cThesis1->cd(1);
	hMeanMC->Draw("E2Same");
	((TH1D*)hMeanMC->Clone("hMeanMCClone"))->Draw("LSame");
	cThesis1->cd(2);
	((TH1D*)hSigmaMC->Clone("hSigmaMCClone"))->Draw("LSame");
	hSigmaMC->Draw("E2LSame");

	// legend
	cThesis1->cd(3);
	legThesis1->AddEntry(hMeanMC,"Monte Carlo","FL");
	legThesis1->Draw();
	cThesis2->cd(2);
	legThesis2->Draw();

	cThesis1   				->SaveAs(Form("%s/cThesis1.pdf"   			 	,outputDirName.Data()));
	cThesis2   				->SaveAs(Form("%s/cThesis2.pdf"   			 	,outputDirName.Data()));

}