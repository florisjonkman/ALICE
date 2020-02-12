
void CorrectedYield(){
	
	// Parameters
	//--------------------------------------------------

	TString crossSectionFile = "../CrossSection/OutputHFPtSpectrum_PfCutsTMVA_FixSigma_1_999.root";
	TString systemicsFile    = "../CrossSection/ResultsSystUncertainties.root";
	TString legendTitle      = "Data, Mult. integrated";
	TString objectTitle      = "hCorrYield_1999";

	Double_t multiplicityCorr = 1.;
	Double_t sigma      	  = 57.8*1e9;  // sigma[pb] // σpp,MB=57.8 mb
	Double_t BR         	  = 0.01010;
	Double_t triggerEff       = 0.92;      // 0.92 for Min Bias

	//--------------------------------------------------

	// Bins: pT
	const Int_t nPtBins = 6;
	Double_t ptBins[nPtBins+1] =  {0, 1, 2,  4,  6,  8, 12}; // pT bin ranges

	// Input
	TFile* inputCrossSec = TFile::Open(crossSectionFile.Data());
	if(!inputCrossSec) {Printf("\'inputCrossSec\' not found!"); return;}

	TFile* inputSyst = TFile::Open(systemicsFile.Data());
	if(!inputSyst) {Printf("\'inputSyst\' not found!"); return;}

	TFile* inputLcpKpi = TFile::Open("FilesCristina/LcpKpiCorrectedYieldPerEvent_MBvspt_ntrkl_1999_19_1029_3059.root");
	if(!inputLcpKpi) {Printf("\'inputLcpKpi\' not found!"); return;}

	TFile* inputD0 = TFile::Open("FilesCristina/finalcrossD0ppMBvspt_ntrklmult0.root");
	if(!inputD0) {Printf("\'inputD0\' not found!"); return;}

	TFile* inputD0Lc = TFile::Open("FilesCristina/LcpKpiOverD0_MBvspt_ntrkl_1999_19_1029_3059.root");
	if(!inputD0Lc) {Printf("\'inputD0Lc\' not found!"); return;}

	TFile* inputSimulation = TFile::Open("SimResults.root");
	if(!inputSimulation) {Printf("\'inputSimulation\' not found!"); return;}

	TFile* inputSimulationLcD0 = TFile::Open("FilesCristina/cLambdaCoverDModes.root");
	if(!inputSimulationLcD0) {Printf("\'inputSimulationLcD0\' not found!"); return;}

	// Open output file 
	TString outputDirName = "CorrectedYield";
	gSystem->Exec(Form("rm -r %s/",outputDirName.Data()));
	gSystem->Exec(Form("mkdir %s",outputDirName.Data()));

	// Get cross-sections and systematics
	TH1D* hSigmaCorr 	  		   = (TH1D*) inputCrossSec->Get("histoSigmaCorr");
	AliHFSystErr* systLc	       = (AliHFSystErr*) inputCrossSec->Get("AliHFSystErr");
	TGraphAsymmErrors* hSystTot    = (TGraphAsymmErrors*) inputSyst->Get("gTotErr");
	TGraphAsymmErrors* hSystNoFD   = (TGraphAsymmErrors*) inputSyst->Get("gTotErrNoFD");
	TGraphAsymmErrors* hSystFD     = (TGraphAsymmErrors*) inputSyst->Get("grErrFeeddown");
	TGraphAsymmErrors* hSystCutVar = (TGraphAsymmErrors*) inputSyst->Get("gCutVariation");

	TH1D* hYieldLcpKpi = (TH1D*) inputLcpKpi->Get("histoSigmaCorr_0");
	TGraphAsymmErrors* hSystTotLcpKpi = (TGraphAsymmErrors*) inputLcpKpi->Get("gr_TotSyst_0");
	TGraphAsymmErrors* hSystFDLcpKpi = (TGraphAsymmErrors*) inputLcpKpi->Get("gr_FDSyst_0");
	TGraphAsymmErrors* hSystNoFDLcpKpi = (TGraphAsymmErrors*) inputLcpKpi->Get("gr_TotSyst_woFD_0");

	TH1D* hLcpKpiLcD0 = (TH1D*) inputD0Lc->Get("histoSigmaCorr_Lc0_ratio");
	TGraphAsymmErrors* hSystLcpKpiLcD0_tot  = (TGraphAsymmErrors*) inputD0Lc->Get("gr_TotSyst_0");
	TGraphAsymmErrors* hSystLcpKpiLcD0_fd   = (TGraphAsymmErrors*) inputD0Lc->Get("gr_FDSyst_0");
	TGraphAsymmErrors* hSystLcpKpiLcD0_nofd = (TGraphAsymmErrors*) inputD0Lc->Get("gr_TotSyst_woFD_0");

	AliHFSystErr* systD0 = (AliHFSystErr*) inputD0->Get("AliHFSystErr");
	TH1D* hSigmaCorrD0   = (TH1D*) inputD0->Get("histoSigmaCorr");
	TGraphAsymmErrors* hSystFDD0   	 = (TGraphAsymmErrors*) inputD0->Get("gFcConservative");

	TH1D* hYieldLc_HardCR   = (TH1D*) inputSimulation->Get("hLambdac_0_999_HardCR");
	TH1D* hYieldLc_HardNoCR = (TH1D*) inputSimulation->Get("hLambdac_0_999_HardNoCR");
	TH1D* hYieldLc_SoftCR   = (TH1D*) inputSimulation->Get("hLambdac_0_999_SoftCR");
	TH1D* hYieldLc_SoftNoCR = (TH1D*) inputSimulation->Get("hLambdac_0_999_SoftNoCR");

	TH1D* hYieldD0_HardCR   = (TH1D*) inputSimulation->Get("hDZero_0_999_HardCR");
	TH1D* hYieldD0_HardNoCR = (TH1D*) inputSimulation->Get("hDZero_0_999_HardNoCR");
	TH1D* hYieldD0_SoftCR   = (TH1D*) inputSimulation->Get("hDZero_0_999_SoftCR");
	TH1D* hYieldD0_SoftNoCR = (TH1D*) inputSimulation->Get("hDZero_0_999_SoftNoCR");

	TH1D* hYieldLc_HardCR_Merged   = new TH1D("hYieldLc_HardCR_Merged","hYieldLc_HardCR_Merged",nPtBins,ptBins);
	TH1D* hYieldLc_HardNoCR_Merged = new TH1D("hYieldLc_HardNoCR_Merged","hYieldLc_HardNoCR_Merged",nPtBins,ptBins);
	TH1D* hYieldLc_SoftCR_Merged   = new TH1D("hYieldLc_SoftCR_Merged","hYieldLc_SoftCR_Merged",nPtBins,ptBins);
	TH1D* hYieldLc_SoftNoCR_Merged = new TH1D("hYieldLc_SoftNoCR_Merged","hYieldLc_SoftNoCR_Merged",nPtBins,ptBins);

	TH1D* hYieldD0_HardCR_Merged   = new TH1D("hYieldD0_HardCR_Merged","hYieldD0_HardCR_Merged",nPtBins,ptBins);
	TH1D* hYieldD0_HardNoCR_Merged = new TH1D("hYieldD0_HardNoCR_Merged","hYieldD0_HardNoCR_Merged",nPtBins,ptBins);
	TH1D* hYieldD0_SoftCR_Merged   = new TH1D("hYieldD0_SoftCR_Merged","hYieldD0_SoftCR_Merged",nPtBins,ptBins);
	TH1D* hYieldD0_SoftNoCR_Merged = new TH1D("hYieldD0_SoftNoCR_Merged","hYieldD0_SoftNoCR_Merged",nPtBins,ptBins);

	TCanvas* canvasSimulations = (TCanvas*) inputSimulationLcD0->Get("cLambdaCoverDModes");

	Double_t nAdditions[nPtBins] = {0.0};
	for (Int_t i = 1; i <= hYieldLc_HardCR->GetNbinsX(); ++i)
	{
		Double_t centre = hYieldLc_HardCR->GetBinCenter(i);	
		Int_t bin = hYieldLc_HardCR_Merged->FindBin(centre);

		hYieldLc_HardCR_Merged  ->SetBinContent(bin,hYieldLc_HardCR_Merged  ->GetBinContent(bin)+hYieldLc_HardCR  ->GetBinContent(i));
		hYieldLc_HardNoCR_Merged->SetBinContent(bin,hYieldLc_HardNoCR_Merged->GetBinContent(bin)+hYieldLc_HardNoCR->GetBinContent(i));
		hYieldLc_SoftCR_Merged  ->SetBinContent(bin,hYieldLc_SoftCR_Merged  ->GetBinContent(bin)+hYieldLc_SoftCR  ->GetBinContent(i));
		hYieldLc_SoftNoCR_Merged->SetBinContent(bin,hYieldLc_SoftNoCR_Merged->GetBinContent(bin)+hYieldLc_SoftNoCR->GetBinContent(i));
	
		for(Int_t iPt = 0; iPt < nPtBins; ++iPt){
			if( centre > ptBins[iPt] &&  centre < ptBins[iPt+1]){
				nAdditions[iPt]+=1.;
				break;
			}
		}
	}
	for(Int_t iPt = 0; iPt < nPtBins; ++iPt){
		hYieldLc_HardCR_Merged  ->SetBinContent(iPt+1,hYieldLc_HardCR_Merged  ->GetBinContent(iPt+1)/nAdditions[iPt]);
		hYieldLc_HardNoCR_Merged->SetBinContent(iPt+1,hYieldLc_HardNoCR_Merged->GetBinContent(iPt+1)/nAdditions[iPt]);
		hYieldLc_SoftCR_Merged  ->SetBinContent(iPt+1,hYieldLc_SoftCR_Merged  ->GetBinContent(iPt+1)/nAdditions[iPt]);
		hYieldLc_SoftNoCR_Merged->SetBinContent(iPt+1,hYieldLc_SoftNoCR_Merged->GetBinContent(iPt+1)/nAdditions[iPt]);
		nAdditions[iPt] = 0;
	}

	for (Int_t i = 1; i <= hYieldD0_HardCR->GetNbinsX(); ++i)
	{
		Double_t centre = hYieldD0_HardCR->GetBinCenter(i);	
		Int_t bin = hYieldD0_HardCR_Merged->FindBin(centre);

		hYieldD0_HardCR_Merged  ->SetBinContent(bin,hYieldD0_HardCR_Merged  ->GetBinContent(bin)+hYieldD0_HardCR  ->GetBinContent(i));
		hYieldD0_HardNoCR_Merged->SetBinContent(bin,hYieldD0_HardNoCR_Merged->GetBinContent(bin)+hYieldD0_HardNoCR->GetBinContent(i));
		hYieldD0_SoftCR_Merged  ->SetBinContent(bin,hYieldD0_SoftCR_Merged  ->GetBinContent(bin)+hYieldD0_SoftCR  ->GetBinContent(i));
		hYieldD0_SoftNoCR_Merged->SetBinContent(bin,hYieldD0_SoftNoCR_Merged->GetBinContent(bin)+hYieldD0_SoftNoCR->GetBinContent(i));
	
		for(Int_t iPt = 0; iPt < nPtBins; ++iPt){
			if( centre > ptBins[iPt] &&  centre < ptBins[iPt+1]){
				nAdditions[iPt]+=1.;
				break;
			}
		}
	}
	for(Int_t iPt = 0; iPt < nPtBins; ++iPt){
		hYieldD0_HardCR_Merged  ->SetBinContent(iPt+1,hYieldD0_HardCR_Merged  ->GetBinContent(iPt+1)/nAdditions[iPt]);
		hYieldD0_HardNoCR_Merged->SetBinContent(iPt+1,hYieldD0_HardNoCR_Merged->GetBinContent(iPt+1)/nAdditions[iPt]);
		hYieldD0_SoftCR_Merged  ->SetBinContent(iPt+1,hYieldD0_SoftCR_Merged  ->GetBinContent(iPt+1)/nAdditions[iPt]);
		hYieldD0_SoftNoCR_Merged->SetBinContent(iPt+1,hYieldD0_SoftNoCR_Merged->GetBinContent(iPt+1)/nAdditions[iPt]);
	}

	//--------------------------------------------------------------------------------
	// First graph
	// My results

	// Corrected yield per event
	TCanvas* cCorrYield = new TCanvas("cCorrYield","cCorrYield",600,600);  
	cCorrYield->cd();
	gPad->SetLogy(1);
	gStyle->SetOptStat(0);
	gROOT->ForceStyle();

	// Legend
	TLegend *legend=new TLegend(0.560201 ,0.745645,0.88796,0.897213);
	legend->SetBorderSize(0);
	legend->SetFillStyle(0);
	legend->SetTextFont(42);
	// legend->SetTextAlign(13);
	legend->SetTextAlign(12);
	legend->SetTextColor(kBlack);

	// Legend
	TLegend *legendSyst=new TLegend(0.560201  ,0.688153  ,0.88796 ,0.801394);
	legendSyst->SetBorderSize(0);
	legendSyst->SetFillStyle(0);
	legendSyst->SetTextFont(42);
	// legend->SetTextAlign(13);
	legendSyst->SetTextAlign(12);
	legendSyst->SetTextColor(kBlack);

	// My results, LcK0Sp
	TH1D* hBasis = new TH1D("hBasis","hBasis",nPtBins,ptBins);
	TH1D* hCorrYield = (TH1D*) hSigmaCorr->Clone("hCorrYield");
	Double_t correction = 1./(sigma*BR*triggerEff*multiplicityCorr);
	hCorrYield->Scale(correction);
	hCorrYield->SetLineColor(kBlack);
	hCorrYield->SetMarkerStyle(21);
	hCorrYield->SetMarkerSize(0.8);
	hCorrYield->SetMarkerColor(kBlack);
	
	// hCorrYield->SetMinimum(0);
	// hCorrYield->SetMaximum(12);
	legend->AddEntry(hCorrYield,legendTitle.Data(),"LP");

	// My results, systematics
	TGraphAsymmErrors* hSystTotScaledNoFD = new TGraphAsymmErrors(0); 
	hSystTotScaledNoFD->SetNameTitle("hSystTotScaledNoFD","hSystTotScaledNoFD");
	TGraphAsymmErrors* hSystScaledFD = new TGraphAsymmErrors(0); 
	hSystScaledFD->SetNameTitle("hSystScaledFD","hSystScaledFD");
	
	Int_t nentriesNoFD = hSystNoFD->GetN();
	Int_t nentriesFD = hSystFD->GetN();
	Int_t nbins = hCorrYield->GetNbinsX();

	Printf("");
	Printf("%4.4s  %4.4s  %4.4s  %4.4s  %10.10s  %10.10s  %10.10s  %10.10s  %10.10s","i","j","k","pt","yield","errlcomb","errhcomb","errl_fd","errh_fd");
	for(Int_t i=1;i<=nbins;i++) {
    	Double_t pt = hCorrYield->GetBinCenter(i);
    	Double_t ptwidth = hCorrYield->GetBinWidth(i);
    	Double_t yield =  hCorrYield->GetBinContent(i);

    	Double_t errylcomb=0., erryhcomb=0, erryl_fd=0,erryh_fd=0, x=0, y=0, errx=0,x2=0, y2=0, errx2=0;

	 	for(Int_t j=0; j<nentriesNoFD; j++) {
	 		hSystNoFD->GetPoint(j,x,y);
	 		errx = hSystNoFD->GetErrorXlow(j);

	 		for(Int_t k=0; k<nentriesFD; k++) {
	 			hSystFD->GetPoint(k,x2,y2);
	 			errx2 = hSystFD->GetErrorXlow(k);

	 			if ( ( ( (x-errx) <= pt) && ( (x+errx) >= pt) ) && ( ( (x2-errx2) <= pt) && ( (x2+errx2) >= pt) ) ){
	 				errylcomb = hSystNoFD->GetErrorYlow(j) * yield ;
				  	erryhcomb = hSystNoFD->GetErrorYhigh(j) * yield ;
				  	hSystTotScaledNoFD->SetPoint(j,x,yield);
				  	hSystTotScaledNoFD->SetPointError(j,errx,errx,errylcomb,erryhcomb);

				  	erryl_fd = hSystFD->GetErrorYlow(k) * yield ;
				  	erryh_fd = hSystFD->GetErrorYhigh(k) * yield ;
				  	hSystScaledFD->SetPoint(j,x,yield);
				  	hSystScaledFD->SetPointError(j,errx,errx,erryl_fd,erryh_fd);
				  	Printf("%4.1d  %4.1d  %4.1d  %4.1f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f",i,j,k,pt,yield,errylcomb,erryhcomb,erryl_fd,erryh_fd);
	 			}
	 		}
	 	}
	}
	Printf("");

	hBasis->SetTitle("");
	hBasis->SetStats(kFALSE);
	hBasis->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	hBasis->GetYaxis()->SetTitle("1/#it{N}_{ev} d#it{N}/d#it{p}_{T}#lbar_{|#it{y}|<0.5} (GeV^{-1}/#it{c})");
	hBasis->GetYaxis()->SetTitleOffset(1.4);
	hBasis->GetYaxis()->SetRangeUser(0.000002,0.1);
	hBasis->GetYaxis()->SetLabelSize(0.028);
	hBasis->Draw("");

	hSystScaledFD->SetLineWidth(1.);
	hSystScaledFD->SetMarkerStyle(22);
	// hSystScaledFD->SetFillStyle(3001);
	hSystScaledFD->SetFillColorAlpha(kBlack,0.4);
	hSystScaledFD->SetLineColorAlpha(kBlack,0.5);
	hSystScaledFD->Draw("E5Same");

	hSystTotScaledNoFD->SetLineColor(kBlack);
	hSystTotScaledNoFD->SetFillStyle(0);
	hSystTotScaledNoFD->Draw("E5Same");

	hCorrYield->Draw("Same");

	// Some text info
	TPaveText *prodInfo=new TPaveText(0.123746,0.792683,0.43311,0.890244,"NDC");
	prodInfo->SetBorderSize(0);
	prodInfo->SetFillStyle(0);
	prodInfo->AddText(Form("pp, #sqrt{s} = 13 TeV"));
	prodInfo->AddText(Form("Prompt #Lambda_{c}^{+}, |#it{y}| < 0.5"));
	prodInfo->SetTextFont(42);
	prodInfo->SetTextAlign(12);
	prodInfo->SetTextColor(kBlack);
	prodInfo->Draw();

	TPaveText *multInfo=new TPaveText(0.560201,0.848432,0.884615,0.912892,"NDC");
	multInfo->SetBorderSize(0);
	multInfo->SetFillStyle(0);
	multInfo->AddText(Form("Multiplicity: |#it{#eta}| < 1.0           "));
	multInfo->SetTextFont(42);
	multInfo->SetTextAlign(13);
	multInfo->SetTextColor(kBlack);
	multInfo->Draw();

	TPaveText *multLumBR=new TPaveText(0.0936455 ,0.102787,0.665552,  0.167247,"NDC");
	multLumBR->SetBorderSize(0);
	multLumBR->SetFillStyle(0);
	multLumBR->AddText(Form("#pm 5.0 %% lumi, #pm 5.2 %% BR uncertainty not shown"));
	multLumBR->SetTextFont(42);
	multLumBR->SetTextAlign(12);
	multLumBR->SetTextColor(kBlack);
	multLumBR->Draw();

	legendSyst->AddEntry(hSystTotScaledNoFD,"Syst. from data","f");
	legendSyst->AddEntry(hSystScaledFD,"Syst. from B feed-down","f");
	legendSyst->Draw();

	legend->Draw();

	cCorrYield->SaveAs(Form("%s/CorrectedYield.pdf",outputDirName.Data()));
	cCorrYield->Close();

	//--------------------------------------------------------------------------------
	// Second graph	
	// LcKppi, with ratio plot

	// Legend
	TLegend *legend2=new TLegend(0.560201 ,0.707532,0.88796 ,0.841747 );
	legend2->SetBorderSize(0);
	legend2->SetFillStyle(0);
	legend2->SetTextFont(42);
	// legend->SetTextAlign(13);
	legend2->SetTextAlign(12);
	legend2->SetTextColor(kBlack);
	legend2->AddEntry(hCorrYield,"#Lambda_{c}^{+} #rightarrow pK_{S}^{0}             ","LP");
	legend2->AddEntry(hYieldLcpKpi,"#Lambda_{c}^{+} #rightarrow pK#pi              ","LP");

	// Legend
	TLegend *legendSyst2=new TLegend(0.560201  ,0.589343   ,0.88796 ,0.699519 );
	legendSyst2->SetBorderSize(0);
	legendSyst2->SetFillStyle(0);
	legendSyst2->SetTextFont(42);
	// legend->SetTextAlign(13);
	legendSyst2->SetTextAlign(12);
	legendSyst2->SetTextColor(kBlack);
	legendSyst2->AddEntry(hSystTotScaledNoFD,"Syst. from data","f");
	legendSyst2->AddEntry(hSystScaledFD,"Syst. from B feed-down","f");

	// Define the Canvas
	TCanvas *cCorrYieldLcpKpi   =new TCanvas("cCorrYieldLcpKpi","cCorrYieldLcpKpi",600,650);
	gPad->SetLogy(0);

	// Upper plot will be in pad1
	TPad *pad1 = new TPad("pad1", "pad1", 0, 0.2, 1, 1.0);
	pad1->SetBottomMargin(0); // Upper and lower plot are joined
	pad1->SetLogy(1);
	// pad1->SetGridx();         // Vertical grid
	pad1->Draw();             // Draw the upper pad: pad1
	pad1->cd();               // pad1 becomes the current pad
	TH1D* hBasisClone = (TH1D*) hBasis->Clone("hBasisClone");
	hBasisClone->GetYaxis()->SetRangeUser(0.000002,0.1);
	hBasisClone->GetYaxis()->SetTitleSize(0.04);
	hBasisClone->GetYaxis()->SetTitleOffset(1.21);
	hBasisClone->GetYaxis()->SetLabelSize(0.032);
	hBasisClone->Draw();

	hSystScaledFD->Draw("E5Same");
	hSystTotScaledNoFD->Draw("E5Same");
	hCorrYield->Draw("Same");

	hSystFDLcpKpi->SetLineWidth(1);
	hSystFDLcpKpi->SetMarkerStyle(22);
	// hSystFDLcpKpi->SetFillStyle(3001);
	hSystFDLcpKpi->SetFillColorAlpha(kRed-7,0.4);
	hSystFDLcpKpi->SetLineColorAlpha(kRed-7,0.5);
	hSystFDLcpKpi->Draw("E5Same");

	hSystNoFDLcpKpi->SetLineColor(kRed-7);
	hSystNoFDLcpKpi->SetFillStyle(0);
	hSystNoFDLcpKpi->Draw("E5Same");

	hYieldLcpKpi->SetLineColor(kRed-7);
	hYieldLcpKpi->SetMarkerStyle(20);
	hYieldLcpKpi->SetMarkerColor(kRed-7);
	hYieldLcpKpi->Draw("Same");

	legend2->Draw();
	legendSyst2->Draw();
	prodInfo->Draw();

	TPaveText *multLumBR2=new TPaveText(0.09699 ,0.0264423,0.643813, 0.0905449,"NDC");
	multLumBR2->SetBorderSize(0);
	multLumBR2->SetFillStyle(0);
	multLumBR2->AddText(Form("#pm 5.0 %% lumi, BR uncertainty not shown"));
	multLumBR2->SetTextFont(42);
	multLumBR2->SetTextAlign(13);
	multLumBR2->SetTextColor(kBlack);
	multLumBR2->Draw();

	TPaveText *multInfo2=new TPaveText(0.550167,0.845753, 0.879599,0.909856,"NDC");
	multInfo2->SetBorderSize(0);
	multInfo2->SetFillStyle(0);
	multInfo2->AddText(Form("Multiplicity: |#it{#eta}| < 1.0         "));
	multInfo2->SetTextFont(42);
	multInfo2->SetTextAlign(13);
	multInfo2->SetTextColor(kBlack);
	multInfo2->Draw();

	cCorrYieldLcpKpi->cd();          // Go back to the main canvas before defining pad2
	TPad *pad2 = new TPad("pad2", "pad2", 0, 0.03, 1, 0.2);
	pad2->SetTopMargin(0);
	pad2->SetBottomMargin(0.2);
	pad2->Draw();
	pad2->cd();

	TH1D* hBasisClone2 = (TH1D*) hBasis->Clone("hBasisClone2");
	hBasisClone2->GetYaxis()->SetRangeUser(0.6,1.5);
	hBasisClone2->GetYaxis()->SetTitle("Ratio");
	hBasisClone2->GetYaxis()->SetLabelSize(0.13);
	hBasisClone2->GetYaxis()->SetTitleSize(0.13);
	hBasisClone2->GetYaxis()->CenterTitle(kTRUE);
	hBasisClone2->GetYaxis()->SetTitleOffset(0.3);
	hBasisClone2->GetYaxis()->SetNdivisions(508);
	hBasisClone2->GetYaxis()->ChangeLabel(1,-1,-1,-1,-1,-1," ");
	// hBasisClone2->GetYaxis()->ChangeLabel(7,-1,-1,-1,-1,-1," ");

	// hBasisClone2->GetXaxis()->SetNdivisions(999);
	hBasisClone2->GetXaxis()->SetTitle("#it{p}_{T}");
	// hBasisClone2->GetXaxis()->SetTitleOffset(0.00);
	hBasisClone2->GetXaxis()->SetLabelSize(0.16);
	hBasisClone2->GetXaxis()->SetTickSize(0.1);
	hBasisClone2->GetXaxis()->SetLabelOffset(0.01);
	hBasisClone2->GetXaxis()->SetNdivisions(510);
	// hBasisClone2->GetXaxis()->SetNdivisions(6);
	hBasisClone2->Draw("");

	TH1D* hCorrYieldRatio = (TH1D*) hCorrYield->Clone("hCorrYieldClone");
	nbins = hCorrYieldRatio->GetNbinsX();

	Printf("");
	Printf("%4.4s  %5.5s  %5.5s  %5.5s  %5.5s","ibin","value","err1","err2","errf");
	// hCorrYieldRatio->Divide(hCorrYieldRatio,hYieldLcpKpi,1,1,"B");	
	for(Int_t ibin=1;ibin<=nbins;ibin++) {
    	Double_t value =  hYieldLcpKpi->GetBinContent(ibin) > 0 ? hCorrYield->GetBinContent(ibin)/hYieldLcpKpi->GetBinContent(ibin) : 0;
    	Double_t err1  =  hYieldLcpKpi->GetBinContent(ibin) > 0 ? hYieldLcpKpi->GetBinError(ibin)/hYieldLcpKpi->GetBinContent(ibin) : 0;
    	Double_t err2  =  hCorrYield->GetBinContent(ibin) > 0 ? hCorrYield->GetBinError(ibin)/hCorrYield->GetBinContent(ibin) : 0;

    	Double_t errf  =  value*TMath::Sqrt(TMath::Abs(err1*err1+err2*err2));
    	hCorrYieldRatio->SetBinContent(ibin, value);
    	hCorrYieldRatio->SetBinError(ibin, errf);
    	Printf("%4.1d  %5.3f  %5.3f  %5.3f  %5.3f",ibin,value,err1,err2,errf);
    }
    Printf("");

	hCorrYieldRatio->Draw("Same");
	TLine* line1 = new TLine(0,1,12,1);
	line1->SetLineColor(kBlack); line1->SetLineStyle(2); line1->SetLineWidth(1);
	line1->Draw();

	TGraphAsymmErrors* hSystRatioNoFD = (TGraphAsymmErrors*) hSystTotScaledNoFD->Clone("hSystRatioNoFD");
	TGraphAsymmErrors* hSystRatioFD   = (TGraphAsymmErrors*) hSystScaledFD->Clone("hSystRatioFD");
	
	Printf("");
	Printf("%4.4s  %4.4s  %2.2s  %6.6s  %2.2s  %6.6s  %6.6s  %6.6s | %6.6s  %6.6s  %6.6s | %6.6s  %6.6s  %6.6s | %10.10s  %10.10s  %10.10s | %10.10s  %10.10s  %10.10s","ibin","pt","i","x","j","x2","value","errx","erryl","erryl2","errylfin","erryh","erryh2","erryhfin","erryl_fd","erryl2_fd","errylfin_fd","erryh_fd","erryh2_fd","erryhfin_fd");
	for(Int_t ibin=1;ibin<=nbins;ibin++) {
    	Double_t pt = hCorrYieldRatio->GetBinCenter(ibin);
    	Double_t ptwidth = hCorrYieldRatio->GetBinWidth(ibin);
    	Double_t value =  hCorrYieldRatio->GetBinContent(ibin);

	 	for(Int_t i=0; i<hSystTotScaledNoFD->GetN(); i++) {
	 		Double_t x=0, y=0, errx=0, erryl=0, erryh=0,erryl_fd=0,erryh_fd=0, erryl2_fd=0,erryh2_fd=0;
	 		Double_t erryfinl = 0, erryfinh = 0, erryfinl_fd = 0, erryfinh_fd = 0;
	 		
	 		hSystTotScaledNoFD->GetPoint(i,x,y);
	 		errx  = hSystTotScaledNoFD->GetErrorXlow(i) ;
	 		erryl = hSystTotScaledNoFD->GetErrorYlow(i)/y ;
	 		erryh = hSystTotScaledNoFD->GetErrorYhigh(i)/y ;

	 		erryl_fd = hSystScaledFD->GetErrorYlow(i)/y ;
	 		erryh_fd = hSystScaledFD->GetErrorYhigh(i)/y ;

	 		for(Int_t j=0; j<hSystTotLcpKpi->GetN(); j++){
	 			Double_t x2=0, y2=0, errx2=0, erryl2=0, erryh2=0;

	 			hSystTotLcpKpi->GetPoint(j,x2,y2);
		 		errx2  = hSystNoFDLcpKpi->GetErrorXlow(j) ;
		 		erryl2 = hSystNoFDLcpKpi->GetErrorYlow(j)/y2 ;
		 		erryh2 = hSystNoFDLcpKpi->GetErrorYhigh(j)/y2 ;

		 		erryl2_fd = hSystFDLcpKpi->GetErrorYlow(j)/y2 ;
		 		erryh2_fd = hSystFDLcpKpi->GetErrorYhigh(j)/y2 ;

		 		if(  (( (x-errx) <= pt) && ( (x+errx) >= pt)) && (( (x2-errx2) <= pt) && ( (x2+errx2) >= pt))   ){
		 			erryfinl = value* TMath::Sqrt(TMath::Abs(erryl*erryl+erryl2*erryl2));
		 			erryfinh = value* TMath::Sqrt(TMath::Abs(erryh*erryh+erryh2*erryh2));

		 			erryfinl_fd = value* TMath::Sqrt(TMath::Abs(erryl_fd*erryl_fd-erryl2_fd*erryl2_fd));
		 			erryfinh_fd = value* TMath::Sqrt(TMath::Abs(erryh_fd*erryh_fd-erryh2_fd*erryh2_fd));

		 			hSystRatioNoFD->SetPoint(i,x,value);
		 			hSystRatioNoFD->SetPointError(i,errx,errx,erryfinl,erryfinh);

		 			hSystRatioFD->SetPoint(i,x,value);
		 			hSystRatioFD->SetPointError(i,errx,errx,erryfinl_fd,erryfinh_fd);

		 			Printf("%4.1d  %4.1f  %2.1d  %6.2f  %2.1d  %6.2f  %6.3f  %6.3f | %6.4f  %6.4f  %6.4f | %6.4f  %6.4f  %6.4f | %10.4f  %10.4f  %10.4f | %10.4f  %10.4f  %10.4f",ibin,pt,i,x,j,x2,value,errx,erryl,erryl2,erryfinl,erryh,erryh2,erryfinh,erryl_fd,erryl2_fd,erryfinl_fd,erryh_fd,erryh2_fd,erryfinh_fd);
		 		}
	 		}
	 	}
 	}
 	hSystRatioFD->SetLineWidth(1);
	hSystRatioFD->SetMarkerStyle(22);
	// hSystRatioFD->SetFillStyle(3001);
	hSystRatioFD->SetFillColorAlpha(kBlack,0.4);
	hSystRatioFD->SetLineColorAlpha(kBlack,0.5);
	hSystRatioFD->Draw("E5Same");


	hSystRatioNoFD->SetLineColor(kBlack);
	hSystRatioNoFD->SetFillStyle(0);
	hSystRatioNoFD->Draw("E5Same");

	hCorrYieldRatio->Draw("Same");


	TGraphAsymmErrors* gPlotForFit = new TGraphAsymmErrors(0); 
	gPlotForFit->SetNameTitle("gPlotForFit","gPlotForFit");

	for(Int_t ibin=1;ibin<=nbins;ibin++) {
    	Double_t pt = hCorrYieldRatio->GetBinCenter(ibin);
    	Double_t ptwidth = hCorrYieldRatio->GetBinWidth(ibin);
    	Double_t value =  hCorrYieldRatio->GetBinContent(ibin);
    	Double_t staterr = hCorrYieldRatio->GetBinError(ibin);

	 	for(Int_t i=0; i<hSystRatioNoFD->GetN(); i++) {
	 		Double_t x=0, y=0, errx=0, erryl=0, erryh=0;

	 		hSystRatioNoFD->GetPoint(i,x,y);
	 		errx  = hSystRatioNoFD->GetErrorXlow(i) ;
	 		erryl = hSystRatioNoFD->GetErrorYlow(i);
	 		erryh = hSystRatioNoFD->GetErrorYhigh(i);
	 		
	 		for(Int_t j=0; j<hSystRatioFD->GetN(); j++){
	 			Double_t x2=0, y2=0, errx2=0, erryl2=0, erryh2=0;

	 			hSystRatioFD->GetPoint(j,x2,y2);
	 			errx2  = hSystRatioFD->GetErrorXlow(j) ;
	 			erryl2 = hSystRatioFD->GetErrorYlow(j);
	 			erryh2 = hSystRatioFD->GetErrorYhigh(j);

		 		if(  (( (x-errx) <= pt) && ( (x+errx) >= pt)) && (( (x2-errx2) <= pt) && ( (x2+errx2) >= pt))   ){
		 			
		 			Double_t erryfinl = TMath::Sqrt(staterr*staterr+erryl*erryl);
		 			Double_t erryfinh = TMath::Sqrt(staterr*staterr+erryh*erryh);
		 			
		 			gPlotForFit->SetPoint(i,x,value);
		 			gPlotForFit->SetPointError(i,errx,errx,erryfinl,erryfinh);

		 		}
	 		}
	 	}
 	}
 	gPlotForFit->SetLineColor(kRed);
 	gPlotForFit->SetMarkerColor(kRed);
 	gPlotForFit->Draw("Same");

 	TF1* fitFunc = new TF1("fitFunc","pol0");
 	gPlotForFit->Fit("fitFunc","0M");
 	TF1 *fitResult = gPlotForFit->GetFunction("fitFunc");
	TF1* fit;
	// if(i == 2 || i ==0) fit = new TF1(Form("fit%d",i), "pol0(0)",1,12);
	// else				
	fit = new TF1("fit", "pol0(0)",1,24);
	fit->SetParameter(0,fitResult->GetParameter(0));
	fit->SetParError(0,fitResult->GetParError(0));
	// fit->SetParameter(1,fitResult->GetParameter(1));
	// fit->SetParError(1,fitResult->GetParError(1));
	fit->SetLineWidth(4); fit->SetLineStyle(3); fit->SetLineColorAlpha(kRed,0.8);
	fit->Draw("Same");
	Printf("Reduced chi-squared = %.2f",fitResult->GetChisquare() / Double_t(fitResult->GetNDF()));

	cCorrYieldLcpKpi->SaveAs(Form("%s/CorrectedYieldComparison.pdf",outputDirName.Data()));
	cCorrYieldLcpKpi->Close();

	//--------------------------------------------------------------------------------
	// Third graph	
	// Results with simulation

	// Corrected yield per event
	TCanvas* cCorrYieldSim = new TCanvas("cCorrYieldSim","cCorrYieldSim",600,700);  
	cCorrYieldSim->cd();
	// gPad->SetLogy(1);
	gStyle->SetOptStat(0);
	gROOT->ForceStyle();

	// Upper plot will be in pad1
	TPad *pad3 = new TPad("pad3", "pad3", 0, 0.33, 1, 1.0);
	pad3->SetBottomMargin(0); // Upper and lower plot are joined
	pad3->SetTopMargin(0.05); // Upper and lower plot are joined
	pad3->SetLogy(1);
	// pad1->SetGridx();         // Vertical grid
	pad3->Draw();             // Draw the upper pad: pad1
	pad3->cd();               // pad1 becomes the current pad

	// Legend
	TLegend *legend3=new TLegend(0.5301 ,0.689977,0.876254, 0.889278);
	legend3->SetBorderSize(0);
	legend3->SetFillStyle(0);
	legend3->SetTextFont(42);
	// legend3->SetTextAlign(13);
	legend3->SetTextAlign(12);
	legend3->SetTextColor(kBlack);

	TH1D* hBasisClone3 = (TH1D*) hBasis->Clone("hBasisClone3");
	hBasisClone3->GetYaxis()->SetRangeUser(0.00000005,0.1);
	hBasisClone3->GetYaxis()->SetLabelSize(0.033);
	hBasisClone3->GetYaxis()->SetTitleSize(0.042);
	hBasisClone3->GetYaxis()->SetTitleOffset(1.1);
	hBasisClone3->Draw("");

	legend3->AddEntry(hCorrYield,legendTitle.Data(),"LP");

	hSystScaledFD->SetLineWidth(1.);
	hSystScaledFD->SetMarkerStyle(22);
	// hSystScaledFD->SetFillStyle(3001);
	hSystScaledFD->SetFillColorAlpha(kBlack,0.4);
	hSystScaledFD->SetLineColorAlpha(kBlack,0.5);
	hSystScaledFD->Draw("E5Same");

 	hSystTotScaledNoFD->SetLineColor(kBlack);
	hSystTotScaledNoFD->SetFillStyle(0);
	hSystTotScaledNoFD->Draw("E5Same");

	hCorrYield->Draw("Same");

	legend3->Draw();

	// Some text info
	TPaveText *prodInfo3=new TPaveText(0.123746,0.844989,0.431438,0.942424 ,"NDC");
	prodInfo3->SetBorderSize(0);
	prodInfo3->SetFillStyle(0);
	prodInfo3->AddText(Form("pp, #sqrt{s} = 13 TeV"));
	prodInfo3->AddText(Form("Prompt #Lambda_{c}^{+}, |#it{y}| < 0.5"));
	prodInfo3->SetTextFont(42);
	prodInfo3->SetTextAlign(12);
	prodInfo3->SetTextColor(kBlack);
	prodInfo3->Draw();

	TPaveText *multInfo3=new TPaveText(0.518395, 0.915851,0.841137,0.96014 ,"NDC");
	multInfo3->SetBorderSize(0);
	multInfo3->SetFillStyle(0);
	multInfo3->AddText(Form("Multiplicity: |#it{#eta}| < 1.0           "));
	multInfo3->SetTextFont(42);
	multInfo3->SetTextAlign(13);
	multInfo3->SetTextColor(kBlack);
	multInfo3->Draw();

	TPaveText *multLumBR3=new TPaveText( 0.0936455 ,0.0101422,0.665552,0.0810045,"NDC");
	multLumBR3->SetBorderSize(0);
	multLumBR3->SetFillStyle(0);
	multLumBR3->AddText(Form("#pm 5.0 %% lumi, #pm 5.2 %% BR uncertainty not shown"));
	multLumBR3->SetTextFont(42);
	multLumBR3->SetTextAlign(12);
	multLumBR3->SetTextColor(kBlack);
	multLumBR3->Draw();


	// Simulation
	hYieldLc_HardNoCR_Merged->SetLineColor(kGreen+3);
	hYieldLc_HardNoCR_Merged->SetMarkerColor(kGreen+3);
	hYieldLc_HardNoCR_Merged->SetLineWidth(3);
	hYieldLc_HardNoCR_Merged->SetLineStyle(1);
	hYieldLc_HardNoCR_Merged->Draw("HistSame");
	legend3->AddEntry(hYieldLc_HardNoCR_Merged,"PYTHIA 8, Hard QCD","LP");

	hYieldLc_HardCR_Merged->SetLineColor(kRed);
	hYieldLc_HardCR_Merged->SetMarkerColor(kRed);
	hYieldLc_HardCR_Merged->SetLineWidth(3);
	hYieldLc_HardCR_Merged->SetLineStyle(8);
	hYieldLc_HardCR_Merged->Draw("HistSame");
	legend3->AddEntry(hYieldLc_HardCR_Merged,"PYTHIA 8, Hard QCD + CR","LP");

	hYieldLc_SoftNoCR_Merged->SetLineColor(kMagenta);
	hYieldLc_SoftNoCR_Merged->SetMarkerColor(kMagenta);
	hYieldLc_SoftNoCR_Merged->SetLineWidth(3);
	hYieldLc_SoftNoCR_Merged->SetLineStyle(7);
	hYieldLc_SoftNoCR_Merged->Draw("HistSame");
	legend3->AddEntry(hYieldLc_SoftNoCR_Merged,"PYTHIA 8, Soft QCD","LP");

	hYieldLc_SoftCR_Merged->SetLineColor(kOrange-2);
	hYieldLc_SoftCR_Merged->SetMarkerColor(kOrange-2);
	hYieldLc_SoftCR_Merged->SetLineWidth(3);
	hYieldLc_SoftCR_Merged->SetLineStyle(9);
	hYieldLc_SoftCR_Merged->Draw("HistSame");
	legend3->AddEntry(hYieldLc_SoftCR_Merged,"PYTHIA 8, Soft QCD + CR","LP");

	// Upper plot will be in pad1
	cCorrYieldSim->cd();
	TPad *pad4 = new TPad("pad4", "pad4", 0, 0.18, 1, 0.33);
	pad4->SetBottomMargin(0); // Upper and lower plot are joined
	pad4->SetTopMargin(0); // Upper and lower plot are joined
	pad4->SetLogy(1);
	// pad1->SetGridx();         // Vertical grid
	pad4->Draw();             // Draw the upper pad: pad1
	pad4->cd();               // pad1 becomes the current pad

	// Legend
	TLegend *legend32=new TLegend( 0.79097 ,0.61998,0.943144,  0.94639);
	legend32->SetBorderSize(0);
	legend32->SetFillStyle(0);
	legend32->SetTextFont(42);
	// legend3->SetTextAlign(13);
	legend32->SetTextAlign(12);
	legend32->SetTextColor(kBlack);

	TH1D* hBasisClone4 = (TH1D*) hBasis->Clone("hBasisClone4");
	hBasisClone4->GetYaxis()->SetRangeUser(20,250);
	hBasisClone4->GetYaxis()->SetTitle("#frac{Data}{Hard QCD}");
	hBasisClone4->GetYaxis()->SetLabelSize(0.14);
	hBasisClone4->GetYaxis()->SetTitleSize(0.13);
	hBasisClone4->GetYaxis()->CenterTitle(kTRUE);
	hBasisClone4->GetYaxis()->SetTitleOffset(0.36);
	hBasisClone4->GetYaxis()->SetLabelOffset(0.0);
	hBasisClone4->GetYaxis()->SetNdivisions(510);

	hBasisClone4->GetXaxis()->SetTitle(" ");
	// hBasisClone5->GetXaxis()->SetTitleOffset(0.00);
	// hBasisClone4->GetXaxis()->SetLabelSize(0.16);
	hBasisClone4->GetXaxis()->SetTickSize(0.25);
	// hBasisClone4->GetXaxis()->SetLabelOffset(0.01);
	hBasisClone4->GetXaxis()->SetNdivisions(510);
	hBasisClone4->Draw("");

	TH1D* hCorrYieldRatioHardCR = (TH1D*) hCorrYield->Clone("hCorrYieldRatioHardCR");
	TH1D* hCorrYieldRatioHardNoCR = (TH1D*) hCorrYield->Clone("hCorrYieldRatioHardNoCR");
	TH1D* hCorrYieldRatioSoftCR = (TH1D*) hCorrYield->Clone("hCorrYieldRatioSoftCR");
	TH1D* hCorrYieldRatioSoftNoCR = (TH1D*) hCorrYield->Clone("hCorrYieldRatioSoftNoCR");

	TGraphAsymmErrors* hSystRatioHardCRNoFD = (TGraphAsymmErrors*) hSystTotScaledNoFD->Clone("hSystRatioHardCRNoFD");
	TGraphAsymmErrors* hSystRatioHardCRFD   = (TGraphAsymmErrors*) hSystScaledFD->Clone("hSystRatioHardCRFD");

	TGraphAsymmErrors* hSystRatioHardNoFD = (TGraphAsymmErrors*) hSystTotScaledNoFD->Clone("hSystRatioHardNoFD");
	TGraphAsymmErrors* hSystRatioHardFD   = (TGraphAsymmErrors*) hSystScaledFD->Clone("hSystRatioHardFD");

	TGraphAsymmErrors* hSystRatioSoftCRNoFD = (TGraphAsymmErrors*) hSystTotScaledNoFD->Clone("hSystRatioSoftCRNoFD");
	TGraphAsymmErrors* hSystRatioSoftCRFD   = (TGraphAsymmErrors*) hSystScaledFD->Clone("hSystRatioSoftCRFD");

	TGraphAsymmErrors* hSystRatioSoftNoFD = (TGraphAsymmErrors*) hSystTotScaledNoFD->Clone("hSystRatioSoftNoFD");
	TGraphAsymmErrors* hSystRatioSoftFD   = (TGraphAsymmErrors*) hSystScaledFD->Clone("hSystRatioSoftFD");

	for (Int_t i = 1; i <= hCorrYieldRatioHardCR->GetNbinsX(); ++i)
	{
		Double_t pt = hCorrYieldRatioHardCR->GetBinCenter(i);
		Int_t bin = hYieldLc_HardCR_Merged->FindBin(pt);

		Double_t ratioHardCR   = hYieldLc_HardCR_Merged->GetBinContent(bin)>0 ? hCorrYieldRatioHardCR->GetBinContent(i)/hYieldLc_HardCR_Merged->GetBinContent(bin) : 0;
		Double_t ratioHardNoCR = hYieldLc_HardNoCR_Merged->GetBinContent(bin)>0 ? hCorrYieldRatioHardNoCR->GetBinContent(i)/hYieldLc_HardNoCR_Merged->GetBinContent(bin) : 0;
		Double_t ratioSoftCR   = hYieldLc_SoftCR_Merged->GetBinContent(bin)>0 ? hCorrYieldRatioSoftCR->GetBinContent(i)/hYieldLc_SoftCR_Merged->GetBinContent(bin) : 0;
		Double_t ratioSoftNoCr = hYieldLc_SoftNoCR_Merged->GetBinContent(bin)>0 ? hCorrYieldRatioSoftNoCR->GetBinContent(i)/hYieldLc_SoftNoCR_Merged->GetBinContent(bin) : 0;


		hCorrYieldRatioHardCR->SetBinContent(i,ratioHardCR);	
		hCorrYieldRatioHardNoCR->SetBinContent(i,ratioHardNoCR);
		hCorrYieldRatioSoftCR->SetBinContent(i,ratioSoftCR);	
		hCorrYieldRatioSoftNoCR->SetBinContent(i,ratioSoftNoCr);

		hCorrYieldRatioHardCR->SetBinError(i,ratioHardCR*hCorrYield->GetBinError(i));	
		hCorrYieldRatioHardNoCR->SetBinError(i,ratioHardNoCR*hCorrYield->GetBinError(i));
		hCorrYieldRatioSoftCR->SetBinError(i,ratioSoftCR*hCorrYield->GetBinError(i));	
		hCorrYieldRatioSoftNoCR->SetBinError(i,ratioSoftNoCr*hCorrYield->GetBinError(i));	

		for (Int_t j = 0; j < hSystTotScaledNoFD->GetN(); ++j)
		{
			Double_t x=0, y=0, errx=0, erryl=0, erryh=0, erryl_fd=0,erryh_fd=0;

	 		hSystTotScaledNoFD->GetPoint(j,x,y);
	 		errx  = hSystTotScaledNoFD->GetErrorXlow(j) ;
	 		erryl = hSystTotScaledNoFD->GetErrorYlow(j)/y ;
	 		erryh = hSystTotScaledNoFD->GetErrorYhigh(j)/y ;

	 		erryl_fd = hSystScaledFD->GetErrorYlow(j)/y ;
	 		erryh_fd = hSystScaledFD->GetErrorYhigh(j)/y ;

	 		if(  (( (x-errx) <= pt) && ( (x+errx) >= pt)) ){
	 			hSystRatioHardCRNoFD	->SetPoint(j,x,ratioHardCR);
	 			hSystRatioHardCRNoFD	->SetPointError(j,errx,errx,ratioHardCR*erryl,ratioHardCR*erryh);
	 			hSystRatioHardCRFD	    ->SetPoint(j,x,ratioHardCR);
	 			hSystRatioHardCRFD      ->SetPointError(j,errx,errx,ratioHardCR*erryl_fd,ratioHardCR*erryh_fd);

	 			hSystRatioHardNoFD		->SetPoint(j,x,ratioHardNoCR);
	 			hSystRatioHardNoFD		->SetPointError(j,errx,errx,ratioHardNoCR*erryl,ratioHardNoCR*erryh);
	 			hSystRatioHardFD	    ->SetPoint(j,x,ratioHardNoCR);
	 			hSystRatioHardFD      	->SetPointError(j,errx,errx,ratioHardNoCR*erryl_fd,ratioHardNoCR*erryh_fd);

	 			hSystRatioSoftCRNoFD	->SetPoint(j,x,ratioSoftCR);
	 			hSystRatioSoftCRNoFD	->SetPointError(j,errx,errx,ratioSoftCR*erryl,ratioSoftCR*erryh);
	 			hSystRatioSoftCRFD	    ->SetPoint(j,x,ratioSoftCR);
	 			hSystRatioSoftCRFD      ->SetPointError(j,errx,errx,ratioSoftCR*erryl_fd,ratioSoftCR*erryh_fd);

	 			hSystRatioSoftNoFD	->SetPoint(j,x,ratioSoftNoCr);
	 			hSystRatioSoftNoFD	->SetPointError(j,errx,errx,ratioSoftNoCr*erryl,ratioSoftNoCr*erryh);
	 			hSystRatioSoftFD	->SetPoint(j,x,ratioSoftNoCr);
	 			hSystRatioSoftFD    ->SetPointError(j,errx,errx,ratioSoftNoCr*erryl_fd,ratioSoftNoCr*erryh_fd);
	 		}

	 	}

	}

	hSystRatioHardCRFD->SetLineWidth(1.);
	hSystRatioHardCRFD->SetMarkerStyle(21);
	// hSystScaledFD->SetFillStyle(3001);
	hSystRatioHardCRFD->SetFillColorAlpha(kRed,0.4);
	hSystRatioHardCRFD->SetLineColorAlpha(kRed,0.5);
	hSystRatioHardCRFD->Draw("E5Same");

 	hSystRatioHardCRNoFD->SetLineColor(kRed);
	hSystRatioHardCRNoFD->SetFillStyle(0);
	hSystRatioHardCRNoFD->Draw("E5Same");

	hCorrYieldRatioHardCR->SetLineColor(kRed);
	hCorrYieldRatioHardCR->SetMarkerStyle(21);
	hCorrYieldRatioHardCR->SetMarkerSize(0.8);
	hCorrYieldRatioHardCR->SetMarkerColor(kRed);
	hCorrYieldRatioHardCR->Draw("Same");

	hSystRatioHardFD->SetLineWidth(1.);
	hSystRatioHardFD->SetMarkerStyle(21);
	// hSystRatioHardFD->SetFillStyle(3001);
	hSystRatioHardFD->SetFillColorAlpha(kGreen+3,0.4);
	hSystRatioHardFD->SetLineColorAlpha(kGreen+3,0.5);
	hSystRatioHardFD->Draw("E5Same");

 	hSystRatioHardNoFD->SetLineColor(kGreen+3);
	hSystRatioHardNoFD->SetFillStyle(0);
	hSystRatioHardNoFD->Draw("E5Same");

	hCorrYieldRatioHardNoCR->SetLineColor(kGreen+3);
	hCorrYieldRatioHardNoCR->SetMarkerStyle(21);
	hCorrYieldRatioHardNoCR->SetMarkerSize(0.8);
	hCorrYieldRatioHardNoCR->SetMarkerColor(kGreen+3);
	hCorrYieldRatioHardNoCR->Draw("Same");

	legend32->AddEntry(hCorrYieldRatioHardNoCR,"w/o CR","LP");
	legend32->AddEntry(hCorrYieldRatioHardCR,"CR","LP");
	legend32->Draw();

	// Upper plot will be in pad1
	cCorrYieldSim->cd();
	TPad *pad5 = new TPad("pad5", "pad5", 0, 0.03, 1, 0.195);
	pad5->SetTopMargin(0); // Upper and lower plot are joined
	pad5->SetBottomMargin(0.2); // Upper and lower plot are joined

	pad5->SetLogy(1);
	// pad1->SetGridx();         // Vertical grid
	pad5->Draw();             // Draw the upper pad: pad1
	pad5->cd();               // pad1 becomes the current pad

	// Legend
	TLegend *legend33=new TLegend( 0.79097 ,0.61998,0.943144,  0.94639);
	legend33->SetBorderSize(0);
	legend33->SetFillStyle(0);
	legend33->SetTextFont(42);
	// legend3->SetTextAlign(13);
	legend33->SetTextAlign(12);
	legend33->SetTextColor(kBlack);

	TH1D* hBasisClone5 = (TH1D*) hBasis->Clone("hBasisClone5");
	hBasisClone5->GetYaxis()->SetRangeUser(1,80);
	hBasisClone5->GetYaxis()->SetTitle("#frac{Data}{Soft QCD}");
	hBasisClone5->GetYaxis()->SetLabelSize(0.13);
	hBasisClone5->GetYaxis()->SetTitleSize(0.13);
	hBasisClone5->GetYaxis()->CenterTitle(kTRUE);
	hBasisClone5->GetYaxis()->SetTitleOffset(0.36);
	hBasisClone5->GetYaxis()->SetLabelOffset(0.0);
	hBasisClone5->GetYaxis()->SetNdivisions(510);
	hBasisClone5->GetYaxis()->ChangeLabel(1,-1,-1,-1,-1,-1," ");

	hBasisClone5->GetXaxis()->SetTitle("#it{p}_{T}");
	// hBasisClone5->GetXaxis()->SetTitleOffset(0.00);
	hBasisClone5->GetXaxis()->SetLabelSize(0.16);
	hBasisClone5->GetXaxis()->SetTickSize(0.1);
	hBasisClone5->GetXaxis()->SetLabelOffset(0.01);
	hBasisClone5->GetXaxis()->SetNdivisions(510);
	hBasisClone5->Draw("");

	hSystRatioSoftCRFD->SetLineWidth(1.);
	hSystRatioSoftCRFD->SetMarkerStyle(21);
	// hSystScaledFD->SetFillStyle(3001);
	hSystRatioSoftCRFD->SetFillColorAlpha(kOrange-2,0.4);
	hSystRatioSoftCRFD->SetLineColorAlpha(kOrange-2,0.5);
	hSystRatioSoftCRFD->Draw("E5Same");

 	hSystRatioSoftCRNoFD->SetLineColor(kOrange-2);
	hSystRatioSoftCRNoFD->SetFillStyle(0);
	hSystRatioSoftCRNoFD->Draw("E5Same");

	hCorrYieldRatioSoftCR->SetLineColor(kOrange-2);
	hCorrYieldRatioSoftCR->SetMarkerStyle(21);
	hCorrYieldRatioSoftCR->SetMarkerSize(0.8);
	hCorrYieldRatioSoftCR->SetMarkerColor(kOrange-2);
	hCorrYieldRatioSoftCR->Draw("Same");

	hSystRatioSoftFD->SetLineWidth(1.);
	hSystRatioSoftFD->SetMarkerStyle(21);
	// hSystRatioSoftFD->SetFillStyle(3001);
	hSystRatioSoftFD->SetFillColorAlpha(kMagenta,0.4);
	hSystRatioSoftFD->SetLineColorAlpha(kMagenta,0.5);
	hSystRatioSoftFD->Draw("E5Same");

 	hSystRatioSoftNoFD->SetLineColor(kMagenta);
	hSystRatioSoftNoFD->SetFillStyle(0);
	hSystRatioSoftNoFD->Draw("E5Same");

	hCorrYieldRatioSoftNoCR->SetLineColor(kMagenta);
	hCorrYieldRatioSoftNoCR->SetMarkerStyle(21);
	hCorrYieldRatioSoftNoCR->SetMarkerSize(0.8);
	hCorrYieldRatioSoftNoCR->SetMarkerColor(kMagenta);
	hCorrYieldRatioSoftNoCR->Draw("Same");

	legend33->AddEntry(hCorrYieldRatioSoftNoCR,"w/o CR","LP");
	legend33->AddEntry(hCorrYieldRatioSoftCR,"CR","LP");
	legend33->Draw();


	cCorrYieldSim->SaveAs(Form("%s/CorrectedYieldSim.pdf",outputDirName.Data()));
	cCorrYieldSim->Close();

	//--------------------------------------------------------------------------------
	// Fourth graph
	// Lc/D0

	// Corrected yield per event
	TCanvas* cLcD0 = new TCanvas("cLcD0","cLcD0",600,600);  
	cLcD0->cd();
	// gPad->SetLogy(1);
	gStyle->SetOptStat(0);
	gROOT->ForceStyle();

	// Legend
	TLegend *legendLcD0=new TLegend( 0.560201 ,0.709059,0.88796 , 0.848432);
	legendLcD0->SetBorderSize(0);
	legendLcD0->SetFillStyle(0);
	legendLcD0->SetTextFont(42);
	// legend->SetTextAlign(13);
	legendLcD0->SetTextAlign(12);
	legendLcD0->SetTextColor(kBlack);

	// Legend
	TLegend *legendSystLcD0=new TLegend(0.560201  ,0.594077   ,0.88796 ,0.707317);
	legendSystLcD0->SetBorderSize(0);
	legendSystLcD0->SetFillStyle(0);
	legendSystLcD0->SetTextFont(42);
	// leglegendSystLcD0end->SetTextAlign(13);
	legendSystLcD0->SetTextAlign(12);
	legendSystLcD0->SetTextColor(kBlack);

	// My results, LcK0Sp
	TH1D* hBasisLcD0 = new TH1D("hBasisLcD0","hBasisLcD0",nPtBins,ptBins);
	TH1D* hLcD0      = (TH1D*) hCorrYield->Clone("hLcD0");
	TH1D* hCorrCross = (TH1D*) hSigmaCorr->Clone("hCorrCross");
	hCorrCross->Scale(1/(BR));
	hSigmaCorrD0->Scale(1/(0.0388*0.94));

	for (Int_t bin = 1; bin < hLcD0->GetNbinsX(); ++bin)
	{
		
		Double_t crossLc = hCorrCross->GetBinContent(bin);
		Double_t errLc   = hCorrCross->GetBinError(bin);
		Double_t pt      = hCorrCross->GetBinCenter(bin);
		Double_t pterr   = hCorrCross->GetBinError(bin);

		Double_t crossD0 = hSigmaCorrD0->GetBinContent(bin);
		Double_t errD0   = hSigmaCorrD0->GetBinError(bin);
		Double_t pt2      = hSigmaCorrD0->GetBinCenter(bin);

		Double_t ratio = crossD0 > 0 ? crossLc/crossD0 : 0.;

		Double_t errRatio = 0.;
		if(crossLc > 0 && crossD0 > 0){
			errRatio = ratio*TMath::Sqrt(TMath::Abs( TMath::Power(errLc/crossLc,2.0) + TMath::Power(errD0/crossD0,2.0)   ));
		}

		if(TMath::Abs(pt-pt2)<0.01){
			hLcD0->SetBinContent(bin,ratio);
			hLcD0->SetBinError(bin,errRatio);
		}

	}

	legendLcD0->AddEntry(hLcD0,legendTitle.Data(),"LP");

	hBasisLcD0->SetTitle("");
	hBasisLcD0->SetStats(kFALSE);
	hBasisLcD0->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
	hBasisLcD0->GetYaxis()->SetTitle("#Lambda_{c}^{+} / D^{0}");
	hBasisLcD0->GetYaxis()->SetTitleOffset(1.3);
	hBasisLcD0->GetYaxis()->SetRangeUser(0.0,1.0);
	hBasisLcD0->GetXaxis()->SetRangeUser(0,12);
	hBasisLcD0->Draw("");

	// Systematics
	// My results, systematics
	TGraphAsymmErrors* hSystLcD0NoFD = new TGraphAsymmErrors(0); 
	hSystLcD0NoFD->SetNameTitle("hSystLcD0NoFD","hSystLcD0NoFD");
	TGraphAsymmErrors* hSystLcD0FD = new TGraphAsymmErrors(0); 
	hSystLcD0FD->SetNameTitle("hSystLcD0FD","hSystLcD0FD");

	Double_t toterr = 0.;
	Double_t temperrLc = 0.;
	Double_t temperrD0 = 0.;
	Double_t cuterrl = 0.;
	Double_t cuterrh = 0.;
	Double_t FDerrl = 0.;
	Double_t FDerrh = 0.;

	Printf("");
	Printf("%4.4s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s  %6.6s","bin","pt","ratio","mcLc","mcD0","pidLc","pidD0","trkLc","trkD0","rawLc","rawD0","fdLc_l","fdLc_h","fdD0_l","fdD0_h","ctLc_l","ctLc_h","cutD0");
	for (Int_t bin = 1; bin < hLcD0->GetNbinsX(); ++bin)
	{
		Double_t pt = hLcD0->GetBinCenter(bin);
		Double_t ratio = hLcD0->GetBinContent(bin);

		printf("%4.1d  %6.1f  %6.2f  ",bin,pt,ratio);

		// MC-Pt-Shape (Uncorrelated)
		temperrLc = systLc->GetMCPtShapeErr(pt);
		temperrD0 = systD0->GetMCPtShapeErr(pt);
		toterr += temperrLc*temperrLc+temperrD0*temperrD0;

		printf("%6.3f  %6.3f  ",temperrLc,temperrD0);

		// PID (Correlated)
		temperrLc = systLc->GetPIDEffErr(pt);
		temperrD0 = systD0->GetPIDEffErr(pt);
		toterr += TMath::Abs(temperrLc*temperrLc-temperrD0*temperrD0);

		printf("%6.3f  %6.3f  ",temperrLc,temperrD0);

		// Tracking (Correlated)
		temperrLc = systLc->GetTrackingEffErr(pt);
		temperrD0 = systD0->GetTrackingEffErr(pt);
		toterr += TMath::Abs(temperrLc*temperrLc-temperrD0*temperrD0);

		printf("%6.3f  %6.3f  ",temperrLc,temperrD0);

		// Raw-yield (Uncorrelated)
		temperrLc = systLc->GetRawYieldErr(pt);
		temperrD0 = systD0->GetRawYieldErr(pt);
		toterr += temperrLc*temperrLc+temperrD0*temperrD0;

		printf("%6.3f  %6.3f  ",temperrLc,temperrD0);

		// ------------

		// Feed-down (Correlated)
		for (Int_t i = 0; i < hSystFDD0->GetN(); ++i)
		{
			Double_t x=0,y=0,errx=0,yerrl=0,yerrh=0;
			hSystFDD0->GetPoint(i,x,y);
			errx = hSystFDD0->GetErrorXlow(i);
			yerrl = y>0 ? hSystFDD0->GetErrorYlow(i)/y : 0;
			yerrh = y>0 ? hSystFDD0->GetErrorYhigh(i)/y : 0;

			for (Int_t j = 0; j < hSystFD->GetN(); ++j)
			{
				Double_t x2=0,y2=0,errx2=0,yerrl2=0,yerrh2=0;
				hSystFD->GetPoint(j,x2,y2);
				errx2 = hSystFD->GetErrorXlow(j);
				yerrl2 = hSystFD->GetErrorYlow(j);
				yerrh2 = hSystFD->GetErrorYhigh(j);

				if(  (( (x-errx) <= pt) && ( (x+errx) >= pt)) && (( (x2-errx2) <= pt) && ( (x2+errx2) >= pt)) ){
					temperrD0 = systD0->GetCutsEffErr(pt);
					FDerrl = ratio*TMath::Sqrt(TMath::Abs(yerrl*yerrl-yerrl2*yerrl2));
					FDerrh = ratio*TMath::Sqrt(TMath::Abs(yerrh*yerrh-yerrh2*yerrh2));

					printf("%6.3f  %6.3f  %6.3f  %6.3f  ",yerrl2,yerrh2,yerrl,yerrh);

					hSystLcD0FD->SetPoint(i,x2,ratio);
					hSystLcD0FD->SetPointError(i,0.3,0.3,FDerrl,FDerrh);
				}

			}
		}

		// Cut variation (Uncorrelated)
		for (Int_t k = 0; k < hSystCutVar->GetN(); ++k)
		{
			Double_t x=0,y=0,errx=0,yerrl=0,yerrh=0;
			hSystCutVar->GetPoint(k,x,y);
			errx = hSystCutVar->GetErrorXlow(k);
			yerrl = hSystCutVar->GetErrorYlow(k);
			yerrh = hSystCutVar->GetErrorYhigh(k);

			if(  (( (x-errx) <= pt) && ( (x+errx) >= pt)) ){
				temperrD0 = systD0->GetCutsEffErr(pt);
				cuterrl = temperrD0*temperrD0+yerrl*yerrl;
				cuterrh = temperrD0*temperrD0+yerrh*yerrh;

				printf("%6.3f  %6.3f  %6.3f ",yerrl,yerrh,temperrD0);
			}
		}

		hSystLcD0NoFD->SetPoint(bin-1,pt,ratio);
		hSystLcD0NoFD->SetPointError(bin-1,0.3,0.3,ratio*TMath::Sqrt(toterr+cuterrl),ratio*TMath::Sqrt(toterr+cuterrh));

		printf("\n");
	}
	Printf("BR unc. uncorr. = ± %.3f", TMath::Sqrt(TMath::Power(0.0005,2.0)+TMath::Power(0.051,2.0)));

	hSystLcD0FD->SetLineWidth(1.);
	hSystLcD0FD->SetMarkerStyle(22);
	// hSystScaledFD->SetFillStyle(3001);
	hSystLcD0FD->SetFillColorAlpha(kBlack,0.4);
	hSystLcD0FD->SetLineColorAlpha(kBlack,0.5);
	hSystLcD0FD->Draw("E5Same");

	hSystLcD0NoFD->SetLineColor(kBlack);
	hSystLcD0NoFD->SetFillStyle(0);
	hSystLcD0NoFD->Draw("E5Same");

	hLcD0->Draw("Same");

	// hSystLcpKpiLcD0_fd->SetLineWidth(1);
	// hSystLcpKpiLcD0_fd->SetMarkerStyle(22);
	// hSystLcpKpiLcD0_fd->SetFillColorAlpha(kRed-7,0.4);
	// hSystLcpKpiLcD0_fd->SetLineColorAlpha(kRed-7,0.5);
	// hSystLcpKpiLcD0_fd->Draw("E5Same");

	// hSystLcpKpiLcD0_nofd->SetLineColor(kRed-7);
	// hSystLcpKpiLcD0_nofd->SetFillStyle(0);
	// hSystLcpKpiLcD0_nofd->Draw("E5Same");

	// hLcpKpiLcD0->SetLineColor(kRed-7);
	// hLcpKpiLcD0->SetMarkerStyle(20);
	// hLcpKpiLcD0->SetMarkerColor(kRed-7);
	// hLcpKpiLcD0->Draw("Same");

	// Simulation ratio plot
	TH1D* hSimLcD0_HardCR   = (TH1D*) hYieldLc_HardCR_Merged->Clone("hSimLcD0_HardCR");
	TH1D* hSimLcD0_HardNoCR = (TH1D*) hYieldLc_HardNoCR_Merged->Clone("hSimLcD0_HardNoCR");
	TH1D* hSimLcD0_SoftCR   = (TH1D*) hYieldLc_SoftCR_Merged->Clone("hSimLcD0_SoftCR");
	TH1D* hSimLcD0_SoftNoCR   = (TH1D*) hYieldLc_SoftNoCR_Merged->Clone("hSimLcD0_SoftNoCR");

	for (Int_t i = 1; i <= hYieldLc_HardCR_Merged->GetNbinsX(); ++i)
	{	
		// HardCR
		Double_t yieldLc = hYieldLc_HardCR_Merged->GetBinContent(i);
		Double_t yieldD0 = hYieldD0_HardCR_Merged->GetBinContent(i);
		Double_t errLc   = hYieldLc_HardCR_Merged->GetBinError(i);
		Double_t errD0   = hYieldD0_HardCR_Merged->GetBinError(i);
		
		Double_t ratio = yieldD0 > 0 ? yieldLc/yieldD0 : 0;
		Double_t errtot = ratio*TMath::Sqrt(errLc*errLc+errD0+errD0);
		
		hSimLcD0_HardCR->SetBinContent(i,ratio);
		// hSimLcD0_HardCR->SetBinError(i,errtot);

		// HardNoCR
		yieldLc = hYieldLc_HardNoCR_Merged->GetBinContent(i);
		yieldD0 = hYieldD0_HardNoCR_Merged->GetBinContent(i);
		errLc   = hYieldLc_HardNoCR_Merged->GetBinError(i);
		errD0   = hYieldD0_HardNoCR_Merged->GetBinError(i);
		
		ratio = yieldD0 > 0 ? yieldLc/yieldD0 : 0;
		errtot = ratio*TMath::Sqrt(errLc*errLc+errD0+errD0);
		
		hSimLcD0_HardNoCR->SetBinContent(i,ratio);
		// hSimLcD0_HardNoCR->SetBinError(i,errtot);

		// SoftCR
		yieldLc = hYieldLc_SoftCR_Merged->GetBinContent(i);
		yieldD0 = hYieldD0_SoftCR_Merged->GetBinContent(i);
		errLc   = hYieldLc_SoftCR_Merged->GetBinError(i);
		errD0   = hYieldD0_SoftCR_Merged->GetBinError(i);
		
		ratio = yieldD0 > 0 ? yieldLc/yieldD0 : 0;
		errtot = ratio*TMath::Sqrt(errLc*errLc+errD0+errD0);
		
		hSimLcD0_SoftCR->SetBinContent(i,ratio);
		// hSimLcD0_SoftCR->SetBinError(i,errtot);

		// SoftNoCR
		yieldLc = hYieldLc_SoftNoCR_Merged->GetBinContent(i);
		yieldD0 = hYieldD0_SoftNoCR_Merged->GetBinContent(i);
		errLc   = hYieldLc_SoftNoCR_Merged->GetBinError(i);
		errD0   = hYieldD0_SoftNoCR_Merged->GetBinError(i);
		
		ratio = yieldD0 > 0 ? yieldLc/yieldD0 : 0;
		errtot = ratio*TMath::Sqrt(errLc*errLc+errD0+errD0);
		
		hSimLcD0_SoftNoCR->SetBinContent(i,ratio);
		// hSimLcD0_SoftNoCR->SetBinError(i,errtot);

	}

	// hSimLcD0_HardCR->Draw("LSame");
	// hSimLcD0_HardNoCR->Draw("LSame");
	hSimLcD0_SoftCR->Draw("LSame");
	// hSimLcD0_SoftNoCR->Draw("LSame");

	legendLcD0->AddEntry(hSimLcD0_SoftCR,"PYTHIA8, Soft QCD + CR","L");

	TPad* padSimulation = (TPad*) canvasSimulations->GetPrimitive("cLambdaCoverD_1");
	// TH1D *h1DtempLc_Mon2013def    =  	(TH1D*) canvasSimulations->GetPrimitive("h1DtempLc_Mon2013def");
	TH1D *h1DtempLcY05_Mon2013def =  	(TH1D*) padSimulation->GetPrimitive("h1DtempLcY05_Mon2013def");
	// TH1D *h1DtempLc_Md2_13tev     = 	(TH1D*) canvasSimulations->GetPrimitive("h1DtempLc_Md2_13tev");
	TH1D *h1DtempLcY05_Md2_13tev  =  	(TH1D*) padSimulation->GetPrimitive("h1DtempLcY05_Md2_13tev");
	// TH1D *h1DtempLc_Md0_13tev     = 	(TH1D*) canvasSimulations->GetPrimitive("h1DtempLc_Md0_13tev");
	TH1D *h1DtempLcY05_Md0_13tev  =  	(TH1D*) padSimulation->GetPrimitive("h1DtempLcY05_Md0_13tev");
	// TH1D *h1DtempLc_Md3_13tev     = 	(TH1D*) canvasSimulations->GetPrimitive("h1DtempLc_Md3_13tev");
	TH1D *h1DtempLcY05_Md3_13tev  =  	(TH1D*) padSimulation->GetPrimitive("h1DtempLcY05_Md3_13tev");


	TH1D* hLcD0_Mon2013def   = new TH1D("hLcD0_Mon2013def","hLcD0_Mon2013def",nPtBins,ptBins);
	TH1D* hLcD0_Md2_13tev    = new TH1D("hLcD0_Md2_13tev","hLcD0_Md2_13tev",nPtBins,ptBins);
	TH1D* hLcD0_Md0_13tev    = new TH1D("hLcD0_Md0_13tev","hLcD0_Md0_13tev",nPtBins,ptBins);
	TH1D* hLcD0_Md3_13tev    = new TH1D("hLcD0_Md3_13tev","hLcD0_Md3_13tev",nPtBins,ptBins);

	Double_t nAdditions2[nPtBins]={0.0};
	for (Int_t i = 1; i <= h1DtempLcY05_Mon2013def->GetNbinsX(); ++i)
	{
		Double_t centre = h1DtempLcY05_Mon2013def->GetBinCenter(i);	
		Int_t bin = hLcD0_Mon2013def->FindBin(centre);

		hLcD0_Mon2013def  ->SetBinContent(bin,hLcD0_Mon2013def  ->GetBinContent(bin)+h1DtempLcY05_Mon2013def  ->GetBinContent(i));
		hLcD0_Md2_13tev->SetBinContent(bin,hLcD0_Md2_13tev->GetBinContent(bin)+h1DtempLcY05_Md2_13tev->GetBinContent(i));
		hLcD0_Md0_13tev  ->SetBinContent(bin,hLcD0_Md0_13tev  ->GetBinContent(bin)+h1DtempLcY05_Md0_13tev  ->GetBinContent(i));
		hLcD0_Md3_13tev->SetBinContent(bin,hLcD0_Md3_13tev->GetBinContent(bin)+h1DtempLcY05_Md3_13tev->GetBinContent(i));
	
		for(Int_t iPt = 0; iPt < nPtBins; ++iPt){
			if( centre > ptBins[iPt] &&  centre < ptBins[iPt+1]){
				nAdditions2[iPt]+=1.;
				break;
			}
		}
	}
	for(Int_t iPt = 0; iPt < nPtBins; ++iPt){
		hLcD0_Mon2013def  ->SetBinContent(iPt+1,hLcD0_Mon2013def  ->GetBinContent(iPt+1)/nAdditions2[iPt]);
		hLcD0_Md2_13tev->SetBinContent(iPt+1,hLcD0_Md2_13tev->GetBinContent(iPt+1)/nAdditions2[iPt]);
		hLcD0_Md0_13tev  ->SetBinContent(iPt+1,hLcD0_Md0_13tev  ->GetBinContent(iPt+1)/nAdditions2[iPt]);
		hLcD0_Md3_13tev->SetBinContent(iPt+1,hLcD0_Md3_13tev->GetBinContent(iPt+1)/nAdditions2[iPt]);
	}


	hLcD0_Mon2013def->SetLineColor(kCyan+2);
	hLcD0_Mon2013def->SetMarkerColor(kCyan+2);
	hLcD0_Mon2013def->SetLineWidth(3);
	hLcD0_Mon2013def->SetMarkerSize(0);
	hLcD0_Mon2013def->SetMarkerStyle(1);
	hLcD0_Mon2013def->SetLineStyle(10);
	hLcD0_Mon2013def->Draw("LSame");

	hLcD0_Md2_13tev->SetLineColor(kMagenta-2);
	hLcD0_Md2_13tev->SetMarkerColor(kMagenta-2);
	hLcD0_Md2_13tev->SetLineWidth(3);
	hLcD0_Md2_13tev->SetMarkerSize(0);
	hLcD0_Md2_13tev->SetMarkerStyle(1);
	hLcD0_Md2_13tev->SetLineStyle(5);
	hLcD0_Md2_13tev->Draw("LSame");

	legendLcD0->AddEntry(hLcD0_Mon2013def,"PYTHIA8, Monash","L");
	legendLcD0->AddEntry(hLcD0_Md2_13tev,"PYTHIA8, Mode2","L");

	// Some text info
	TPaveText *prodInfoLcD0=new TPaveText(0.113712,0.80662,0.503344,0.891986,"NDC");
	prodInfoLcD0->SetBorderSize(0);
	prodInfoLcD0->SetFillStyle(0);
	prodInfoLcD0->AddText(Form("pp, #sqrt{s} = 13 TeV, |#it{y}| < 0.5"));
	prodInfoLcD0->AddText(Form("BR #pm 5.2 %% uncertainty not shown"));
	prodInfoLcD0->SetTextFont(42);
	prodInfoLcD0->SetTextAlign(12);
	prodInfoLcD0->SetTextColor(kBlack);
	prodInfoLcD0->Draw();

	TPaveText *multInfoLcD0=new TPaveText(0.560201,0.848432,0.884615,0.912892,"NDC");
	multInfoLcD0->SetBorderSize(0);
	multInfoLcD0->SetFillStyle(0);
	multInfoLcD0->AddText(Form("Multiplicity: |#it{#eta}| < 1.0           "));
	multInfoLcD0->SetTextFont(42);
	multInfoLcD0->SetTextAlign(13);
	multInfoLcD0->SetTextColor(kBlack);
	multInfoLcD0->Draw();

	// TPaveText *multLumBRLcD0=new TPaveText(0.11204,0.114983,0.520067,0.181185,"NDC");
	// multLumBRLcD0->SetBorderSize(0);
	// multLumBRLcD0->SetFillStyle(0);
	
	// multLumBRLcD0->SetTextFont(42);
	// multLumBRLcD0->SetTextAlign(12);
	// multLumBRLcD0->SetTextColor(kBlack);
	// multLumBRLcD0->Draw();

	legendSystLcD0->AddEntry(hSystTotScaledNoFD,"Syst. from data","f");
	legendSystLcD0->AddEntry(hSystScaledFD,"Syst. from B feed-down","f");
	legendSystLcD0->Draw();

	legendLcD0->Draw();

	cLcD0->SaveAs(Form("%s/LcD0.pdf",outputDirName.Data()));
	// cLcD0->Close();

	//--------------------------------------------------------------------------------
	// Saving output

	TFile* fout = TFile::Open(Form("%s/Output.root",outputDirName.Data()),"RECREATE");
	hCorrYield->Write();
	hSystTot->Write();
	hSystNoFD->Write();
	hSystFD->Write();
	hSystTotScaledNoFD->Write();
	hSystScaledFD->Write();

	hYieldLcpKpi->Write();
	hSystTotLcpKpi->Write();
	hSystFDLcpKpi->Write();
	hSystNoFDLcpKpi->Write();

	hSystRatioNoFD->Write();
	hSystRatioFD->Write();

	fout->Close();


	return NULL;

}