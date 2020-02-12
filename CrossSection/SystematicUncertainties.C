#include "AliHFSystErr.h"

Double_t GetTotalSystErr(AliHFSystErr* systemtic, Double_t pt, Double_t feeddownErr = 0, Double_t cutvarErr = 0){
  //
  // Get total syst error (except norm. error)
  //

  Double_t err=0.;

  err += systemtic->GetRawYieldErr(pt)*systemtic->GetRawYieldErr(pt);
  err += systemtic->GetTrackingEffErr(pt)*systemtic->GetTrackingEffErr(pt);
  //  if(fBR) err += GetBRErr()*GetBRErr();
  err += systemtic->GetCutsEffErr(pt)*systemtic->GetCutsEffErr(pt);
  err += systemtic->GetPIDEffErr(pt)*systemtic->GetPIDEffErr(pt);
  err += systemtic->GetMCPtShapeErr(pt)*systemtic->GetMCPtShapeErr(pt);
  // err += systemtic->GetPartAntipartErr(pt)*systemtic->GetPartAntipartErr(pt);

  if(feeddownErr>0.0) err += feeddownErr*feeddownErr;
  if(cutvarErr>0.0) err += cutvarErr*cutvarErr;

  return TMath::Sqrt(err);
}

TH1F* ReflectHisto(TH1F *hin) {
  //
  // Clones and reflects histogram
  //
  TH1F *hout=(TH1F*)hin->Clone("hout");
  hout->Scale(-1.);

  return hout;
}

void SystematicUncertainties(){

	// Bins: pT
	const Int_t nPtBins = 5;
	Double_t ptBins[nPtBins+1] =  {  1,  2,  4,  6,  8, 12}; // pT bin ranges

	// Input
	TFile* input = TFile::Open("OutputHFPtSpectrum_PfCutsTMVA_FixSigma_1_999.root");
	if(!input) {Printf("input not found!"); return;}
	AliHFSystErr* systemtic = (AliHFSystErr*) input->Get("AliHFSystErr");
	if(!systemtic) {Printf("systemtic not found!"); return;}

	// Total systemtic output
	TH1F *hTotErr= new TH1F("hTotErr","hTotErr",nPtBins,ptBins);
	TGraphAsymmErrors *gTotErr = new TGraphAsymmErrors(0);
	gTotErr->SetNameTitle("gTotErr","gTotErr");

  TH1F *hTotErrNoFD= new TH1F("hTotErr_noFD","hTotErr_noFD",nPtBins,ptBins);
  TGraphAsymmErrors *gTotErrNoFD = new TGraphAsymmErrors(0);
  gTotErrNoFD->SetNameTitle("gTotErrNoFD","gTotErrNoFD");

	// Input different histograms
	TGraphAsymmErrors *gFcConservative =  (TGraphAsymmErrors*) input->Get("gFcConservative");
	gFcConservative->SetNameTitle("gFcConservative","gFcConservative");

	TH1F *fNorm	        = new TH1F("fNorm","fNorm",nPtBins,ptBins);
	TH1F *fTrackingEff	= new TH1F("fTrackingEff","fTrackingEff",nPtBins,ptBins);     /// tracking efficiency
	TH1F *fRawYield	    = new TH1F("fRawYield","fRawYield",nPtBins,ptBins);        /// raw yield
	TH1F *fBR	        = new TH1F("fBR","fBR",nPtBins,ptBins);              /// branching ratio
	// TH1F *fCutsEff	    = new TH1F("fCutsEff","fCutsEff",nPtBins,ptBins);         /// cuts efficiency
	TH1F *fPIDEff	    = new TH1F("fPIDEff","fPIDEff",nPtBins,ptBins);          /// PID efficiency
	TH1F *fMCPtShape	= new TH1F("fMCPtShape","fMCPtShape",nPtBins,ptBins);       /// MC dNdpt

	TGraphAsymmErrors *grErrFeeddown = new TGraphAsymmErrors(0);
	grErrFeeddown->SetNameTitle("grErrFeeddown","grErrFeeddown");

	for (Int_t bin = 1; bin <= nPtBins; ++bin)
	{
		Double_t pt = (ptBins[bin-1]+ptBins[bin])/2.;
		// Printf("bin = %d, norm = %.3f",bin,systemtic->GetNormErr());
		fNorm->SetBinContent(bin,0.05);
		fTrackingEff->SetBinContent(bin,systemtic->GetTrackingEffErr(pt));
		fRawYield->SetBinContent(bin,systemtic->GetRawYieldErr(pt));
		fBR->SetBinContent(bin,0.051);
		fPIDEff->SetBinContent(bin,systemtic->GetPIDEffErr(pt));
		fMCPtShape->SetBinContent(bin,systemtic->GetMCPtShapeErr(pt));
	}


	// Asymmetric Cutvariation!
	Double_t cutVarX       [6] = {  1.5,    3,    5,    7,   10,   18};
	Double_t cutVarErrX    [6] = {  0.5,    1,    1,    1,    2,    6};

	Double_t cutVarY       [6] = {    0,    0,    0,    0,    0,   0};
	Double_t cutVarErrYLow [6] = {  .08,  .08,  .08,  .08,  .08,  .00};
	Double_t cutVarErrYHigh[6] = {  .04,  .04,  .04,  .04,  .04,  .00};
	TGraphAsymmErrors *gCutVariation =  new TGraphAsymmErrors(nPtBins,cutVarX,cutVarY,cutVarErrX,cutVarErrX,cutVarErrYLow,cutVarErrYHigh);
	gCutVariation->SetNameTitle("gCutVariation","gCutVariation");

	Int_t nbins = hTotErr->GetNbinsX();
	Printf("%2.2s  %5.4s  %5.4s  %5.4s  %9.9s  %9.9s  %9.9s  %9.9s  %9.9s  %9.9s","i","pTLow","pTUp","Middle","erryl","erryh","errylCut","erryhCut","toterryl","toterryh");
	for(Int_t i=1;i<=nbins;i++) {
    Double_t pt = hTotErr->GetBinCenter(i);
    Double_t ptwidth = hTotErr->GetBinWidth(i);
    
    if(gFcConservative) {
      Double_t toterryl=0., toterryh=0.;
      Double_t toterrylNoFD=0., toterryhNoFD=0.;

      Double_t x=0., y=0., errxl=0., errxh=0., erryl=0., erryh=0.;
      for(Int_t j=0; j<gFcConservative->GetN(); j++) {
        gFcConservative->GetPoint(j,x,y);
        errxh = gFcConservative->GetErrorXhigh(j);
        errxl = gFcConservative->GetErrorXlow(j);
        if ( ( (x-errxl) <= pt) && ( (x+errxh) >= pt) && y>0. ) {
          erryh = gFcConservative->GetErrorYhigh(j) / y;
          erryl = gFcConservative->GetErrorYlow(j) / y;
          grErrFeeddown->SetPoint(j,x,0.);
          grErrFeeddown->SetPointError(i,errxl,errxh,erryl,erryh);
          Printf("pt = %.2f, erryl = %.3f, erryh = %.3f",pt,erryl,erryh);
        }
        
      }
      Double_t xCut=0., yCut=0., errxlCut=0., errxhCut=0., errylCut=0., erryhCut=0.;
      for(Int_t j=0; j<gCutVariation->GetN(); j++) {
        gCutVariation->GetPoint(j,xCut,yCut);
        errxhCut = gCutVariation->GetErrorXhigh(j);
        errxlCut = gCutVariation->GetErrorXlow(j);
        if ( ( (xCut-errxlCut) <= pt) && ( (xCut+errxhCut) >= pt) ) {
          erryhCut = gCutVariation->GetErrorYhigh(j) ;
          errylCut = gCutVariation->GetErrorYlow(j);
        }
      }

      if (erryl>=1e-3) toterryl = GetTotalSystErr(systemtic,pt,erryl,errylCut);
      else toterryl = GetTotalSystErr(systemtic,pt);
      if (erryh>=1e-3) toterryh = GetTotalSystErr(systemtic,pt,erryh,erryhCut);
      else toterryh = GetTotalSystErr(systemtic,pt);

      toterrylNoFD = GetTotalSystErr(systemtic,pt,0,errylCut);
      toterryhNoFD = GetTotalSystErr(systemtic,pt,0,erryhCut);

      hTotErr->SetBinContent(i,toterryh);
      gTotErr->SetPoint(i-1,pt,0.);
      gTotErr->SetPointError(i-1,0.3/*ptwidth/2.*/,0.3  /*ptwidth/2.*/,toterryl,toterryh); // i, exl, exh, eyl, eyh
      gTotErrNoFD->SetPoint(i-1,pt,0.);
      gTotErrNoFD->SetPointError(i-1,0.3/*ptwidth/2.*/,0.3  /*ptwidth/2.*/,toterrylNoFD,toterryhNoFD); // i, exl, exh, eyl, eyh

      Printf("%2.1d  %5.2f  %5.2f  %5.2f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f  %9.3f",i,pt - ptwidth/2.,pt + ptwidth/2.,pt,erryl,erryh,errylCut,erryhCut,toterryl,toterryh);
    }
  }

  gStyle->SetOptStat(0);

  TCanvas *cSystErr = new TCanvas("cSystErr","Systematic Errors",300,80,1000,600);
  cSystErr->Range(0.20,-0.5,18.4,0.34);
  cSystErr->SetRightMargin(0.318);
  cSystErr->SetFillColor(0);

  TH2F *hFrame = new TH2F("hFrame","Systematic errors; p_{T} (GeV/c); Relative Syst. Unc.",50,0,50,100,-1,+1);
  hFrame->SetAxisRange(1.,11.9,"X");
  hFrame->SetAxisRange(-0.23,0.23,"Y");
  hFrame->Draw();

  TLegend *leg = new TLegend(0.69,0.44,0.98,0.86,NULL,"brNDC");
  leg->SetTextSize(0.03601695);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);

  gTotErr->SetLineColor(kBlack);
  gTotErr->SetFillColor(kRed);
  gTotErr->SetFillStyle(3002);
  gTotErr->Draw("2 same");
  leg->AddEntry(gTotErr,"Total (excl. norm. and BR)","f");

  fNorm->SetFillColor(1);
  fNorm->SetFillStyle(3002);
  //fNorm->Draw("same");
  //TH1F *hNormRefl = ReflectHisto(fNorm);
  //hNormRefl->Draw("same");
  Double_t norm = fNorm->GetBinContent(1)*100;
  leg->AddEntry(fNorm,Form("Normalization (#pm %.1f%s)",norm,"%"),"");
  if(fBR) {
    Double_t brsys=fBR->GetBinContent(1)*100;
    leg->AddEntry(fBR,Form("Branching ratio(#pm %.1f%s)",0.051*100,"%"),"");
  }
  if(fRawYield) {
    Int_t ci;   // for color index setting
    ci = TColor::GetColor("#00cc00");
    fRawYield->SetLineColor(ci);
    //    fRawYield->SetLineColor(3);
    fRawYield->SetLineWidth(3);
    fRawYield->Draw("same hist ");
    TH1F *hRawYieldRefl = ReflectHisto(fRawYield);
    hRawYieldRefl->Draw("same hist ");
    leg->AddEntry(fRawYield,"Yield extraction","l");
  }
  if(fTrackingEff) {
    fTrackingEff->SetFillColor(4);
    fTrackingEff->SetFillStyle(3005);
    fTrackingEff->Draw("same hist ");
    TH1F *hTrackingEffRefl = ReflectHisto(fTrackingEff);
    hTrackingEffRefl->Draw("same hist ");
    leg->AddEntry(fTrackingEff,"Tracking efficiency","f");
  }
  if(gCutVariation) {
  	gCutVariation->SetFillColorAlpha(kGray+1,0.4);
    // gCutVariation->SetFillStyle(3006);
    gCutVariation->Draw("2 same ");
    leg->AddEntry(gCutVariation,"Cut efficiency","f");

    // fCutsEff->SetLineColor(4);
    // fCutsEff->SetLineWidth(3);
    // fCutsEff->Draw("same hist ");
    // TH1F *hCutsEffRefl = ReflectHisto(fCutsEff);
    // hCutsEffRefl->Draw("same hist ");
    // leg->AddEntry(fCutsEff,"Cut efficiency","l");
  }
  if(fPIDEff) {
    fPIDEff->SetLineColor(7);
    fPIDEff->SetLineWidth(3);
    fPIDEff->Draw("same hist ");
    TH1F *hPIDEffRefl = ReflectHisto(fPIDEff);
    hPIDEffRefl->Draw("same hist ");
    leg->AddEntry(fPIDEff,"PID efficiency","l");
  }
  if(fMCPtShape) {
    Int_t ci = TColor::GetColor("#9933ff");
    fMCPtShape->SetLineColor(ci);
    //    fMCPtShape->SetLineColor(8);
    fMCPtShape->SetLineWidth(3);
    fMCPtShape->Draw("same hist ");
    TH1F *hMCPtShapeRefl = ReflectHisto(fMCPtShape);
    hMCPtShapeRefl->Draw("same hist ");
    leg->AddEntry(fMCPtShape,"MC p_{t} shape","l");
  }
  if(grErrFeeddown) {
    grErrFeeddown->SetFillColor(kTeal-8);
    grErrFeeddown->SetFillStyle(3001);
    grErrFeeddown->Draw("2 same ");
    leg->AddEntry(grErrFeeddown,"Feed-down from B","f");
  }

  leg->Draw();

  cSystErr->SaveAs("SystUncertainties.pdf");

  Printf("");
  for(Int_t j=0; j<gTotErr->GetN(); j++) {
        Double_t x = 0, y = 0, erryh = 0, erryl = 0;
        gTotErr->GetPoint(j,x,y);
        erryh = gTotErr->GetErrorYhigh(j);
        erryl = gTotErr->GetErrorYlow(j);
        Printf("x = %.1f, y = %.1f, erryl = %.3f, erryh = %.3f",x,y,erryl,erryh);
  }

  TFile* fout = new TFile("ResultsSystUncertainties.root","RECREATE");
  fNorm->Write();
  fTrackingEff->Write();
  gCutVariation->Write();
  fRawYield->Write();
  fBR->Write();
  fPIDEff->Write();
  fMCPtShape->Write();
  grErrFeeddown->Write();
  gTotErr->Write();
  gTotErrNoFD->Write();
  fout->Close();

}