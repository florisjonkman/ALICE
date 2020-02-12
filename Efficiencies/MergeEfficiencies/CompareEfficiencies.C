#include "TFile.h"
#include "TString.h"
#include "TH1D.h"

// Int_t   colors[4] = {kBlack,kBlue+1,kRed+1,kOrange+1};// {kBlue+1,kBlue+1,kRed+1,kRed+1}; // {kBlack,kBlue+1,kGreen+2,kRed+1};
// Int_t   styles[4] = {20,21,22,23}; //{20,24,21,25}; // {20,21,22,34};

// Int_t   colors[4] = {kBlack,kBlack,kRed+1,kRed+1};// {kBlue+1,kBlue+1,kRed+1,kRed+1}; // {kBlack,kBlue+1,kGreen+2,kRed+1};
// Int_t   styles[4] = {20,24,21,25}; //{20,24,21,25}; // {20,21,22,34};

Int_t   colors[5] = {kBlack,kBlue+1,kCyan,kRed-9,kRed+1};
Int_t   styles[5]= {20,22,22,23,23};
TString titles[10];
TString filesNames[10];
TFile* files[10];

TH1D* effPrompt[10];
TH1D* effBfd[10];

void CompareEfficienciesList(TString suffix) {
    
    // gROOT->SetStyle("Plain");
    // gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    // // gStyle->SetPalette(1);
    // gStyle->SetCanvasColor(0);
    // gStyle->SetFrameFillColor(0);

    // const Int_t nFiles = 2;

    // files[0] = TFile::Open("NoPCts/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999.root");
    // files[1] = TFile::Open("NoPCts/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999_FONLL.root");
   
    // titles[0] = "PYTHIA";
    // titles[1] = "FONLL";
    
    const Int_t nFiles = 5;

    files[0] = TFile::Open("NoPCts/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999.root");
    files[1] = TFile::Open("NoPCts/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999_loose2.root");
    files[2] = TFile::Open("NoPCts/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999_loose1.root");
    files[3] = TFile::Open("NoPCts/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999_tight1.root");
    files[4] = TFile::Open("NoPCts/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999_tight2.root");

    titles[0] = "Pf. + MVA (Central), 1 - 999";
    titles[1] = "Loose2";
    titles[2] = "Loose1";
    titles[3] = "Tight1";
    titles[4] = "Tight2";

    // const Int_t nFiles = 3;

    // files[0] = TFile::Open("NoPCts/LcAccEff_9_over_0_stdCuts_161718_wNtrklWeights_1_999.root");
    // files[1] = TFile::Open("NoPCts/LcAccEff_9_over_0_stdCutsTMVA_161718_wNtrklWeights_1_999.root");
    // files[2] = TFile::Open("NoPCts/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999.root");
  
    // titles[0] = "Std., 1 - 999";
    // titles[1] = "Std. + MVA, 1 - 999";
    // titles[2] = "Pf. + MVA, 1 - 999";

    TCanvas* cEffPrompt = new TCanvas("cEffPrompt","cEffPrompt",1200,600);
    cEffPrompt->Divide(2,1);

    TCanvas* cEffBfd = new TCanvas("cEffBfd","cEffBfd",1200,600);
    cEffBfd->Divide(2,1);

    TCanvas* cEffPromptBfd = new TCanvas("cEffPromptBfd","cEffPromptBfd",1200,600);
    cEffPromptBfd->Divide(2,1,0.00001,0.00001);

    TCanvas* cEffPromptBfdRatio = new TCanvas("cEffPromptBfdRatio","cEffPromptBfdRatio",1200,400);
    cEffPromptBfdRatio->Divide(2,1,0.00001,0.00001);

    TH1D* effPromptDen;
    TH1D* effPromptRatio;
    TH1D* effBfdDen;
    TH1D* effBfdRatio;

    TLegend* legEffPrompt = new TLegend(0.128721 ,  0.734658, 0.441743  , 0.880431);

    for (Int_t i = 0; i < nFiles; ++i)
    {
    	if(!files[i]){
    		Printf("File %d not found",i);
    		return;
    	}

    	effPrompt[i] = (TH1D*) files[i]->Get("hEffD");
    	effBfd[i]    = (TH1D*) files[i]->Get("hEffB");
    
    	cEffPrompt->cd(1);
    	effPrompt[i]->GetYaxis()->SetTitle("Acc #times Eff");
    	effPrompt[i]->GetXaxis()->SetTitle("#it{p}_{T}");
    	effPrompt[i]->SetLineColor(colors[i]);
    	effPrompt[i]->SetMarkerColor(colors[i]);
    	effPrompt[i]->SetMarkerStyle(styles[i]);
        effPrompt[i]->GetYaxis()->SetRangeUser(0,0.9);
    	effPrompt[i]->SetTitle("#Lambda_{c}-prompt");
    	if(i == 0 ) effPrompt[i]->Draw("PL");
    	else        effPrompt[i]->Draw("SAME");
        cEffPromptBfd->cd(1);
        if(i == 0 ) effPrompt[i]->Draw("PL");
        else        effPrompt[i]->Draw("SAME");

    	legEffPrompt->AddEntry(effPrompt[i],titles[i],"PL");

    	cEffPrompt->cd(2);
    	if(i==0) effPromptDen = (TH1D*) effPrompt[i]->Clone(Form("effPromptDen%d",i));
    	effPromptRatio = (TH1D*) effPrompt[i]->Clone(Form("effPromptRatio%d",i));
        // effPromptRatio->Divide(effPromptDen); // Uncorrelated
        effPromptRatio->Divide(effPromptRatio,effPromptDen,1,1,"B"); // Correlated
        effPromptRatio->SetTitle(Form("Acc #times Eff / ( %s )",titles[0].Data()));
    	effPromptRatio->GetYaxis()->SetTitle(Form("Acc #times Eff / ( %s )",titles[0].Data()));
        // effPromptRatio->GetYaxis()->SetRangeUser(0.8,1.2);
        effPromptRatio->GetYaxis()->SetRangeUser(0.0,2.0);
    	effPromptRatio->GetXaxis()->SetTitle("#it{p}_{T}");
    	effPromptRatio->SetLineColor(colors[i]);
    	effPromptRatio->SetMarkerColor(colors[i]);
    	effPromptRatio->SetMarkerStyle(styles[i]);
    	if(i == 0 ) effPromptRatio->Draw("PL");
    	else        effPromptRatio->Draw("SAME");
        cEffPromptBfdRatio->cd(1);
       if(i == 0 ){
            effPromptRatio->GetYaxis()->SetTitle("Ratio");
            effPromptRatio->Draw("PL");
        }
        else        effPromptRatio->Draw("SAME");

        cEffBfd->cd(1);
        effBfd[i]->GetYaxis()->SetTitle("Acc #times Eff");
        effBfd[i]->GetXaxis()->SetTitle("#it{p}_{T}");
        effBfd[i]->SetLineColor(colors[i]);
        effBfd[i]->SetMarkerColor(colors[i]);
        effBfd[i]->SetMarkerStyle(styles[i]);
        effBfd[i]->GetYaxis()->SetRangeUser(0,0.9);
        effBfd[i]->SetTitle("#Lambda_{c}-Bfd");
        if(i == 0 ) effBfd[i]->Draw("PL");
        else        effBfd[i]->Draw("SAME");
        cEffPromptBfd->cd(2);
        if(i == 0 ) effBfd[i]->Draw("PL");
        else        effBfd[i]->Draw("SAME");

        cEffBfd->cd(2);
        if(i==0) effBfdDen = (TH1D*) effBfd[i]->Clone(Form("effBfdDen%d",i));
        effBfdRatio = (TH1D*) effBfd[i]->Clone(Form("effBfdRatio%d",i));
        // effBfdRatio->Divide(effBfdDen); // Uncorrelated
        effBfdRatio->Divide(effBfdRatio,effBfdDen,1,1,"B"); // Correlated
        effBfdRatio->SetTitle(Form("Acc #times Eff / ( %s )",titles[0].Data()));
        effBfdRatio->GetYaxis()->SetTitle(Form("Acc #times Acc #times Eff / ( %s )",titles[0].Data()));
        effBfdRatio->GetXaxis()->SetTitle("#it{p}_{T}");
        // effBfdRatio->GetYaxis()->SetRangeUser(0.8,1.2);
        effBfdRatio->GetYaxis()->SetRangeUser(0.0,2.0);
        effBfdRatio->SetLineColor(colors[i]);
        effBfdRatio->SetMarkerColor(colors[i]);
        effBfdRatio->SetMarkerStyle(styles[i]);
        if(i == 0 ) effBfdRatio->Draw("PL");
        else        effBfdRatio->Draw("SAME");
        cEffPromptBfdRatio->cd(2);
        if(i == 0 ){
            effBfdRatio->GetYaxis()->SetTitle("Ratio");
            effBfdRatio->Draw("PL");
        }
        else  effBfdRatio->Draw("SAME");
    }

    cEffPrompt->cd(1);
    legEffPrompt->Draw();

    cEffBfd->cd(1);
    legEffPrompt->Draw();

    cEffPromptBfd->cd(1);
    legEffPrompt->Draw();

    TString outputDirName = Form("CompareEff_%s",suffix.Data());
    // Open output file 
    gSystem->Exec(Form("rm -r %s/",outputDirName.Data()));
    gSystem->Exec(Form("mkdir %s",outputDirName.Data()));

    cEffPrompt->SaveAs(Form("%s/EffPrompt.pdf",outputDirName.Data()));
    cEffBfd->SaveAs(Form("%s/EffBfd.pdf",outputDirName.Data()));
    cEffPromptBfd->SaveAs(Form("%s/EffPromptBfd.pdf",outputDirName.Data()));
    cEffPromptBfdRatio->SaveAs(Form("%s/EffPromptBfdRatio.pdf",outputDirName.Data()));
}

void MultiplicityPerPt(TString suffix) {
    
    // gROOT->SetStyle("Plain");
    // gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    // // gStyle->SetPalette(1);
    // gStyle->SetCanvasColor(0);
    // gStyle->SetFrameFillColor(0);

    
    const Int_t nFiles = 4;

    files[0] = TFile::Open("CorrectTPCCuts/LcAccEff_9_over_0_stdCuts_161718_wNtrklWeights_1_9.root");
    files[1] = TFile::Open("CorrectTPCCuts/LcAccEff_9_over_0_stdCuts_161718_wNtrklWeights_10_29.root");
    files[2] = TFile::Open("CorrectTPCCuts/LcAccEff_9_over_0_stdCuts_161718_wNtrklWeights_30_59.root");
    files[3] = TFile::Open("CorrectTPCCuts/LcAccEff_9_over_0_stdCuts_161718_wNtrklWeights_1_999.root");
  
    titles[0] = "Std., 1 - 9";
    titles[1] = "Std., 10 - 29";
    titles[2] = "Std., 30 - 59";
    titles[3] = "Std., 1 - 999";

    // Bins: pT
    const Int_t nPtBins = 6;
    Double_t ptBins[nPtBins+1] =  {  1,      2,      4,      6,     8,    12,    24}; // pT bin ranges
    Int_t colorsPt[nPtBins]    =  {  kRed+1,kBlue+1,kGreen+2,kOrange+2,kViolet+2,kGray+2};

    const Int_t nMultBins = 4;
    // Double_t multBins[nMultBins+1] =  {1,10,30,60,100}; 
    TString label[nMultBins] = {"1 - 9", "10 - 29", " 30 - 59", "Mult. Int."};

    TH1D* hEffCVsMult[nPtBins];

    TLegend* legEffMult=new TLegend(0.125418,0.662021, 0.438127,0.890244 );
    legEffMult->SetFillStyle(0);
    legEffMult->SetTextFont(42);
    legEffMult->SetBorderSize(0);

    TCanvas* cEffMultPrompt=new TCanvas("cEffMultPrompt","cEffMultPrompt",600,600);
    cEffMultPrompt->cd();

    for (Int_t i = 0; i < nFiles; ++i)
    {
        if(!files[i]){
            Printf("File %d not found",i);
            return;
        }

        effPrompt[i] = (TH1D*) files[i]->Get("hEffD");
        effBfd[i]    = (TH1D*) files[i]->Get("hEffB");
    }
    
    for(Int_t iPtBin = 0; iPtBin < nPtBins; iPtBin++) {
                
        hEffCVsMult[iPtBin] = new TH1D(Form("hEffPromptVsMult_pT_%0.f_%0.f",ptBins[iPtBin],ptBins[iPtBin+1]),Form("hEffPromptVsMult_pT_%0.f_%0.f",ptBins[iPtBin],ptBins[iPtBin+1]),4,-0.5,3.5);
        for (Int_t ibin = 1; ibin <= hEffCVsMult[iPtBin]->GetNbinsX(); ibin++){
            hEffCVsMult[iPtBin]->GetXaxis()->SetBinLabel(ibin, label[ibin-1].Data());
        }

        hEffCVsMult[iPtBin]->SetBinContent(1,effPrompt[0]->GetBinContent(iPtBin+1));
        hEffCVsMult[iPtBin]->SetBinContent(2,effPrompt[1]->GetBinContent(iPtBin+1));
        hEffCVsMult[iPtBin]->SetBinContent(3,effPrompt[2]->GetBinContent(iPtBin+1));
        hEffCVsMult[iPtBin]->SetBinContent(4,effPrompt[3]->GetBinContent(iPtBin+1));

        hEffCVsMult[iPtBin]->SetBinError(1,effPrompt[0]->GetBinError(iPtBin+1));
        hEffCVsMult[iPtBin]->SetBinError(2,effPrompt[1]->GetBinError(iPtBin+1));
        hEffCVsMult[iPtBin]->SetBinError(3,effPrompt[2]->GetBinError(iPtBin+1));
        hEffCVsMult[iPtBin]->SetBinError(4,effPrompt[2]->GetBinError(iPtBin+1));

        hEffCVsMult[iPtBin]->SetMarkerColor(colorsPt[iPtBin]);
        hEffCVsMult[iPtBin]->SetLineColor(colorsPt[iPtBin]);
        hEffCVsMult[iPtBin]->SetMarkerStyle(39+iPtBin*2);
        // hEffCVsMult[iPtBin]->GetXaxis()->SetTitle("Multiplicity");
        hEffCVsMult[iPtBin]->SetTitle("");
        hEffCVsMult[iPtBin]->GetYaxis()->SetTitle("Prompt #Lambda_{c} (Eff #times Acc)");
        hEffCVsMult[iPtBin]->GetYaxis()->SetRangeUser(0.1,1.0);

        if(iPtBin > 0) hEffCVsMult[iPtBin]->Draw("Same");
        else           hEffCVsMult[iPtBin]->Draw();

        legEffMult->AddEntry(hEffCVsMult[iPtBin],Form("%.0f < #it{p}_{T} < %0.f",ptBins[iPtBin],ptBins[iPtBin+1]),"P");
    }
    TLine* line = new TLine(2.5,0.1,2.5,1.0);
    line->SetLineColorAlpha(kBlack,1.0); line->SetLineWidth(1); line->SetLineStyle(2);
    line->Draw("Same");

    legEffMult->Draw();
    cEffMultPrompt->SaveAs(Form("Output/EffMultPrompt_%s.pdf",suffix.Data()));
}

void DrawOneEfficiency(TString filename){

    // gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);

    TFile* file = TFile::Open("NoPCts/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999.root");

    TCanvas* cEfficiency = new TCanvas("cEfficiency","cEfficiency",600,600);
    TLegend* leg = new TLegend(0.583612 , 0.137631, 0.874582,  0.278746);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    TH1D* effP = (TH1D*) file->Get("hEffD");
    TH1D* effB = (TH1D*) file->Get("hEffB");

    cEfficiency->cd();
    effP->GetYaxis()->SetTitle("Acceptance #times Efficiency");
    effP->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    effP->SetLineColor(kRed);
    effP->SetMarkerColor(kRed);
    effP->SetMarkerStyle(20);

    // effP->GetYaxis()->SetRangeUser(0.001,1);
    // gPad->SetLogy(1);
    effP->GetYaxis()->SetRangeUser(0.001,0.8);

    effP->SetTitle("");
    effP->Draw("HIST");
    TH1D* effP2 = (TH1D*) effP->Clone("effP2");
    effP2->Draw("PSAME");

    leg->AddEntry(effP,"Prompt #Lambda_{c}","PL");

    cEfficiency->cd();
    effB->GetYaxis()->SetTitle("Acceptance #times Efficiency");
    effB->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    effB->SetLineColor(kBlue);
    effB->SetLineStyle(3);
    effB->SetMarkerColor(kBlue);
    effB->SetMarkerStyle(25);

    // effB->GetYaxis()->SetRangeUser(0.001,1);
    // gPad->SetLogy(1);
    effB->GetYaxis()->SetRangeUser(0.001,0.8);

    effB->SetTitle("");
    effB->Draw("HISTSAME");
    TH1D* effB2 = (TH1D*) effB->Clone("effB2");
    effB2->Draw("PSAME");

    leg->AddEntry(effB,"Feed-down #Lambda_{c}","PL");
    
    leg->Draw();

    TPaveText *prodInfo=new TPaveText(0.588629, 0.280488 ,0.879599 ,0.343206,"NDC");
    prodInfo->SetBorderSize(0);
    prodInfo->SetFillStyle(0);
    prodInfo->AddText(Form("#Lambda_{c}#rightarrowpK_{S}^{0}, Prefilt."));
    prodInfo->SetTextFont(42);
    prodInfo->SetTextAlign(12);
    prodInfo->SetTextColor(kBlack);
    prodInfo->Draw();

    TPaveText *prodInfo2=new TPaveText(0.123746, 0.829268 ,0.386288,0.890244,"NDC");
    prodInfo2->SetBorderSize(0);
    prodInfo2->SetFillStyle(0);
    prodInfo2->AddText(Form("pp, #sqrt{#it{s}} = 13 TeV"));
    prodInfo2->SetTextFont(42);
    prodInfo2->SetTextAlign(12);
    prodInfo2->SetTextColor(kBlack);
    prodInfo2->Draw();

    cEfficiency->SaveAs(Form("%s.pdf",filename.Data()));

    for (int ibin = 1; ibin <= effB->GetNbinsX(); ++ibin)
    {
        Double_t valueP = effP->GetBinContent(ibin);
        Double_t valueB = effB->GetBinContent(ibin); 
        Double_t valuePerr = effP->GetBinError(ibin);
        Double_t valueBerr = effB->GetBinError(ibin); 

        Printf("ibin = %d, eff prompt = %.3f,rel err prompt = %.3f, eff feed-down = %.3f, rel err feed-down = %.3f",ibin,valueP,valuePerr/valueP,valueB,valueBerr/valueB);
    }


    return;

}

void DrawThreeEfficiencies(TString filename){

    // gROOT->SetStyle("Plain");
    gStyle->SetPalette(1);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(1);
    gStyle->SetCanvasColor(0);
    gStyle->SetFrameFillColor(0);

    TFile* files[3];
    TString titles[3];

    files[0] = TFile::Open("NoPCts/LcAccEff_9_over_0_stdCuts_161718_wNtrklWeights_1_999.root");
    files[1] = TFile::Open("NoPCts/LcAccEff_9_over_0_stdCutsTMVA_161718_wNtrklWeights_1_999.root");
    files[2] = TFile::Open("NoPCts/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999.root");
    titles[0] = "Standard";
    titles[1] = "Standard + MVA";
    titles[2] = "Prefiltering + MVA";

    TCanvas* cEfficiency = new TCanvas("cEfficiency","cEfficiency",1200,400);
    cEfficiency->Divide(3,1,0.001,0);
    
    TLegend* leg = new TLegend(0.689812 , 0.187766, 0.981623,  0.329075);
    leg->SetTextFont(42);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    for (Int_t i = 0; i < 3; ++i)
    {
        TH1D* effP = (TH1D*) files[i]->Get("hEffD");
        TH1D* effB = (TH1D*) files[i]->Get("hEffB");

        cEfficiency->cd(i+1);
        // gPad->SetLogy(1);
        effP->GetYaxis()->SetTitle("Acceptance #times Efficiency");
        effP->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        effP->SetLineColor(kRed);
        effP->SetMarkerColor(kRed);
        effP->SetMarkerStyle(20);

        effP->GetYaxis()->SetRangeUser(0.06,0.62);
        effP->GetXaxis()->SetNdivisions(22);
        effP->GetXaxis()->ChangeLabel(12,-1,-1,4,-1,-1);

        effP->SetTitle("");
        effP->Draw("HIST");
        TH1D* effP2 = (TH1D*) effP->Clone("effP2");
        effP2->Draw("PSAME");

        if(i==0) leg->AddEntry(effP,"Prompt #Lambda_{c}","PL");  

        effB->GetYaxis()->SetTitle("Acceptance #times Efficiency");
        effB->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
        effB->SetLineColor(kBlue);
        effB->SetLineStyle(3);
        effB->SetMarkerColor(kBlue);
        effB->SetMarkerStyle(25);

        // effB->GetYaxis()->SetRangeUser(0.001,1);
        
        effB->GetYaxis()->SetRangeUser(0.06,0.62);
        // effB->GetXaxis()->SetRangeUser(1,24.1);

        effB->SetTitle("");
        effB->Draw("HISTSAME");
        TH1D* effB2 = (TH1D*) effB->Clone("effB2");
        effB2->Draw("PSAME");

        if(i==0) leg->AddEntry(effB,"Feed-down #Lambda_{c}","PL");

        TPaveText *prodInfo=new TPaveText(0.69161, 0.123887 ,0.982124 ,0.187766,"NDC");
        prodInfo->SetBorderSize(0);
        prodInfo->SetFillStyle(0);
        prodInfo->AddText(titles[i].Data());
        prodInfo->SetTextFont(42);
        prodInfo->SetTextAlign(12);
        prodInfo->SetTextColor(kBlack);
        prodInfo->Draw();

    }

    cEfficiency->cd(1);
    
    leg->Draw();


    TPaveText *prodInfo2=new TPaveText(0.145042, 0.917538, 0.407499,0.977545,"NDC");
    prodInfo2->SetBorderSize(0);
    prodInfo2->SetFillStyle(0);
    prodInfo2->AddText(Form("pp, #sqrt{#it{s}} = 13 TeV"));
    prodInfo2->SetTextFont(42);
    prodInfo2->SetTextAlign(31);
    prodInfo2->SetTextColor(kBlack);
    prodInfo2->Draw();

    TPaveText *prodInfo3=new TPaveText(0.0483473, 0.859466 , 0.302171 ,0.913666,"NDC");
    prodInfo3->SetBorderSize(0);
    prodInfo3->SetFillStyle(0);
    prodInfo3->AddText(Form("#Lambda_{c}#rightarrowpK_{S}^{0}"));
    prodInfo3->SetTextFont(42);
    prodInfo3->SetTextAlign(31);
    prodInfo3->SetTextColor(kBlack);
    prodInfo3->Draw();

    cEfficiency->SaveAs(Form("%s.pdf",filename.Data()));

    // for (int ibin = 1; ibin <= effB->GetNbinsX(); ++ibin)
    // {
    //     Double_t valueP = effP->GetBinContent(ibin);
    //     Double_t valueB = effB->GetBinContent(ibin); 
    //     Double_t valuePerr = effP->GetBinError(ibin);
    //     Double_t valueBerr = effB->GetBinError(ibin); 

    //     Printf("ibin = %d, eff prompt = %.3f,rel err prompt = %.3f, eff feed-down = %.3f, rel err feed-down = %.3f",ibin,valueP,valuePerr/valueP,valueB,valueBerr/valueB);
    // }


    return;

}

    