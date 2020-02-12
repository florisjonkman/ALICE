#include "TFile.h"
#include "TString.h"
#include "TH1D.h"

const int nmult = 1;
TString mult[nmult] = {"1_999"};
TString fileNameEffTMVA = "../TMVAEfficiencies/TMVAEffOutput_PfCuts_RandomName/TMVAEff_PfCuts_RandomName.root";

void doLcEffMergingBunchesTMVA(TString cuts,Int_t num, Int_t den) {
       
    // TString fileNameAcc = "../../../../../../../MC_6Maggio/D0Acc_FONLL5mio.root";
    // TFile *fileacc = TFile::Open(fileNameAcc.Data());
    // TH1D *hAcc = (TH1D*)fileacc->Get("hAcc");

    Printf("TMVAEff file: %s",fileNameEffTMVA.Data());
    Printf("Mult : %s",mult[nmult].Data());

    TFile *fileEffTMVA = TFile::Open(fileNameEffTMVA.Data());
    TH1D *hEffTMVAPrompt = (TH1D*)fileEffTMVA->Get("hPtLcEffPrompt");
    TH1D *hEffTMVABfd = (TH1D*)fileEffTMVA->Get("hPtLcEffBfd");

    TFile *_file0;
    TH1D *hnum;
    TH1D *hnumB;
    TH1D *hden;
    TH1D *hdenB;
    
    TH1D *hNumer[nmult];
    TH1D *hNumerB[nmult];
    TH1D *hDenom[nmult];
    TH1D *hDenomB[nmult];
    TH1D *hEff[nmult];
    TH1D *hEffB[nmult];
    
    for(int ibunch=16; ibunch<=18; ibunch++) {
        _file0 = TFile::Open(Form("../CutEfficiencies/LcEff_%d_over_%d_%s_20%d.root",num,den,cuts.Data(),ibunch));
        _file0->ls(); 
        for(int jmult=0; jmult<nmult; jmult++) {
          cout<<  " nome mult " <<mult[jmult] <<endl;
            hnum  = (TH1D*)_file0->Get(Form("hNum_Mult_%s",mult[jmult].Data()));
            hnumB = (TH1D*)_file0->Get(Form("hNumB_Mult_%s",mult[jmult].Data()));
            hden  = (TH1D*)_file0->Get(Form("hDen_Mult_%s",mult[jmult].Data()));
            hdenB = (TH1D*)_file0->Get(Form("hDenB_Mult_%s",mult[jmult].Data()));
            
            hnum->Sumw2();
            hnumB->Sumw2();
            hden->Sumw2();
            hdenB->Sumw2();
            
            if(ibunch==16) {
                hNumer[jmult] = (TH1D*)hnum->Clone(Form("hNumMult%s",mult[jmult].Data()));
                hNumerB[jmult] = (TH1D*)hnumB->Clone(Form("hNumB%s",mult[jmult].Data()));
                hDenom[jmult] = (TH1D*)hden->Clone(Form("hDen%s",mult[jmult].Data()));
                hDenomB[jmult] = (TH1D*)hdenB->Clone(Form("hDenB%s",mult[jmult].Data()));
                hNumer[jmult]->Sumw2();
                hNumerB[jmult]->Sumw2();
                hDenom[jmult]->Sumw2();
                hDenomB[jmult]->Sumw2();
            }
            else {
                hNumer[jmult]->Add(hnum);
                hNumerB[jmult]->Add(hnumB);
                hDenom[jmult]->Add(hden);
                hDenomB[jmult]->Add(hdenB);
            }
            
        }
        
    }
    
    
    for(int jmult=0; jmult<nmult; jmult++) {
        hEff[jmult] = 0x0;
        hEffB[jmult] = 0x0;
        
        hEff[jmult]  = (TH1D*)hNumer[jmult]->Clone(Form("hEff_C_Mult%s",mult[jmult].Data()));
        hEffB[jmult] = (TH1D*)hNumerB[jmult]->Clone(Form("hEff_B_Mult%s",mult[jmult].Data()));
   
        hEff[jmult]->Divide(hNumer[jmult],hDenom[jmult],1,1,"B");
        hEffB[jmult]->Divide(hNumerB[jmult],hDenomB[jmult],1,1,"B");
        hEff[jmult]->SetName("hEffD");//Form("hEff_C_Mult%s",mult[jmult].Data()));
        hEffB[jmult]->SetName("hEffB");//Form("hEff_B_Mult%s",mult[jmult].Data()));

        for(int ibin=1;ibin<=hEff[jmult]->GetNbinsX();ibin++){
            double ptcenter = hEff[jmult]->GetBinCenter(ibin);
            // int accbin = hAcc->FindBin(ptcenter);

            int effTMVAbin = hEffTMVAPrompt->FindBin(ptcenter);

            // Prompt cut efficiency
            Double_t effPrompt        = hEff[jmult]->GetBinContent(ibin);
            Double_t effPromptErr     = hEff[jmult]->GetBinError(ibin);
            
            // Feed-down cut efficiency
            Double_t effBfd           = hEffB[jmult]->GetBinContent(ibin);
            Double_t effBfdErr        = hEffB[jmult]->GetBinError(ibin);
            
            // Prompt tmva efficiency
            Double_t effTMVAPrompt    = hEffTMVAPrompt->GetBinContent(ibin);
            Double_t effTMVAPromptErr = hEffTMVAPrompt->GetBinError(ibin);

            // Feed-down tmva efficiency
            Double_t effTMVABfd       = hEffTMVABfd->GetBinContent(ibin);
            Double_t effTMVABfdErr    = hEffTMVABfd->GetBinError(ibin);

            hEff[jmult]->SetBinContent(ibin,effPrompt*effTMVAPrompt*1./* hAcc->GetBinContent(accbin) */);
            hEff[jmult]->SetBinError  (ibin,TMath::Sqrt(effPromptErr*effPromptErr*effTMVAPrompt*effTMVAPrompt+effPrompt*effPrompt*effTMVAPromptErr*effTMVAPromptErr)*1./* hAcc->GetBinContent(accbin) */);

            hEffB[jmult]->SetBinContent(ibin,effBfd*effTMVABfd*1./* hAcc->GetBinContent(accbin) */);
            hEffB[jmult]->SetBinError  (ibin,TMath::Sqrt(effBfdErr*effBfdErr*effTMVABfd*effTMVABfd+effBfd*effBfd*effTMVABfdErr*effTMVABfdErr)*1./* hAcc->GetBinContent(accbin) */);
            
        }
        
        TFile* outfil=new TFile(Form("LcAccEff_%d_over_%d_%sTMVA_161718_wNtrklWeights_%s.root",num,den,cuts.Data(),mult[jmult].Data()),"recreate");
        hEff[jmult]->Write();
        hEffB[jmult]->Write();
        outfil->Close();

    }
    
    
}
