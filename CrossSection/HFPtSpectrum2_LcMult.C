#if !defined(__CINT__) || defined(__MAKECINT__)
#include <Riostream.h>
#include "TH1D.h"
#include "TH1.h"
#include "TH2F.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TSystem.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "AliHFSystErr.h"
#endif

//
// Macro to use mimic AliHFPtSpectrum class using HFPtSpectrum file as input
// used for high multiplicity analyses where fPrompt calculation with Nb 
// method fails due to multiplicity-independent FONLL prediction.
//
// Besides using the fPrompt fraction from the input file, also the syst.
// errors (AliHFSystErr class) is copied.
//
// NB: So far only works for pp, using Nb method.
//

enum decay { kD0Kpi, kDplusKpipi, kDstarD0pi, kDsKKpi, kLctopKpi, kLcK0Sp};
enum centrality{ kpp8, kpp7, kpp5, kpp276, k07half, kpPb0100, k010, k1020, k020, k1030, k2040, k2030, k3040, k4050, k3050, k5060, k4060, k6080, k4080, k5080, k80100, kpPb010, kpPb020, kpPb1020, kpPb2040, kpPb4060, kpPb60100, kpp13 };
enum centestimator{ kV0M, kV0A, kZNA, kCL1 };
enum energy{ k276, k5dot023, k55,k13 };
enum datayear{k2010, k2011, k2012, k2013, k2015, k2016, k2017, k2018};
enum RaavsEP {kPhiIntegrated, kInPlane, kOutOfPlane};
enum rapidity{ kdefault, k08to04, k07to04, k04to01, k01to01, k01to04, k04to07, k04to08, k01to05 };
enum particularity{ kTopological, kLowPt, kPP7TeVPass4, kBDT };
//../../HFPtSpectrum_D0_multInt_16Oct.root
//../../HFPtSpectrum_LcpKpi_STD_pp_13TeV_NbNbx2_corrected.root
//HFPtSpectrum2("/mnt/downloads/HFPtSpectrum_D09nov.root", kD0Kpi, "/mnt/downloads/D0AccEff_161718_wNtrklWeights_19_EvWithD_18SPD.root", "hEff_C", "hEff_B", "/mnt/downloads/D0_InvMass_withR_TOTint.root", "hRawYields_Mult0", "", "/mnt/downloads/HFPtSpectrum_D0_19_testForCristina.root", 862335392.000000, 57.8e9, kTRUE, kpp7, k2010, k276, kV0M, kPhiIntegrated, kdefault, kTopological, kTRUE)

// Double_t nEvents;

// Standard Analysis
// nEvents = 1691226245; // stdCuts (0-999)
// nEvents = 1587704199; // stdCuts (1-999)

// Std + MVA
// nEvents = 1574197260; // stdCuts + MVA (1-999)

// Pf + MVA
// nEvents = 1574421742; // PfCuts + MVA (1-999)

void HFPtSpectrum2_LcMult(const char *inputCrossSection="OutputHFPtSpectrum_PfCutsTMVA_FixSigma_1_999.root",
                    Int_t decay=kLcK0Sp,
//                    const char *efffilename="../../LcAccEff_final_cut_0998_mb.root",
                    const char *efffilename="../Efficiencies/MergeEfficiencies/LcAccEff_9_over_0_PfCutsTMVA_161718_wNtrklWeights_1_999.root",
                    const char *nameeffprompt= "hEffD",
                    const char *nameefffeed = "hEffB",
//                    const char *recofilename="../../outputfileLcAB_mb_sigmafixtoMB_AB.root",
                    const char *recofilename="../InvariantMassFit/OutputHF/FitResults.root",
                    const char *recohistoname="hRawyield",
                    const char *nevhistoname="",
                    const char *outfilename="OuputHFPtSpectrumLcMult_PfCutsTMVA_FixSigma_1_999.root",
                    Double_t nevents= 1574421742, // overriden by nevhistoname
                    Double_t sigma=57.8e9, // sigma[pb]
                    Bool_t isParticlePlusAntiParticleYield=true,
                    Int_t cc=kpp13,
                    Int_t year=k2018,
                    Int_t Energy=k13,
                    Int_t ccestimator = kV0M,
                    Int_t isRaavsEP=kPhiIntegrated,
                    Int_t rapiditySlice=kdefault,
                    Int_t analysisSpeciality=kTopological,
                    Bool_t setUsePtDependentEffUncertainty=true){

  //
  // Get the histograms from the files
  //
  TH1D *hDirectMCpt=0;           // Input MC c-->D spectra
  TH1D *hFeedDownMCpt=0;         // Input MC b-->D spectra
  TH1D *hDirectMCptMax=0;        // Input MC maximum c-->D spectra
  TH1D *hDirectMCptMin=0;        // Input MC minimum c-->D spectra
  TH1D *hFeedDownMCptMax=0;      // Input MC maximum b-->D spectra
  TH1D *hFeedDownMCptMin=0;      // Input MC minimum b-->D spectra
  TGraphAsymmErrors * gFcConservative = 0; // Input fPrompt fraction

  TH1D *hDirectEffpt=0;          // c-->D Acceptance and efficiency correction
  TH1D *hFeedDownEffpt=0;        // b-->D Acceptance and efficiency correction
  TH1D *hRECpt=0;                // all reconstructed D
  //
  // Get theory predictions from cross section file for f_prompt
  //
  if(gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",inputCrossSection)) !=0){
    printf("File %s with input fPrompt does not exist -> exiting\n",inputCrossSection);
    return;
  }
  TFile* inputcrossfile = new TFile(inputCrossSection,"read");
  if(!inputcrossfile){
    printf("File %s with input fPrompt not opened -> exiting\n",inputCrossSection);
    return;
  }

  hDirectMCpt = (TH1D*)inputcrossfile->Get("hDirectMCpt");
  hFeedDownMCpt = (TH1D*)inputcrossfile->Get("hFeedDownMCpt");
  hDirectMCptMax = (TH1D*)inputcrossfile->Get("hDirectMCptMax");
  hDirectMCptMin = (TH1D*)inputcrossfile->Get("hDirectMCptMin");
  hFeedDownMCptMax = (TH1D*)inputcrossfile->Get("hFeedDownMCptMax");
  hFeedDownMCptMin = (TH1D*)inputcrossfile->Get("hFeedDownMCptMin");
  gFcConservative = (TGraphAsymmErrors*)inputcrossfile->Get("gFcConservative");//gFcCorrConservative

  //
  // Get efficiencies for cross section calculation
  //
  if(gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",efffilename)) !=0){
    printf("File %s with efficiencies does not exist -> exiting\n",efffilename);
    return;
  }
  TFile * efffile = new TFile(efffilename,"read");
  if(!efffile){
    printf("File %s with efficiencies not opened -> exiting\n",efffilename);
    return;
  }
  hDirectEffpt = (TH1D*)efffile->Get(nameeffprompt);
  hDirectEffpt->SetNameTitle("hDirectEffpt","direct acc x eff");
  hFeedDownEffpt = (TH1D*)efffile->Get(nameefffeed);
  hFeedDownEffpt->SetNameTitle("hFeedDownEffpt","feed-down acc x eff");

  //
  // Get raw yield for cross section calculation
  //
  if(gSystem->Exec(Form("ls -l %s > /dev/null 2>&1",recofilename)) !=0){
    printf("File %s with raw yield does not exist -> exiting\n",recofilename);
    return;
  }
  TFile * recofile = new TFile(recofilename,"read");
  if(!recofile){
    printf("File %s with raw yields not opened -> exiting\n",recofilename);
    return;
  }
  hRECpt = (TH1D*)recofile->Get(recohistoname);
  hRECpt->SetNameTitle("hRECpt","Reconstructed spectra");

  //
  // Get (corrected) number of analysd events for cross section calculation
  //
  TH1F* hNorm=(TH1F*)recofile->Get(nevhistoname);
  if(hNorm){
    nevents=hNorm->GetBinContent(1);
    cout<< " number of events " << nevents << " it should be 1683928970.619030 " << endl;
  }else{
    printf("Histogram with number of events for norm not found in raw yiled file\n");
    printf("  nevents = %.0f will be used\n",nevents);
  }

  Int_t fnPtBins = hRECpt->GetNbinsX();
  Double_t *fPtBinLimits = new Double_t[fnPtBins+1];
  Double_t *fPtBinWidths = new Double_t[fnPtBins];
  Double_t xlow=0., binwidth=0.;
  for(Int_t i=1; i<=fnPtBins; i++){
    binwidth = hRECpt->GetBinWidth(i);
    xlow = hRECpt->GetBinLowEdge(i);
    fPtBinLimits[i-1] = xlow;
    fPtBinWidths[i-1] = binwidth;
  }
  fPtBinLimits[fnPtBins] = xlow + binwidth;

  //
  // Define the remaining output histograms to be calculated here
  //
  TH1D *histoYieldCorr = new TH1D("histoYieldCorr","corrected yield",fnPtBins,fPtBinLimits);
  TH1D *histoYieldCorrMax = new TH1D("histoYieldCorrMax","max corrected yield (no feed-down corr)",fnPtBins,fPtBinLimits);
  TH1D *histoYieldCorrMin = new TH1D("histoYieldCorrMin","min corrected yield",fnPtBins,fPtBinLimits);
  TGraphAsymmErrors * gYieldCorr = new TGraphAsymmErrors(fnPtBins+1);
  TGraphAsymmErrors * gYieldCorrExtreme = new TGraphAsymmErrors(fnPtBins+1);
  TGraphAsymmErrors * gYieldCorrConservative = new TGraphAsymmErrors(fnPtBins+1);
  gYieldCorr->SetNameTitle("gYieldCorr","gYieldCorr (by Nb)");
  gYieldCorrExtreme->SetNameTitle("gYieldCorrExtreme","Extreme gYieldCorr (by Nb)");
  gYieldCorrConservative->SetNameTitle("gYieldCorrConservative","Conservative gYieldCorr (by Nb)");

  TH1D *histoSigmaCorr = new TH1D("histoSigmaCorr","corrected invariant cross-section",fnPtBins,fPtBinLimits);
  TH1D *histoSigmaCorrMax = new TH1D("histoSigmaCorrMax","max corrected invariant cross-section",fnPtBins,fPtBinLimits);
  TH1D *histoSigmaCorrMin = new TH1D("histoSigmaCorrMin","min corrected invariant cross-section",fnPtBins,fPtBinLimits);
  TGraphAsymmErrors * gSigmaCorr = new TGraphAsymmErrors(fnPtBins+1);
  TGraphAsymmErrors * gSigmaCorrExtreme = new TGraphAsymmErrors(fnPtBins+1);
  TGraphAsymmErrors * gSigmaCorrConservative = new TGraphAsymmErrors(fnPtBins+1);
  gSigmaCorr->SetNameTitle("gSigmaCorr","gSigmaCorr (by Nb)");
  gSigmaCorrExtreme->SetNameTitle("gSigmaCorrExtreme","Extreme gSigmaCorr (by Nb)");
  gSigmaCorrConservative->SetNameTitle("gSigmaCorrConservative","Conservative gSigmaCorr (by Nb)");

  // NB: Don't care for the moment about fhStatUncEffc/bSigma fhStatUncEffc/bFD
  TH1D *hStatUncEffcSigma = new TH1D("fhStatUncEffcSigma","direct charm stat unc on the cross section",fnPtBins,fPtBinLimits);
  TH1D *hStatUncEffbSigma = new TH1D("fhStatUncEffbSigma","secondary charm stat unc on the cross section",fnPtBins,fPtBinLimits);
  TH1D *hStatUncEffcFD = new TH1D("fhStatUncEffcFD","direct charm stat unc on the feed-down correction",fnPtBins,fPtBinLimits);
  TH1D *hStatUncEffbFD = new TH1D("fhStatUncEffbFD","secondary charm stat unc on the feed-down correction",fnPtBins,fPtBinLimits);

  //
  // Do the corrected yield calculation.
  // NB: Don't care for the moment about histoYieldCorrMin/Max gYieldCorrExtreme/Conservative
  // 
  Double_t value = 0., errvalue = 0., errvalueMax = 0., errvalueMin = 0.;
  for (Int_t ibin=1; ibin<=fnPtBins; ibin++) {
    // Calculate the value
    //    physics =  [ reco  - (lumi * delta_y * BR_b * eff_trig * eff_b * Nb_th) ] / bin-width
    //            =    reco * fprompt_NB
    Double_t frac = 1.0, errfrac =0.;

    // Variables initialization
    value = 0.; errvalue = 0.; errvalueMax = 0.; errvalueMin = 0.;

    // Get fPrompt from input value
    Double_t x = 0., correction = 0;
    gFcConservative->GetPoint(ibin,x,correction);

    // Calculate corrected yield (= raw yield * fprompt)
    if( hRECpt->GetBinContent(ibin)>0. && hRECpt->GetBinContent(ibin)!=0. /*&& hFeedDownMCpt->GetBinContent(ibin)>0.*/ && hFeedDownEffpt->GetBinContent(ibin)>0. ) value = hRECpt->GetBinContent(ibin) * correction;
    value /= hRECpt->GetBinWidth(ibin);
    if (value<0.) value =0.;

    //  Statistical uncertainty:   delta_physics = sqrt ( (delta_reco)^2 )  / bin-width
    if (value!=0. && hRECpt->GetBinError(ibin) && hRECpt->GetBinError(ibin)!=0.) errvalue = hRECpt->GetBinError(ibin);
    errvalue /= hRECpt->GetBinWidth(ibin);

    histoYieldCorr->SetBinContent(ibin,value);
    histoYieldCorr->SetBinError(ibin,errvalue);
    gYieldCorr->SetPoint(ibin,x,value);
    gYieldCorr->SetPointError(ibin,(fPtBinWidths[ibin-1]/2.),(fPtBinWidths[ibin-1]/2.),errvalueMin,errvalueMax);

  }

  //
  // Do the corrected sigma calculation.
  // NB: Don't care for the moment about histoSigmaCorrMin/Max and gSigmaCorrExtreme/Conservative
  // 
  Double_t fLuminosity[2] = {nevents / sigma, 0.04 * nevents / sigma};
  //Double_t fTrigEfficiency[2] = {0.862, 0.018};
  //0.862 +/- 0.018
  Double_t fTrigEfficiency[2] = {1.0, 0};
  Double_t fGlobalEfficiencyUncertainties[2] = {0.05, 0.05};
  Double_t deltaY = 1.0;
  Double_t branchingRatioC = 1.0;
  Double_t branchingRatioBintoFinalDecay = 1.0;
  Bool_t fGlobalEfficiencyPtDependent = setUsePtDependentEffUncertainty;
  Int_t fParticleAntiParticle = 1;
  if( isParticlePlusAntiParticleYield ) fParticleAntiParticle = 2;
  printf("\n\n     Correcting the spectra with : \n   luminosity = %2.2e +- %2.2e, trigger efficiency = %2.2e +- %2.2e, \n    delta_y = %2.2f, BR_c = %2.2e, BR_b_decay = %2.2e \n    %2.2f percent uncertainty on the efficiencies, and %2.2f percent uncertainty on the b/c efficiencies ratio \n    usage of pt-dependent efficiency uncertainty for Nb uncertainy calculation? %1.0d \n\n",fLuminosity[0],fLuminosity[1],fTrigEfficiency[0],fTrigEfficiency[1],deltaY,branchingRatioC,branchingRatioBintoFinalDecay,fGlobalEfficiencyUncertainties[0],fGlobalEfficiencyUncertainties[1],fGlobalEfficiencyPtDependent);

  // protect against null denominator
  if (deltaY==0. || fLuminosity[0]==0. || fTrigEfficiency[0]==0. || branchingRatioC==0.) {
    printf(" Hey you ! Why luminosity or trigger-efficiency or the c-BR or delta_y are set to zero ?! ");
    return;
  }

  for (Int_t ibin=1; ibin<=fnPtBins; ibin++) {

    // Variables initialization
    value=0.; errvalue=0.;

    Double_t x = histoYieldCorr->GetBinCenter(ibin);

    // Sigma calculation
    //   Sigma = ( 1. / (lumi * delta_y * BR_c * ParticleAntiPartFactor * eff_trig * eff_c ) ) * spectra (corrected for feed-down)
    if (hDirectEffpt->GetBinContent(ibin) && hDirectEffpt->GetBinContent(ibin)!=0. && hRECpt->GetBinContent(ibin)>0.) { 
      value = histoYieldCorr->GetBinContent(ibin) / ( deltaY * branchingRatioC * fParticleAntiParticle * fLuminosity[0] * fTrigEfficiency[0] * hDirectEffpt->GetBinContent(ibin) );
    }

    // Sigma statistical uncertainty:
    //   delta_sigma = sigma * sqrt ( (delta_spectra/spectra)^2 )
    if (value!=0.) {
      errvalue = value * (histoYieldCorr->GetBinError(ibin)/histoYieldCorr->GetBinContent(ibin));
    }

    histoSigmaCorr->SetBinContent(ibin,value);
    histoSigmaCorr->SetBinError(ibin,errvalue);
    gSigmaCorr->SetPoint(ibin,x,value); // i,x,y
    gSigmaCorr->SetPointError(ibin,(fPtBinWidths[ibin-1]/2.),(fPtBinWidths[ibin-1]/2.),errvalueMin,errvalueMax); // i,xl,xh,yl,yh
  }
  
  //Set the systematics externally
  Bool_t combineFeedDown = true;
  AliHFSystErr *systematics = new AliHFSystErr();
  if(year==k2010) systematics->SetRunNumber(10);
  else if(year==k2011) systematics->SetRunNumber(11);
  else if(year==k2012) systematics->SetRunNumber(12);
  else if(year==k2013) systematics->SetRunNumber(13);
  else if(year==k2015) systematics->SetRunNumber(15);
  else if(year==k2016) systematics->SetRunNumber(16);
  else if(year==k2017) systematics->SetRunNumber(17);
  else if(year==k2018) systematics->SetRunNumber(18);
  
  if( cc==kpp276 ) {
    systematics->SetIsLowEnergy(true);
  }
  else if (cc==kpp8){
    systematics->SetRunNumber(12);
  }
  else if (cc==kpp5){
    systematics->SetIs5TeVAnalysis(true);
    systematics->SetCollisionType(0);
  }
  else if ( cc == kpPb0100 || cc == kpPb010 || cc == kpPb020 || cc == kpPb1020 || cc == kpPb2040 || cc == kpPb4060 || cc == kpPb60100 ) {
    systematics->SetCollisionType(2);
    
    // Rapidity slices
    if(rapiditySlice!=kdefault){
      systematics->SetIspPb2011RapidityScan(true);
      TString rapidity="";
      switch(rapiditySlice) {
        case k08to04: rapidity="0804"; break;
        case k07to04: rapidity="0804"; break;
        case k04to01: rapidity="0401"; break;
        case k01to01: rapidity="0101"; break;
        case k01to04: rapidity="0104"; break;
        case k04to07: rapidity="0408"; break;
        case k04to08: rapidity="0408"; break;
        case k01to05: rapidity="0401"; break;
      }
      systematics->SetRapidity(rapidity);
    }
    // Centrality slices
    if(ccestimator==kV0A) {
      if(cc == kpPb020) systematics->SetCentrality("020V0A");
      else if(cc == kpPb2040) systematics->SetCentrality("2040V0A");
      else if(cc == kpPb4060) systematics->SetCentrality("4060V0A");
      else if(cc == kpPb60100) systematics->SetCentrality("60100V0A");
    } else if (ccestimator==kZNA) {
      if(cc == kpPb010) {systematics->SetCentrality("010ZNA");}
      else if(cc == kpPb020) systematics->SetCentrality("020ZNA");
      else if(cc == kpPb1020) systematics->SetCentrality("1020ZNA");
      else if(cc == kpPb2040) systematics->SetCentrality("2040ZNA");
      else if(cc == kpPb4060) systematics->SetCentrality("4060ZNA");
      else if(cc == kpPb60100) systematics->SetCentrality("60100ZNA");
    } else if (ccestimator==kCL1) {
      if(cc == kpPb020) systematics->SetCentrality("020CL1");
      else if(cc == kpPb2040) systematics->SetCentrality("2040CL1");
      else if(cc == kpPb4060) systematics->SetCentrality("4060CL1");
      else if(cc == kpPb60100) systematics->SetCentrality("60100CL1");
    } else {
      if(!(cc == kpPb0100)) {
        cout <<" Error on the pPb options"<<endl;
        return;
      }
    }
  }
  //
  else if( cc!=kpp7 )  {
    systematics->SetCollisionType(1);
    if(Energy==k276){
      if ( cc == k07half ) systematics->SetCentrality("07half");
      else if ( cc == k010 )  systematics->SetCentrality("010");
      else if ( cc == k1020 )  systematics->SetCentrality("1020");
      else if ( cc == k020 )  systematics->SetCentrality("020");
      else if ( cc == k2040 || cc == k2030 || cc == k3040 ) {
        systematics->SetCentrality("2040");
        systematics->SetIsPbPb2010EnergyScan(true);
      }
      else if ( cc == k3050 ) {
        if (isRaavsEP == kPhiIntegrated) systematics->SetCentrality("4080");
        else if (isRaavsEP == kInPlane) systematics->SetCentrality("3050InPlane");
        else if (isRaavsEP == kOutOfPlane) systematics->SetCentrality("3050OutOfPlane");
      }
      else if ( cc == k4060 || cc == k4050 || cc == k5060 )  systematics->SetCentrality("4060");
      else if ( cc == k6080 )  systematics->SetCentrality("6080");
      else if ( cc == k4080 ) systematics->SetCentrality("4080");
    } else if (Energy==k5dot023){
      if ( cc == k010 ){
        systematics->SetCentrality("010");
      } else if ( cc == k1030 ) {
        systematics->SetCentrality("3050"); //no systematics available for 10--30
      } else if ( cc == k3050 ) {
        systematics->SetCentrality("3050");
      } else if ( cc == k6080 ) {
        systematics->SetCentrality("6080");
      }
    }else if (Energy == k13){
        systematics->SetCentrality("0100");
        systematics->SetCollisionType(0);
        systematics->SetRunNumber(18);
    }else {
      cout << " Systematics not yet implemented " << endl;
      return;
    }
  } else { systematics->SetCollisionType(0); }
  if(analysisSpeciality==kLowPt){
    systematics->SetIsLowPtAnalysis(true);
  }
  else if(analysisSpeciality==kPP7TeVPass4){
    systematics->SetIsPass4Analysis(kTRUE);
  }
  else if(analysisSpeciality==kBDT){
    systematics->SetIsBDTAnalysis(kTRUE);
  }
  //
  systematics->Init(decay+1);
  systematics->DrawErrors(gFcConservative);

  //
  // Define output file, in same style as original HFPtSpectrum macro
  //
  TFile *out = new TFile(outfilename,"recreate");
  out->cd();

//  hDirectMCpt->Write();
//  hFeedDownMCpt->Write();
//  hDirectMCptMax->Write();
//  hDirectMCptMin->Write();
//  hFeedDownMCptMax->Write();
//  hFeedDownMCptMin->Write();

  hDirectEffpt->Write();
  hFeedDownEffpt->Write();
  hRECpt->Write();

  histoYieldCorr->Write();
  histoYieldCorrMax->Write();
  histoYieldCorrMin->Write();

  histoSigmaCorr->Write();
  histoSigmaCorrMax->Write();
  histoSigmaCorrMin->Write();

  gYieldCorr->Write();
  gSigmaCorr->Write();
  gYieldCorrExtreme->Write();
  gSigmaCorrExtreme->Write();
  gYieldCorrConservative->Write();
  gSigmaCorrConservative->Write();

  gFcConservative->Write();

  hStatUncEffcSigma->Write();
  hStatUncEffbSigma->Write();
  hStatUncEffcFD->Write();
  hStatUncEffbFD->Write();

  systematics->Write();

  out->Close();
  recofile->Close();
  efffile->Close();
  inputcrossfile->Close();
}
