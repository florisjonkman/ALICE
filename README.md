# ALICE

Dear scientist,

I have setup the current repository to document, store and preserve all the macro's, files and codes that I have created in my analysis. The main goal of my analysis was to analyze 'charm baryon production at central rapidity as a function of multiplicity in proton-proton collisions at √13 TeV'. For a full overview of my analysis, I would like to refer you to my [Master Thesis](https://www.dropbox.com/s/i491e4obhvqkvv5/Thesis_Jonkman.pdf?dl=0) (Currently still under construction and not approved).

In the upcoming weeks until the end of my thesis I will keep updating this repository, and try to update and debug my code even more. If you encounter any bugs or faults, please contact me, so that I can update this.

Thanks,

Floris Jonkman

## Analysis Chain
A full overview of my analysis chain can be found in the thesis. I have two different analysis chains, one for the analysis with the use of BDT, one for the analysis without BDT. 

### Analysis with BDT
The analysis with BDT has 7 consecutive steps. After each steps is listed the referred directory.
1. Prepare signal and background candidates (BoostedDecisionTree/PrepareSgnBkgTree.C)
2. Train the BDT (BoostedDecisionTree/BDTTraining.C)
3. Apply model on full data statistics (PWGHF/vertexingHF/AliAnalysisTaskSELc2V0bachelorTMVAApp.cxx .h)
    - Different cut files in CutFiles
    - Please see for example D2H_pp trains, 2653 - 2656
4. Determine optimal BDT cut, using BDTToolkit (SignalExtraction/BDTToolkit)
5. Extract rawyield (SignalExtraction/HFMassFitter.C)
    - Compare results (SignalExtraction/CompareMassFits.C)
6. Compute Efficiencies
    - Cut efficiencies (Efficiencies/CutEfficiencies)
    - TMVA efficiencies (Efficiencies/TMVAEfficiencies)
    - Combine Efficiencies (Efficiencies/doLcEffMergingBunchesTMVA.C)
7. Compute corrected yield per event
    - Compute cross-section with all the ingredients (CorrectedYield/HFPtSpectrum.C) 
    - Compute corrected yield per event (CorrectedYield/ComputeTotalYieldPerNEvents.C) 
    
### Analysis w/o BDT
The analysis without BDT has 5 consecutive steps. After each steps is listed the referred directory.
1. Run on full stastics (PWGHF/vertexingHF/DvsMultiplicityTask.cxx .h)
    - Please see for example D2H_pp trains, 2579 - 2582
2. Generate invariant mass histograms (SignalExtraction/GetInvariantMassFromTH3F.C)
3. Extract rawyield (SignalExtraction/HFMassFitter.C)
    - Compare results (SignalExtraction/CompareMassFits.C)
4. Compute Efficiencies
    - Cut efficiencies (Efficiencies/CutEfficiencies)
    - Combine Efficiencies (Efficiencies/doLcEffMergingBunches.C)
5. Compute corrected yield per event
    - Compute cross-section with all the ingredients (CorrectedYield/HFPtSpectrum.C) 
    - Compute corrected yield per event (CorrectedYield/ComputeTotalYieldPerNEvents.C) 

