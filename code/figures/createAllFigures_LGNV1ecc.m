% This m-file reproduces the figures presented in "Eccentricity variation 
% in temporal sensitivity of primary visual cortex related to non-selective 
% midget wiring"

% MRI pre-processing
% The directory submitGears contains CSV files that provide the parameters
% used to process the raw fMRI data

% Figure 2 - example MRI data
run averageTimeSeries.m % Downloads Flywheel Results files, creates plots of average V1 time series for each subject

% Figure 3 - full field LGN and V1 TSFs, peak response, and peak frequency
run ExtractTTFdataDifExpFitArea.m % downloads MRI data from flywheel, and calculates bootstrapped fits for TSFs in fullfield LGN and V1, and saves them
run TTF_figures_areaExp.m % creates TSF figure, and peak response and peak frequency

% Figures 4 and 5 - TSFs in V1 across eccentricity
run ExtractTTFdataDifExpFitEccentricity.m % downloads MRI data from flywheel, and calculates bootstrapped fits for TSFs across V1 eccentricity, and saves them
run TTF_figures_eccExp.m % creates TSF figure, and peak response and peak frequency

% Figure 6 - midget amplification model
run midgetMixingSimulation.m % Runs Wool 2018 DoG code, then create Figure 6