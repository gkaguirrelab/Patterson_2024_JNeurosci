% This m-file reproduces the figures presented in "Eccentricity variation 
% in temporal sensitivity of primary visual cortex related to non-selective 
% midget wiring"

% MRI pre-processing


% Figure 2 - example MRI data


% Figure 3 - full field LGN and V1 TSFs, peak response, and peak frequency
run ExtractTTFdataDifExpFitArea.m % downloads MRI data from flywheel, and calculates bootstrapped fits for TSFs in fullfield LGN and V1, and saves them
run TTF_figures_areaExp.m % creates TSF figure, and peak response and peak frequency

% Figures 4 and 5 - TSFs in V1 across eccentricity
run ExtractTTFdataDifExpFitEccentricity.m % downloads MRI data from flywheel, and calculates bootstrapped fits for TSFs across V1 eccentricity, and saves them
run TTF_figures_eccExp.m % creates TSF figure, and peak response and peak frequency

% Figure 6 - midget amplification model