%% ====== DEPENDENCIES ======
% First we'll need to set up dependencies (SST) with TbTb
tbUse('SilentSubstitutionToolbox');

%% ====== FILES ======
% We need to locate the files. We'll just use the spot check
% measurements, which we did before and after the experiment.
%
% Light flux
% ==========
% Pre-experiment:
%thePathLightFluxPre = "cache/stimuli/Cache-LightFluxXEccentricity/BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm/11-Apr-2016_17_07_28/validation/15-Apr-2016_08_32_51/Cache-LightFluxXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat";
thePathLightFluxPre = "data/Cache-LightFluxXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat";

% Post-experiment (same day in the afternoon):
%thePathLightFluxPost = "cache/stimuli/Cache-LightFluxXEccentricity/BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm/11-Apr-2016_17_07_28/validation/15-Apr-2016_15_56_12/Cache-LightFluxXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat";

% L-M
% ==========
% Pre-experiment:
%thePathLMinusMPre = "cache/stimuli/Cache-LMinusMDirectedXEccentricity/BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm/11-Apr-2016_17_07_28/validation/15-Apr-2016_08_34_29/Cache-LMinusMDirectedXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat";
thePathLMinusMPre = "data/Cache-LMinusMDirectedXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat";

% Post-experiment (same day in the afternoon):
%thePathLMinusMPost = "cache/stimuli/Cache-LMinusMDirectedXEccentricity/BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm/11-Apr-2016_17_07_28/validation/15-Apr-2016_15_57_48/Cache-LMinusMDirectedXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat";

% S
% ==========
% Pre-experiment:
%thePathSPre = "cache/stimuli/Cache-SDirectedXEccentricity/BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm/11-Apr-2016_17_07_28/validation/15-Apr-2016_08_35_09/Cache-SDirectedXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat";
thePathSPre = "data/Cache-SDirectedXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat";

% Post-experiment (same day in the afternoon):
%thePathSPost = "cache/stimuli/Cache-SDirectedXEccentricity/BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm/11-Apr-2016_17_07_28/validation/15-Apr-2016_15_58_29/Cache-SDirectedXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat";

%% Define which direction to analyse
whichDirection = 'LMinusM';
switch whichDirection
    case 'LightFlux'
        thePathToLoad = thePathLightFluxPre;
        axisLimits = [80 100 ; 80 100 ; 80 100 ; 80 100 ; -13 13 ; -13 13];
        nominalContrasts = [90 90 90 90 0 0];
    case 'LMinusM'
        thePathToLoad = thePathLMinusMPre;
        axisLimits = [-20 20 ; -20 20 ; -10 10 ; -10 10 ; -20 20 ; -10 10];
        nominalContrasts = [9 -9 0 0 9 0 0];
    case 'S'
        thePathToLoad = thePathSPre;
        axisLimits = [-2.5 2.5 ; -2.5 2.5 ; 40 60 ; -2.5 2.5 ; -2.5 2.5 ; 40 60];
        nominalContrasts = [0 0 50 0 0 50];
end


%% Loading
% Now, we'll load in the files
tmp = load(thePathToLoad)

bgSpd = tmp.cals{1}.modulationAllMeas(2).meas.pr650.spectrum;
modPosSpd = tmp.cals{1}.modulationAllMeas(3).meas.pr650.spectrum;
modNegSpd = tmp.cals{1}.modulationAllMeas(1).meas.pr650.spectrum;

%% Calculate and plot the chromaticities
load T_xyz1931
S = tmp.cals{1}.modulationAllMeas(1).meas.pr650.S;
T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
bg_photopicLuminanceCdM2 = T_xyz(2,:)*bgSpd;
bg_chromaticityXY = T_xyz(1:2,:)*bgSpd/sum(T_xyz*bgSpd);
modPos_chromaticityXY = T_xyz(1:2,:)*modPosSpd/sum(T_xyz*modPosSpd);
modNeg_chromaticityXY = T_xyz(1:2,:)*modNegSpd/sum(T_xyz*modNegSpd);


%% Set up receptor object
% We'll create multiple variants, one for each field size
fieldSizeDeg = 2:60;
for ii = 1:length(fieldSizeDeg)
    % Generate a receptor with standard age (32a), assumed pupil diameter (7mm) and varying field size
    tmpReceptor = SSTReceptorHuman('obsPupilDiameterMm', 7, 'fieldSizeDeg', fieldSizeDeg(ii));
    
    % Just pull out the spectral sensitivities
    T = tmpReceptor.T.T_energyNormalized;
    
    % Calculate the contrast
    receptorContrast(:, ii) = T*(modPosSpd-bgSpd) ./ (T*bgSpd);
    postreceptoralContrast(:, ii) = ComputePostreceptoralContrastsFromLMSContrasts(receptorContrast(1:3, ii));
    
   % Resample
   NSamples = 200;
   tmpReceptor.makeSpectralSensitivitiesStochastic('NSamples', NSamples);
   for ij = 1:NSamples
       Ts = tmpReceptor.Ts{ij}.T_energyNormalized;
           receptorContrastStochastic(:, ij) = Ts*(modPosSpd-bgSpd) ./ (T*bgSpd);
           postreceptoralContrastStochastic(:, ij) = ComputePostreceptoralContrastsFromLMSContrasts(receptorContrastStochastic(1:3, ij));
   end
   receptorContrastLowCI(ii, :) = prctile(receptorContrastStochastic', 2.5)';
   receptorContrastHighCI(ii, :) = prctile(receptorContrastStochastic', 97.5)';
   postreceptoralContrastLowCI(ii, :) = prctile(postreceptoralContrastStochastic', 2.5)';
   postreceptoralContrastHighCI(ii, :) = prctile(postreceptoralContrastStochastic', 97.5)';
end

%% Plot
% Make the contrast plot for the backdrop
fieldSizeFine = 2:0.1:60;
contrastFine = -150:0.1:150;
for ij = 1:length(contrastFine)
    for ii = 1:length(fieldSizeFine)
        theContrastHere(ij, ii) = contrastFine(ij);
    end
end

% Set up some colours for plotting
theRGB = SSTDefaultReceptorColors;

subplot(2, 3, 1); % L cone
imagesc(fieldSizeFine, contrastFine, theContrastHere); hold on
plot(fieldSizeDeg, 100*receptorContrast(1, :), '-', 'Color', theRGB(1, :), 'LineWidth', 2)
plot(fieldSizeDeg, 100*receptorContrastLowCI(:, 1), '-', 'Color', theRGB(1, :));
plot(fieldSizeDeg, 100*receptorContrastHighCI(:, 1), '-', 'Color', theRGB(1, :));
title('L cone contrast');
xlabel('Field size [deg]');
ylabel('Contrast [%]');

subplot(2, 3, 2); % M cone
imagesc(fieldSizeFine, contrastFine, theContrastHere); hold on
plot(fieldSizeDeg, 100*receptorContrast(2, :), '-', 'Color', theRGB(2, :), 'LineWidth', 2)
plot(fieldSizeDeg, 100*receptorContrastLowCI(:, 2), '-', 'Color', theRGB(2, :));
plot(fieldSizeDeg, 100*receptorContrastHighCI(:, 2), '-', 'Color', theRGB(2, :));
title('M cone contrast');
xlabel('Field size [deg]');
ylabel('Contrast [%]');

subplot(2, 3, 3); % S cone
imagesc(fieldSizeFine, contrastFine, theContrastHere); hold on
plot(fieldSizeDeg, 100*receptorContrast(3, :), '-', 'Color', theRGB(3, :), 'LineWidth', 2)
plot(fieldSizeDeg, 100*receptorContrastLowCI(:, 3), '-', 'Color', theRGB(3, :));
plot(fieldSizeDeg, 100*receptorContrastHighCI(:, 3), '-', 'Color', theRGB(3, :));
title('S cone contrast');
xlabel('Field size [deg]');
ylabel('Contrast [%]');

subplot(2, 3, 4); % LMS
imagesc(fieldSizeFine, contrastFine, theContrastHere); hold on
plot(fieldSizeDeg, 100*postreceptoralContrast(1, :), 'Color', [150 80 0]/255, 'LineWidth', 2)
plot(fieldSizeDeg, 100*postreceptoralContrastLowCI(:, 1), '-', 'Color', [150 80 0]/255);
plot(fieldSizeDeg, 100*postreceptoralContrastHighCI(:, 1), '-', 'Color', [150 80 0]/255);
title('LMS contrast');
xlabel('Field size [deg]');
ylabel('Contrast [%]');

subplot(2, 3, 5); % LMinusM
imagesc(fieldSizeFine, contrastFine, theContrastHere); hold on
plot(fieldSizeDeg, 100*postreceptoralContrast(2, :), '-k', 'LineWidth', 2)
plot(fieldSizeDeg, 100*postreceptoralContrastLowCI(:, 2), '-', 'Color', 'k');
plot(fieldSizeDeg, 100*postreceptoralContrastHighCI(:, 2), '-', 'Color', 'k');
title('L-M contrast');
xlabel('Field size [deg]');
ylabel('Contrast [%]');

subplot(2, 3, 6); % S-[L+M]
imagesc(fieldSizeFine, contrastFine, theContrastHere); hold on
plot(fieldSizeDeg, 100*postreceptoralContrast(3, :), '-k', 'LineWidth', 2)
plot(fieldSizeDeg, 100*postreceptoralContrastLowCI(:, 3), '-', 'Color', 'k');
plot(fieldSizeDeg, 100*postreceptoralContrastHighCI(:, 3), '-', 'Color', 'k');
title('S-[L+M] contrast');
xlabel('Field size [deg]');
ylabel('Contrast [%]');



% Update the axis limits
for ii = 1:6
    subplot(2, 3, ii);
    
    % Set the colormap
    colormap(lbmap(1000, 'BrownBlue'))
    
    % Set axis limits
    ylim(axisLimits(ii, :))
    xlim([fieldSizeDeg(1) fieldSizeDeg(end)]);
    
    % Draw the nominal contrast line
    hold on;
    plot([fieldSizeDeg(1) fieldSizeDeg(end)], [nominalContrasts(ii) nominalContrasts(ii)], ':w');
    
    % Some aesthetics
    box off;
    set(gca, 'TickDir', 'out');
    set(gca, 'YDir', 'normal');
    pbaspect([1 0.5 1]);
    colorbar;
end

% Contrast plots
set(gcf, 'PaperPosition', [0 0 30 10]);
set(gcf, 'PaperSize', [30 10]);
saveas(gcf, ['figures/splatter_' whichDirection], 'pdf');
