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
whichDirections = {'LightFlux','LMinusM','S'};

for dd=1:length(whichDirections)

    whichDirection = whichDirections{dd};

    switch whichDirection
        case 'LightFlux'
            thePathToLoad = thePathLightFluxPre;
            nominalContrasts = [90 90 90 90 0 0];
        case 'LMinusM'
            thePathToLoad = thePathLMinusMPre;
            nominalContrasts = [9 -9 0 0 9 0 0];
        case 'S'
            thePathToLoad = thePathSPre;
            nominalContrasts = [0 0 50 0 0 50];
    end


    %% Loading
    % Now, we'll load in the files
    tmp = load(thePathToLoad);

    bgSpd = tmp.cals{1}.modulationAllMeas(2).meas.pr650.spectrum;
    modPosSpd = tmp.cals{1}.modulationAllMeas(3).meas.pr650.spectrum;
    modNegSpd = tmp.cals{1}.modulationAllMeas(1).meas.pr650.spectrum;

    %% Calculate the chromaticities
    load T_xyz1931
    S = tmp.cals{1}.modulationAllMeas(1).meas.pr650.S;
    T_xyz = SplineCmf(S_xyz1931,683*T_xyz1931,S);
    bg_photopicLuminanceCdM2 = T_xyz(2,:)*bgSpd;
    bg_chromaticityXY = T_xyz(1:2,:)*bgSpd/sum(T_xyz*bgSpd);
    modPos_chromaticityXY = T_xyz(1:2,:)*modPosSpd/sum(T_xyz*modPosSpd);
    modNeg_chromaticityXY = T_xyz(1:2,:)*modNegSpd/sum(T_xyz*modNegSpd);


    %% Set up receptor object
    % We'll create multiple variants, one for each field size
    fieldSizeDeg = [2 30];
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

    % Store the values
    Tline = [];
    for ff=1:2
        for pp=1:3
            Tline = [Tline postreceptoralContrast(pp,ff)*100, postreceptoralContrastLowCI(ff,pp)*100, postreceptoralContrastHighCI(ff,pp)*100];
        end
    end

    Tmat(dd,:)=Tline;
end

Tmat