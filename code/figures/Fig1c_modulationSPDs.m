%% spdModulationFigure
%
% Load the validation files and plot the SPDs


% Get the localSaveDir pref
localSaveDir = getpref('Patterson_2024_JNeurosci','localSaveDir');

% Labels for directions
modDirections = {'Light flux','L-M','S'};

% Path to validation files (LF, L-M, S)
validationFilePath{1} = fullfile(fileparts(localSaveDir),'MtSinaiCalibrationFiles/cache/stimuli/Cache-LightFluxXEccentricity/BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm/11-Apr-2016_17_07_28/validation/14-Apr-2016_15_29_55/Cache-LightFluxXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat');
validationFilePath{2} = fullfile(fileparts(localSaveDir),'MtSinaiCalibrationFiles/cache/stimuli/Cache-LMinusMDirectedXEccentricity/BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm/11-Apr-2016_17_07_28/validation/14-Apr-2016_15_31_33/Cache-LMinusMDirectedXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat');
validationFilePath{3} = fullfile(fileparts(localSaveDir),'MtSinaiCalibrationFiles/cache/stimuli/Cache-SDirectedXEccentricity/BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm/11-Apr-2016_17_07_28/validation/14-Apr-2016_15_32_14/Cache-SDirectedXEccentricity-BoxCRandomizedLongCableCStubby1NoLens_ND10_ContactLens_0_5mm-SpotCheck.mat');

% Define where we want to save these figures
resultsSaveDir = fullfile(localSaveDir,'Fig 1 - experimental design');
mkdir(resultsSaveDir);

% Create a figure
figHandle = figure();

% Load a validation file
load(validationFilePath{1});

% Create the SPD x-axis vector
S = cals{1}.modulationMaxMeas.meas.pr650.S;
wavelengths = S(1):S(2):S(3)*S(2)+S(1)-1;

% Grab the photoreceptors
T_receptors = cals{1}.describe.cache.data(end).describe.T_receptors;
photoreceptors = cals{1}.describe.cache.data(end).describe.photoreceptors;

% Plot the cone sensitivity functions
subplot(2,2,1);
for ii = 1:3
    switch photoreceptors{ii}
        case 'LCone2DegTabulatedSS'
            lineColor = 'r';
        case 'MCone2DegTabulatedSS'
            lineColor = 'g';
        case 'SCone2DegTabulatedSS'
            lineColor = 'b';
    end
    plot(wavelengths,T_receptors(ii,:),[':' lineColor],'LineWidth',2)
    hold on
    plot(wavelengths,T_receptors(ii+3,:),['-' lineColor],'LineWidth',2)
    p=patch([wavelengths fliplr(wavelengths)], [T_receptors(ii,:) fliplr(T_receptors(ii+3,:))], lineColor, 'FaceAlpha',0.1);
end
ylabel('Sensitivity [au]');
xlabel('Wavelength [nm]');
xlim([380 780]);
ylim([-0.05 1]);
legend({'2°','10°','','','','','','',''});
    box off
set(gca,'TickDir','out'); 

% Plot the three modulations
for ii = 1:3
    subplot(2,2,ii+1);
    load(validationFilePath{ii});

    % Plot the predicted SPD for background, positive, and negative arms Do
    % we need to adjust the radiance by S(2) to express the value per nm,
    % as opposed to per measurement bin?
    plot(wavelengths,cals{1}.modulationBGMeas.predictedSpd,'--','Color',[0.5 0.5 0.5],'LineWidth',2);
    hold on
    plot(wavelengths,cals{1}.modulationMaxMeas.predictedSpd,'Color','r','LineWidth',1);
    plot(wavelengths,cals{1}.modulationMinMeas.predictedSpd,'Color','k','LineWidth',1);

    % Add some labels
    xlim([380 780]);
    ylim([-0.001 0.025])
    ylabel('Radiance [W/m^2/sr/nm]');
    xlabel('Wavelength [nm]');
    title(modDirections{ii});
    box off
    set(gca,'TickDir','out'); 
end

% Save the figure
fileName = fullfile(resultsSaveDir,'spdModulations.pdf');
print(fileName,'-dpdf');
