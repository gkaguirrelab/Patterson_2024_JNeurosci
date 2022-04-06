% This simulation examines how changes in migdet chromatic sensitivity
% might impact the temporal sensitivity function of visual cortex across
% eccentricity.
%
% The function begins by calling code from Wool 2018 J Neuroscience to
% estimate the response of midget retinal ganglion cells to low-spatial
% frequency L-M and L+M contrast. The response of modeled midgets to
% achromatic (~spatially uniform) contrast remains constant across
% eccentricity. Due to random wiring, however, the response to chromatic
% contrast drops into the periphery.
%
% The simulation then uses the proportional L+M and L-M responses at each
% eccentricity to mix together the temporal sensitivity functions measured
% at the fovea in V1 cortex for luminance and L-M modulations. The idea is
% that the response to a luminance modulation is conveyed by parasol cells,
% and by an increasing proportion of midgets further into the visual
% periphery.
%

% Housekeeping
clear
close all

% Set a flag that forces the routine to re-run the cell simulation, even if
% the simulation results are present.
forceRecalc = true;

% Define a set of eccentricities in degrees, log spaced and centered within
% each of the bins of the Mt Sinai data
eD = logspace(log10(1),log10(64),13); % eccentricities in degrees,

% The Wool model operates in units of mm in the macaque retina, so convert
% ecentricity to macaque mm.
eM = (eD.*223)./10000; % convert to macaque mm (M. mulatta, Perry & Cowey 1985)

% The cellSim results filename
cellSimResultsFile = fullfile(fileparts(mfilename('fullpath')),'midgetMixingSimResults.mat');

% Check if we already have results

if isfile(cellSimResultsFile) && ~forceRecalc

    % Load the previous simulation
    load(cellSimResultsFile,'Data')

    % How many cells?
    tmp = fieldnames(Data);
    nCells = length(fieldnames(Data.(tmp{1})));

else

    % The number of cells in the simulation
    nCells = 1000;

    % Loop through the eccentricities
    for ee=1:length(eM)
        EccentricityLabel=strcat('Ecc',num2str(eM(ee)),'mm');
        for cc=1:nCells
            CellLabel=strcat('Cell',num2str(cc));
            Cell=DoG(eM(ee));
            Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel))=Cell;
            disp(strcat(CellLabel,' @ ',EccentricityLabel,' complete.'));
        end
        disp(strcat(EccentricityLabel,' complete.'));
    end

    % Save
    save(cellSimResultsFile,'Data');
    
end

% Extract the L+M and L-M responses of the cells
for ee=1:length(eM)
    EccentricityLabel=matlab.lang.makeValidName(strcat('Ecc',num2str(eM(ee)),'mm'));
    thisEccLMSumResponse = []; thisEccLMDiffResponse = [];
    for cc = 1: nCells
        CellLabel=strcat('Cell',num2str(cc));
        thisEccLMSumResponse(cc) = Data.(EccentricityLabel).(CellLabel).ResponseFunctions.LMSumResponse.Amplitude;
        thisEccLMDiffResponse(cc) = Data.(EccentricityLabel).(CellLabel).ResponseFunctions.LMDiffResponse.Amplitude;
    end
    LMSumResponse(ee) = mean(thisEccLMSumResponse);
    LMDiffResponse(ee) = mean(thisEccLMDiffResponse);
end

% Calculate the relative L+M vs L-M response
lumIntrusionRelToFovea = LMSumResponse./LMDiffResponse;
lumIntrusionRelToFovea = (lumIntrusionRelToFovea./lumIntrusionRelToFovea(1));

% Subject the luminance intrusion to an expansive non-linearity. This is a
% "boost" of the chromatic component.
epsilon = 4;
lumIntrusionRelToFoveaExpansiveNonLin = lumIntrusionRelToFovea.^epsilon;
lumIntrusionRelToFoveaExpansiveNonLin = (lumIntrusionRelToFoveaExpansiveNonLin-1)./max(lumIntrusionRelToFoveaExpansiveNonLin-1);
lumIntrusionRelToFoveaExpansiveNonLin = lumIntrusionRelToFoveaExpansiveNonLin.*max(lumIntrusionRelToFovea-1)+1;
chromaticBoost = lumIntrusionRelToFovea./lumIntrusionRelToFoveaExpansiveNonLin;

% Plot the response to chromatic contrast as a function of eccentricity
figure
subplot(2,1,1)
tmp = LMDiffResponse ./ max(LMDiffResponse);
loglog(eD,tmp,'-k');
hold on
tmp = tmp.*chromaticBoost;
loglog(eD,tmp,'-r');
ylabel({'Midget L-M sensitivity','[relative to fovea]'})
xlabel('Eccentricity [deg]')
ylim([0.1 1.1])
xlim([0.9 70])
xticks([1 2 4 8 16 32 64])
xticklabels({'1','2','4','8','16','32','64'})
yticks([0.125 0.25 0.5 1])
yticklabels({'0.125','0.25','0.5','1'})
title('Midget L-M sensitivity as a function of eccentricity')

% Plot the intrusion of luminance signals into the midget pathway
subplot(2,1,2)
semilogx(eD,lumIntrusionRelToFovea,'-k');
hold on
semilogx(eD,lumIntrusionRelToFoveaExpansiveNonLin,'-r');
semilogx(eD,ones(size(eD)),':k');
semilogx([21.5 21.5],[0 7],':k');
ylabel('L+M intrusion relative to fovea')
xlim([0.9 70])
xlabel('Eccentricity')
xticks([1 2 4 8 16 32 64])
xticklabels({'1','2','4','8','16','32','64'})

% Load the TSFs for P1 for the fovea
v1DataFile = fullfile(fileparts(mfilename('fullpath')),'mtSinaiData','V1fovea_Mean_and_Fits_fMRIflicker.mat');
load(v1DataFile,'P1_fovea_fit');

% Pull out the luminance and L-M fits
fitLum = P1_fovea_fit(3,:)';
fitRG = P1_fovea_fit(1,:)';

% Define the frequency support for the TSFs
freqSupport = logspace(log10(1),log10(128),3500);

% Plot the TSFs
figure
semilogx(freqSupport,fitLum,'-k')
hold on
semilogx(freqSupport,fitRG,'-r')

% Create a mix of the lum and RG TSFs
figure
for mm=1:length(eM)
    k = fitLum + (lumIntrusionRelToFovea(mm)-1).*fitRG;
    k = k./(max(k));
    [~,idx] = max(k);
    peakTTF(mm) = freqSupport(idx);

    k = fitLum + (lumIntrusionRelToFoveaExpansiveNonLin(mm)-1).*fitRG;
    k = k./(max(k));
    [~,idx] = max(k);
    peakTTFExpanded(mm) = freqSupport(idx);
end

% Plot these
figure
semilogx(eD,peakTTF,'-k');
hold on
semilogx(eD,peakTTFExpanded,'-r');
ylabel('Peak temporal sensitivity [Hz]')
xlim([0.9 70])
xlabel('Eccentricity')
xticks([1 2 4 8 16 32 64])
xticklabels({'1','2','4','8','16','32','64'})
