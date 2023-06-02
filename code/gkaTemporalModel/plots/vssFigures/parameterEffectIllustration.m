% scriptCreatePlots

% Housekeeping
clear

% Properties of which model to plot
freqsForPlotting = logspace(0,2,50);

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% Define the eccentricity location of the data
studiedEccentricites = 8;

stimulusDirections = {'LminusM','S','LMS'};
nStims = length(stimulusDirections);

% Fixed features of the model
nCells = 3; nParams = 3; nEccs = 6;

figHandles = figure('Renderer','painters');
figuresize(600,400,'pt');
tiledlayout(2,2,'TileSpacing','tight','Padding','tight')

freqSet = [60,15];

for ff=1:length(freqSet)

    % Create a p vector with a 60 Hz corner frequency
    pMRI = [1 repmat([freqSet(ff) 1 0.15],1,nCells*nEccs)];

    % Get the modeled response
    response = returnResponse(pMRI,stimulusDirections,studiedEccentricites,freqsForPlotting);

    % Get the complex filter
    cellEquation = stageSecondOrderLP(freqSet(ff),1);

    nexttile
    vec = squeeze(response(3,1,:));
    vec = vec./max(vec);
    plot(log10(freqsForPlotting),vec,'-','Color',[0.5 0.5 0.5],'LineWidth',2);
    xlabel('log frequency')

    nexttile
    irfWindowSecs = 0.1;
    myFreqs = linspace(0,1000,201);
    ttfComplex = double(subs(cellEquation,myFreqs));
    [irf, sampleRate] = simpleIFFT( myFreqs, abs(ttfComplex), angle(ttfComplex));
    myTime = 0:sampleRate:(length(irf)-1)*sampleRate;
    [~,windowIdx] = min(abs(myTime-irfWindowSecs));
    irf = irf(1:windowIdx); myTime = myTime(1:windowIdx);
    irf = irf ./ max(irf);
    plot(myTime,irf,'-','Color',[0.5 0.5 0.5],'LineWidth',2);
    xlabel('time [msecs]')
end

% Save the plots
plotNamesPDF = 'parameterEffectIllo.pdf';
saveas(figHandles,fullfile(savePath,plotNamesPDF));



function [response,rfsAtEcc] = returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs)
% Assemble the response across eccentricity locations

nCells = 3;
nParams = 3;
blockLength = nParams*nCells;

for ee = 1:length(studiedEccentricites)

    % Assemble the sub parameters
    startIdx = (ee-1)*blockLength + 1 + 1;
    subP = [p(1) p(startIdx:startIdx+blockLength-1)];

    % Obtain the response at this eccentricity
    [ttfAtEcc{ee},rfsAtEcc{ee}] = returnTTFAtEcc(subP,stimulusDirections,studiedEccentricites(ee),studiedFreqs);

end

% Reshape the responses into the dimension stim x ecc x freqs
for ee = 1:length(studiedEccentricites)
    response(:,ee,:) = ttfAtEcc{ee};
end

end
