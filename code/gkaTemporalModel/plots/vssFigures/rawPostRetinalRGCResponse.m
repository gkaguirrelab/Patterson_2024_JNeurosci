% scriptCreatePlots

% Housekeeping
clear

% Properties of which model to plot
freqsForPlotting = logspace(0,2,50);

% Load the RGC temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the V1 cortical bins
nEccs = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
studiedEccentricites = eccDegBinEdges(4:2:14);

stimulusDirections = {'LminusM','S','LMS'};
nStims = length(stimulusDirections);


% Fixed features of the model
nCells = 3; nParams = 3;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

% Params that allows the plots to appear in the order LMS, L-M, S
stimOrder = [2 3 1];

% The colors used for the plots
plotColor={[0.75 0.75 0.75],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k','r','b'};
faceAlpha = 0.4; % Transparency of the shaded error region
shift_ttf = [0 3 6 9 11 13]; % shifts each ttf down so they can be presented tightly on the same figure

% Create a p vector that produces no post-retinal modification

pMRI = [1 repmat([200 1 0.15],1,nCells*nEccs)];

% Get the modeled response
[response, rfsAtEcc] = returnResponse(pMRI,stimulusDirections,studiedEccentricites,freqsForPlotting,rgcTemporalModel);

% Prepare the figures
figHandles = figure('Renderer','painters');
figuresize(600,400,'pt');
tiledlayout(1,3,'TileSpacing','tight','Padding','tight')

% Loop over stimuli and plot
for whichStim = 1:nStims

    % Get the average response across eccentricity
    for ee=1:nEccs

        % Select the plot of the correct stimulus direction
        nexttile(stimOrder(whichStim));

        % Add the model fit
        plot(log10(freqsForPlotting),squeeze(response(whichStim,ee,:))-shift_ttf(ee),...
            ['-' lineColor{stimOrder(whichStim)}],...
            'LineWidth',2);
        hold on

        % If this is LMS, show the separate midget and parasol components
        lineSpecs = {'-','--'};
        if whichStim == 3
            for cc=1:2
                ttfComplex = double(subs(rfsAtEcc{ee}{3}(cc),freqsForPlotting));
                vec = abs(ttfComplex);
                plot(log10(freqsForPlotting),vec-shift_ttf(ee),...
                    lineSpecs{cc},'Color',[0.5 0.5 0.5],...
                    'LineWidth',2);
            end
        end

        % Add reference lines
        if ee==1 && whichStim == 3
            plot(log10([1 1]),[0 2],'-k');
        end
        plot(log10([1 2]),[0 0]-shift_ttf(ee),':k');
            plot(log10([50 100]),[0 0]-shift_ttf(ee),':k');

    end

end

% Clean up
for ss=1:nStims
    for ee = 1:nEccs
        nexttile(ss);
        xlim(log10([0.5 150]))
        ylim([-14 5])
        a=gca;
        a.YTick = [0,2,4,6];
        a.YTickLabel = {'0','2','4','6'};
        a.XTick = log10([2,4,8,16,32,64]);
        a.XTickLabel = {'2','4','8','16','32','64'};
        a.XTickLabelRotation = 0;
        a.XMinorTick = 'off';
        a.YAxis.Visible = 'off';
        box off
        if ss>1
            a.XAxis.Visible = 'off';
        end
    end
end

% Save the plots
plotNamesPDF = 'postRetinalResponseModel.pdf';
saveas(figHandles,fullfile(savePath,plotNamesPDF));



function [response,rfsAtEcc] = returnResponse(p,stimulusDirections,studiedEccentricites,studiedFreqs,rgcTemporalModel)
% Assemble the response across eccentricity locations

nCells = 3;
nParams = 3;
blockLength = nParams*nCells;

for ee = 1:length(studiedEccentricites)

    % Assemble the sub parameters
    startIdx = (ee-1)*blockLength + 1 + 1;
    subP = [p(1) p(startIdx:startIdx+blockLength-1)];

    % Obtain the response at this eccentricity
    [ttfAtEcc{ee},rfsAtEcc{ee}] = returnTTFAtEcc(subP,stimulusDirections,studiedEccentricites(ee),studiedFreqs,rgcTemporalModel);

end

% Reshape the responses into the dimension stim x ecc x freqs
for ee = 1:length(studiedEccentricites)
    response(:,ee,:) = ttfAtEcc{ee};
end

end
