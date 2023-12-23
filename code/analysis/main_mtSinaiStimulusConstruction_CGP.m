function main_mtSinaiStimulusConstruction_CGP()
% This function constructs stimuli for mtSinai forwardModel Flywheel runs.
% Modified to work for the data collected in 2023 for HEROcgp1
%
% Syntax:
%  mtSinaiStimulusConstruction()
%
% Description:
%   This function constructs stimuli for mtSinai forwardModel Flywheel runs
%
% Inputs:
%   This function does not require any inputs.
%
% Optional key/value pairs:
%   NA
%
% Outputs:
%   NA
%

% These values account for the slight timing error in the stimulus result
% files created for the CGP scanning sessions
trialDur = 11.945;
timeDivs = 12;
deltaT = trialDur/timeDivs;
attentionEventTimeScale = 0.9955;
modelDur = 336;
nTimeSamples = ceil(modelDur/deltaT);

% Path to the directory of stimulus files
stimFileDir = '/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/_Papers/Patterson_2023_EccentricityFlicker/HERO_cpg1_2023/scanResultFiles';

% List of stimFiles
stimFileList = dir([stimFileDir '/*/*mat']);

stimulus = {};
stimTime = {};
stimLabels = '(stimLabels),{';
dirLabels = {'LMS','LminusM','S'};
condLabels = {'f0Hz','f2Hz','f4Hz','f8Hz','f16Hz','f32Hz','f64Hz','attention'};
freqVals = [0,2,4,8,16,32,64];
dirVals = {'LightFlux','LminusM_wide','S_wide'};
nConds = length(condLabels);
nRows = nConds*length(stimFileList);
offset = 0;
for ii=1:length(stimFileList)
    load(fullfile(stimFileList(ii).folder,stimFileList(ii).name),'results');
    dirOrderSet{ii} = results.thisDir;
end
for gg = [2 3 1]
    stimFileIdx = find(strcmp(dirOrderSet,dirVals{gg}));
    for ii=1:length(stimFileIdx)
        load(fullfile(stimFileList(stimFileIdx(ii)).folder,stimFileList(stimFileIdx(ii)).name),'results');
        stimMat = zeros(nRows,nTimeSamples);
        offset = (ii-1)*nConds;
        freqEvents = [results.trialEvents.stimFreqHz];
        for tt = 1:length(freqEvents)
            rowIdx = (tt-1)*timeDivs+1;
            freqIdx = find(freqVals==freqEvents(tt));
            stimMat(offset+freqIdx,rowIdx:rowIdx+timeDivs) = 1;
        end
        stimMat(offset+1,337:end) = 1;
        attentionEventTimes = [results.attentionEvents.eventTimeSecs];
        for tt = 1:length(attentionEventTimes)
            rowIdx = round(attentionEventTimes*attentionEventTimeScale/deltaT);
            stimMat(offset+nConds,rowIdx) = 1;
        end
        stimulus{end+1}=stimMat;
        thisDir = dirLabels{strcmp(results.thisDir,dirVals)};
        theseLabels = cellfun(@(x) ['(' x '_' thisDir '_' sprintf('%0.2d',ii) '),'],condLabels,'UniformOutput',false);
        theseLabels = strcat(theseLabels{:});
        stimLabels = [stimLabels theseLabels];

        % Create the stimTime vector
        stimTime{end+1} = 0:deltaT:(nTimeSamples-1)*deltaT;
    end
end

% Remove trailing comma from the stimLabels and cap with bracket
stimLabels = [stimLabels(1:end-1) '}\n' ];
fprintf(stimLabels);

% This command outputs the TR indices to be used for averaging
fprintf('{ '); for ii=1:6; fprintf(sprintf('[%d:%d,%d:%d,%d:%d,%d:%d,%d:%d,%d:%d],',[(ii-1)*336+1,ii*336,(ii+5)*336+1,(ii+6)*336], [(ii-1+12)*336+1,(ii+12)*336,(ii+5+12)*336+1,(ii+6+12)*336], [(ii-1+24)*336+1,(ii+24)*336,(ii+5+24)*336+1,(ii+6+24)*336] ));end; fprintf(' }\n');

% Save the results
save('stimulus_HEROcgp1_all.mat', 'stimulus','stimTime')

end
