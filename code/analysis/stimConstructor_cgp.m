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
seqA = [5 2 4 5 7 2 5 4 2 7 3 6 5 1 2 3 7 6 1 7 1 3 5 6 4 1 1 1];
seqB = [1 6 6 2 1 5 3 2 2 6 7 4 6 3 3 4 4 3 1 1 4 7 7 5 5 1 1 1];
freqVals = [0,2,4,8,16,32,64];
dirOrder = [2 3 1];
seqVals = {'A','B'};
dirVals = {'LightFlux','LminusM_wide','S_wide'};
nConds = length(condLabels);
nRows = nConds*length(stimFileList);
for ii=1:length(stimFileList)
    load(fullfile(stimFileList(ii).folder,stimFileList(ii).name),'results');
    dirOrderSet{ii} = results.thisDir;
    switch char([results.trialEvents.stimFreqHz])
        case char(arrayfun(@(x) freqVals(x),seqA))
            seqOrderSet{ii} = 'A';
        case char(arrayfun(@(x) freqVals(x),seqB))
            seqOrderSet{ii} = 'B';
        otherwise
            error('I do not recognize this sequence of stimuli')
    end
end
for dd = 1:length(dirOrder)
    for ss = 1:length(seqVals)
        stimFileIdx = find(...
            strcmp(dirOrderSet,dirVals{dirOrder(dd)})...
            .* ...
            strcmp(seqOrderSet,seqVals{ss}) ...
            );
        for ii=1:length(stimFileIdx)
            load(fullfile(stimFileList(stimFileIdx(ii)).folder,stimFileList(stimFileIdx(ii)).name),'results');
            stimMat = zeros(nRows,nTimeSamples);
            offset = (ii-1)*nConds+(ss-1)*nConds*6+(dd-1)*nConds*12;
            freqEvents = [results.trialEvents.stimFreqHz];
            for tt = 1:length(freqEvents)
                rowIdx = (tt-1)*timeDivs+1;
                freqIdx = find(freqVals==freqEvents(tt));
                stimMat(offset+freqIdx,rowIdx:rowIdx+timeDivs-1) = 1;
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
end

% Remove trailing comma from the stimLabels and cap with bracket
stimLabels = [stimLabels(1:end-1) '}\n' ];
fprintf(stimLabels);

% This command outputs the TR indices to be used for averaging
fprintf('{ '); for ii=1:6; fprintf(sprintf('[%d:%d,%d:%d,%d:%d,%d:%d,%d:%d,%d:%d],',[(ii-1)*336+1,ii*336,(ii+5)*336+1,(ii+6)*336], [(ii-1+12)*336+1,(ii+12)*336,(ii+5+12)*336+1,(ii+6+12)*336], [(ii-1+24)*336+1,(ii+24)*336,(ii+5+24)*336+1,(ii+6+24)*336] ));end; fprintf(' }\n');

% Save the results
save('stimulus_HERO_cgp1_allxAcq.mat', 'stimulus','stimTime')

end
