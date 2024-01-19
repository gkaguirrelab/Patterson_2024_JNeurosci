
% These values account for the slight timing error in the stimulus result
% files created for the CGP scanning sessions
trialDur = 11.945;
timeDivs = 12;
deltaT = trialDur/timeDivs;
attentionEventTimeScale = 0.9955;
modelDur = 336;
nTimeSamples = ceil(modelDur/deltaT);

% Path to the directory of stimulus files
stimFileDir = '/Users/aguirre/Aguirre-Brainard Lab Dropbox/Geoffrey Aguirre/_Papers/Patterson_2024_EccentricityFlicker/HERO_cpg1_2023/scanResultFiles';

% List of stimFiles
stimFileList = dir([stimFileDir '/*/*mat']);

% Tally up the correct detections
nTrials = 0;
nResponses = 0;
for ii=1:length(stimFileList)
    load(fullfile(stimFileList(ii).folder,stimFileList(ii).name),'results');
    nTrials = nTrials + length([results.attentionEvents.responseTimeSecs]);
    nResponses = nResponses + sum(~isnan([results.attentionEvents.responseTimeSecs]));
end

% Report the values
fprintf('Attention performance for CGP was %2.1f % correct\d',100*(nResponses/nTrials));