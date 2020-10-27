function [stimStructCellArray] = loadStimStructCellArray(userName)
% function [stimStructCellAray] = fmriBFDM_LoadStimStructCellArray()
%
%

% Set Dropbox path used to load the stimulus specification files
if ismac
dropboxDataDir = fullfile('/Users', userName, '/Dropbox (Aguirre-Brainard Lab)/MELA_data_Through061317');
elseif isunix
dropboxDataDir = fullfile('/home', userName, '/Dropbox (Aguirre-Brainard Lab)/MELA_data_Through061317');    
elseif ispc
dropboxDataDir = fullfile('C:\Users', userName, '/Dropbox (Aguirre-Brainard Lab)/MELA_data_Through061317');     
end

% Define the stimulus file session directories
stimulusSessDirs = {...
    'fmriBlockFreqDirection7T/HERO_asb1/041416' ...
    'fmriBlockFreqDirection7T/HERO_asb1/041516' ...
    'fmriBlockFreqDirection7T/HERO_gka1/041416' ...
    'fmriBlockFreqDirection7T/HERO_gka1/041516' ...
    };

%stimulusSessDirs = {'fmriBlockFreqDirection7T/HERO_asb1/041416'};

% Define which sessions we'd like to merge
whichSessionsToMerge = {[1 2], [3 4]};
%whichSessionsToMerge = {[1]};

for ss = 1:length(stimulusSessDirs)
    
    % Extract some information about this session and put it into the
    % params variable that will be passed to MakeStrimStruct
    tmp = strsplit(stimulusSessDirs{ss}, '/');
    makeStimStructParams.sessionType = tmp{1};

    % Hack the name of the sessionType here because of some inconsistent
    % names on the DropBox directory between the top level of the directory
    % and the names of individual files
    makeStimStructParams.sessionType = tmp{1};
    makeStimStructParams.sessionType = 'HCLV_Photo_7T';

    makeStimStructParams.sessionObserver = tmp{2};
    makeStimStructParams.sessionDate = tmp{3};
    makeStimStructParams.packetType       = 'fMRI';
    makeStimStructParams.stimulusDir       = fullfile(dropboxDataDir, stimulusSessDirs{ss});
    nRuns = length(listdir(fullfile(makeStimStructParams.stimulusDir, 'MatFiles', '*.mat'), 'files'));
    
    % Display some useful information
    fprintf('>> Processing <strong>%s</strong> | <strong>%s</strong> | <strong>%s</strong>\n', makeStimStructParams.sessionType, makeStimStructParams.sessionObserver, makeStimStructParams.sessionDate);

    % Iterate over runs
    for ii = 1:nRuns
        fprintf('\t* Run <strong>%g</strong> / <strong>%g</strong>\n', ii, nRuns);
        % Further define the params
        makeStimStructParams.runNum           = ii;
        makeStimStructParams.stimulusFile     = fullfile(makeStimStructParams.stimulusDir, 'MatFiles', [makeStimStructParams.sessionObserver '-' makeStimStructParams.sessionType '-' num2str(ii, '%02.f') '.mat']);

        % Make the stimulus structure
        [preMergeStimStructCellArray{ss, ii}.values, ...
         preMergeStimStructCellArray{ss, ii}.timebase, ...
         preMergeStimStructCellArray{ss, ii}.metaData] = fmriBFDM_MakeStimStruct(makeStimStructParams);
    end
    fprintf('\n');
end

%% Merge sessions
NSessionsMerged = length(whichSessionsToMerge);
for mm = 1:NSessionsMerged
    mergeIdx = whichSessionsToMerge{mm};
    tempMerge = {preMergeStimStructCellArray{mergeIdx, :}};
    tempMerge = tempMerge(~cellfun('isempty', tempMerge));
    stimStructCellArray{mm} = tempMerge;
end
