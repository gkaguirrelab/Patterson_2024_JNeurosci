function main_mtSinaiStimulusConstruction()
% This function constructs stimuli for mtSinai forwardModel Flywheel runs.
%
% Syntax:
%  mtSinaiStimulusConstruction()
%
% Description:
%   This function constructs stimuli for mtSinai forwardModel Flywheel runs
%
% Inputs:
%   This function does not require any inputs. It gets all from the local
%   hook.
%
% Optional key/value pairs:
%   NA
%
% Outputs:
%   NA
%

% Get the username
userName = getpref('mriSinaiAnalysis', 'userName');

% Get the flywheel key
flywheelAPIkey = getpref('flywheelMRSupport', 'flywheelAPIKey');

% Get the stimulus file from Mela Data
stimStructCellArray = loadStimStructCellArray(userName);

% Path to the dataDir within dropbox
subjectNames = {'HERO_asb1' 'HERO_gka1'};

for ss = 1:length(subjectNames)
    subject = subjectNames{ss};
    
    % Get the subjects
    if strcmp(subject, 'HERO_asb1')
        sub = stimStructCellArray{1};
    elseif strcmp(subject, 'HERO_gka1')
        sub = stimStructCellArray{2};
    else
        error('Subject is not recognized. Enter gka or asb')
    end
    
    %% Construct subject
    % Create empty cells for dividing 3 modulations
    LightFluxSessions = {};
    LMinusMSessions = {};
    SSessions = {};
    
    % Loop through sessions and group them together
    lengthSes = length(sub);
    
    % First append the condition A session 1
    for ii = 1:lengthSes
        if sub{ii}.metaData.modulationDirection == "LightFlux" && sub{ii}.metaData.stimulusOrderAorB == 'A' && strcmp(sub{ii}.metaData.sessionDate, '041416')
            LightFluxSessions{end+1} = sub{ii};
        elseif sub{ii}.metaData.modulationDirection == "L-M" && sub{ii}.metaData.stimulusOrderAorB == 'A' && strcmp(sub{ii}.metaData.sessionDate, '041416')
            LMinusMSessions{end+1} = sub{ii};
        elseif sub{ii}.metaData.modulationDirection == "S" && sub{ii}.metaData.stimulusOrderAorB == 'A' && strcmp(sub{ii}.metaData.sessionDate, '041416')
            SSessions{end+1} = sub{ii};
        end
    end
    % condition A session 2
    for ii = 1:lengthSes
        if sub{ii}.metaData.modulationDirection == "LightFlux" && sub{ii}.metaData.stimulusOrderAorB == 'A' && strcmp(sub{ii}.metaData.sessionDate, '041516')
            LightFluxSessions{end+1} = sub{ii};
        elseif sub{ii}.metaData.modulationDirection == "L-M" && sub{ii}.metaData.stimulusOrderAorB == 'A' && strcmp(sub{ii}.metaData.sessionDate, '041516')
            LMinusMSessions{end+1} = sub{ii};
        elseif sub{ii}.metaData.modulationDirection == "S" && sub{ii}.metaData.stimulusOrderAorB == 'A' && strcmp(sub{ii}.metaData.sessionDate, '041516')
            SSessions{end+1} = sub{ii};
        end
    end
    % Now append the condition B session 1
    for ii = 1:lengthSes
        if sub{ii}.metaData.modulationDirection == "LightFlux" && sub{ii}.metaData.stimulusOrderAorB == 'B' && strcmp(sub{ii}.metaData.sessionDate, '041416')
            LightFluxSessions{end+1} = sub{ii};
        elseif sub{ii}.metaData.modulationDirection == "L-M" && sub{ii}.metaData.stimulusOrderAorB == 'B' && strcmp(sub{ii}.metaData.sessionDate, '041416')
            LMinusMSessions{end+1} = sub{ii};
        elseif sub{ii}.metaData.modulationDirection == "S" && sub{ii}.metaData.stimulusOrderAorB == 'B' && strcmp(sub{ii}.metaData.sessionDate, '041416')
            SSessions{end+1} = sub{ii};
        end
    end
    % condition B session 2
    for ii = 1:lengthSes
        if sub{ii}.metaData.modulationDirection == "LightFlux" && sub{ii}.metaData.stimulusOrderAorB == 'B' && strcmp(sub{ii}.metaData.sessionDate, '041516')
            LightFluxSessions{end+1} = sub{ii};
        elseif sub{ii}.metaData.modulationDirection == "L-M" && sub{ii}.metaData.stimulusOrderAorB == 'B' && strcmp(sub{ii}.metaData.sessionDate, '041516')
            LMinusMSessions{end+1} = sub{ii};
        elseif sub{ii}.metaData.modulationDirection == "S" && sub{ii}.metaData.stimulusOrderAorB == 'B' && strcmp(sub{ii}.metaData.sessionDate, '041516')
            SSessions{end+1} = sub{ii};
        end
    end
    
    % For asb1, replace the original LigtFLuxB_run1 and SA_run2 with the 
    % extra ones that were acquired in the second session.
    if strcmp(subject, 'HERO_asb1')
        LightFluxSessions(7) = LightFluxSessions(13);
        LightFluxSessions(13) = [];
        
        SSessions(2) = SSessions(7);
        SSessions(7) = [];
    end    
    
    % Combine the input matrices
    inputMatrix = {LightFluxSessions LMinusMSessions SSessions};
    
    % Construct the modulation matrices and combine
    lightFlux = {};
    lMinusM = {};
    s = {};
    outputMatrix = {lightFlux lMinusM s};
    
    % Loop through input matrices
    for modality = 1:length(inputMatrix)
        
        % Loop through acquisitions
        for ii = 1:length(inputMatrix{modality})
            
            % Create an empty stimulus matrix 8x336
            emptyMat = zeros(8,336);
            
            % Get the order of different Hz in the current acquisition
            zeroHz = find(inputMatrix{modality}{1, ii}.metaData.stimTypes == 1);
            TwoHz = find(inputMatrix{modality}{1, ii}.metaData.stimTypes == 2);
            FourHz = find(inputMatrix{modality}{1, ii}.metaData.stimTypes == 3);
            EightHz = find(inputMatrix{modality}{1, ii}.metaData.stimTypes == 4);
            SixteenHz = find(inputMatrix{modality}{1, ii}.metaData.stimTypes == 5);
            ThirtytwoHz = find(inputMatrix{modality}{1, ii}.metaData.stimTypes == 6);
            SixtyfourHz = find(inputMatrix{modality}{1, ii}.metaData.stimTypes == 7);
            
            % Combine the frequencies in cell
            frequencies = {zeroHz' TwoHz' FourHz' EightHz' SixteenHz' ThirtytwoHz' SixtyfourHz'};
            
            % Loop through frequencies in the current acquisition and determine
            % the starting and stopping point
            for freq = 1:length(frequencies)
                for item = frequencies{freq}
                    start = (item-1)*11 + item; % Here we multiply by 11 as each freq is 12 sec long
                    if start == 0
                        start = 1;
                    end
                    stop = start+11;
                    emptyMat(freq, start:stop) = 1;
                end
            end
            
            % Get the attention events in each acquisition and place them in
            % the matrix
            for event = inputMatrix{modality}{1, ii}.metaData.eventTimesArray
                attentionVal = round(event/1000);
                emptyMat(length(frequencies)+1,attentionVal) = 1;
            end
            
            % Put the final acquisition matrix in the output matrix
            outputMatrix{modality}{1,ii} = emptyMat;
        end
    end
    
    % Assign the cell values to stimulus and save them separately
    userName = getenv('username');
    if ismac || isunix
        outputFolder = '/tmp/outputStimSinai';
    elseif ispc
        outputFolder = fullfile('C:\Users', userName, 'outputStimSinai');
    end
    makeFolder = ['mkdir' ' ' outputFolder];
    system(makeFolder)
    
    % Save lightflux
    stimulus = outputMatrix{1,1};
    lfmod = fullfile(outputFolder, ['stimulus_' subject '_LightFlux.mat']);
    save(lfmod, 'stimulus')
    
    % Save LM
    stimulus = outputMatrix{1,2};
    lmmod = fullfile(outputFolder, ['stimulus_' subject '_L-M.mat']);
    save(lmmod, 'stimulus')
    
    % S
    stimulus = outputMatrix{1,3};
    smod = fullfile(outputFolder, ['stimulus_' subject '_S.mat']);
    save(smod, 'stimulus')
    
    % Save using the CLI
    system(['fw login' ' ' flywheelAPIkey ';' ' ' 'fw upload' ' ' lfmod ' ' 'fw://gkaguirrelab/mtSinaiFlicker/'])
    system(['fw login' ' ' flywheelAPIkey ';' ' ' 'fw upload' ' ' lmmod ' ' 'fw://gkaguirrelab/mtSinaiFlicker/'])
    system(['fw login' ' ' flywheelAPIkey ';' ' ' 'fw upload' ' ' smod ' ' 'fw://gkaguirrelab/mtSinaiFlicker/'])

    % Create a combined stimulus matrix by concatenating the L-M, S, and LF
    % stimulus sequences, adjusting the index values for the stimulus
    % trials, and combining into a single vector the baseline (0Hz) and
    % attention trials
    stimulus = {};
    for ii=1:length(outputMatrix{1,2})
        stimMat = zeros(20,336);
        stimMat(1,:) = outputMatrix{1,2}{ii}(1,:);
        stimMat(20,:) = outputMatrix{1,2}{ii}(8,:);
        stimMat(2:7,:) = outputMatrix{1,2}{ii}(2:7,:);
        stimulus{end+1}=stimMat;
    end
    for ii=1:length(outputMatrix{1,3})
        stimMat = zeros(20,336);
        stimMat(1,:) = outputMatrix{1,3}{ii}(1,:);
        stimMat(20,:) = outputMatrix{1,3}{ii}(8,:);
        stimMat(8:13,:) = outputMatrix{1,3}{ii}(2:7,:);
        stimulus{end+1}=stimMat;
    end
    for ii=1:length(outputMatrix{1,1})
        stimMat = zeros(20,336);
        stimMat(1,:) = outputMatrix{1,1}{ii}(1,:);
        stimMat(20,:) = outputMatrix{1,1}{ii}(8,:);
        stimMat(14:19,:) = outputMatrix{1,1}{ii}(2:7,:);
        stimulus{end+1}=stimMat;
    end
    allmod = fullfile(outputFolder, ['stimulus_' subject '_all.mat']);
    save(allmod, 'stimulus')
    system(['fw login' ' ' flywheelAPIkey ';' ' ' 'fw upload' ' ' allmod ' ' 'fw://gkaguirrelab/mtSinaiFlicker/'])
    

end

end % main