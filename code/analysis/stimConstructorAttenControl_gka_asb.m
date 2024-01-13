function stimConstructorAttenControl_gka_asb()
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
userName = getpref('Patterson_2024_JNeurosci', 'userName');

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

            % Get the vector of attention event times
            attentionEventTimes = inputMatrix{modality}{1, ii}.metaData.eventTimesArray/1000;

            % Create an empty stimulus matrix
            emptyMat = zeros(2*length(frequencies)+1,336);

            % Loop through frequencies in the current acquisition and determine
            % the starting and stopping point
            for freq = 1:length(frequencies)
                for item = frequencies{freq}
                    trialStart = (item-1)*11 + item; % Here we multiply by 11 as each freq is 12 sec long
                    if trialStart == 0
                        trialStart = 1;
                    end
                    trialEnd = trialStart+11;
                    attenFlag = any(~isnan(discretize(attentionEventTimes,[trialStart trialEnd])));
                    emptyMat(freq+attenFlag*length(frequencies), trialStart:trialEnd) = 1;
                end
            end
            
            % Get the attention events in each acquisition and place them in
            % the matrix
            for event = inputMatrix{modality}{1, ii}.metaData.eventTimesArray
                attentionVal = round(event/1000);
                emptyMat(2*length(frequencies)+1,attentionVal) = 1;
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
    
    % Create a combined stimulus matrix by concatenating the L-M, S, and LF
    % stimulus sequences, adjusting the index values for the stimulus
    % trials, and combining into a single vector the baseline (0Hz) and
    % attention trials
    stimulus = {};
    stimLabels = '(stimLabels),{';
    dirLabels = {'LMS','LminusM','S'};
    condLabels = {'f0Hz','f2Hz','f4Hz','f8Hz','f16Hz','f32Hz','f64Hz','Af0Hz','Af2Hz','Af4Hz','Af8Hz','Af16Hz','Af32Hz','Af64Hz','attention'};
    nConds = size(outputMatrix{1,1}{1},1);
    nCols = nConds * sum(cellfun(@(x) length(x),outputMatrix));
    offset = 0;
    for gg = [2 3 1]
        for ii=1:length(outputMatrix{1,gg})
            stimMat = zeros(nCols,336);
            stimMat(1+offset:nConds+offset,:) = outputMatrix{1,gg}{ii};
            offset = offset + nConds;
            stimulus{end+1}=stimMat;
            theseLabels = cellfun(@(x) ['(' x '_' dirLabels{gg} '_' sprintf('%0.2d',ii) '),'],condLabels,'UniformOutput',false);
            theseLabels = strcat(theseLabels{:});
            stimLabels = [stimLabels theseLabels];
        end
    end
    
    % Remove trailing comma from the stimLabels and cap with bracket
    stimLabels = [stimLabels(1:end-1) '}\n' ];
	fprintf(stimLabels);
        
    allmod = fullfile(outputFolder, ['stimulus_' subject '_all.mat']);
    save(allmod, 'stimulus')
    
    % This command outputs the indices to be used for averaging
    %{
        fprintf('{ '); for ii=1:6; fprintf(sprintf('[%d:%d,%d:%d,%d:%d,%d:%d,%d:%d,%d:%d],',[(ii-1)*336+1,ii*336,(ii+5)*336+1,(ii+6)*336], [(ii-1+12)*336+1,(ii+12)*336,(ii+5+12)*336+1,(ii+6+12)*336], [(ii-1+24)*336+1,(ii+24)*336,(ii+5+24)*336+1,(ii+6+24)*336] ));end; fprintf(' }\n');
    %}

end

end % main