function signal = neuralForward(obj, x)
% Neural forward model. Provides the predicted neural signal for a set of
% parameters
%
% Syntax:
%   signal = obj.neuralForward(x)
%
% Description:
%   Returns a time-series vector that is the predicted neural response to
%   the stimulus, based upon the parameters provided in x.
%
% Inputs:
%   x                     - 1xnParams vector.
%
% Optional key/value pairs:
%   none
%
% Outputs:
%   signal                   - 1xtime vector.
%
% Obj variables
stimulus = obj.stimulus;
stimDeltaT = obj.stimDeltaT;
nGainParams = obj.nGainParams;
nAdaptParams = obj.nAdaptParams;
stimDelatTsPerTrial = obj.stimDelatTsPerTrial;
stimClassSet = obj.stimClassSet;

% Temporal support for a single trial
stimTimeTrial = 0:stimDeltaT:(stimDelatTsPerTrial-1)*stimDeltaT;

% A vector that provides an index for each separate trial across the entire
% concatenated time series
trialAcqGroups = cell2mat(arrayfun(@(x) repmat(x,1,stimDelatTsPerTrial),1:size(stimulus,1)/stimDelatTsPerTrial,'UniformOutput',false));

% The indices within the x parameter vector that hold the tau and asymptote
% parameters
tauIdx = nGainParams+1:1:nGainParams+nAdaptParams;

% Initialize an empty neural signal, and then loop over the stimulus
% classes for which we have an assigned tau exponential decay parameter.
neuralSignal = zeros(size(stimulus,1),1);
for ss = 1:length(stimClassSet)

    % Get the modeled amplitudes for this stimClass
    stimParamIdx = find(contains(obj.stimLabels,stimClassSet{ss}));
    thisSignal = stimulus(:,stimParamIdx)*x(stimParamIdx)';

    % Create an exponential kernel and apply this to each trial
    tau = x(tauIdx(ss));
    exponentialKernel = exp(-1/tau*stimTimeTrial);
    exponentialKernel = (exponentialKernel-mean(exponentialKernel))+1;
    thisSignal = conv2run(thisSignal,exponentialKernel,trialAcqGroups);

    % Add this signal to the neuralSignal
    neuralSignal = neuralSignal + thisSignal;
end

% Now add in the stimulus components for which we did not have an assigned
% tau parameter.
stimParamIdx = find(~contains(obj.stimLabels,horzcat(stimClassSet{:})));
thisSignal = stimulus(:,stimParamIdx)*x(stimParamIdx)';
signal = neuralSignal + thisSignal;

end

