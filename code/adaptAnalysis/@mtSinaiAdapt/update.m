function x0 = update(obj,x,x0,floatSet,signal)
% Update x0 with the result of the search
%
% Syntax:
%   x0 = obj.update(x,x0,floatSet)
%
% Description:
%   Update the x0 parameters based upon the results of a stage of search.
%   Here we use regression to obtain the ampliude parameters of the model,
%   saving on time-consuming non-linear search.
%
% Inputs:
%   x                     - Vector of parameters that were the subject of
%                           the prior stage of search.
%   x0                    - Vector of the entire set of parameters of the
%                           model at the time of the start of the
%                           just-completed search
%   floatSet              - Vector of index values that indicate the
%                           mapping of vector x into vector x0 
%   signal                - 1 x time vector. The data to be fit.
%
% Outputs:
%   x0                    - Vector of the entire set of parameters of the
%                           model to be used in the next stage of search
%


% Obj variables
stimulus = obj.stimulus;
stimCols = size(stimulus,2);
stimAcqGroups = obj.stimAcqGroups;
stimTime = obj.stimTime;

% Update x0 with the passed x values
x0(floatSet) = x;

% Create the HRF
hrf = obj.flobsbasis*x(end-2:end)';

% Normalize the kernel to have unit area
hrf = hrf/sum(abs(hrf));

% Create a regression matrix for the HRF implied by the FLOBS params
X = zeros(size(obj.dataTime,1),stimCols);

for ss = 1:size(stimulus,2)
    
    % Zero all parameters except for this column of the stimulus matrix
    subx = x0;
    subx(1:obj.nGainParams) = 0;
    subx(ss) = 1;

    % Get the neural signal for this column, at the stimulus temporal
    % resolution
    neuralSignal = obj.neuralForward(subx);

    % Convolve the neuralSignal by the hrf, respecting acquisition boundaries
    fit = conv2run(neuralSignal,hrf,stimAcqGroups);
    
    % If the stimTime variable is not empty, resample the fit to match
    % the temporal support of the data.
    if ~isempty(stimTime)
        dataAcqGroups = obj.dataAcqGroups;
        dataTime = obj.dataTime;
        fit = resamp2run(fit,stimAcqGroups,stimTime,dataAcqGroups,dataTime);
    end
    
    % Apply the cleaning step
    fit = obj.clean(fit);
    
    % Add the vector to the design matrix
    X(:,ss) = fit;
    
end

% Obtain the beta values
b = X\signal;

% Assemble the updated x0
x0(1:obj.nGainParams) = b';

end