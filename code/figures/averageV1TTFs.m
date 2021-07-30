%% averageV1TTFs
%


% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% Define where we want to save these figures
resultsSaveDir = fullfile(localSaveDir,'Fig X - average V1 TTFs');
mkdir(resultsSaveDir);

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'L-M','S','LF'};
plotColors = {'r','b',[0.75 0.75 0.75]};

% Load the retino maps
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Create a "subcortical" map
subcorticalMap = zeros(size(vArea));
subcorticalMap(1:26298)=1;

% This is the threshold for the goodness of fit to the fMRI time-series
% data
r2Thresh = 0.25;
attentionThresh = 0;

% This is the visual area and eccentricity range to grab. The visual areas
% are: V1 = 1, V2 = 2, V3 = 3, hV4/LO = [4 5], MT/MST = [8 9]
area = 1;

% Some values used for fitting peak frequency
freqs = [2 4 8 16 32 64];
nFreqs = length(freqs);
deltaF10 = min(diff(log10(freqs)));
fitScaleUp = 10;
freqsFit = 10.^(log10(min(freqs))-deltaF10+deltaF10/fitScaleUp:deltaF10/fitScaleUp:log10(max(freqs))+deltaF10);

% Which TTF model to use
ttfModel = 'watson';

% x0 and bounds
switch ttfModel
    case 'watson'
        p0 = [1.72, 0.87, 0.013];
        lb = [0 0 0];
        ub = [5 1 0.025];
    case 'beta'
        p0 = [1.8518    5.2050    6.7683   11.7180];
        lb = [0 0 0 1];
        ub = [100 8 100 100];
end

% define some search options
options = optimoptions(@fmincon,...
    'Diagnostics','off',...
    'Display','off');

% Set up a variable to hold the model fit results
pVals = [];

% Create a figure
figHandle0 = figure();

% Loop through the subjects and directions
for ss = 1:2
    
    subplot(1,2,ss);
    
    allResults = cell(3,1);
    
    % Load the results files
    for dd = 1:3
        
        % Load the results file for this subject / direction
        filePath = fullfile(localSaveDir,'resultsFiles',[analysisLabels{dd} '_' subjectNames{ss} '_agtcOL_results.mat']);
        load(filePath,'results')
        
        % Store the results file in a cell arrray
        allResults{dd}=results;
        
    end
    
    % All vertices in V1 and
    % for which there was a response to the attention task above
    % threshold.
    goodIdxA = ( (vArea>=1) .* (vArea<=1) .* (eccenMap > eccenDivs(ee)) .* (eccenMap < eccenDivs(ee+1)) ) .* ...
        ( ((allResults{1}.attention-allResults{1}.f0Hz)>attentionThresh) .* ((allResults{2}.attention-allResults{2}.f0Hz)>attentionThresh) .* ((allResults{3}.attention-allResults{3}.f0Hz)>attentionThresh) );
    
    % Also require that at least one of the three modulation directions
    % was associated with an fMRI time-series fit above the R2
    % threshold.
    goodIdxB = (allResults{1}.R2 > r2Thresh) + (allResults{2}.R2 > r2Thresh) + (allResults{3}.R2 > r2Thresh);
    goodIdxB(goodIdxB>0)=1;
    goodIdx = logical(goodIdxA .* goodIdxB);
    
    % Plot
    for dd = 1:3
        
        % Get the beta values for these indices
        yVals = allResults{dd}.params(goodIdx,1:nFreqs+1);
        
        % Get the mean betas across vertices / voxels
        yVals = nanmean(yVals);
        
        % Obtain the set of beta valyes, relative to the baseline condition
        yVals = yVals(2:end)-yVals(1);
        
        % Define the objective and non-linear constraint
        switch ttfModel
            case 'watson'
                myObj = @(p) norm(yVals - watsonTTF(p,freqs));
                myNonlcon = [];
            case 'beta'
                myObj = @(p) norm(yVals - betaTTF(p,freqs));
                myNonlcon = @(p) betaTTF(p,freqsFit,yVals);
        end
        
        % Search
        p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
        
        % Get the fitted response
        switch ttfModel
            case 'watson'
                myFit = watsonTTF(p,freqsFit);
            case 'beta'
                myFit = betaTTF(p,freqsFit);
        end
        myFit(~isfinite(myFit))=nan;
        
        % Get the interpolated max
        yMax=max(myFit);
        
        % Add this fit to the plot
        semilogx(freqs,yVals./yMax,'o','Color',plotColors{dd})
        hold on
        semilogx(freqsFit,myFit./yMax,'-','Color',plotColors{dd})
        semilogx([1 64],[0 0],':k','LineWidth',1)
        ylim([-0.25 1]);
        
    end
%    axis off
    
    
end

% Save the TTF plot
outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_TTFbyEccenbyChannel.pdf']);
print(figHandle0,outFile,'-dpdf','-bestfit');

