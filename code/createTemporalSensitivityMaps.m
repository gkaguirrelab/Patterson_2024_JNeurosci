% Script that downloads the forwardModel results from the Mt Sinai flicker
% frequency experiments, and then fits a difference-of-exponentials model
% to the data. The resulting fits are saved back as maps.

% To find the analysis IDs, get the ID for the session (which is in the URL
% of the web GUI, and then use this command to get a list of the analyses
% associated with that session, and then find the analysis ID the we want.
%
%{
    toolboxName = 'flywheelMRSupport';
    fw = flywheel.Flywheel(getpref(toolboxName,'flywheelAPIKey'));
    sessionID = '5a1fa4b33b71e50019fd55dd';
    analysisList = fw.getSessionAnalyses(sessionID);
%}


% Save location for the maps
subjectNames = {'HEROgka1','HEROasb1'};
analysisIDs = { {'60ca690869059f3228c9a883','60ca68f6ba295e18031aaa35','60ca68e3f90bf6d5775e9e2b'} , ...
    {'60ca692f00f93080b0e14068', '60ca692379649f8717b0f6ec', '60ca69155f24eb23f2e14335'} };
retinoMapIDs = {'5dc88aaee74aa3005e169380','5dc88aaee74aa3005e169380' };
retinoFileNames = {'TOME_3021_cifti_maps.zip','TOME_3021_cifti_maps.zip'};

% Analysis parameters
scratchSaveDir = getpref('flywheelMRSupport','flywheelScratchDir');
scratchSaveDir = '/Users/aguirre/Desktop/tempFiles';
analysisLabels = {'LF','L-M','S'};

% Create the functional tmp save dir if it does not exist
saveDir = fullfile(scratchSaveDir,'v0','output');
if ~exist(saveDir,'dir')
    mkdir(saveDir);
end

% Create a flywheel object
fw = flywheel.Flywheel(getpref('forwardModelWrapper','flywheelAPIKey'));

% Loop over subjects
for ss = 1: length(subjectNames)
    
    % Set up the paths for this subject
    fileStem = [subjectNames{ss} '_agtcOL_'];
    resultsSaveDir = ['/Users/aguirre/Desktop/' subjectNames{ss}];
    mkdir(resultsSaveDir);
    
    % Download and unzip the retino maps
    fileName = retinoFileNames{ss};
    tmpPath = fullfile(saveDir,fileName);
    fw.downloadOutputFromAnalysis(retinoMapIDs{ss},fileName,tmpPath);
    command = ['unzip -q -n ' tmpPath ' -d ' saveDir];
    system(command);
    
    % Load the retino maps
    tmpPath = fullfile(saveDir,strrep(fileName,'_cifti_maps.zip','_inferred_varea.dtseries.nii'));
    vArea = cifti_read(tmpPath);
    vArea = vArea.cdata;
    tmpPath = fullfile(saveDir,strrep(fileName,'_cifti_maps.zip','_inferred_eccen.dtseries.nii'));
    eccenMap = cifti_read(tmpPath);
    eccenMap = eccenMap.cdata;
    tmpPath = fullfile(saveDir,strrep(fileName,'_cifti_maps.zip','_inferred_angle.dtseries.nii'));
    polarMap = cifti_read(tmpPath);
    polarMap = polarMap.cdata;
    tmpPath = fullfile(saveDir,strrep(fileName,'_cifti_maps.zip','_inferred_sigma.dtseries.nii'));
    sigmaMap = cifti_read(tmpPath);
    sigmaMap = sigmaMap.cdata;
    
    % Loop over analysis IDs
    for aa = 1:length(analysisIDs{ss})
        
        aa
        
        % Download the results file
        fileName = [fileStem 'results.mat'];
        tmpPath = fullfile(saveDir,[analysisLabels{aa} '_' fileName]);
        fw.downloadOutputFromAnalysis(analysisIDs{ss}{aa},fileName,tmpPath);
        
        % Load the result file into memory and delete the downloaded file
        clear results
        load(tmpPath,'results')
        %        delete(tmpPath)
        
        % Download the templateImage file
        fileName = [fileStem 'templateImage.mat'];
        tmpPath = fullfile(saveDir,[analysisLabels{aa} '_' fileName]);
        fw.downloadOutputFromAnalysis(analysisIDs{ss}{aa},fileName,tmpPath);
        
        % Load the result file into memory and delete the downloaded file
        clear templateImage
        load(tmpPath,'templateImage')
        %        delete(tmpPath)
        
        % Fit the TTF model
        [resultsFit,fieldNames] = fitTTFModel(results);
        
        % Save the TTF model fits
        outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_' analysisLabels{aa} '_ResultsFit.mat']);
        save(outFile,'resultsFit');
        
        % Add the original time-series R2 fit
        fieldNames = [fieldNames 'R2'];
        
        % Save the map results into images
        for ff = 1:length(fieldNames)
            % The initial, CIFTI space image
            outCIFTIFile = fullfile(resultsSaveDir, [subjectNames{ss} '_' analysisLabels{aa} '_' fieldNames{ff} '.dtseries.nii']);
            outData = templateImage;
            outData.cdata = single(resultsFit.(fieldNames{ff}));
            outData.diminfo{1,2}.length = 1;
            cifti_write(outData, outCIFTIFile)
        end
        
        % generate visual field maps
        outFileStem = fullfile(resultsSaveDir, [subjectNames{ss} '_' analysisLabels{aa} ]);
        makeVisualFieldMap(resultsFit,eccenMap,polarMap,vArea,sigmaMap,outFileStem);
        
        foo=1;
    end
    
end



function [results,fieldNames] = fitTTFModel(results)

%% Fit the difference-of-exponentials model
freqs = [2 4 8 16 32 64];
nFreqs = length(freqs);
freqsIdx = 1:nFreqs;
deltaF10 = min(diff(log10(freqs)));
fitScaleUp = 10;
freqsFit = 10.^(log10(min(freqs))-deltaF10+deltaF10/fitScaleUp:deltaF10/fitScaleUp:log10(max(freqs))+deltaF10);
freqsFitIdx = 1/fitScaleUp:1/fitScaleUp:nFreqs+1;

nV = size(results.params,1);

% Variables to hold the results
fitPeakAmp = nan(nV,1);
fitPeakFreq = nan(nV,1);
fitMaxFreq = nan(nV,1);
fitWidthFreqDB = nan(nV,1);
fit64HzResid = nan(nV,1);
fitR2 = nan(nV,1);
fitStuff = struct('p',{},'yVals',{},'myFit',{});
fitStuff_p = nan(nV,4);
fitStuff_yVals = nan(nV,nFreqs);
fitStuff_myFit = nan(nV,length(freqsFit));


% The fitting function is a non-central beta, that is further modified to
% allow adjustment of the bounded interval and scaling of the overall
% amplitude. The function is constrained to hold the first degree parameter
% to an arbitrarily small value, and to apply non-linear constraints in the
% fitting (described below).
%
% There is no particular theoretical motivation for using this fitting
% form. The fit does reflect the following expectations:
% - The amplitude of response will approach zero as the stimulus
%   frequency approaches 1 Hz.
% - The amplitude of the response will return to zero at higher frequencies
% - There is an enforced degree of smoothness in the change in the function
%   with frequency

myFunc = @(f,A,B,C,D)  C.*ncbeta(f./D, 1e-6, A, B );

% x0 and bounds
x0 = [1 1 1 8];
lb = [0 1 0 0];
ub = [100 100 100 8];

% define some search options
options = optimoptions(@fmincon,...
    'Diagnostics','off',...
    'Display','off');

% Loop through the vertices / voxels
parfor vv = 1:nV
    if results.R2(vv)>0.15
        
        % save the current warning status and silence anticipated warnings
        warningState = warning;
        warning('off','MATLAB:singularMatrix');
        warning('off','MATLAB:nearlySingularMatrix');
        warning('off','MATLAB:illConditionedMatrix');
                
        % Get the beta values
        yVals = results.params(vv,1:nFreqs+1);
        
        % The params have an explicit coding for the blank screen, so we
        % adjust for this
        yVals = yVals(2:end) - yVals(1);
        
        % Define the objective and non-linear constraint
        myObj = @(p) norm(yVals - myFunc(freqsIdx,p(1),p(2),p(3),p(4)));
        myNonlcon = @(p) betaNonlcon(p,yVals,freqsFitIdx);
        
        % Fit
        p = fmincon(myObj,x0,[],[],[],[],lb,ub,myNonlcon,options);
        R2 = corr(yVals',myFunc(freqsIdx,p(1),p(2),p(3),p(4))').^2;
        myFit = myFunc(freqsFitIdx,p(1),p(2),p(3),p(4));
        myFit(~isfinite(myFit))=nan;
        
        % Plot
        %{
        plot(1:nFreqs,yVals,'*k');
        hold on
        plot(freqsFitIdx,myFit,'-r');
        hold off
        pause
        %}
        %{
        plot(freqs,yVals,'*k');
        hold on
        plot(freqsFit,myFit,'-r');
        hold off
        pause
        %}
        
        % Store the values from this vertex
        [a,idx] = max(myFit);
        fitPeakAmp(vv) = a;
        fitPeakFreq(vv) = freqsFit(idx);
        [~,idx]=min(abs(freqsFitIdx-p(4)));
        fitMaxFreq(vv) = freqsFit(idx);
        [~,idx]=mink(abs(myFit-a/2),2);
        fitWidthFreqDB(vv) = abs(diff(idx));
        fit64HzResid(vv) = yVals(end)-myFunc(max(freqs),p(1),p(2),p(3),p(4));
        fitR2(vv) = R2;
        fitStuff_p(vv,:) = p;
        fitStuff_yVals(vv,:) = yVals;
        fitStuff_myFit(vv,:) = myFit;
        
        % Restore the warning state
        warning(warningState);
        
    end
end

% Place the analysis in the results variable
results.fitPeakAmp = fitPeakAmp;
results.fitPeakFreq = fitPeakFreq;
results.fitMaxFreq = fitMaxFreq;
results.fitWidthFreqDB = fitWidthFreqDB;
results.fit64HzResid = fit64HzResid;
results.fitR2 = fitR2;
results.fitSupport.freqs = freqs;
results.fitSupport.freqsFit = freqsFit;
results.fitStuff.p = fitStuff_p;
results.fitStuff.yVals = fitStuff_yVals;
results.fitStuff.myFit = fitStuff_myFit;

fieldNames = {'fitPeakAmp','fitPeakFreq','fitMaxFreq','fitWidthFreqDB','fit64HzResid','fitR2'};

end


function [c, ceq] = betaNonlcon(p,y,f)

% Evaluate the function
yFit = p(3).*ncbeta(f./p(4), 1e-6, p(1), p(2) );
yFit = yFit(isfinite(yFit));

% The rate of change should not exceed 1 unit
d = diff(yFit);
d = d(isfinite(d));
c = max(abs(d))-1;

% The interpolated peak should not be more than 25% the size of the
% largest y value
peak = (max(yFit)-max(y))./max(y) - 0.5;
if peak > 0
    ceq = peak;
else
    ceq = 0;
end

end


function [ pz ] = ncbeta (x, a, b, lam )
%This function computes the probability density function for the
%noncentral beta distribution using a transformation of variables to put
%the desired density function in terms of a noncentral F, which is included
%in Matlabs statsitics toolbox already.
%
% USAGES
% [ pz ] = ncbeta (x, a, b, lambda )
%
% INPUT
% x:    Vector of possible argument values for the pdf.
% a:    The first degree of freedom/shape parameter.
% b:    The second degree of freedom/shape parameter.
% lam:  The noncentrality parameter.
%
% OUTPUT
% pz:   The probability density function for x, given the input parameters.
%
%-----------------------------------------------------------------------
% Latest Edit: 18.Feb.2014
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------
const   = a./b;
pz      = @(r) ncfpdf( r./(const*(1-r)), a, b, lam).* 1./(const*(1-r).^2);
pz      = pz(x);

end



function makeVisualFieldMap(resultsFit,eccenMap,polarMap,vArea,sigmaMap,outFileStem)


% Identify the vertices with fits above the threshold
goodIdx = logical( (resultsFit.R2 > 0.15) .* (resultsFit.fitR2 > 0.8) .* (eccenMap>00).*(vArea==1).*(vArea==1) );

% Separate the polar angles for the left and right hemispheres
%polarMap(32492:end)=-polarMap(32492:end);

fieldNames = {'R2','fitPeakFreq','fitMaxFreq','fit64HzResid','fitWidthFreqDB'};
ranges = { [0.15 0.5], [0 20], [0 64], [-1 1], [0 40] };

for ff = 1:length(fieldNames)
    
    % Which result to plot
    thisField = fieldNames{ff};
    
    % Make the fieldMap
    figHandle = createFieldMap(resultsFit.(thisField)(goodIdx),polarMap(goodIdx),eccenMap(goodIdx),sigmaMap(goodIdx),30,true);
    
    % Clean up
    title(thisField);
    caxis(ranges{ff})
    outFile = [outFileStem '_fieldMap_ ' thisField '.pdf'];
    saveas(figHandle,outFile);
    close(figHandle)
end

end


