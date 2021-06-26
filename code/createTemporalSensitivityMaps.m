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
    {'6048d45f2fc5506ee3c84f3e', '6048d45349868fea27c850ab', '6048d447171bd2f8468932a8'} };
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
    for aa = 2:length(analysisIDs{ss})
        
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
        
        % Fit the DoE model
        [resultsFit,fieldNames] = fitDoEModel(results);
        
        % Add the original time-series R2 fit
        resultsFit.initialR2 = results.R2;
        fieldNames = [fieldNames 'initialR2'];
        
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
        makeVisualFieldMap(resultsFit,eccenMap,polarMap,vArea,sigmaMap);
        
        foo=1;
    end
    
end



function [results,fieldNames] = fitDoEModel(results)

%% Fit the difference-of-exponentials model
nFreqs = 6;
freqs = [2 4 8 16 32 64];
freqsFit = logspace(log10(1),log10(128),1000);

nV = size(results.params,1);

% Variables to hold the results
fitPeakAmp = nan(nV,1);
fitPeakFreq = nan(nV,1);
fitOffset = nan(nV,1);
maxF = 7;

%myFunc =  @(f,A,B,C,D,E) C.*evpdf( ((E+f)./D),A,B);
%myFunc =  @(f,A,B,C,D,E) C.*A.*B.*((f+D)/E).^(A-1).*(1-((f+D)/E).^(A)).^(B-1);
myFunc = @(f,A,B,C,D)  D.*ncbeta(f./maxF, A, B, C );

[~,idxSet]=maxk(results.R2,10);

figure
for ii=1:9
subplot(3,3,ii)
idx = idxSet(ii);
y = results.params(idx,2:7)-min(results.params(idx,1));
myObj = @(p) norm(y - myFunc(1:nFreqs,p(1),p(2),p(3),p(4)));
p=fmincon(myObj,[1,1,1,1]);
myObj(p)
plot(1:nFreqs,y,'*k');
hold on
plot(0:0.1:nFreqs+1,myFunc(0:0.1:nFreqs+1,p(1),p(2),p(3),p(4)),'-r');
end

% Loop through the vertices / voxels
for vv = 1:nV
    if results.R2(vv)>0.25
        
        % Get the beta values
        yVals = results.params(vv,1:nFreqs+1);
        
        % The params have an explicit coding for the blank screen, so we
        % adjust for this
        yVals = yVals(2:end) - yVals(1);
        
        % Handle a negative offset
        if min(yVals)<0
            offset = min(yVals);
            yVals = yVals-offset;
        else
            offset = 0;
        end
                        
        myFit = spline(0:nFreqs-1,yVals,linspace(0,nFreqs-1,1000));
        
        [a,idx] = max(myFit);
        fitPeakAmp(vv) = a+offset;
        fitPeakFreq(vv) = freqsFit(idx);
        fitOffset(vv) = offset;
        fitStuff(vv).yVals = yVals;
        fitStuff(vv).myFit = myFit;
    end
end

% Place the analysis in the results variable
results.fitPeakAmp = fitPeakAmp;
results.fitPeakFreq = fitPeakFreq;
results.fitOffset = fitOffset;
results.fitSupport.freqs = freqs;
results.fitSupport.freqsFit = freqsFit;
results.fitStuff = fitStuff;
fieldNames = {'fitPeakAmp','fitPeakFreq','fitOffset'};

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


function [results,fieldNames] = fitSplineModel(results)

%% Fit the difference-of-exponentials model
nFreqs = 6;
freqs = [2 4 8 16 32 64];
freqsFit = logspace(log10(freqs(1)),log10(freqs(end)),1000);

nV = size(results.params,1);

% Variables to hold the results
fitPeakAmp = nan(nV,1);
fitPeakFreq = nan(nV,1);
fitOffset = nan(nV,1);

% Loop through the vertices / voxels
for vv = 1:nV
    if results.R2(vv)>0.25
        
        % Get the beta values
        yVals = results.params(vv,1:nFreqs+1);
        
        % The params have an explicit coding for the blank screen, so we
        % adjust for this
        yVals = yVals(2:end) - yVals(1);
        
        % Handle a negative offset
        if min(yVals)<0
            offset = min(yVals);
            yVals = yVals-offset;
        else
            offset = 0;
        end
                        
        myFit = spline(0:nFreqs-1,yVals,linspace(0,nFreqs-1,1000));
        
        [a,idx] = max(myFit);
        fitPeakAmp(vv) = a+offset;
        fitPeakFreq(vv) = freqsFit(idx);
        fitOffset(vv) = offset;
        fitStuff(vv).yVals = yVals;
        fitStuff(vv).myFit = myFit;
    end
end

% Place the analysis in the results variable
results.fitPeakAmp = fitPeakAmp;
results.fitPeakFreq = fitPeakFreq;
results.fitOffset = fitOffset;
results.fitSupport.freqs = freqs;
results.fitSupport.freqsFit = freqsFit;
results.fitStuff = fitStuff;
fieldNames = {'fitPeakAmp','fitPeakFreq','fitOffset'};

end




function figHandle = makeVisualFieldMap(resultsFit,eccenMap,polarMap,vArea,sigmaMap)


%figHandle = figure('visible','on');
%set(figHandle,'PaperOrientation','landscape');
%set(figHandle,'PaperUnits','normalized');
%set(gcf,'Units','points','Position',[100 100 400 400]);

% Identify the vertices with fits above the threshold
goodIdx = logical((resultsFit.R2 > 0.25).*( (resultsFit.fitOffset./resultsFit.fitPeakAmp) < 0.1).*(eccenMap>0.1).*(vArea==1).*(vArea==1));

% Left and right hemisphere
%polarMap(32492:end)=-polarMap(32492:end);

% % Map R2 value to a red-green plot symbol color
nColors = 200;
[mycolormap] = make_colormap(nColors);

valMin = 0;
valMax = 16;
val = resultsFit.fitPeakFreq(goodIdx);
ind = floor((nColors-1).* (val-valMin)./(valMax-valMin))+1;

% valMin = 0;
% valMax = 10;
% val = resultsDoE.fitPeakAmpDoE(goodIdx) ./ abs(resultsDoE.fitOffset(goodIdx));
% ind = floor((nColors-1).* (val-valMin)./(valMax-valMin))+1;

ind(ind>nColors)=nColors;
ind(ind<1)=1;
colorTriple = mycolormap(ind,:);


markSize = sigmaMap(goodIdx).*100;
markSize(markSize==0) = 25;

% The polar angle map has the dorsal V1 border at +180, and the ventral V1
% border at 0. By adding 90 degrees to the polar angle, this places the
% dorsal V1 border at 270, which on the matlab polar scatter corresponds to
% the 6 o'clock position on the plot. Therefore, the polarscatter plot is
% in visual field coordinates(i.e., 6 o'clock represents the inferior
% vertical meridian of the visual field).
%h = polarscatter(deg2rad(polarMap(goodIdx)+90),eccenMap(goodIdx),markSize, ...
%    colorTriple,'filled','o','MarkerFaceAlpha',1/8);

%rlim([0 60])
% 
 fieldMap = createFieldMap(resultsFit.fitPeakFreq(goodIdx),polarMap(goodIdx),eccenMap(goodIdx),sigmaMap(goodIdx));
% fieldMap = createFieldMap(resultsFit.fitPeakAmp(goodIdx),polarMap(goodIdx),eccenMap(goodIdx),sigmaMap(goodIdx));

% subplot(2,1,1)
% plot(eccenMap(goodIdx),resultsFit.fitPeakAmp(goodIdx),'.r');
% ylim([0 50]);
% xlim([0 80]);
% subplot(2,1,2)
% plot(eccenMap(goodIdx),resultsFit.fitPeakFreq(goodIdx),'.r');
% ylim([0 64]);
% xlim([0 80]);

end

function [mycolormap] = make_colormap(mapres)



%% Create colormap
mycolormap = zeros(mapres,3);
mycolormap(:,1)=1;

% red to yellow
mycolormap(:,2) = linspace(0,1,mapres);

end

