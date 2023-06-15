% close all

% Place to save figures
savePath = '~/Desktop/VSS 2023/';

% These variables define the subject names, stimulus directions.
subjectNames = {'HEROgka1','HEROasb1'};
subjects = {'gka','asb'};
stimulusDirections = {'LminusM','S','LMS'};
stimPlotColors = {'r','b','k'};
nSubs = length(subjects);
nStims = length(stimulusDirections);

% Fixed features of the model
nAcqs = 12; nCells = 3; nParams = 3;

% The frequencies studied. We also define a set of interpolated frequencies
% so that the model fit does not wiggle too much in between the studied
% frequencies. Finally, we define a high-resolution set of the frequencies
% for plotting.
allFreqs = [0,2,4,8,16,32,64];
studiedFreqs = [2 4 8 16 32 64];
nFreqs = length(studiedFreqs);
initialInterpFreqs = logspace(0,2,50);
fineInterpFreqs = logspace(0,2,500);



% Download Mt Sinai results

% Define the localSaveDir
localDataDir = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data');

% Load the retino maps
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localDataDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Save a template map variable so we can create new maps below
templateImage = cifti_read(tmpPath);

% This is the threshold for the goodness of fit to the fMRI time-series
% data. We only analyze those voxels with this quality fit or better
r2Thresh = 0.2;


% Loop through subjects and fit each vertex
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the stimLabels
    stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

    % Initialize or load the fitResults
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_fit_results.mat']);
    load(filePath,'fitResults')

    % Create a map of the fVal fit
    newMap = templateImage;
    newMap.cdata = single(fitResults.fVal);
    fileOut = fullfile(savePath,[subjectNames{ss} '_fVal.dtseries.nii']);
    cifti_write(newMap, fileOut);

    figure

    % Loop over stimulus directions
    for whichStim = 1:3

        % Find valid V1 voxels
        posRespIdx = zeros(size(results.R2));
        posRespIdx(fitResults.fVal > 0) = cellfun(@(x) sum(x(whichStim,:))>0, fitResults.Y(fitResults.fVal > 0));

        goodIdx = find(logical( (results.R2 > r2Thresh) .* (vArea==1) .* (fitResults.fVal<5) .* posRespIdx  ));
        nGood = length(goodIdx);

        peakFreq = [];
        for gg = 1:nGood

            k=interp1(freqsForPlotting,fitResults.yFit{goodIdx(gg)}(whichStim,:),fineInterpFreqs,'spline');
            [~,idx]=max(k);
            peakFreq(gg) = fineInterpFreqs(idx);

            if whichStim == 3 && peakFreq(gg) > 25
                foo = 1;
                %{
                fitResults.fVal(goodIdx(gg))
                hold off               
                plot(log10(studiedFreqs),fitResults.Y{goodIdx(gg)}','*')
                hold on
                 plot(log10(initialInterpFreqs),fitResults.yFit{goodIdx(gg)}','-')
                %}
            end

        end

        x = log10(fitResults.eccDeg(goodIdx));
        [x, sortedIdx] = sort(x);
        xq = 0:0.01:1.8;
        x(x<0) = 0;
        v = peakFreq(sortedIdx);
        plot(x,v,['.',stimPlotColors{whichStim}]);
        hold on
        sp = spaps(x,v,-150);
        vq = fnval(sp,xq);
        plot(xq,vq,['-' stimPlotColors{whichStim}],'LineWidth',3)

            % save a peakFreq map
    newMap = templateImage;
    newMap.cdata = single(zeros(size(fitResults.fVal)));
    newMap.cdata(goodIdx) = single(peakFreq);
    fileOut = fullfile(savePath,[subjectNames{ss} '_' stimulusDirections{whichStim} '_peakFreq.dtseries.nii']);
    cifti_write(newMap, fileOut);


    end

    %
    %
    % % Get the peak temporal freququency for the interpolated fit at each
    % % vertex
    % k = cell2mat(cellfun(@(x) interpFreqs(find(((x' == max(x'))))-[0; 50; 100]),fitResults.yFit(goodIdx),'UniformOutput',false));
    %
    % figure
    % x = log10(fitResults.eccDeg(goodIdx));
    % [x, sortedIdx] = sort(x);
    % xq = 0:0.01:1.8;
    % x(x<0) = 0;
    %      v = k(:,whichStim);
    %     v = v(sortedIdx);
    %     plot(x,v,['.',stimPlotColors{whichStim}]);
    %     hold on
    %         sp = spaps(x,v,-150);
    %         vq = fnval(sp,xq);
    %         plot(xq,vq,['-' stimPlotColors{whichStim}],'LineWidth',3)
    % end
    %
    % % Plot params
    % figure
    % paramNames = {'corner Freq','exponent','gain'};
    % sTol = [-100 -100 -1500];
    % for pp = 1:nParams
    %     subplot(1,nParams,pp)
    %
    %     for cc = 1:nCells
    %         v = fitResults.p(goodIdx,3*(cc-1)+pp+1);
    %         v = v(sortedIdx);
    %         plot(x,v,['.' stimPlotColors{cc}]);
    %         hold on
    %         sp = spaps(x,v,sTol(pp));
    %         vq = fnval(sp,xq);
    %         plot(xq,vq,['-' stimPlotColors{cc}],'LineWidth',3)
    %
    %     end
    %     title(paramNames{pp});
    %     if pp == 2
    %         refline(0,1);
    %     end
    %     if pp == 3
    %         a = gca();
    %         a.YScale = 'log';
    %     end
    % end


end
