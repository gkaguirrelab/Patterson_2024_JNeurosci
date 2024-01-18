
clear
close all

% Place to save figures
savePath = '~/Desktop/Patterson_2024_EccentricityFlicker/';

% Define the localDataDir
localDataDir = fullfile(tbLocateProjectSilent('Patterson_2024_JNeurosci'),'data');

% These variables define the subject names and stimulus directions
subjectNames = {'HEROgka1','HEROasb1','HEROcgp1'};
subjects = {'gka','asb','cgp'};

% The modelClass for the analysis
modelClass = 'mtSinai';

% The stimulus frequencies
allFreqs = [0,2,4,8,16,32,64];

% The stimulus directions
directions = {'chromatic','achromatic'};

% Color map
cmap = [ linspace(0,1,255);[linspace(0,0.5,127) linspace(0.5,0,128)];[linspace(0,0.5,127) linspace(0.5,0,128)]]';

% Create the diret and carry-over vector for a given stimulus direction
directSeqA = [5 2 4 5 7 2 5 4 2 7 3 6 5 1 2 3 7 6 1 7 1 3 5 6 4 1 1 1];
directSeqB = [1 6 6 2 1 5 3 2 2 6 7 4 6 3 3 4 4 3 1 1 4 7 7 5 5 1 1 1];
directSeq = [directSeqA directSeqB];
directSeq = repmat(directSeq,1,6);
carryVec = [];
for ii=1:length(directSeq)
    if ii==1
        carryVec(ii) = sub2ind([7 7],directSeq(ii),1);
    else
        carryVec(ii) = sub2ind([7 7],directSeq(ii),directSeq(ii-1));
    end
end

% Create the cary-over labels
for dd=1:2
    for ii=1:49
        [rr,cc] = ind2sub([7 7],ii);
        thisLabel = sprintf(['co_f%d->f%d_' directions{dd}],allFreqs(cc),allFreqs(rr));
        carryOverLabels{(dd-1)*49+ii} = thisLabel;
    end
end

% Loop through the subjects
for ss = 1:length(subjectNames)

    % Load the results file for this subject
    filePath = fullfile(localDataDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_avgV1_mtSinai_results.mat']);
    load(filePath,'results')

    % Grab the tr
    tr = results.meta.tr;

    % Obtain the residual V1 response
    residts = results.data.datats' - results.data.modelts';

    % Create the carry-over covariates
    stimulus = results.model.inputs{2};
    nTimePoints = size(stimulus{1},2);
    carryMat = zeros(49*2,nTimePoints*36);
    for dd = 1:2
        % Loop through the carry-over crossings. Create a vector that
        % is positive for the particular crossing, and negative for the
        % other crossings for that stimulus type
        for cc=1:49
            [currStim,priorStim] = ind2sub([7 7],cc);
            posIdx = carryVec == cc;
            negIdx = double((directSeq == currStim) .* ~posIdx);
            modVals = posIdx / sum(posIdx) - (negIdx / sum(negIdx));
            modVals = modVals ./ range(modVals);
            for tt=1:length(negIdx)
                rowIdx = (dd-1)*49+cc;
                switch dd
                    case 1
                        colIdx = (tt-1)*12 + 1;
                        carryMat(rowIdx,colIdx:colIdx+11) = modVals(tt);
                        colIdx = nTimePoints*12 + (tt-1)*12 + 1;
                        carryMat(rowIdx,colIdx:colIdx+11) = modVals(tt);
                    case 2
                        colIdx = (dd-1)*nTimePoints*12*2 + (tt-1)*12 + 1;
                        carryMat(rowIdx,colIdx:colIdx+11) = modVals(tt);
                end
            end
        end
    end

    % Apply the HRF to the carry-over matrix
    hrf = results.data.hrf;
    stimAcqGroups = sort(repmat(1:36,1,nTimePoints))';
    carryMatCov = conv2run(carryMat',hrf,stimAcqGroups)';

    % Resample the carryMatCov to match the data
    modelOpts = results.model.opts;
    if ~isempty(modelOpts{10})
        stimTime = cell2mat(modelOpts(10:end))';
        dataTime = repmat(0:tr:tr*335,1,36)';
        dataAcqGroups = sort(repmat(1:36,1,336))';
        for ii=1:size(carryMatCov,1)
            vec = carryMatCov(ii,:);
            resampVec = resamp2run(vec,stimAcqGroups,stimTime,dataAcqGroups,dataTime);
            X(ii,:) = resampVec;
        end
    else
        X = carryMatCov;
    end

    % Obtain the beta values
    b = X'\residts';

    % Get the fit
    fit = b'*X;

    % Report the R2
    fprintf([subjects{ss} ' - R2 carry-over model: %2.2f\n'],corr(residts',fit')^2)

    % reshape the carry-over beta values into a matrix
    for dd = 1:2
        subB = b((dd-1)*49+1:dd*49);
        for ii=1:49
            [rr,cc] = ind2sub([7 7],ii);
            coMatrix{ss,dd}(rr,cc)=subB(ii);
        end
    end

end

figHandleA = figure();
figuresize(400,600,'pt');
t = tiledlayout(3,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';

figHandleB = figure();
figuresize(400,200,'pt');
t = tiledlayout(1,2);
t.TileSpacing = 'compact';
t.Padding = 'compact';


for dd=2:-1:1
    for ss = 1:3
        figure(figHandleA)
        nexttile((ss-1)*2+dd);
        im = round(coMatrix{ss,dd}*2*128+128);
        image(im);
        colormap(cmap)
        axis square
        title([subjects{ss} ' - ' directions{dd} ]);
        drawnow
        if dd == 1
            plotCleanUp(allFreqs);
        else
            axis off
        end
        allMat(ss,:,:) = coMatrix{ss,dd};
    end

    figure(figHandleB)
    nexttile();
    meanMat = squeeze(mean(allMat));
    semMat = squeeze(std(allMat))/sqrt(3);
    im = round(meanMat*3*128+128);
    image(im);
    colormap(cmap)
    colormap(cmap)
    axis square
    starMat = abs(meanMat./semMat)>2;
    hold on
    for rr=1:7
        for cc=1:7
            if starMat(rr,cc)
                plot(cc,rr,'*w')
            end
        end
    end
    title(['Avg Sub - ' directions{dd} ]);
    if dd == 2
        plotCleanUp(allFreqs);
    else
        axis off
    end

end

% Save the figures
plotNamesPDF = 'Fig7-1_individSubjectCarryOver.pdf';
saveas(figHandleA,fullfile(savePath,plotNamesPDF));

plotNamesPDF = 'Fig7_avgSubjectCarryOver.pdf';
saveas(figHandleB,fullfile(savePath,plotNamesPDF));


%% LOCAL FUNCTION
function plotCleanUp(allFreqs)
a = gca();
a.XTick = 1:length(allFreqs);
a.YTick = 1:length(allFreqs);
a.XTickLabels = arrayfun(@(x) {num2str(x)},allFreqs);
a.YTickLabels = arrayfun(@(x) {num2str(x)},allFreqs);
a.XAxis.TickLength = [0 0];
a.YAxis.TickLength = [0 0];
xlabel('current stimulus [Hz]');
ylabel('prior stimulus [Hz]');
box off
end


