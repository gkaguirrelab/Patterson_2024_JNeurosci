% scriptCreatePlots

% Housekeeping
clear

modelType = 'stimulus';
paramSearch = 'full';

% Load the empirical RGC data
rcgData = loadRGCResponseData();

% Load the RGC model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','rgcTemporalModel.mat');
load(loadPath,'rgcTemporalModel');

% Load the MRI temporal model
loadPath = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))),'data','temporalModelResults','v1',modelType);
load(fullfile(loadPath,['mriFullResultSet_' paramSearch '.mat']),'mriFullResultSet');

savePath = fullfile('~','Desktop','mtSinaiTemporalModelPlots','MRIData_FullModel',modelType);

% Extract some meta info from the mriTemporalModel
studiedFreqs = mriFullResultSet.meta.studiedFreqs;
studiedEccentricites = mriFullResultSet.meta.studiedEccentricites;
subjects = mriFullResultSet.meta.subjects;
stimulusDirections = mriFullResultSet.meta.stimulusDirections;
paramCounts = mriFullResultSet.meta.paramCounts;
cellClasses = {'midget','bistratified','parasol'};
nEccs = length(studiedEccentricites);
nFreqs = length(studiedFreqs);
myFreqs = logspace(log10(1),log10(100),101);
nFreqsForPlotting = length(myFreqs);
nCells = length(cellClasses);
subjectLineSpec = {'-','-'};


% Params that control the plot appearance
spacing = 1;
stimOrder = [2 3 1];

plotColor={[0.85 0.85 0.85],[0.85 0.55 0.55],[0.75 0.75 1]};
lineColor={'k',[.5 0.25 0.25],[0.25 0.25 0.5]};
faceAlpha = [0,1,1];


for whichSub = 1:length(subjects)


    % Get the model params and data
    pMRI = mean(mriFullResultSet.(subjects{whichSub}).pMRI,1);
    pMRISEM = std(mriFullResultSet.(subjects{whichSub}).pMRI,0,1);
    v1Y = mean(mriFullResultSet.(subjects{whichSub}).v1Y,1);
    v1YSEM = std(mriFullResultSet.(subjects{whichSub}).v1Y,0,1);
    lgnY = mean(mriFullResultSet.(subjects{whichSub}).lgnY,1);
    lgnYSEM = std(mriFullResultSet.(subjects{whichSub}).lgnY,0,1);

    % Get the model fits
    [~,v1YFitMatrix,v1RFMatrix] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,myFreqs,rgcTemporalModel,paramCounts,modelType);

    % Zero-out the delayed surround effect
    pMRI([11:16,24:29,37:42])=0;
    [~,~,v1RFMatrixNoSurround] = assembleV1Response(pMRI,cellClasses,stimulusDirections,studiedEccentricites,myFreqs,rgcTemporalModel,paramCounts,modelType);

    figure('Name',subjects{whichSub});
    tiledlayout(nEccs,length(stimulusDirections))

    % Loop over eccentricities
    for ee=1:nEccs

        % Loop over stimuli and plot
        for whichStim = 1:length(stimulusDirections)

            % Initialize variables to hold the average V1 response
            v1AvgY = zeros(1,length(studiedFreqs));
            v1YAvglow = zeros(1,length(studiedFreqs));
            v1YAvghigh = zeros(1,length(studiedFreqs));
            v1AvgFit = zeros(1,length(myFreqs));
            v1AvgIRF = zeros(1,401);
            v1AvgIRFNoSurround = zeros(1,401);


            % The indices of the data to be plotted in the big vector
            v1DataIndices = 1+(whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nFreqs): ...
                (whichStim-1)*(nEccs*nFreqs)+(ee-1)*(nEccs)+nFreqs;

            % Assemble the average V1 response for plotting later
            v1AvgY = v1AvgY + v1Y(v1DataIndices)./length(studiedEccentricites);
            v1YAvglow = v1YAvglow + (v1Y(v1DataIndices)-v1YSEM(v1DataIndices))./length(studiedEccentricites);
            v1YAvghigh = v1YAvghigh + (v1Y(v1DataIndices)+v1YSEM(v1DataIndices))./length(studiedEccentricites);
            v1AvgFit = v1AvgFit + squeeze(v1YFitMatrix(whichStim,ee,:))'/length(studiedEccentricites);

            myIRFFreqs = linspace(0,1000,201);
            irfWindowSecs = 0.2;
            ttfComplex = double(subs(v1RFMatrix(whichStim,ee),myIRFFreqs));
            [irf, sampleRate] = simpleIFFT( myIRFFreqs, abs(ttfComplex), angle(ttfComplex));
            myTime = 0:sampleRate:(length(irf)-1)*sampleRate;
            [~,windowIdx] = min(abs(myTime-irfWindowSecs));
            irf = irf(1:windowIdx); myTime = myTime(1:windowIdx);
            v1AvgIRF = v1AvgIRF + irf;

            myIRFFreqs = linspace(0,1000,201);
            irfWindowSecs = 0.2;
            ttfComplex = double(subs(v1RFMatrixNoSurround(whichStim,ee),myIRFFreqs));
            [irf, sampleRate] = simpleIFFT( myIRFFreqs, abs(ttfComplex), angle(ttfComplex));
            myTime = 0:sampleRate:(length(irf)-1)*sampleRate;
            [~,windowIdx] = min(abs(myTime-irfWindowSecs));
            irf = irf(1:windowIdx); myTime = myTime(1:windowIdx);
            v1AvgIRFNoSurround = v1AvgIRFNoSurround + irf;


            nexttile

            % V1
            plot(v1AvgIRF./max(v1AvgIRF),'-k');
            hold on
            plot(v1AvgIRFNoSurround./max(v1AvgIRFNoSurround),'-r');
            if ee==1 && whichStim == 1
            xlabel('time [secs]');
            ylabel('relative response');
            end
            if ee == 1
                            title([stimulusDirections{whichStim} ' ' num2str(studiedEccentricites(ee)) ' deg']);
            else
                title([num2str(studiedEccentricites(ee)) ' deg']);
            end
            box off

            % semilogx(myFreqs,v1AvgFit,'-k');
            % hold on

            % Add the data symbols
            % semilogx(studiedFreqs,v1AvgY,'*');
        end

    end


end
