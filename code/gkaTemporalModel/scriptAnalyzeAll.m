
%% Housekeeping
clear
close all

%% Create the RGC temporal sensitivity model
fitRGCFResponseData


%% Load the Mt. Sinai data
% This should be the V1 and LGN area bold fMRI signal mean, and 95% CI. The
% matrix is subject (GKA 1, ASB 2) x channel (L-M 1, S 2, LMS 3) x area
% (LGN 1, V1 2) x flicker freqency x bootstrap (1st value is 2.5%tile, 2nd
% value is 50%tile, 3rd value is 97.5%tile)
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','amplitudeResults','gka_asb_lgn_V1_BOLD.mat');
load(loadPath,'LGN_V1mri');

% Load the Mt. Sinai data
% This should be the V1 across eccentricity bold fMRI signal mean, and 95%
% CI. The matrix is subject (GKA 1, ASB 2) x channel (L-M 1, S 2, LMS 3) x
% eccentricity x flicker freqency x bootstrap (1st value is 2.5%tile, 2nd
% value is 50%tile, 3rd value is 97.5%tile)
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','amplitudeResults','gka_asb_V1_ecc_BOLD.mat');
load(loadPath,'V1ecc_mri');

% Do we want to conduct a search, or just use the p0 values and make plots?
searchFlag = true;

% Do we wish to use the monotonic constraint upon surround index in the
% search?
useMonotonicConstraint = false;

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the bins for the cortical
nEcc = 6;
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
eccDegVals = eccDegBinEdges(4:2:14);

% The identities of the stims and subjects
subjects = {'gka','asb'};
stimuli = {'L-M','S','LMS'};
plotColor = {'r','b','k'};
modelType = {'chromatic','chromatic','luminance'};

% The number of parameters in the model that are fixed across eccentricity.
% We need this information later when we make a plot of the surround index
% values across ecccentricity.
nFixed = 4;

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

% Loop over stimulus directions and subjects
for whichStim = 1:3
    for whichSub = 1:2

        % Extract the relevant LGN data
        lgnFreqX = studiedFreqs;
        lgnY = squeeze(LGN_V1mri(whichSub,whichStim,1,:,2))';
        lgnW = 1./(squeeze(LGN_V1mri(whichSub,whichStim,1,:,3))'-squeeze(LGN_V1mri(whichSub,whichStim,1,:,1))');

        % Extract the relevant V1 data acros eccentricities
        v1Eccentricity = []; v1FreqX = []; v1Y = []; v1W = [];
        for ee = 1:nEcc
            v1Eccentricity = [v1Eccentricity repmat(eccDegVals(ee),1,6)];
            v1FreqX = [v1FreqX studiedFreqs];
            v1Y = [v1Y squeeze(V1ecc_mri(whichSub,whichStim,ee,:,2))'];
            v1W = [v1W 1./(squeeze(V1ecc_mri(whichSub,whichStim,ee,:,3))'-squeeze(V1ecc_mri(whichSub,whichStim,ee,:,1))')];
        end

        % Stored p0 vectors
        switch whichStim
            case 1 % L-M chromatic
                switch whichSub
                    case 1
                        p0 = [0.0041   17.2462    0.9075   30.9342    0.7069    0.3565    0.2228    0.1369    0.1368    0.1111    0.0053    0.0110    0.0317    0.0320    0.0383    0.1294];
                    case 2
                        p0 = [0.0120   15.8006    0.5299   26.1496    0.4394    0.3332    0.1586    0.0069    0.0058    0.0043    0.0106    0.0202    0.0572    0.0700    0.1925    0.7719];
                end

            case 2 % S chromatic
                switch whichSub
                    case 1
                        p0 = [0.0030   15.8839    0.3990   28.3153    0.9261    0.4802    0.3507    0.1678    0.1367    0.1222    0.0016    0.0080    0.0323    0.0460    0.0663    0.1226];
                    case 2
                        p0 = [0.0026   12.9277    0.3248   19.9888    0.3326    0.3316    0.2216    0.0000    0.0000    0.0000    0.0028    0.0142    0.0652    0.0896    0.2406    0.8036];
                end

            case 3 % LMS luminance
                switch whichSub
                    case 1
                        p0 = [0.0460   25.8102    0.1533   11.3143    0.5645    0.5645    0.5166    0.4001    0.0822    0.0000    0.1506    0.2075    0.2673    0.2485    0.1404    0.2098];
                    case 2
                        p0 = [0.0445   19.9141    0.1284   15.1347    1.0000    0.7509    0.6211    0.4745    0.2943    0.1984    0.1363    0.2714    0.3897    0.3690    0.3764    0.7561];
                end
        end

        % Perform the search
        if searchFlag
            p = fitMRIResponseData(p0,v1FreqX, v1Eccentricity, v1Y, v1W, ...
                lgnFreqX, lgnY, lgnW, ...
                modelType{whichStim}, useMonotonicConstraint  );
        else
            p = p0;
        end

        % plot
        figure

        freqsForPlotting = logspace(0,2,50);
        v1Y = reshape(v1Y,6,1,nEcc);

        % Loop over eccentricities
        for ee=1:nEcc
            subplot(2,4,ee+(ee>3))
            semilogx(studiedFreqs,squeeze(v1Y(:,1,ee)),['o' plotColor{whichStim}]);
            hold on
            pBlock = [p(2:4) p(nFixed+ee) p(nFixed+nEcc+ee)];
            switch modelType{whichStim}
                case 'chromatic'
                    yFit = returnV1ChromEccTTFFit(pBlock,freqsForPlotting,eccDegVals(ee));
                case 'luminance'
                    yFit = returnV1LumEccTTFFit(pBlock,freqsForPlotting,eccDegVals(ee));
            end
            semilogx(freqsForPlotting,yFit,['-' plotColor{whichStim}]);
            refline(0,0);
            title([stimuli{whichStim} ', ' subjects{whichSub} ', ecc = ' num2str(eccDegVals(ee),2) '°']);
            ylim([-1 7]);
        end

        % Add the LGN response
        switch modelType{whichStim}
            case 'chromatic'
                lgnTTFFit = returnlgnLumTTFFit(p,freqsForPlotting,v1Eccentricity);
            case 'luminance'
                lgnTTFFit = returnlgnLumTTFFit(p,freqsForPlotting,v1Eccentricity);
        end
        subplot(2,4,8)
        semilogx(lgnFreqX,lgnY,['o' plotColor{whichStim}]);
        hold on
        semilogx(freqsForPlotting,lgnTTFFit,['-' plotColor{whichStim}]);
        refline(0,0);
        title([stimuli{whichStim} ', ' subjects{whichSub} ', LGN']);
        ylim([-0.5 4]);

        % Plot the surround suppression index vs. eccentricity
        subplot(2,4,4)
        plot(log10(eccDegVals),p(:,nFixed+1:nFixed+nEcc),['*' plotColor{whichStim}]);
        xlabel('Eccentricity [log deg]');
        ylabel('Suppression index');
        ylim([0 1]);

        % Save the plot
        plotName = [stimuli{whichStim} '_' subjects{whichSub} '_ModelFit.pdf' ];
        saveas(gcf,fullfile('~/Desktop',plotName));

    end

end
