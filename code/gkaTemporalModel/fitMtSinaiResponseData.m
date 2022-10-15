% Fit the Mt Sinai fMRI data with a model that starts with the RGC temporal
% sensitivity functions

% Plotting mode. Copy and paste this to the console
%{
    searchFlag = false;
    for whichStim = 1:3; for whichSub = 1:2; fitMtSinaiResponseData; end; end
%}

% Search mode. Copy and paste this to the console
%{
    searchFlag = true;
    whichStim = 1; whichSub = 1;
    fitMtSinaiResponseData
%}

% The identities of the stims and subjects
stimuli = {'L-M','S','LMS'};
subjects = {'gka','asb'};

% The number of parameters in the model that are fixed across eccentricity
nFixed = 4;

% Load the Mt. Sinai data
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

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

% Extract the relevant LGN data
lgnFreqX = studiedFreqs;% studiedFreqs];
lgnRedGreenY = squeeze(LGN_V1mri(whichSub,1,1,:,2))';
lgnRedGreenW = 1./(squeeze(LGN_V1mri(whichSub,1,1,:,3))'-squeeze(LGN_V1mri(whichSub,1,1,:,1))');
lgnBlueYellowY = squeeze(LGN_V1mri(whichSub,2,1,:,2))';
lgnBlueYellowW = 1./(squeeze(LGN_V1mri(whichSub,2,1,:,3))'-squeeze(LGN_V1mri(whichSub,2,1,:,1))');
lgnLumY = squeeze(LGN_V1mri(whichSub,3,1,:,2))';
lgnLumW = 1./(squeeze(LGN_V1mri(whichSub,3,1,:,3))'-squeeze(LGN_V1mri(whichSub,3,1,:,1))');

% Define the eccentricity locations of the data. We use the log-mid point
% within each of the bins for the cortical
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
eccDegVals = eccDegBinEdges(4:2:14);
v1Eccentricity = []; v1FreqX = []; 
v1RedGreenY = []; v1RedGreenW = [];
v1BlueYellowY = []; v1BlueYellowW = [];
v1LumY = []; v1LumW = []; 
nEcc = 6;
for ee = 1:nEcc
    v1Eccentricity = [v1Eccentricity repmat(eccDegVals(ee),1,6)];
    v1FreqX = [v1FreqX studiedFreqs];
    v1RedGreenY = [v1RedGreenY squeeze(V1ecc_mri(whichSub,1,ee,:,2))'];
    v1RedGreenW = [v1RedGreenW 1./(squeeze(V1ecc_mri(whichSub,1,ee,:,3))'-squeeze(V1ecc_mri(whichSub,1,ee,:,1))')];
    v1BlueYellowY = [v1BlueYellowY squeeze(V1ecc_mri(whichSub,2,ee,:,2))'];
    v1BlueYellowW = [v1BlueYellowW 1./(squeeze(V1ecc_mri(whichSub,2,ee,:,3))'-squeeze(V1ecc_mri(whichSub,2,ee,:,1))')];
    v1LumY = [v1LumY squeeze(V1ecc_mri(whichSub,3,ee,:,2))'];
    v1LumW = [v1LumW 1./(squeeze(V1ecc_mri(whichSub,3,ee,:,3))'-squeeze(V1ecc_mri(whichSub,3,ee,:,1))')];
end

% Define the objective, plausible bound, and p0, depending upon subject and
% stimulus
switch whichStim
    case 1 % L-M chromatic

        % Returns the TTF, and handles reshaping into a linear vector
        myV1ChromTTF = @(p) assembleV1ChromResponseAcrossEcc(p,v1FreqX,v1Eccentricity);

        % The weighted objective
        myObj = @(p) norm(v1RedGreenW.*(v1RedGreenY - myV1ChromTTF(p))) + ...
            norm(lgnRedGreenW.*(lgnRedGreenY - returnlgnChromTTFFit(p,lgnFreqX,v1Eccentricity)));
        plb = [ 0.01 10 0.5 20 repmat(0.2,1,nEcc) zeros(1,nEcc)];
        pub = [ 0.10 20 1.0 30 repmat(0.8,1,nEcc) ones(1,nEcc)];
        switch whichSub
            case 1
                p0 = [0.0041   17.2462    0.9075   30.9342    0.7069    0.3565    0.2228    0.1369    0.1368    0.1111    0.0053    0.0110    0.0317    0.0320    0.0383    0.1294];
            case 2
                p0 = [0.0120   15.8006    0.5299   26.1496    0.4394    0.3332    0.1586    0.0069    0.0058    0.0043    0.0106    0.0202    0.0572    0.0700    0.1925    0.7719];
        end

    case 2 % S chromatic

        % Returns the TTF, and handles reshaping into a linear vector
        myV1ChromTTF = @(p) assembleV1ChromResponseAcrossEcc(p,v1FreqX,v1Eccentricity);

        % The weighted objective
        myObj = @(p) norm(v1BlueYellowW.*(v1BlueYellowY - myV1ChromTTF(p))) + ...
            norm(lgnBlueYellowW.*(lgnBlueYellowY - returnlgnChromTTFFit(p,lgnFreqX,v1Eccentricity)));
        plb = [ 0.01 10 0.5 20 repmat(0.2,1,nEcc) zeros(1,nEcc)];
        pub = [ 0.10 20 1.0 30 repmat(0.8,1,nEcc) ones(1,nEcc)];
        switch whichSub
            case 1
                p0 = [0.0030   15.8839    0.3990   28.3153    0.9261    0.4802    0.3507    0.1678    0.1367    0.1222    0.0016    0.0080    0.0323    0.0460    0.0663    0.1226];
            case 2
                p0 = [0.0026   12.9277    0.3248   19.9888    0.3326    0.3316    0.2216    0.0000    0.0000    0.0000    0.0028    0.0142    0.0652    0.0896    0.2406    0.8036];
        end

    case 3 % LMS luminance

        % Returns the TTF, and handles reshaping into a linear vector
        myV1LumTTF = @(p) assembleV1LumResponseAcrossEcc(p,v1FreqX,v1Eccentricity);

        myObj = @(p) norm(v1LumW.*(v1LumY - myV1LumTTF(p))) + ...
            norm(lgnLumW.*(lgnLumY - returnlgnLumTTFFit(p,lgnFreqX,v1Eccentricity)));
        plb = [ 0.01 20 0.1 10 repmat(0.2,1,nEcc) zeros(1,nEcc)];
        pub = [ 0.10 30 0.3 20 repmat(0.8,1,nEcc) ones(1,nEcc)];
        switch whichSub
            case 1
                p0 = [0.0460   25.8102    0.1533   11.3143    0.5645    0.5645    0.5166    0.4001    0.0822    0.0000    0.1506    0.2075    0.2673    0.2485    0.1404    0.2098];
            case 2
                p0 = [0.0445   19.9141    0.1284   15.1347    1.0000    0.7509    0.6211    0.4745    0.2943    0.1984    0.1363    0.2714    0.3897    0.3690    0.3764    0.7561];
        end
end

% hard bounds
lb = [ 0  10 0.01 05 zeros(1,nEcc) zeros(1,nEcc)];
ub = [ 1  50 2.00 40 ones(1,nEcc) ones(1,nEcc)];

% Non-linear constraint that surround index decrease with eccentricity
myNonbcon = @(p) nonbcon(p);

% Options - the objective function is deterministic
optionsBADS.UncertaintyHandling = 0;

% search
if searchFlag
    p = bads(myObj,p0,lb,ub,plb,pub,myNonbcon,optionsBADS);
else
    p = p0;
end

% plot
freqsForPlotting = logspace(0,2,50);

switch whichStim
    case 1
        figure
        v1RedGreenY = reshape(v1RedGreenY,6,1,nEcc);
        for ee=1:6
            subplot(2,4,ee+(ee>3))
            semilogx(studiedFreqs,squeeze(v1RedGreenY(:,1,ee)),'ok');
            hold on
            pBlock = [p(2:4) p(nFixed+ee) p(nFixed+nEcc+ee)];
            yFit = returnV1ChromEccTTFFit(pBlock,freqsForPlotting,eccDegVals(ee));
            semilogx(freqsForPlotting,yFit,'-r');
            refline(0,0);
            title([stimuli{whichStim} ', ' subjects{whichSub} ', ecc = ' num2str(eccDegVals(ee),2) '°']);
            ylim([-1 7]);
        end

        subplot(2,4,8)
        lgnChromTTFFit = returnlgnChromTTFFit(p,freqsForPlotting,v1Eccentricity);
        semilogx(lgnFreqX,lgnRedGreenY,'ok');
        hold on
        semilogx(freqsForPlotting,lgnChromTTFFit,'-r');
        title([stimuli{whichStim} ', ' subjects{whichSub} ', LGN']);
        refline(0,0);
        ylim([-0.5 4]);

    case 2
        figure
        v1BlueYellowY = reshape(v1BlueYellowY,6,1,nEcc);
        for ee=1:6
            subplot(2,4,ee+(ee>3))
            semilogx(studiedFreqs,squeeze(v1BlueYellowY(:,1,ee)),'ok');
            hold on
            pBlock = [p(2:4) p(nFixed+ee) p(nFixed+nEcc+ee)];
            yFit = returnV1ChromEccTTFFit(pBlock,freqsForPlotting,eccDegVals(ee));
            semilogx(freqsForPlotting,yFit,'-b');
            refline(0,0);
            ylim([-1 7]);
            title([stimuli{whichStim} ', ' subjects{whichSub} ', ecc = ' num2str(eccDegVals(ee),2) '°']);
        end

        subplot(2,4,8)
        lgnChromTTFFit = returnlgnChromTTFFit(p,freqsForPlotting,v1Eccentricity);
        semilogx(lgnFreqX,lgnBlueYellowY,'ok');
        hold on
        semilogx(freqsForPlotting,lgnChromTTFFit,'-b');
        refline(0,0);
        title([stimuli{whichStim} ', ' subjects{whichSub} ', LGN']);
        ylim([-0.5 4]);

    case 3
        figure
        v1LumY = reshape(v1LumY,6,1,nEcc);
        for ee=1:6
            subplot(2,4,ee+(ee>3))
            semilogx(studiedFreqs,squeeze(v1LumY(:,1,ee)),'ok');
            hold on
            pBlock = [p(2:4) p(nFixed+ee) p(nFixed+nEcc+ee)];
            yFit = returnV1LumEccTTFFit(pBlock,freqsForPlotting,eccDegVals(ee));
            semilogx(freqsForPlotting,yFit,'-k');
            refline(0,0);
            ylim([-1 7]);
            title([stimuli{whichStim} ', ' subjects{whichSub} ', ecc = ' num2str(eccDegVals(ee),2) '°']);
        end

        subplot(2,4,8)
        lgnLumTTFFit = returnlgnLumTTFFit(p,freqsForPlotting,v1Eccentricity);
        semilogx(lgnFreqX,lgnLumY,'ok');
        hold on
        semilogx(freqsForPlotting,lgnLumTTFFit,'-k');
        refline(0,0);
        title([stimuli{whichStim} ', ' subjects{whichSub} ', LGN']);
        ylim([-0.5 4]);

end

% Plot the surround suppression index vs. eccentricity
subplot(2,4,4)
plot(log10(eccDegVals),p(:,nFixed+1:nFixed+nEcc),'*k');
xlabel('Eccentricity [log deg]');
ylabel('Suppression index');
ylim([0 1]);

plotName = [stimuli{whichStim} '_' subjects{whichSub} '_ModelFit.pdf' ];
saveas(gcf,fullfile('~/Desktop',plotName));

%% LOCAL FUNCTIONS

%% nonbcon

% Enforce constraint of declining surround index with eccentricity
function c = nonbcon(p)
nEcc = 6; nFixed = 4;
surroundIndex = p(:,nFixed+1:nFixed+nEcc);
c = sum(diff(surroundIndex,1,2)>0,2);
end

function response = assembleV1ChromResponseAcrossEcc(p,v1FreqX,v1Eccentricity)
% Loop through eccentricities and obtain modeled responses
eccDegVals = unique(v1Eccentricity);
studiedFreqs = unique(v1FreqX);
% Info needed to unpack the param vector
nFixed = 4;
nEcc = length(eccDegVals);
% Build the response vector
response = [];
parfor ee=1:length(eccDegVals)
    pBlock = [p(2:4) p(nFixed+ee) p(nFixed+nEcc+ee)];
    response(ee,:) = returnV1ChromEccTTFFit(pBlock,studiedFreqs,eccDegVals(ee));
end
response = reshape(response',1,length(v1FreqX));
end


function response = assembleV1LumResponseAcrossEcc(p,v1FreqX,v1Eccentricity)
% Loop through eccentricities and obtain modeled responses
eccDegVals = unique(v1Eccentricity);
studiedFreqs = unique(v1FreqX);
% Info needed to unpack the param vector
nFixed = 4;
nEcc = length(eccDegVals);
% Build the response vector
response = [];
parfor ee=1:length(eccDegVals)
    pBlock = [p(2:4) p(nFixed+ee) p(nFixed+nEcc+ee)];
    response(ee,:) = returnV1LumEccTTFFit(pBlock,studiedFreqs,eccDegVals(ee));
end
response = reshape(response',1,length(v1FreqX));
end
