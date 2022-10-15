% Load the Mt Sinai TTF results

searchFlag = true ;

whichStim = 1; % 1 = L-M; 3 = LMS
whichSub = 2;  % 1 = gka; 2 = asb

% This should be the V1 and LGN area bold fMRI signal mean, and 95% CI. The
% matrix is subject (GKA 1, ASB 2) x channel (L-M 1, S 2, LMS 3) x area
% (LGN 1, V1 2) x flicker freqency x bootstrap (1st value is 2.5%tile, 2nd
% value is 50%tile, 3rd value is 97.5%tile)
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','amplitudeResults','gka_asb_lgn_V1_BOLD.mat');
load(loadPath,'LGN_V1mri');

% This should be the V1 across eccentricity bold fMRI signal mean, and 95%
% CI. The matrix is subject (GKA 1, ASB 2) x channel (L-M 1, S 2, LMS 3) x
% eccentricity x flicker freqency x bootstrap (1st value is 2.5%tile, 2nd
% value is 50%tile, 3rd value is 97.5%tile)
loadPath = fullfile(fileparts(fileparts(fileparts(mfilename('fullpath')))),'data','amplitudeResults','gka_asb_V1_ecc_BOLD.mat');
load(loadPath,'V1ecc_mri');

% The frequencies studied
studiedFreqs = [2 4 8 16 32 64];

% Extract the relevant LGN data. Concatenate gka then asb
lgnFreqX = studiedFreqs;% studiedFreqs];
lgnLumY = squeeze(LGN_V1mri(whichSub,3,1,:,2))';
lgnLumW = 1./(squeeze(LGN_V1mri(whichSub,3,1,:,3))'-squeeze(LGN_V1mri(whichSub,3,1,:,1))');
lgnChromY = squeeze(LGN_V1mri(whichSub,1,1,:,2))';
lgnChromW = 1./(squeeze(LGN_V1mri(whichSub,1,1,:,3))'-squeeze(LGN_V1mri(whichSub,1,1,:,1))');

% Extract the relevant V1 data across eccentricity

% This returns the edges of each bin, along with the log-positioned
% mid-point within each bin
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
eccDegVals = eccDegBinEdges(4:2:14);
v1Eccentricity = []; v1FreqX = []; v1LumY = []; v1ChromY = []; v1LumW = []; v1ChromW = [];
nEcc = 6;
for ee = 1:nEcc
    v1Eccentricity = [v1Eccentricity repmat(eccDegVals(ee),1,6)];
    v1FreqX = [v1FreqX studiedFreqs];
    v1LumY = [v1LumY squeeze(V1ecc_mri(whichSub,3,ee,:,2))'];
    v1LumW = [v1LumW 1./(squeeze(V1ecc_mri(whichSub,3,ee,:,3))'-squeeze(V1ecc_mri(whichSub,3,ee,:,1))')];
    v1ChromY = [v1ChromY squeeze(V1ecc_mri(whichSub,1,ee,:,2))'];
    v1ChromW = [v1ChromW 1./(squeeze(V1ecc_mri(whichSub,1,ee,:,3))'-squeeze(V1ecc_mri(whichSub,1,ee,:,1))')];
end

% Define the objective and p0, depending upon subject and stimulus
switch whichStim
    case 1 % L-M chromatic
        myObj = @(p) norm(v1ChromW.*(v1ChromY - returnV1ChromEccTTFFit(p,v1FreqX,v1Eccentricity))) + ...
            norm(lgnChromW.*(lgnChromY - returnlgnChromTTFFit(p,lgnFreqX,v1Eccentricity)));
        switch whichSub
            case 1
                p0 = [0.0038   17.3918    0.9091   31.2138    0.7094    0.3494    0.2154    0.1251    0.1201    0.1078    0.0052    0.0111    0.0298    0.0318    0.0376    0.0808];
            case 2
                p0 = [0.0109   15.3284    0.5645   26.8723    0.5004    0.3325    0.1850    0.0048    0.0043    0.0000    0.0098    0.0192    0.0577    0.0696    0.1809    0.7404];
        end
    case 3 % LMS luminance
        myObj = @(p) norm(v1LumW.*(v1LumY - returnV1LumEccTTFFit(p,v1FreqX,v1Eccentricity))) + ...
            norm(lgnLumW.*(lgnLumY - returnlgnLumTTFFit(p,lgnFreqX,v1Eccentricity)));
        switch whichSub
            case 1
                p0 = [0.0462   25.8171    0.1522   10.9342    0.5731    0.5731    0.5309    0.4144    0.0891    0.0000    0.1586    0.2156    0.2784    0.2571    0.1427    0.2091];
            case 2
                p0 = [0.0453   18.7557    0.1920   14.6532    1.0000    0.7155    0.6070    0.4503    0.2440    0.1536    0.1010    0.2044    0.2876    0.2751    0.2830    0.5669];
        end
end

% bounds
lb = [ 0  05 0.01 05 zeros(1,nEcc) zeros(1,nEcc)];
ub = [ 1  50 2.00 40 ones(1,nEcc) ones(1,nEcc)];

plb = [ 0.01 10 0.1 10 repmat(0.2,1,nEcc) zeros(1,nEcc)];
pub = [ 0.10 30 1.5 30 repmat(0.8,1,nEcc) ones(1,nEcc)];

% Non-linear constraint that surround index decrease with eccentricity
myNonbcon = @(p) nonbcon(p);

% Options - the objective function is deterministic
optionsBADS.UncertaintyHandling = 0;

% search
if searchFlag
    p = bads(myObj,p0,lb,ub,plb,pub,myNonbcon,optionsBADS);
else
    p=p0;
end

% plot
switch whichStim
    case 1
        figure
        v1ChromTTFFit = returnV1ChromEccTTFFit(p,v1FreqX,v1Eccentricity);
        v1ChromTTFFit = reshape(v1ChromTTFFit,6,1,nEcc);
        v1ChromY = reshape(v1ChromY,6,1,nEcc);
        for ee=1:6
            subplot(3,2,ee)
            semilogx(studiedFreqs,squeeze(v1ChromY(:,1,ee)),'.k');
            hold on
            semilogx(studiedFreqs,squeeze(v1ChromTTFFit(:,1,ee)),'-r');
        end

        figure
        lgnChromTTFFit = returnlgnChromTTFFit(p,lgnFreqX,v1Eccentricity);
        semilogx(lgnFreqX,lgnChromY,'.k');
        hold on
        semilogx(lgnFreqX,lgnChromTTFFit,'-r');

    case 3
        figure
        v1LumTTFFit = returnV1LumEccTTFFit(p,v1FreqX,v1Eccentricity);
        v1LumTTFFit = reshape(v1LumTTFFit,6,1,nEcc);
        v1LumY = reshape(v1LumY,6,1,nEcc);
        for ee=1:6
            subplot(3,2,ee)
            semilogx(studiedFreqs,squeeze(v1LumY(:,1,ee)),'*k');
            hold on
            semilogx(studiedFreqs,squeeze(v1LumTTFFit(:,1,ee)),'-r');
            title(num2str(eccDegVals(ee),2));
        end

        figure
        lgnLumTTFFit = returnlgnLumTTFFit(p,lgnFreqX,v1Eccentricity);
        semilogx(lgnFreqX,lgnLumY,'.k');
        hold on
        semilogx(lgnFreqX,lgnLumTTFFit,'-r');
end


foo=1;



%% LOCAL FUNCTIONS

%% nonbcon

% Enforce constraint of declining surround index with eccentricity
function c = nonbcon(p)
nEcc = 6; nFixed = 4;
surroundIndex = p(:,nFixed+1:nFixed+nEcc);
c = sum(diff(surroundIndex,1,2)>0,2);
end
