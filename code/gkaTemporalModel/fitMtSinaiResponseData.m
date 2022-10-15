% Load the Mt Sinai TTF results

searchFlag = true ;

whichStim = 3; % 1 = L-M; 3 = LMS
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
lgnChromY = squeeze(LGN_V1mri(whichSub,1,1,:,2))';

% Extract the relevant V1 data across eccentricity

% This returns the edges of each bin, along with the log-positioned
% mid-point within each bin
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
eccDegVals = eccDegBinEdges(4:2:14);
v1Eccentricity = []; v1FreqX = []; v1LumY = []; v1ChromY = [];
nEcc = 6;
for ee = 1:nEcc
    v1Eccentricity = [v1Eccentricity repmat(eccDegVals(ee),1,6)];
    v1FreqX = [v1FreqX studiedFreqs];
    v1LumY = [v1LumY squeeze(V1ecc_mri(whichSub,3,ee,:,2))'];
    v1ChromY = [v1ChromY squeeze(V1ecc_mri(whichSub,1,ee,:,2))'];
end

% Define the objective and p0, depending upon subject and stimulus
switch whichStim
    case 1 % L-M chromatic
        myObj = @(p) norm(v1ChromY - returnV1ChromEccTTFFit(p,v1FreqX,v1Eccentricity)) + ...
            norm(lgnChromY - returnlgnChromTTFFit(p,lgnFreqX,v1Eccentricity));
        switch whichSub
            case 1
                p0 = [0.0039   17.7473    0.9175   30.8000    0.6467    0.3519    0.2082    0.1368    0.1270    0.0562    0.0056    0.0110    0.0298    0.0311    0.0358    0.0658];
            case 2
                p0 = [0.0109   15.3284    0.5645   26.8723    0.5004    0.3325    0.1850    0.0048    0.0043    0.0000    0.0098    0.0192    0.0577    0.0696    0.1809    0.7404];
        end
    case 3 % LMS luminance
        myObj = @(p) norm(v1LumY - returnV1LumEccTTFFit(p,v1FreqX,v1Eccentricity)) + ...
            norm(lgnLumY - returnlgnLumTTFFit(p,lgnFreqX,v1Eccentricity));
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

function [v1LumAmplitude,v1LumPhase] = returnV1LumEccTTFFit(p,v1FreqX,v1Eccentricity)

nEcc = 6;
nFixed = 4;

% Unpack model parameters
lgnGain = p(1);
secondOrderFc = p(2);
secondOrderQ = p(3);
surroundDelay = p(4);
surroundIndex = p(nFixed+1:nFixed+nEcc);
nSubtractions = 2; % Two delayed surround stages: LGN and then V1
v1LumGain = p(nFixed+nEcc+1:end);

% Loop through eccentricities and obtain modeled responses
eccDegVals = unique(v1Eccentricity);
studiedFreqs = unique(v1FreqX);
v1LumAmplitude=[];v1LumPhase=[];

parfor ee=1:length(eccDegVals)
    [~,rfLumV1Ecc] = returnPostRetinalResponses(eccDegVals(ee),secondOrderFc,secondOrderQ,surroundIndex(ee),surroundDelay,nSubtractions);
    gainVals = zeros(size(studiedFreqs)); angleVals = zeros(size(studiedFreqs));
    for ii=1:length(rfLumV1Ecc)
        ttfComplex = double(subs(rfLumV1Ecc{ii},studiedFreqs));
        gainVals = gainVals + abs(ttfComplex).*(1/length(rfLumV1Ecc));
        angleVals = angleVals + unwrap(angle(ttfComplex)).*(1/length(rfLumV1Ecc));
    end
    v1LumAmplitude(ee,:) = v1LumGain(ee).*gainVals;
    v1LumPhase(ee,:) = angleVals;
end

v1LumAmplitude = reshape(v1LumAmplitude',1,nEcc*length(studiedFreqs));
v1LumPhase = reshape(v1LumPhase',1,nEcc*length(studiedFreqs));

end


function [v1ChromAmplitude,v1ChromPhase] = returnV1ChromEccTTFFit(p,v1FreqX,v1Eccentricity)

nEcc = 6;
nFixed = 4;

% Unpack model parameters
lgnGain = p(1);
secondOrderFc = p(2);
secondOrderQ = p(3);
surroundDelay = p(4);
surroundIndex = p(nFixed+1:nFixed+nEcc);
nSubtractions = 2; % Two delayed surround stages: LGN and then V1
v1LumGain = p(nFixed+nEcc+1:end);
v1ChromGain = v1LumGain;

% Loop through eccentricities and obtain modeled responses
eccDegVals = unique(v1Eccentricity);
studiedFreqs = unique(v1FreqX);
v1ChromAmplitude=[];v1ChromPhase=[];

parfor ee=1:length(eccDegVals)
    [rfChromV1Ecc] = returnPostRetinalResponses(eccDegVals(ee),secondOrderFc,secondOrderQ,surroundIndex(ee),surroundDelay,nSubtractions);
    ttfComplex = double(subs(rfChromV1Ecc,studiedFreqs));
    v1ChromAmplitude(ee,:) = v1ChromGain(ee).*abs(ttfComplex);
    v1ChromPhase(ee,:) = unwrap(angle(ttfComplex));
end

v1ChromAmplitude = reshape(v1ChromAmplitude',1,nEcc*length(studiedFreqs));
v1ChromPhase = reshape(v1ChromPhase',1,nEcc*length(studiedFreqs));

end



function [lgnLumAmplitude,lgnLumPhase] = returnlgnLumTTFFit(p,lgnFreqX,v1Eccentricity)

nEcc = 6;
nFixed = 4;

% Unpack model parameters
lgnGain = p(1);
secondOrderFc = p(2);
secondOrderQ = p(3);
surroundDelay = p(4);
surroundIndex = p(nFixed+1:nFixed+nEcc);
nSubtractions = 1; % One delayed surround stage at the LGN

% Retinal eccentricities to be averaged. Note that these are linearly
% spaced, as we wish to give equal weight to each annular area of the
% retina in contributing its population of RGCs to the average response.
% The effective weight of an annular area may be low due to a small number
% of RGCs, but each area is treated equivalently.
eccDegVals = 1:20:61;

% Interpolate the surround index values for the modeled eccentricities
surroundIndex = interp1(log10(unique(v1Eccentricity)),surroundIndex,log10(eccDegVals),'linear','extrap');

% Loop through eccentricities and obtain modeled responses
rfLumLGN = {};
parfor ee=1:length(eccDegVals)
    [~,rfLumLGN{ee}] = returnPostRetinalResponses(eccDegVals(ee),secondOrderFc,secondOrderQ,surroundIndex(ee),surroundDelay,nSubtractions);
end
rfLumLGN = reshape(vertcat(rfLumLGN{:}),length(eccDegVals)*2,1);

% Obtain the average luminance TTF across these eccentricities
lgnLumAmplitude = zeros(size(lgnFreqX)); lgnLumPhase = zeros(size(lgnFreqX));
for ii=1:length(rfLumLGN)
    ttfComplex = double(subs(rfLumLGN{ii},lgnFreqX));
    lgnLumAmplitude = lgnLumAmplitude + abs(ttfComplex).*(1/length(rfLumLGN));
    lgnLumPhase = lgnLumPhase + unwrap(angle(ttfComplex)).*(1/length(rfLumLGN));
end

lgnLumAmplitude = lgnLumAmplitude.*lgnGain;

end


function [lgnChromAmplitude,lgnChromPhase] = returnlgnChromTTFFit(p,lgnFreqX,v1Eccentricity)

nEcc = 6;
nFixed = 4;

% Unpack model parameters
lgnGain = p(1);
secondOrderFc = p(2);
secondOrderQ = p(3);
surroundDelay = p(4);
surroundIndex = p(nFixed+1:nFixed+nEcc);
nSubtractions = 1; % One delayed surround stage at the LGN

% Retinal eccentricities to be averaged. Note that these are linearly
% spaced, as we wish to give equal weight to each annular area of the
% retina in contributing its population of RGCs to the average response.
% The effective weight of an annular area may be low due to a small number
% of RGCs, but each area is treated equivalently.
eccDegVals = 1:20:61;

% Interpolate the surround index values for the modeled eccentricities
surroundIndex = interp1(log10(unique(v1Eccentricity)),surroundIndex,log10(eccDegVals),'linear','extrap');

% Loop through eccentricities and obtain modeled responses
rfChromLGN = {};
parfor ee=1:length(eccDegVals)
    [rfChromLGN{ee}] = returnPostRetinalResponses(eccDegVals(ee),secondOrderFc,secondOrderQ,surroundIndex(ee),surroundDelay,nSubtractions);
end

% Obtain the average chromatic TTF across these eccentricities
lgnChromAmplitude = zeros(size(lgnFreqX)); lgnChromPhase = zeros(size(lgnFreqX));
for ii=1:length(rfChromLGN)
    ttfComplex = double(subs(rfChromLGN{ii},lgnFreqX));
    lgnChromAmplitude = lgnChromAmplitude + abs(ttfComplex).*(1/length(rfChromLGN));
    lgnChromPhase = lgnChromPhase + unwrap(angle(ttfComplex)).*(1/length(rfChromLGN));
end

lgnChromAmplitude = lgnChromAmplitude.*lgnGain;

end