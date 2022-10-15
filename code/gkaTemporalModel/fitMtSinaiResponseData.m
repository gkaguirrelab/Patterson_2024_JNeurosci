% Load the Mt Sinai TTF results

whichSub = 2;

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
lgnFreqX = [studiedFreqs];% studiedFreqs];
lgnLumY = [squeeze(LGN_V1mri(whichSub,3,1,:,2))'];

% Extract the relevant V1 data across eccentricity

% This returns the edges of each bin, along with the log-positioned
% mid-point within each bin
eccDegBinEdges = logspace(log10(0.7031),log10(90),15);
eccDegVals = eccDegBinEdges(4:2:14);
v1Eccentricity = []; v1FreqX = []; v1LumY = []; v1ChromY = [];
nEcc = 6;
for ee = 1:nEcc
%     v1Eccentricity = [v1Eccentricity repmat(eccDegVals(ee),1,6) repmat(eccDegVals(ee),1,6)];
%     v1FreqX = [v1FreqX studiedFreqs studiedFreqs];
%     v1LumY = [v1LumY squeeze(V1ecc_mri(1,3,ee,:,2))' squeeze(V1ecc_mri(2,3,ee,:,2))'];
%     v1ChromY = [v1ChromY squeeze(V1ecc_mri(1,1,ee,:,2))' squeeze(V1ecc_mri(2,1,ee,:,2))'];
    v1Eccentricity = [v1Eccentricity repmat(eccDegVals(ee),1,6)];
    v1FreqX = [v1FreqX studiedFreqs];
    v1LumY = [v1LumY squeeze(V1ecc_mri(whichSub,3,ee,:,2))'];
    v1ChromY = [v1ChromY squeeze(V1ecc_mri(whichSub,1,ee,:,2))'];
end

% Define the objective
myObj = @(p) norm(v1LumY - returnV1LumEccTTFFit(p,v1FreqX,v1Eccentricity)) + ...
             norm(lgnLumY - returnlgnLumTTFFit(p,lgnFreqX,v1Eccentricity));

%myObj = @(p) norm(v1ChromY - returnV1ChromEccTTFFit(p,v1FreqX,v1Eccentricity));

options = optimoptions(@fmincon,'Display','iter');

% p0 and bounds
%p0 = [0.1 09 0.29 10 repmat(0.5,1,nEcc) repmat(0.3,1,nEcc)];
p0 = [0.0462   25.7673    0.1516   10.9039    0.5773    0.5777    0.5304    0.4170    0.0914    0.0019    0.1602    0.2178    0.2807    0.2594    0.1439    0.2107];
lb = [ 00 05 0.01 05 repmat(0,1,nEcc) zeros(1,nEcc)];
ub = [inf 30 2.00 20 repmat(1,1,nEcc) inf(1,nEcc)];

% search
p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

% plot
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


foo=1;



%% LOCAL FUNCTIONS

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

% Fit the surround index values with a falling Weibull CDF as a function of
% log eccentricity.
weibullFunc = @(x,g,lambda,k) g.*exp(-(x./lambda).^k);
myObj = @(pWeib) norm(surroundIndex - weibullFunc(log10(unique(v1Eccentricity)),pWeib(1),pWeib(2),pWeib(3)));
options = optimoptions(@fmincon,'Display','off');
pWeib = fmincon(myObj,[1 1 1],[],[],[],[],[],[],[],options);

% Get the interpolated v1SurroundIndex values at the eccentricity values
surroundIndex = weibullFunc(log10(eccDegVals),pWeib(1),pWeib(2),pWeib(3));

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