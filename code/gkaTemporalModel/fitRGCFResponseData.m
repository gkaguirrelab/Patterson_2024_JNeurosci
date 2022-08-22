% fitRGCFResponseData
%



%% Load the flicker response data
[midgetData, parasolData] = rgcFlickerResponseData();


%% Set the p0 values
g = 5; % Overall gain
k = 0.75; % relative strength of the "lead compensators" (feedback stages)
cfLowPass = 17; % Corner frequency of the "cone" low-pass stage
cfInhibit = 28; % Corner frequency of the inhibitory stage
cf2ndStage = 60; % Corner frequency of 2nd order filter
Q = 1.5; % The "quality" parameter of the 2nd order filter
surroundWeight = 0.9; % Weight of the surround relative to center
surroundDelay = 3; % Delay (in msecs) of the surround relative to center
LMRatio = 1.5;


%% Define the bounds
p0Block = [g, k, cfLowPass, cfInhibit, cf2ndStage, Q, surroundWeight, surroundDelay];
lbBlock = [0, 0.6, 5, 5, 20, 1, 0.5, 1];
plbBlock = [3, 0.65, 10, 10, 50, 0.75, 0.75, 2];
pubBlock = [6, 0.85, 20, 25, 90, 1.25, 1, 15];
ubBlock = [10, 0.9, 25, 30, 100, 1.5, 1, 20];
shrinkParams = [false true false true true true true false];

p0 = [p0Block p0Block p0Block LMRatio];
lb = [lbBlock lbBlock lbBlock 0.1];
plb = [plbBlock plbBlock plbBlock 0.33];
pub = [pubBlock pubBlock pubBlock 3];
ub = [ubBlock ubBlock ubBlock 10];


%% Define the objective
myObj = @(p) rgcFitObjective(p,midgetData,shrinkParams);


%% Search
p = bads(myObj,p0,lb,ub,plb,pub);


%% Plot the results
eccFields = {'e0','e20','e30'};
eccDegs = [5, 25, mean([30 47])];

LMRatio = p(end);
p = reshape(p(1:24),[8,3]);

for ee = 1:3

    for cc = 1:1000
        [tmpCenterWeight(cc),tmpSurroundWeight(cc)] = ...
            returnChromaticWeights(eccDegs(ee),LMRatio);
    end
    chromaticCenterWeight = mean(abs(tmpCenterWeight));
    chromaticSurroundWeight = mean(abs(tmpSurroundWeight));

    pBlock = p(:,ee);

    g = pBlock(1); k = pBlock(2);
    cfLowPass = pBlock(3); cfInhibit = pBlock(4); 
    cf2ndStage = pBlock(5); Q = pBlock(6);
    surroundWeight = pBlock(7); surroundDelay = pBlock(8);


    [rfMidgetChrom, rfMidgetLum] = assembleMidgetRFs(...
        g, k, cfLowPass, cfInhibit, cf2ndStage, Q, ...
        surroundWeight, surroundDelay, ...
        chromaticCenterWeight, chromaticSurroundWeight);

    figHandle = figure();
    plotRF(rfMidgetLum,figHandle,'-k');
    plotRF(rfMidgetChrom,figHandle,'-r');
    subplot(3,1,1);
    loglog(midgetData.(eccFields{ee}).chromatic.f,midgetData.(eccFields{ee}).chromatic.g,'*r');
    loglog(midgetData.(eccFields{ee}).luminance.f,midgetData.(eccFields{ee}).luminance.g,'*k');
    subplot(3,1,2);
    semilogx(midgetData.(eccFields{ee}).chromatic.f,midgetData.(eccFields{ee}).chromatic.p,'*r');
    semilogx(midgetData.(eccFields{ee}).luminance.f,midgetData.(eccFields{ee}).luminance.p,'*k');

end




% Local functions

function fVal = rgcFitObjective(p,midgetData,shrinkParams)

eccFields = {'e0','e20','e30'};
eccDegs = [5, 25, mean([30 47])];

LMRatio = p(end);
p = reshape(p(1:24),[8,3]);
fVal = [];

for ee = 1:3

    for cc = 1:1000
        [tmpCenterWeight(cc),tmpSurroundWeight(cc)] = ...
                returnChromaticWeights(eccDegs(ee),LMRatio);
    end
    chromaticCenterWeight = mean(abs(tmpCenterWeight));
    chromaticSurroundWeight = mean(abs(tmpSurroundWeight));

    pBlock = p(:,ee);
    g = pBlock(1); k = pBlock(2);
    cfLowPass = pBlock(3); cfInhibit = pBlock(4); 
    cf2ndStage = pBlock(5); Q = pBlock(6);
    surroundWeight = pBlock(7); surroundDelay = pBlock(8);

    [rfMidgetChrom, rfMidgetLum] = assembleMidgetRFs(...
        g, k, cfLowPass, cfInhibit, cf2ndStage, Q, ...
        surroundWeight, surroundDelay, ...
        chromaticCenterWeight, chromaticSurroundWeight);

    fVal = [ fVal, ...
        norm(midgetData.(eccFields{ee}).chromatic.g - abs(eval(subs(rfMidgetChrom,midgetData.(eccFields{ee}).chromatic.f)))), ...
        norm(midgetData.(eccFields{ee}).luminance.g - abs(eval(subs(rfMidgetLum,midgetData.(eccFields{ee}).luminance.f)))) ...
        ];

end

% Take the L2 norm of all of the model fits
fVal = norm(fVal);

% Add a regularization that will attempt to get the parameter sets to match
% across eccentricity
shrink = norm(std(p(shrinkParams,:),[],2)./mean(p(shrinkParams,:),2));
fprintf('fVal: %2.2f, shrink: %2.2f \n',fVal,shrink);
fVal = fVal + shrink;

if isnan(fVal)
    foo=1;
end

end
