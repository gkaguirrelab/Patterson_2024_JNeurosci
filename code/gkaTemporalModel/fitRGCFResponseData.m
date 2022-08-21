% fitRGCFResponseData
%

g = 5; % Overall gain
k = 0.85; % relative strength of the "lead compensators" (feedback stages)
cfLowPass = 12; % Corner frequency of the "cone" low-pass stage
cfInhibit = 28; % Corner frequency of the inhibitory stage
cf2ndStage = 60; % Corner frequency of 2nd order filter
Q = 1.5; % The "quality" parameter of the 2nd order filter
surroundWeight = 0.9; % Weight of the surround relative to center
surroundDelay = 3; % Delay (in msecs) of the surround relative to center

% Create the set of chromatic weights by eccentricity
%{
eccDegs = [5, 25, mean([30 47])];
for ee = 1:3
    for rr = 1:100
        LogLMDistribution=makedist('Normal','mu',0.47,'sigma',0.74); %Mu and sigma computed from a fit of Dacey et al. (2000) JOSA 17:589-596, Fig 5A
        LMratio=exp(random(LogLMDistribution));
        for cc = 1:1000
            [tmpCenterWeight(rr,cc),tmpSurroundWeight(rr,cc)] = ...
                returnChromaticWeights(eccDegs(ee),LMratio);
        end
    end
    chromaticCenterWeight(ee) = mean(abs(tmpCenterWeight(:)));
    chromaticSurroundWeight(ee) = mean(abs(tmpSurroundWeight(:)));
end
%}
chromaticCenterWeight = [0.9974    0.3943    0.3354];
chromaticSurroundWeight = [0.3259    0.3292    0.3037];


%% Cell
[midgetData, parasolData] = rgcFlickerResponseData();


p0Block = [g, k, cfLowPass, cfInhibit, cf2ndStage, Q, surroundWeight, surroundDelay];
lbBlock = [0, 0.75, 10, 10, 50, 1, 0.5, 2];
ubBlock = [inf, 0.85, 20, 30, 70, 2, 1.5, 15];
shrinkParams = [false true false true true true true false];


p0 = [p0Block p0Block p0Block];
lb = [lbBlock lbBlock lbBlock];
ub = [ubBlock ubBlock ubBlock];

myObj = @(p) rgcFitObjective(p,chromaticCenterWeight,chromaticSurroundWeight,midgetData,shrinkParams);

p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],[]);

% Plot the results
eccFields = {'e0','e20','e30'};
p = reshape(p,[8,3]);

for ee = 1:3

    pBlock = p(:,ee);

    g = pBlock(1); k = pBlock(2);
    cfLowPass = pBlock(3); cfInhibit = pBlock(4); cf2ndStage = p(5); Q = p(6);
    surroundWeight = p(7); surroundDelay = p(8);

    [rfMidgetChrom, rfMidgetLum] = assembleMidgetRFs(...
        g, k, cfLowPass, cfInhibit, cf2ndStage, Q, ...
        surroundWeight, surroundDelay, ...
        chromaticCenterWeight(ee), chromaticSurroundWeight(ee));

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

function fVal = rgcFitObjective(p,chromaticCenterWeight,chromaticSurroundWeight,midgetData,shrinkParams)

fVal = [];
eccFields = {'e0','e20','e30'};

p = reshape(p,[8,3]);

for ee = 1:3

    pBlock = p(:,ee);
    g = pBlock(1); k = pBlock(2);
    cfLowPass = pBlock(3); cfInhibit = pBlock(4); 
    cf2ndStage = pBlock(5); Q = pBlock(6);
    surroundWeight = pBlock(7); surroundDelay = pBlock(8);

    [rfMidgetChrom, rfMidgetLum] = assembleMidgetRFs(...
        g, k, cfLowPass, cfInhibit, cf2ndStage, Q, ...
        surroundWeight, surroundDelay, ...
        chromaticCenterWeight(ee), chromaticSurroundWeight(ee));

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
