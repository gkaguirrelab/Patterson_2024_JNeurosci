function [chromaticCenterWeight,chromaticSurroundWeight] = returnMidgetChromaticWeights(eccentricityDeg,LMratio)

%{
    LogLMDistribution=makedist('Normal','mu',0.47,'sigma',0.74); %Mu and sigma computed from a fit of Dacey et al. (2000) JOSA 17:589-596, Fig 5A
    LMratio=exp(random(LogLMDistribution));
    figure
    plot(0,0,'.')
    hold on
    for ii=0.1:3:51.1; for xx=1:5000; [chromaticCenterWeight(xx),chromaticSurroundWeight(xx)] = returnChromaticWeights(ii,LMratio); end; plot(ii,mean(abs(chromaticCenterWeight)),'xk'); plot(ii,mean(abs(chromaticSurroundWeight)),'or'); end; title(num2str(LMratio));
%}

% Convert from degrees of visual angle to mm in the macaque retina
Em = (eccentricityDeg.*223)./1000;

%Determine some parameters for centers and surrounds
CSRadRatio=6; %Let's say that a surround radius is ~6x greater than a center radius (Croner & Kaplan, Vis. Res. 1995)

% Simulate 1000 RGCs with random sampling from the cone mosaic at the
% specified eccentricity and with the specified LMRatio. For each RGC,
% obtain the chromatic weight on the center and surround. We then take the
% absolute value of these, as we assume that the behavior of an L-center
% RGC is the same as the behavior of an M-center RGC.
nSims = 1000;

for cc = 1:nSims


    %% GKA modification
    % How many cones comprise an RF center at any eccentricity?
    % The cones to center was under a "ceil" operation. The trouble with this
    % is that, instead of a smooth function across eccentricity, the output of
    % the routine takes discrete steps at certain eccentricites. To handle
    % this, we set the number of cones in the center to be a probabilistic
    % outcome of the fractional component of the value
    ConesToCenter=(0.29*(Em)^2+0.83*Em-0.28); %From Crook et al., New Vis. Neurosci., 2014, Fig 2B
    ConesToCenter=floor(ConesToCenter)+fix(rem(ConesToCenter,1)>rand);
    ConesToCenter=max([1 ConesToCenter]);

    %Multiply that value by the C:S area ratio (if Rs=CSRadRatio*Rc, then
    %As=(CSRadRatio^2)*Ac)
    ConesToSurround=(CSRadRatio^2)*ConesToCenter;

    %Determine a 1-mm^2 patch of cones for a particular retinal eccentricity
    ConesPerMM=ceil(19890*Em^(-0.6331)); %Curcio 1990

    % Assign a retinal L:M:S ratio from which any one cell comes
    % (L:M:S ratios vary across retinas, lognormally)
    pL=LMratio/(LMratio+1); %Likelihood of L
    pM=1-pL; %Likelihood of M

    %Assign the patch (a vector of length ConesPerMM) to L, M, and S
    %'zones' (for assigning discrete cones later)
    LUIndex=ceil(ConesPerMM*pL);
    MLIndex=LUIndex+1;
    MUIndex=MLIndex+ceil(ConesPerMM*pM);

    %Determine the dimensions of the NxM square matrix needed to place all
    %the cells in a 2D arrangement
    ConesToCenterDim=ceil(sqrt(ConesToCenter));
    ConesToSurroundDim=ceil(sqrt(ConesToSurround));

    %Compute a responsivity Gaussian for one dimension of the 2D grid

    %Choose a value for the sigma (cone contribution falloff)
    ApertureSigma=.34; %%%%Can't remember why we chose this value...

    %Symmetric Gaussian correction
    if rem(ConesToCenter,2)==1
        CenterStrength=normpdf(1:ConesToCenterDim,ceil(ConesToCenterDim/2),ApertureSigma*ConesToCenterDim);
    else
        CenterStrength=normpdf(1:ConesToCenterDim,(ConesToCenterDim+1)/2,ApertureSigma*ConesToCenterDim);
    end

    if rem(ConesToSurround,2)==1
        SurroundStrength=normpdf(1:ConesToSurroundDim,ceil(ConesToSurroundDim/2),ApertureSigma*ConesToSurroundDim);
    else
        SurroundStrength=normpdf(1:ConesToSurroundDim,(ConesToSurroundDim+1)/2,ApertureSigma*ConesToSurroundDim);
    end

    %Compute the outer product of 2 1D Gaussians to create a 2D Gaussian of the
    %receptive field "center" and "surround" responsivity (there are no cones placed in the
    %patch yet...)
    CenterSurface=CenterStrength'*CenterStrength;
    SurroundSurface=SurroundStrength'*SurroundStrength;

    CenterWeights=sort(reshape(CenterSurface,[1,ConesToCenterDim^2]),'descend');
    SurroundWeights=sort(reshape(SurroundSurface,[1,ConesToSurroundDim^2]),'descend');

    %Now randomly assign L, M, and S cones to the surround
    SurroundAssign=randi(ConesPerMM,1,ConesToSurround);

    %Select a subsection of surround cones to serve as the center
    %(the retina double dips)
    CenterAssignIndex0=1;
    CenterAssignIndexN=ConesToCenter;
    CenterAssign=SurroundAssign(CenterAssignIndex0:CenterAssignIndexN);

    %Determine which center values are Ls, Ms, and Ss...
    CenterAssignL=CenterAssign<=LUIndex;
    CenterAssignM=(CenterAssign>MLIndex & CenterAssign<=MUIndex);

    %Weight the cone contributions given the responsivity Gaussian computed
    %above
    cpL=(sum(CenterAssignL.*CenterWeights(1:ConesToCenter)))/(sum(CenterWeights(1:ConesToCenter)));
    cpM=(sum(CenterAssignM.*CenterWeights(1:ConesToCenter)))/(sum(CenterWeights(1:ConesToCenter)));

    %Determine which surround values are Ls, Ms, and Ss...
    SurroundAssignL=SurroundAssign<=LUIndex;
    SurroundAssignM=(SurroundAssign>MLIndex & SurroundAssign<=MUIndex);

    %Weight the cone contributions given the responsivity Gaussian computed
    %above
    spL=(sum(SurroundAssignL.*SurroundWeights(1:ConesToSurround)))/(sum(SurroundWeights(1:ConesToSurround)));
    spM=(sum(SurroundAssignM.*SurroundWeights(1:ConesToSurround)))/(sum(SurroundWeights(1:ConesToSurround)));


    chromaticCenterWeight(cc) = cpL - cpM;
    chromaticSurroundWeight(cc) = spL - spM;
end

chromaticCenterWeight = mean(abs(chromaticCenterWeight));
chromaticSurroundWeight = mean(abs(chromaticSurroundWeight));


end