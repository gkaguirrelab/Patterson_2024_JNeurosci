function [pL,pM,pS,cnAll,snAll,cpAll,spAll]=RandomizeConeInputs_Gauss1(Em)

CSRadRatio=6; %Let's also say that a surround radius is ~7x greater than a center radius (Croner & Kaplan, Vis. Res. 1995)
CSStrengthRatio=.75; %Finally, let's say that surround gain is ~1/2 that of the center (Croner & Kaplan, Vis. Res. 1995)

ConesToCenter=ceil(0.29*(Em)^2+0.83*Em-0.28); %From Crook et al., New Vis. Neurosci., 2014, Fig 2B
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

ConesPerMM=ceil(19890*Em^(-0.6331)); %Curcio 1990

% %L:M:S ratio
LogLMDistribution=makedist('Normal','mu',0.47,'sigma',0.74); %Mu and sigma computed from a fit of Dacey et al. (2000) JOSA 17:589-596, Fig 5A
RandLM=exp(random(LogLMDistribution));
pL=RandLM/(RandLM+1);
pM=1-pL;
pS=0;

% pL=RetinaWeights/(RetinaWeights+1);
% pM=1-pL;
% pS=0;

LLIndex=0;
LUIndex=ceil(ConesPerMM*pL);
MLIndex=LUIndex+1;
MUIndex=MLIndex+ceil(ConesPerMM*pM);
SLIndex=MUIndex+1;
SUIndex=ConesPerMM;

if rem(ConesToCenter,2)==1
    CenterStrength=normpdf(1:ConesToCenter,ceil(ConesToCenter/2),0.34*ConesToCenter);
else
    CenterStrength=normpdf(1:ConesToCenter,(ConesToCenter+1)/2,0.34*ConesToCenter);
end

if rem(ConesToSurround,2)==1
    SurroundStrength=normpdf(1:ConesToSurround,ceil(ConesToSurround/2),0.34*ConesToSurround);
else
    SurroundStrength=normpdf(1:ConesToSurround,(ConesToSurround+1)/2,0.34*ConesToSurround);
end

SurroundAssign=randi(ConesPerMM,1,ConesToSurround);
CenterAssignIndex0=ceil(ConesToSurround/2)-floor(ConesToCenter/2);
CenterAssignIndexN=CenterAssignIndex0+ConesToCenter-1;
CenterAssign=SurroundAssign(CenterAssignIndex0:CenterAssignIndexN);

CenterAssignL=CenterAssign<=LUIndex;
CenterAssignM=(CenterAssign>MLIndex & CenterAssign<=MUIndex);
CenterAssignS=(CenterAssign>SLIndex & CenterAssign<=SUIndex);

cpL=(sum(CenterAssignL.*CenterStrength))/(sum(CenterStrength));
cpM=(sum(CenterAssignM.*CenterStrength))/(sum(CenterStrength));
cpS=(sum(CenterAssignS.*CenterStrength))/(sum(CenterStrength));

cnL=sum(CenterAssignL);
cnM=sum(CenterAssignM);
cnS=sum(CenterAssignS);

SurroundAssignL=SurroundAssign<=LUIndex;
SurroundAssignM=(SurroundAssign>MLIndex & SurroundAssign<=MUIndex);
SurroundAssignS=(SurroundAssign>SLIndex & SurroundAssign<=SUIndex);

spL=(sum(SurroundAssignL.*SurroundStrength))/(sum(SurroundStrength));
spM=(sum(SurroundAssignM.*SurroundStrength))/(sum(SurroundStrength));
spS=(sum(SurroundAssignS.*SurroundStrength))/(sum(SurroundStrength));

snL=sum(SurroundAssignL);
snM=sum(SurroundAssignM);
snS=sum(SurroundAssignS);

% disp(strcat(['Cones to center: ',num2str(cnL+cnM+cnS),' (',num2str(cnL),' L, ',num2str(cnM),' M, ' ,num2str(cnS),' S)']));
% disp(strcat(['Cones to surround: ',num2str(snL+snM+snS),' (',num2str(snL),' L, ',num2str(snM),' M, ' ,num2str(snS),' S)']));
% disp(strcat(['Center weighting: ',num2str(cnL/(cnL+cnM+cnS)),' L, ',num2str(cnM/(cnL+cnM+cnS)),' M, ' ,num2str(cnS/(cnL+cnM+cnS)),' S']));
% disp(strcat(['Surround weighting: ',num2str(snL/(snL+snM+snS)),' L, ',num2str(snM/(snL+snM+snS)),' M, ' ,num2str(snS/(snL+snM+snS)),' S']));

cnAll=[cnL cnM cnS];
snAll=[snL snM snS];

cpAll=[cpL cpM cpS];
spAll=[spL spM spS];