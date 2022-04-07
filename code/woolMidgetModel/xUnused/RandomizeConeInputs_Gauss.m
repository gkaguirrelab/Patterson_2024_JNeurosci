function [cnAll,snAll,cpAll,spAll]=RandomizeConeInputs_Gauss(Em)
CellType=1;
%Em=5; %Eccentricity in mm
Ed=(Em*1000)/223; %Eccentricity in degrees (M. mulatta, Perry & Cowey 1985)


%Center cone count
SigmasToRFRadius=2; %Let's say that 2 SDs make up a RF radius
CSRadRatio=6; %Let's also say that a surround radius is ~7x greater than a center radius (Croner & Kaplan, Vis. Res. 1995)
CSStrengthRatio=.75; %Finally, let's say that surround gain is ~1/2 that of the center (Croner & Kaplan, Vis. Res. 1995)

ConesToCenter=ceil(0.29*(Em)^2+0.83*Em-0.28); %From Crook et al., New Vis. Neurosci., 2014, Fig 2B
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

DFradius=((0.00056*(Em)^2+0.0014*Em+0.0013)*1000)/223; %From Crook et al. Excel data. This is in degrees.
CenterSigma=DFradius/SigmasToRFRadius; %Here we assign a sigma for our cell given a particular eccentricity.
SurroundSigma=CSRadRatio*CenterSigma; %Ditto.

%L:M:S ratio
ConesPerMM=ceil(19890*Em^(-0.6331)); %Curcio 1990
PopSubstrate=ceil(ConesPerMM/2);

pL=0.58;
pM=0.42;
pS=0;

LLIndex=0;
LUIndex=ceil(PopSubstrate*pL);
MLIndex=LUIndex+1;
MUIndex=MLIndex+ceil(PopSubstrate*pM);
SLIndex=MUIndex+1;
SUIndex=PopSubstrate;

CenterStrength=normpdf([1:ceil(ConesToCenter/2)],1,0.34*(ceil(ConesToCenter/2)));
SurroundStrength=normpdf([1:ceil(ConesToSurround/2)],1,0.34*(ceil(ConesToSurround/2)));

% CenterAssign1=randi(PopSubstrate,1,ceil(ConesToCenter/2)); %THIS ADDS EXTRA CONES
% CenterAssign2=randi(PopSubstrate,1,ceil(ConesToCenter/2)); %THIS ADDS EXTRA CONES
% SurroundAssign1=randi(PopSubstrate,1,ceil(ConesToSurround/2)); %THIS ADDS EXTRA CONES
% SurroundAssign2=randi(PopSubstrate,1,ceil(ConesToSurround/2)); %THIS ADDS EXTRA CONES

SurroundAssign1=randi(PopSubstrate,1,ceil(ConesToSurround/2)); %THIS ADDS EXTRA CONES
SurroundAssign2=randi(PopSubstrate,1,ceil(ConesToSurround/2)); %THIS ADDS EXTRA CONES
CenterAssign1=SurroundAssign1(1:ceil(ConesToCenter/2)); %THIS ADDS EXTRA CONES
CenterAssign2=SurroundAssign2(1:ceil(ConesToCenter/2)); %THIS ADDS EXTRA CONES

CenterAssignL1=CenterAssign1<=LUIndex;
CenterAssignM1=(CenterAssign1>MLIndex & CenterAssign1<=MUIndex);
CenterAssignS1=(CenterAssign1>SLIndex & CenterAssign1<=SUIndex);

CenterAssignL2=CenterAssign2<=LUIndex;
CenterAssignM2=(CenterAssign2>MLIndex & CenterAssign2<=MUIndex);
CenterAssignS2=(CenterAssign2>SLIndex & CenterAssign2<=SUIndex);

cpL=(sum(CenterAssignL1.*CenterStrength)+sum(CenterAssignL2.*CenterStrength))/(2*sum(CenterStrength));
cpM=(sum(CenterAssignM1.*CenterStrength)+sum(CenterAssignM2.*CenterStrength))/(2*sum(CenterStrength));
cpS=(sum(CenterAssignS1.*CenterStrength)+sum(CenterAssignS2.*CenterStrength))/(2*sum(CenterStrength));

cnL=sum(CenterAssignL1)+sum(CenterAssignL2);
cnM=sum(CenterAssignM1)+sum(CenterAssignM2);
cnS=sum(CenterAssignS1)+sum(CenterAssignS2);

SurroundAssignL1=SurroundAssign1<=LUIndex;
SurroundAssignM1=(SurroundAssign1>MLIndex & SurroundAssign1<=MUIndex);
SurroundAssignS1=(SurroundAssign1>SLIndex & SurroundAssign1<=SUIndex);

SurroundAssignL2=SurroundAssign2<=LUIndex;
SurroundAssignM2=(SurroundAssign2>MLIndex & SurroundAssign2<=MUIndex);
SurroundAssignS2=(SurroundAssign2>SLIndex & SurroundAssign2<=SUIndex);

spL=(sum(SurroundAssignL1.*SurroundStrength)+sum(SurroundAssignL2.*SurroundStrength))/(2*sum(SurroundStrength));
spM=(sum(SurroundAssignM1.*SurroundStrength)+sum(SurroundAssignM2.*SurroundStrength))/(2*sum(SurroundStrength));
spS=(sum(SurroundAssignS1.*SurroundStrength)+sum(SurroundAssignS2.*SurroundStrength))/(2*sum(SurroundStrength));

snL=sum(SurroundAssignL1)+sum(SurroundAssignL2);
snM=sum(SurroundAssignM1)+sum(SurroundAssignM2);
snS=sum(SurroundAssignS1)+sum(SurroundAssignS2);

% disp(strcat(['Cones to center: ',num2str(cnL+cnM+cnS),' (',num2str(cnL),' L, ',num2str(cnM),' M, ' ,num2str(cnS),' S)']));
% disp(strcat(['Cones to surround: ',num2str(snL+snM+snS),' (',num2str(snL),' L, ',num2str(snM),' M, ' ,num2str(snS),' S)']));
% disp(strcat(['Center weighting: ',num2str(cnL/(cnL+cnM+cnS)),' L, ',num2str(cnM/(cnL+cnM+cnS)),' M, ' ,num2str(cnS/(cnL+cnM+cnS)),' S']));
% disp(strcat(['Surround weighting: ',num2str(snL/(snL+snM+snS)),' L, ',num2str(snM/(snL+snM+snS)),' M, ' ,num2str(snS/(snL+snM+snS)),' S']));

cnAll=[cnL cnM cnS];
snAll=[snL snM snS];

cpAll=[cpL cpM cpS];
spAll=[spL spM spS];