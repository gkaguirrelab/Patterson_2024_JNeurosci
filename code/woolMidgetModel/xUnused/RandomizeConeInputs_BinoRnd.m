function [cpAll,spAll]=RandomizeConeInputs(Em)
CellType=1;
%Em=5; %Eccentricity in mm
Ed=(Em*1000)/223; %Eccentricity in degrees (M. mulatta, Perry & Cowey 1985)

%Center cone count
SigmasToRFRadius=2; %Let's say that 2 SDs make up a RF radius
CSRadRatio=7; %Let's also say that a surround radius is ~7x greater than a center radius (Croner & Kaplan, Vis. Res. 1995)
CSStrengthRatio=1/2; %Finally, let's say that surround gain is ~1/2 that of the center (Croner & Kaplan, Vis. Res. 1995)

ConesToCenter=ceil(0.29*(Em)^2+0.83*Em-0.28); %From Crook et al., New Vis. Neurosci., 2014, Fig 2B
ConesToSurround=(CSRadRatio^2)*ConesToCenter;

DFradius=((0.00056*(Em)^2+0.0014*Em+0.0013)*1000)/223; %From Crook et al. Excel data. This is in degrees.
CenterSigma=DFradius/SigmasToRFRadius; %Here we assign a sigma for our cell given a particular eccentricity.
SurroundSigma=CSRadRatio*CenterSigma; %Ditto.

%L:M:S ratio
pL=0.58;
pM=0.42;
pS=0;

cnL=0;
cnM=0;
cnS=0;
while (cnL+cnM+cnS)==0
    cnL=CellType*sum(binornd(1:1:1,pL,1,ConesToCenter));
    cnM=CellType*sum(binornd(1:1:1,pM,1,ConesToCenter));
    cnS=CellType*sum(binornd(1:1:1,pS,1,ConesToCenter));
end    

snL=-CellType*sum(binornd(1:1:1,pL,1,ConesToSurround));
snM=-CellType*sum(binornd(1:1:1,pM,1,ConesToSurround));
snS=-CellType*sum(binornd(1:1:1,pS,1,ConesToSurround));

% disp(strcat(['Cones to center: ',num2str(cnL+cnM+cnS),' (',num2str(cnL),' L, ',num2str(cnM),' M, ' ,num2str(cnS),' S)']));
% disp(strcat(['Cones to surround: ',num2str(snL+snM+snS),' (',num2str(snL),' L, ',num2str(snM),' M, ' ,num2str(snS),' S)']));

if cnL+cnM+cnS>1
    ConesToCenter=cnL+cnM+cnS;
    ConesToSurround=(CSRadRatio^2)*ConesToCenter;
    Em=(Ed/1000)*233;
    snL=-CellType*sum(binornd(1:1:1,pL,1,ConesToSurround));
    snM=-CellType*sum(binornd(1:1:1,pM,1,ConesToSurround));
    snS=-CellType*sum(binornd(1:1:1,pS,1,ConesToSurround));
    disp(strcat(['Cones to center: ',num2str(cnL+cnM+cnS),' (',num2str(cnL),' L, ',num2str(cnM),' M, ' ,num2str(cnS),' S)']));
    disp(strcat(['Cones to surround: ',num2str(snL+snM+snS),' (',num2str(snL),' L, ',num2str(snM),' M, ' ,num2str(snS),' S)']));
    disp(strcat(['NOTE: Eccentricity corrected to ',num2str(Em),' mm to reflect additional cones']));
else
    disp(strcat(['Cones to center: ',num2str(cnL+cnM+cnS),' (',num2str(cnL),' L, ',num2str(cnM),' M, ' ,num2str(cnS),' S)']));
    disp(strcat(['Cones to surround: ',num2str(snL+snM+snS),' (',num2str(snL),' L, ',num2str(snM),' M, ' ,num2str(snS),' S)']));
end

disp(strcat(['Center weighting: ',num2str(cnL/(cnL+cnM+cnS)),' L, ',num2str(cnM/(cnL+cnM+cnS)),' M, ' ,num2str(cnS/(cnL+cnM+cnS)),' S']));
disp(strcat(['Surround weighting: ',num2str(snL/(snL+snM+snS)),' L, ',num2str(snM/(snL+snM+snS)),' M, ' ,num2str(snS/(snL+snM+snS)),' S']));

cpL=cnL/(cnL+cnM+cnS);
cpM=cnM/(cnL+cnM+cnS);
cpS=cnS/(cnL+cnM+cnS);

spL=snL/(snL+snM+snS);
spM=snM/(snL+snM+snS);
spS=snS/(snL+snM+snS);

cpAll=[cpL cpM cpS];
spAll=[spL spM spS];