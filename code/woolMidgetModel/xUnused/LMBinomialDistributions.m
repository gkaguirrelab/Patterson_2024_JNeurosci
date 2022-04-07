Sampling=[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20]; %mm
%Sampling=1;
CSRadRatio=7; %Let's also say that a surround radius is ~7x greater than a center radius (Croner & Kaplan, Vis. Res. 1995)
CSStrengthRatio=.75; %Finally, let's say that surround gain is ~1/2 that of the center (Croner & Kaplan, Vis. Res. 1995)

for f=1:length(Sampling)
    Em=Sampling(f);
    Ed=(Em*1000)/223; %Eccentricity in degrees (M. mulatta, Perry & Cowey 1985)

    C=ceil(0.29*(Em)^2+0.83*Em-0.28); %From Crook et al., New Vis. Neurosci., 2014, Fig 2B
    S=(CSRadRatio^2)*C;

    p=pdf('bino',0:C,C,0.5833)'*pdf('bino',0:S,S,0.5833);

    [m,i]=max(p(:));

    sL=floor(i/(length(p(:,1))))-1;
    sM=S-sL;
    
    cL=rem(i,(length(p(:,1))))-1;
    if cL<0
        cL=length(p(:,1))-1;
    end
    cM=C-cL;
    
    PlotCPurity(:,f)=cL/C;
    PlotSPurity(:,f)=sL/S;
end

plot(PlotCPurity);
hold on
plot(PlotSPurity);

% WORKSHOP

% C=8:38;
% S=(CSRadRatio^2).*C;
% index=1;
% 
% p=pdf('bino',0:C(index),C(index),0.5833)'*pdf('bino',0:S(index),S(index),0.5833);
% 
