

figure
Em = logspace(log10(1.5/2/2),log10(18),10);
Ed=(Em*1000)/223; %Eccentricity in degrees (M. mulatta, Perry & Cowey 1985)
eccens = fieldnames(Data);
vals = [];
figure
for ee = 1:length(eccens)
    LplusM = []; LminusM = []; pT = [];
    cells = fieldnames(Data.(eccens{ee}));
    for ii = 1:length(cells)
        LplusM(ii) = Data.(eccens{ee}).(cells{ii}).ResponseFunctions.LMSumResponse.Amplitude;
        LminusM(ii) = Data.(eccens{ee}).(cells{ii}).ResponseFunctions.LMDiffResponse.Amplitude;
        pT(ii)=sum(Data.(eccens{ee}).(cells{ii}).CenterCones);
    end
    vals(ee)=mean(LminusM); % ./LplusM
    coneCenterCount(ee) = mean(pT);
    subplot(4,5,ee)
    histogram(LminusM./LplusM);
end
figure
semilogx(Ed,vals,'*')