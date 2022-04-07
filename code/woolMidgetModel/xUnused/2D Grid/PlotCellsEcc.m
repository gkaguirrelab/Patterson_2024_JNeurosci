%%
Cells=fieldnames(Data);

for i=1:length(Cells)
    ImportCell=Data.(matlab.lang.makeValidName(Cells{i}));

    Em=ImportCell.Eccentricity;
    CenterCones=ImportCell.CenterConeCount;
    ChromTag=ImportCell.ChromTag;
    TotalLData=ImportCell.ResponseFunctions.LIsoResponse;
    TotalMData=ImportCell.ResponseFunctions.MIsoResponse;
    LMSumData=ImportCell.ResponseFunctions.LMSumResponse;
    LMDiffData=ImportCell.ResponseFunctions.LMDiffResponse;
    
    FirstAmpLPhase=TotalLData.Theta0deg.Phase(1);
    FirstAmpMPhase=TotalMData.Theta0deg.Phase(1);
    MaxAmpL=find(TotalLData.Theta0deg.Amplitude==max(TotalLData.Theta0deg.Amplitude));
    MaxAmpLPhase=TotalLData.Theta0deg.Phase(MaxAmpL);
    MaxAmpM=find(TotalMData.Theta0deg.Amplitude==max(TotalMData.Theta0deg.Amplitude));
    MaxAmpMPhase=TotalMData.Theta0deg.Phase(MaxAmpM);
       
    if strcmp(ChromTag,'Achromatic')==1
        NotchTag='N/A';
    elseif strcmp(ChromTag,'L-dominated')==1
        PhaseIndex=abs(FirstAmpMPhase)-abs(MaxAmpMPhase);
        if PhaseIndex>100
            NotchTag='M notch';
        elseif PhaseIndex<100
            NotchTag='No notch';
        end
    elseif strcmp(ChromTag,'M-dominated')==1
        PhaseIndex=abs(FirstAmpLPhase)-abs(MaxAmpLPhase);
        if PhaseIndex>100
            NotchTag='L notch';
        elseif PhaseIndex<100
            NotchTag='No notch';
        end
    end

    IGOR{i,1}=Data.(matlab.lang.makeValidName(Cells{i})).Eccentricity;
    Eccentricities(i,1)=IGOR{i,1};
    IGOR{i,2}=Data.(matlab.lang.makeValidName(Cells{i})).CenterWeights(1);
    IGOR{i,3}=Data.(matlab.lang.makeValidName(Cells{i})).SurroundWeights(1);
    IGOR{i,4}=Data.(matlab.lang.makeValidName(Cells{i})).ChromTag;
    IGOR{i,5}=NotchTag;
    IGOR{i,6}=Data.(matlab.lang.makeValidName(Cells{i})).RetinaWeights(1);
    IGOR{i,8}=Data.(matlab.lang.makeValidName(Cells{i})).CSStrengthRatio(1);
    %IGOR{i,7}=LMDiffData.Theta0deg.Amplitude(1)/LMSumData.Theta0deg.Amplitude(1);
    IGOR{i,7}=LMDiffData.Theta0deg.Amplitude(1);
    %IGOR{i,9}=Data.(matlab.lang.makeValidName(Cells{i})).Polarity;
    LumGainIndex=find(LMSumData.Theta0deg.Amplitude==max(LMSumData.Theta0deg.Amplitude));
    %IGOR{i,8}=(LMSumData.Theta0deg.Amplitude(LumGainIndex));%/(LMDiffData.Theta0deg.Amplitude(LumGainIndex));
    %IGOR{i,8}=(LMSumData.Theta0deg.Amplitude(1));
    IGOR{i,9}=i;
    AllRetinaWeights(i,1)=IGOR{i,6};
    RGGains(i,1)=IGOR{i,7};
    LumGains(i,:)=IGOR{i,8};
    RGPhases(i,1)=IGOR{i,9};
    
    SWeights(i,1)=Data.(matlab.lang.makeValidName(Cells{i})).CenterWeights(3);
    CSStrengthRatio(i,1)=Data.(matlab.lang.makeValidName(Cells{i})).CSStrengthRatio(1);
    CenterLOverTotalL(i,1)=(Data.(matlab.lang.makeValidName(Cells{i})).CenterWeights(1))/(Data.(matlab.lang.makeValidName(Cells{i})).CenterWeights(1)+(Data.(matlab.lang.makeValidName(Cells{i})).CSStrengthRatio*Data.(matlab.lang.makeValidName(Cells{i})).SurroundWeights(1)));
    CenterMOverTotalM(i,1)=(Data.(matlab.lang.makeValidName(Cells{i})).CenterWeights(2))/(Data.(matlab.lang.makeValidName(Cells{i})).CenterWeights(2)+(Data.(matlab.lang.makeValidName(Cells{i})).CSStrengthRatio*Data.(matlab.lang.makeValidName(Cells{i})).SurroundWeights(2)));
end

SortedIGOR=sortrows(IGOR,1);
SortedIGOR_Retina=sortrows(IGOR,6);
SortedEccentricities=sort(Eccentricities);
SortedRetinaWeights=sort(AllRetinaWeights);

%% BIN ALL THE CELLS
bins=linspace(.25,10,10);
bins=[.25 3.5 6.75 10];
bins=[.25 2 4 6 8 10];
IDX=uint32(1:size(SortedEccentricities,1))';
for g=2:length(bins)
    BinLabel=strcat('Bin',num2str(bins(g)));
%     EccsUB=find(SortedEccentricities<bins(g));
%     EccsLB=find(SortedEccentricities>bins(g-1));
%     Eccs=EccsLB(1:EccsUB(end));
    Eccs=IDX(SortedEccentricities>=bins(g-1) & SortedEccentricities<bins(g));
    if isempty(Eccs)
        BinnedCells{1,1}=0;
        BinnedCells{1,2}=0;
        BinnedCells{1,3}=0;
        BinnedCells{1,4}=0;
        BinnedCells{1,5}=0;
        BinnedCells{1,6}=0;
        BinnedCells{1,7}=0;
        BinnedCells{1,8}=0;
        BinnedCells{1,9}=0;
    else
        for f=1:size(Eccs,1)
        BinnedCells{f,1}=SortedIGOR{Eccs(f),1};
        BinnedCells{f,2}=SortedIGOR{Eccs(f),2};
        BinnedCells{f,3}=SortedIGOR{Eccs(f),3};
        BinnedCells{f,4}=SortedIGOR{Eccs(f),4};
        BinnedCells{f,5}=SortedIGOR{Eccs(f),5};
        BinnedCells{f,6}=SortedIGOR{Eccs(f),6};
        BinnedCells{f,7}=SortedIGOR{Eccs(f),7};
        BinnedCells{f,8}=SortedIGOR{Eccs(f),8};
        BinnedCells{f,9}=SortedIGOR{Eccs(f),9};
        end
    end
    Histograms.(matlab.lang.makeValidName(BinLabel)).AllData=BinnedCells;
    clearvars BinnedCells
end

for f=2:length(bins)
    LDoms=0;
    MDoms=0;
    AChroms=0;
    NotchCells=0;
    ChromTag=[];
    Center=[];
    Surround=[];
    CSRatio=[];
    BinLabel=strcat('Bin',num2str(bins(f)));
    BinnedCells=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData;
    for g=1:size(BinnedCells,1)
        Center(g,:)=BinnedCells{g,2};
        Surround(g,:)=BinnedCells{g,3};
        CSRatio(g,:)=BinnedCells{g,8};
        RetinaWeights(g,:)=BinnedCells{g,6};
        RGGain(g,:)=BinnedCells{g,7};
        Polarity(g,:)=BinnedCells{g,9};
        %RGPhase(g,:)=BinnedCells{g,9};
        
        if strcmp(BinnedCells{g,4},'Achromatic')==1
            ChromTag(g,:)=0;
            AChroms=AChroms+1;
        elseif strcmp(BinnedCells{g,4},'L-dominated')==1
            ChromTag(g,:)=1;
            LDoms=LDoms+1;
            if strcmp(BinnedCells{g,5},'M notch')
                NotchCells=NotchCells+1;
            end
        elseif strcmp(BinnedCells{g,4},'M-dominated')==1
            ChromTag(g,:)=-1;
            MDoms=MDoms+1;
            if strcmp(BinnedCells{g,5},'L notch')
                NotchCells=NotchCells+1;
            end
        end
        
    end
    Histograms.(matlab.lang.makeValidName(BinLabel)).CenterWeights=Center;
    Histograms.(matlab.lang.makeValidName(BinLabel)).SurroundWeights=Surround;
    Histograms.(matlab.lang.makeValidName(BinLabel)).CSRatio=CSRatio;
    Histograms.(matlab.lang.makeValidName(BinLabel)).LDoms=LDoms;
    Histograms.(matlab.lang.makeValidName(BinLabel)).MDoms=MDoms;
    Histograms.(matlab.lang.makeValidName(BinLabel)).NotchCells=NotchCells;
    Histograms.(matlab.lang.makeValidName(BinLabel)).AChroms=AChroms;
    Histograms.(matlab.lang.makeValidName(BinLabel)).RetinaWeights=RetinaWeights;
    Histograms.(matlab.lang.makeValidName(BinLabel)).RGGain=RGGain;
    %Histograms.(matlab.lang.makeValidName(BinLabel)).RGPhase=RGPhase;
    Histograms.(matlab.lang.makeValidName(BinLabel)).Polarity=Polarity;
end

% for f=2:length(bins)
%     BinLabel=strcat('Bin',num2str(bins(f)));
%     semilogy(Histograms.(matlab.lang.makeValidName(BinLabel)).CSRatio,'o')
%     hold on;
% end

%% BIN CELLS FROM A CERTAIN SET OF RETINAS

% RetinaLB=find(SortedRetinaWeights>.4, 1 );
% RetinaUB=find(SortedRetinaWeights>.4, 1, 'last' );

RetinaLB=find(SortedRetinaWeights>.65 & SortedRetinaWeights<=.67, 1);
RetinaUB=find(SortedRetinaWeights>.65 & SortedRetinaWeights<=.67, 1,'last');

for x=1:RetinaUB-RetinaLB;
    SubIGOR{x,1}=SortedIGOR_Retina{x+RetinaLB-1,1};
    SubIGOR{x,2}=SortedIGOR_Retina{x+RetinaLB-1,2};
    SubIGOR{x,3}=SortedIGOR_Retina{x+RetinaLB-1,3};
    SubIGOR{x,4}=SortedIGOR_Retina{x+RetinaLB-1,4};
    SubIGOR{x,5}=SortedIGOR_Retina{x+RetinaLB-1,5};
    SubIGOR{x,6}=SortedIGOR_Retina{x+RetinaLB-1,6};
    SubIGOR{x,7}=SortedIGOR_Retina{x+RetinaLB-1,7};
    SubIGOR{x,8}=SortedIGOR_Retina{x+RetinaLB-1,8};
    SubIGOR{x,9}=SortedIGOR_Retina{x+RetinaLB-1,9};
    Eccentricities_Retina(x,1)=SortedIGOR_Retina{x+RetinaLB-1,1};
end
SortedSubIGOR=sortrows(SubIGOR,1);
SortedEccentricities_Retina=sort(Eccentricities_Retina);

bins=linspace(.25,10,40);
IDX=uint32(1:size(SortedEccentricities_Retina,1))';
for g=2:length(bins)
    BinLabel=strcat('Bin',num2str(bins(g)));
%     EccsUB=find(SortedEccentricities<bins(g));
%     EccsLB=find(SortedEccentricities>bins(g-1));
%     Eccs=EccsLB(1:EccsUB(end));
    Eccs=IDX(SortedEccentricities_Retina>=bins(g-1) & SortedEccentricities_Retina<bins(g));
    if isempty(Eccs)
        BinnedCells{1,1}=0;
        BinnedCells{1,2}=0;
        BinnedCells{1,3}=0;
        BinnedCells{1,4}=0;
        BinnedCells{1,5}=0;
        BinnedCells{1,6}=0;
        BinnedCells{1,7}=0;
        BinnedCells{1,8}=0;
        BinnedCells{1,9}=0;
    else
        for f=1:size(Eccs,1)
        BinnedCells{f,1}=SortedSubIGOR{Eccs(f),1};
        BinnedCells{f,2}=SortedSubIGOR{Eccs(f),2};
        BinnedCells{f,3}=SortedSubIGOR{Eccs(f),3};
        BinnedCells{f,4}=SortedSubIGOR{Eccs(f),4};
        BinnedCells{f,5}=SortedSubIGOR{Eccs(f),5};
        BinnedCells{f,6}=SortedSubIGOR{Eccs(f),6};
        BinnedCells{f,7}=SortedSubIGOR{Eccs(f),7};
        BinnedCells{f,8}=SortedSubIGOR{Eccs(f),8};
        BinnedCells{f,9}=SortedSubIGOR{Eccs(f),9};
        end
    end
    Histograms.(matlab.lang.makeValidName(BinLabel)).AllData=BinnedCells;
    clearvars BinnedCells
end

for f=2:length(bins)
    LDoms=0;
    MDoms=0;
    AChroms=0;
    NotchCells=0;
    ChromTag=[];
    Center=[];
    Surround=[];
    CSRatio=[];
    BinLabel=strcat('Bin',num2str(bins(f)));
    BinnedCells=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData;
    for g=1:size(BinnedCells,1)
        Center(g,:)=BinnedCells{g,2};
        Surround(g,:)=BinnedCells{g,3};
        CSRatio(g,:)=BinnedCells{g,2}/BinnedCells{g,3};
        RetinaWeights(g,:)=BinnedCells{g,6};
        RGGain(g,:)=BinnedCells{g,7};
        LumGain(g,:)=BinnedCells{g,8};
        
        if strcmp(BinnedCells{g,4},'Achromatic')==1
            ChromTag(g,:)=0;
            AChroms=AChroms+1;
        elseif strcmp(BinnedCells{g,4},'L-dominated')==1
            ChromTag(g,:)=1;
            LDoms=LDoms+1;
            if strcmp(BinnedCells{g,5},'M notch')
                NotchCells=NotchCells+1;
            end
        elseif strcmp(BinnedCells{g,4},'M-dominated')==1
            ChromTag(g,:)=-1;
            MDoms=MDoms+1;
            if strcmp(BinnedCells{g,5},'L notch')
                NotchCells=NotchCells+1;
            end
        end
    end
    Histograms.(matlab.lang.makeValidName(BinLabel)).CenterWeights=Center;
    Histograms.(matlab.lang.makeValidName(BinLabel)).SurroundWeights=Surround;
    Histograms.(matlab.lang.makeValidName(BinLabel)).CSRatio=CSRatio;
    Histograms.(matlab.lang.makeValidName(BinLabel)).LDoms=LDoms;
    Histograms.(matlab.lang.makeValidName(BinLabel)).MDoms=MDoms;
    Histograms.(matlab.lang.makeValidName(BinLabel)).NotchCells=NotchCells;
    Histograms.(matlab.lang.makeValidName(BinLabel)).AChroms=AChroms;
    Histograms.(matlab.lang.makeValidName(BinLabel)).RetinaWeights=RetinaWeights;
    Histograms.(matlab.lang.makeValidName(BinLabel)).RGGain=RGGain;
    Histograms.(matlab.lang.makeValidName(BinLabel)).LumGain=LumGain;
end

%% ECC-BINNED SUBPLOTS
figure;
hold on;
for e=2:length(bins)
    subplot(1,3,e-1);
    BinLabel=strcat('Bin',num2str(bins(e)));
    BinnedCells=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData;
    for c=1:length(BinnedCells)
        PlotNotch=Histograms.(matlab.lang.makeValidName(BinLabel)).CenterWeights(c);
        PlotSurroundWeight=Histograms.(matlab.lang.makeValidName(BinLabel)).SurroundWeights(c);
        PlotBin=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData{c,4};
        if strcmp(PlotBin,'L-dominated')==1
            PlotL=plot(PlotNotch,PlotSurroundWeight,'ro');
            hold on;
            set(PlotL,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',[0.64 0.08 0.18],...
                'MarkerSize',8,...
                'Marker','o')
        elseif strcmp(PlotBin,'M-dominated')==1
            PlotM=plot(PlotNotch,PlotSurroundWeight,'go');
            hold on;
            set(PlotM,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',[0.47 0.67 0.19],...
                'MarkerSize',8,...
                'Marker','o')
        elseif strcmp(PlotBin,'Achromatic')==1
            PlotA=plot(PlotNotch,PlotSurroundWeight,'ko');
            hold on;
            set(PlotA,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',[.5 .5 .5],...
                'MarkerSize',8,...
                'Marker','o')
        end
        %title(BinLabel);
        axis([0 1 0 1]);
        axis square
        Unity=line([0 1],[0 1]);
        set(Unity,...
            'LineStyle','--',...
            'Color','k')
    end
    clearvars BinnedCells
end


%% ONE BIG PLOT

figure;
%subplot(1,4,4);
hold on;
axis([0 1 0 1]);
axis square
Unity=line([0 1],[0 1]);
set(Unity,...
    'LineStyle','--',...
    'Color','k')
for e=2:length(bins)
    BinLabel=strcat('Bin',num2str(bins(e)));
    BinnedCells=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData;
    for c=1:size(BinnedCells,1)
        PlotNotch=Histograms.(matlab.lang.makeValidName(BinLabel)).CenterWeights(c);
        PlotSurroundWeight=Histograms.(matlab.lang.makeValidName(BinLabel)).SurroundWeights(c);
        PlotBin=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData{c,4};
        if strcmp(PlotBin,'L-dominated')==1
            PlotL=plot(PlotNotch,PlotSurroundWeight,'ro');
            hold on;
            set(PlotL,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',[0.64 0.08 0.18],...
                'MarkerSize',8,...
                'Marker','o')
        elseif strcmp(PlotBin,'M-dominated')==1
            PlotM=plot(PlotNotch,PlotSurroundWeight,'go');
            hold on;
            set(PlotM,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',[0.47 0.67 0.19],...
                'MarkerSize',8,...
                'Marker','o')
        elseif strcmp(PlotBin,'Achromatic')==1
            PlotA=plot(PlotNotch,PlotSurroundWeight,'ko');
            hold on;
            set(PlotA,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',[1 0.84 0],...
                'MarkerSize',8,...
                'Marker','o')
        end
        %title(BinLabel);
        
    end
    clearvars BinnedCells
end

title('Center-Surround Weights for 5000 Model Cells');
xlabel('Center L/(L+M)');
ylabel('Surround L/(L+M)');

%% ONE BIG PLOT - BETTER

LDoms=zeros(1,2);
MDoms=zeros(1,2);
Achroms=zeros(1,2);
    
l=1;
m=1;
a=1; 
for e=1:length(Eccentricities)
    CellLabel=strcat(Cells{e});
    %CellLabel=strcat('Cell',num2str(CellNumber(e)));
    PlotEccentricity=Eccentricities(e);
    PlotCenterWeight=Data.(matlab.lang.makeValidName(CellLabel)).CenterWeights(1);
    PlotSurroundWeight=Data.(matlab.lang.makeValidName(CellLabel)).SurroundWeights(1);
    PlotChromTag=Data.(matlab.lang.makeValidName(CellLabel)).ChromTag;
   
    if strcmp(PlotChromTag,'L-dominated')==1
        LDoms(l,:)=[PlotCenterWeight,PlotSurroundWeight];
        l=l+1;
    elseif strcmp(PlotChromTag,'M-dominated')==1
        MDoms(m,:)=[PlotCenterWeight,PlotSurroundWeight];
        m=m+1;
    elseif strcmp(PlotChromTag,'Achromatic')==1
        Achroms(a,:)=[PlotCenterWeight,PlotSurroundWeight];
        a=a+1;
    end
end
figure;
hold on;
axis([0 1 0 1]);
axis square
Unity=line([0 1],[0 1]);
set(Unity,...
    'LineStyle','--',...
    'Color','k')

PlotL=plot(LDoms(:,1),LDoms(:,2),'ro');
    set(PlotL,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[0.64 0.08 0.18],...
        'MarkerSize',8,...
        'Marker','o')
PlotM=plot(MDoms(:,1),MDoms(:,2),'go');
    set(PlotM,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[0.47 0.67 0.19],...
        'MarkerSize',8,...
        'Marker','o')
PlotA=plot(Achroms(:,1),Achroms(:,2),'go');
    set(PlotA,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[1 0.84 0],...
        'MarkerSize',8,...
        'Marker','o')

title('Center-Surround Weights for 5000 Model Cells');
xlabel('Center L/(L+M)');
ylabel('Surround L/(L+M)');

%% 3D PLOT

LDoms=zeros(1,3);
MDoms=zeros(1,3);
Achroms=zeros(1,3);
    
l=1;
m=1;
a=1;
for e=1:length(Eccentricities)
    CellLabel=strcat(Cells{e});
    %CellLabel=strcat('Cell',num2str(CellNumber(e)));
    PlotEccentricity=Eccentricities(e);
    PlotCenterWeight=Data.(matlab.lang.makeValidName(CellLabel)).CenterWeights(1);
    PlotSurroundWeight=Data.(matlab.lang.makeValidName(CellLabel)).SurroundWeights(1);
    PlotChromTag=Data.(matlab.lang.makeValidName(CellLabel)).ChromTag;
   
    if strcmp(PlotChromTag,'L-dominated')==1
        LDoms(l,:)=[PlotCenterWeight,PlotEccentricity,PlotSurroundWeight];
        l=l+1;
    elseif strcmp(PlotChromTag,'M-dominated')==1
        MDoms(m,:)=[PlotCenterWeight,PlotEccentricity,PlotSurroundWeight];
        m=m+1;
    elseif strcmp(PlotChromTag,'Achromatic')==1
        Achroms(a,:)=[PlotCenterWeight,PlotEccentricity,PlotSurroundWeight];
        a=a+1;
    end
end

figure;
xlabel('Center L/(L+M)');
ylabel('Eccentricity'); 
zlabel('Surround L/(L+M)'); 
axis([0 1 1 10 0 1]);
view([45 22]);
grid on;
hold on;

Unity=fill3([0,1,1,0],[0,0,10,10],[0,1,1,0],'k');
set(Unity,...
    'LineStyle','none',...
    'FaceAlpha',0.125)

PlotL=plot3(LDoms(:,1),LDoms(:,2),LDoms(:,3),'ro');
    set(PlotL,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[0.64 0.08 0.18],...
        'MarkerSize',8,...
        'Marker','o')
PlotM=plot3(MDoms(:,1),MDoms(:,2),MDoms(:,3),'go');
    set(PlotM,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[0.47 0.67 0.19],...
        'MarkerSize',8,...
        'Marker','o')
PlotA=plot3(Achroms(:,1),Achroms(:,2),Achroms(:,3),'go');
    set(PlotA,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[1 0.84 0],...
        'MarkerSize',8,...
        'Marker','o')

%% PERCENT CHROMATIC PER ECC

figure;
hold on;
for e=2:length(bins)
    BinLabel=strcat('Bin',num2str(bins(e)));
    ChromaticCells(e-1,:)=((Histograms.(matlab.lang.makeValidName(BinLabel)).LDoms)+(Histograms.(matlab.lang.makeValidName(BinLabel)).MDoms))/length(Histograms.(matlab.lang.makeValidName(BinLabel)).AllData);
    AchromaticCells(e-1,:)=((Histograms.(matlab.lang.makeValidName(BinLabel)).AChroms))/length(Histograms.(matlab.lang.makeValidName(BinLabel)).AllData);
end
PlotChroms=plot(bins(2:end),ChromaticCells);
hold on;
PlotAchroms=plot(bins(2:end),AchromaticCells); 
set(PlotChroms,...
    'LineStyle','none',...
    'Color','k',...
    'MarkerFaceColor',[0.64 0.08 0.18],...
    'MarkerSize',8,...
    'Marker','o')
set(PlotAchroms,...
    'LineStyle','none',...
    'Color','k',...
    'MarkerFaceColor',[1 0.84 0],...
    'MarkerSize',8,...
    'Marker','o')

%% PERCENT THREE TYPES
figure;
hold on;
for e=2:length(bins)
    BinLabel=strcat('Bin',num2str(bins(e)));
    NotchCells(e-1,:)=(Histograms.(matlab.lang.makeValidName(BinLabel)).NotchCells)/(length(Histograms.(matlab.lang.makeValidName(BinLabel)).AllData));
    ChromaticCells(e-1,:)=((Histograms.(matlab.lang.makeValidName(BinLabel)).LDoms)+(Histograms.(matlab.lang.makeValidName(BinLabel)).MDoms)-(Histograms.(matlab.lang.makeValidName(BinLabel)).NotchCells))/length(Histograms.(matlab.lang.makeValidName(BinLabel)).AllData);
    AchromaticCells(e-1,:)=((Histograms.(matlab.lang.makeValidName(BinLabel)).AChroms))/length(Histograms.(matlab.lang.makeValidName(BinLabel)).AllData);
end

PlotChroms=plot(bins(2:end),ChromaticCells);
hold on;
PlotAchroms=plot(bins(2:end),AchromaticCells); 
PlotNotch=plot(bins(2:end),NotchCells);

set(PlotChroms,...
    'LineStyle','none',...
    'Color','k',...
    'MarkerFaceColor',[0.64 0.08 0.18],...
    'MarkerSize',8,...
    'Marker','o')
set(PlotAchroms,...
    'LineStyle','none',...
    'Color','k',...
    'MarkerFaceColor',[1 0.84 0],...
    'MarkerSize',8,...
    'Marker','o')
set(PlotNotch,...
    'LineStyle','none',...
    'Color','k',...
    'MarkerFaceColor',[0.64 0.08 0.18],...
    'MarkerSize',8,...
    'Marker','o')

%%
figure;
hold on;

for e=2:length(bins)
    BinLabel=strcat('Bin',num2str(bins(e)));
    NotchCells(e-1,:)=(Histograms.(matlab.lang.makeValidName(BinLabel)).NotchCells)/(length(Histograms.(matlab.lang.makeValidName(BinLabel)).AllData));
    ChromaticCells(e-1,:)=((Histograms.(matlab.lang.makeValidName(BinLabel)).LDoms)+(Histograms.(matlab.lang.makeValidName(BinLabel)).MDoms)-(Histograms.(matlab.lang.makeValidName(BinLabel)).NotchCells))/length(Histograms.(matlab.lang.makeValidName(BinLabel)).AllData);
    AchromaticCells(e-1,:)=((Histograms.(matlab.lang.makeValidName(BinLabel)).AChroms))/length(Histograms.(matlab.lang.makeValidName(BinLabel)).AllData);
end

PlotChroms=plot(bins(2:end),ChromaticCells);
hold on;
PlotAchroms=plot(bins(2:end),AchromaticCells); 
PlotNotch=plot(bins(2:end),NotchCells);

set(PlotChroms,...
    'LineStyle','none',...
    'Color','k',...
    'MarkerFaceColor',[0.64 0.08 0.18],...
    'MarkerSize',8,...
    'Marker','o')
set(PlotAchroms,...
    'LineStyle','none',...
    'Color','k',...
    'MarkerFaceColor',[1 0.84 0],...
    'MarkerSize',8,...
    'Marker','o')
set(PlotNotch,...
    'LineStyle','none',...
    'Color','k',...
    'MarkerFaceColor',[0.64 0.08 0.18],...
    'MarkerSize',8,...
    'Marker','o')

%% RG GAIN PLOT

LDoms=zeros(1,4);
MDoms=zeros(1,4);
Achroms=zeros(1,4);
    
l=1;
m=1;
a=1;
for e=1:length(Eccentricities)
    CellLabel=strcat(Cells{e});
    %CellLabel=strcat('Cell',num2str(CellNumber(e)));
    PlotEccentricity=Eccentricities(e);
    PlotCenterWeight=Data.(matlab.lang.makeValidName(CellLabel)).CenterWeights(1);
    PlotSurroundWeight=Data.(matlab.lang.makeValidName(CellLabel)).SurroundWeights(1);
    PlotChromTag=Data.(matlab.lang.makeValidName(CellLabel)).ChromTag;
    PlotGain=LumGains(e);
   
    if strcmp(PlotChromTag,'L-dominated')==1
        LDoms(l,:)=[PlotEccentricity,PlotCenterWeight,PlotSurroundWeight,PlotGain];
        l=l+1;
    elseif strcmp(PlotChromTag,'M-dominated')==1
        MDoms(m,:)=[PlotEccentricity,PlotCenterWeight,PlotSurroundWeight,PlotGain];
        m=m+1;
    elseif strcmp(PlotChromTag,'Achromatic')==1
        Achroms(a,:)=[PlotEccentricity,PlotCenterWeight,PlotSurroundWeight,PlotGain];
        a=a+1;
    end
end
figure;
hold on;
axis([0 10 0 15]);
axis normal
Unity=line([0 10],[1 1]);
set(Unity,...
    'LineStyle','--',...
    'Color','k')

PlotL=plot(LDoms(:,1),LDoms(:,4),'ro');
    set(PlotL,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[0.64 0.08 0.18],...
        'MarkerSize',8,...
        'Marker','o')
PlotM=plot(MDoms(:,1),MDoms(:,4),'go');
    set(PlotM,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[0.47 0.67 0.19],...
        'MarkerSize',8,...
        'Marker','o')
PlotA=plot(Achroms(:,1),Achroms(:,4),'go');
    set(PlotA,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[1 0.84 0],...
        'MarkerSize',8,...
        'Marker','o')

title('Center-Surround Weights for 5000 Model Cells');
xlabel('Eccentricity (mm)');
ylabel('(L-M)/(L+M)');

%% PLOT PERCENT GAIN
figure;
hold on;
a=1;
c=1;
for e=2:length(bins)
    BinLabel=strcat('Bin',num2str(bins(e)));
%     for d=1:length(Histograms.(matlab.lang.makeValidName(BinLabel)).AllData)
%     if strcmp((Histograms.(matlab.lang.makeValidName(BinLabel)).AllData{d,4}),'Achromatic')==1;
%         AchRGGain(a,:)=Histograms.(matlab.lang.makeValidName(BinLabel)).RGGain(d);
%         AchLumGain(a,:)=Histograms.(matlab.lang.makeValidName(BinLabel)).LumGain(d);
%         a=a+1;
%     else
%         ChRGGain(c,:)=Histograms.(matlab.lang.makeValidName(BinLabel)).RGGain(d);
%         ChLumGain(c,:)=Histograms.(matlab.lang.makeValidName(BinLabel)).LumGain(d);
%         c=c+1;
%     end
%     MeanChRGGain(e-1,:)=mean(ChRGGain);
%     SEMChRGGain(e-1,:)=std(ChRGGain)/sqrt(length(ChRGGain));
%     MeanChLumGain(e-1,:)=mean(ChLumGain);
%     SEMChLumGain(e-1,:)=std(ChLumGain)/sqrt(length(ChLumGain));
%     MeanAchRGGain(e-1,:)=mean(AchRGGain);
%     
%     SEMAchRGGain(e-1,:)=std(AchRGGain)/sqrt(length(AchRGGain));
%     MeanAchLumGain(e-1,:)=mean(AchLumGain);
%     SEMAchLumGain(e-1,:)=std(AchLumGain)/sqrt(length(AchLumGain));
%     end
    
    MeanRGGain(e-1,:)=mean(Histograms.(matlab.lang.makeValidName(BinLabel)).RGGain);
    SEMRGGain(e-1,:)=(std(Histograms.(matlab.lang.makeValidName(BinLabel)).RGGain))/sqrt((length(Histograms.(matlab.lang.makeValidName(BinLabel)).RGGain)));
    MeanLumGain(e-1,:)=mean(Histograms.(matlab.lang.makeValidName(BinLabel)).LumGain);
    SEMLumGain(e-1,:)=(std(Histograms.(matlab.lang.makeValidName(BinLabel)).LumGain))/sqrt((length(Histograms.(matlab.lang.makeValidName(BinLabel)).LumGain)));
end
plotbins=bins(2:end);
for j=1:length(plotbins)
    PlotRGSEM=line([plotbins(j) plotbins(j)],[MeanRGGain(j)-SEMRGGain(j) MeanRGGain(j)+SEMRGGain(j)]);
    PlotLumSEM=line([plotbins(j) plotbins(j)],[MeanLumGain(j)-SEMLumGain(j) MeanLumGain(j)+SEMLumGain(j)]);
    set([PlotRGSEM PlotLumSEM],...
    'LineStyle','-',...
    'Color','k')
end
PlotRGGain=plot(bins(2:end),MeanRGGain);
hold on;
PlotLumGain=plot(bins(2:end),MeanLumGain); 


set(PlotRGGain,...
    'LineStyle','none',...
    'Color','k',...
    'MarkerFaceColor',[0.64 0.08 0.18],...
    'MarkerSize',8,...
    'Marker','o')
set(PlotLumGain,...
    'LineStyle','none',...
    'Color','k',...
    'MarkerFaceColor',[1 0.84 0],...
    'MarkerSize',8,...
    'Marker','o')

%% ONE BIG PLOT - WITH S

LDoms=zeros(1,2);
MDoms=zeros(1,2);
Achroms=zeros(1,2);
    
l=1;
m=1;
a=1;
for e=1:length(Eccentricities)
    CellLabel=strcat(Cells{e});
    %CellLabel=strcat('Cell',num2str(CellNumber(e)));
    PlotEccentricity=Eccentricities(e);
    PlotCenterWeight=Data.(matlab.lang.makeValidName(CellLabel)).CenterWeights(1);
    PlotSurroundWeight=Data.(matlab.lang.makeValidName(CellLabel)).SurroundWeights(1);
    PlotChromTag=Data.(matlab.lang.makeValidName(CellLabel)).ChromTag;
   
    if strcmp(PlotChromTag,'L-dominated')==1
        LDoms(l,:)=[PlotCenterWeight,PlotSurroundWeight];
        l=l+1;
    elseif strcmp(PlotChromTag,'M-dominated')==1
        MDoms(m,:)=[PlotCenterWeight,PlotSurroundWeight];
        m=m+1;
    elseif strcmp(PlotChromTag,'Achromatic')==1
        Achroms(a,:)=[PlotCenterWeight,PlotSurroundWeight];
        a=a+1;
    end
    
    if SWeights(e,1)<0.1
        MarkerWeight=0;
    elseif SWeights(e,1)>=0.1 && SWeights(e,1)<0.2
        MarkerWeight=1;
    elseif SWeights(e,1)>=0.3 && SWeights(e,1)<0.3
        MarkerWeight=2;
    elseif SWeights(e,1)>=0.4 && SWeights(e,1)<0.5
        MarkerWeight=3;
    elseif SWeights(e,1)>=0.5 && SWeights(e,1)<0.6
        MarkerWeight=4;
    elseif SWeights(e,1)>=0.6 && SWeights(e,1)<0.7
        MarkerWeight=5;
    end
end
figure;
hold on;
axis([0 1 0 1]);
axis square
Unity=line([0 1],[0 1]);
set(Unity,...
    'LineStyle','--',...
    'Color','k')

PlotL=plot(LDoms(:,1),LDoms(:,2),'ro');
    set(PlotL,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[0.64 0.08 0.18],...
        'MarkerSize',8,...
        'Marker','o')
PlotM=plot(MDoms(:,1),MDoms(:,2),'go');
    set(PlotM,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[0.47 0.67 0.19],...
        'MarkerSize',8,...
        'Marker','o')
PlotA=plot(Achroms(:,1),Achroms(:,2),'go');
    set(PlotA,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[1 0.84 0],...
        'MarkerSize',8,...
        'Marker','o')

title('Center-Surround Weights for 5000 Model Cells');
xlabel('Center L/(L+M)');
ylabel('Surround L/(L+M)');

%% S GAIN PLOT

LDoms=zeros(1,4);
MDoms=zeros(1,4);
Achroms=zeros(1,4);
    
l=1;
m=1;
a=1;

figure;
hold on;
axis([0 10 0 100]);
axis normal
for e=1:length(Eccentricities)
    CellLabel=strcat(Cells{e});
    %CellLabel=strcat('Cell',num2str(CellNumber(e)));
    PlotEccentricity=Eccentricities(e);
    PlotCenterWeight=Data.(matlab.lang.makeValidName(CellLabel)).CenterWeights(1);
    PlotSurroundWeight=Data.(matlab.lang.makeValidName(CellLabel)).SurroundWeights(1);
    PlotChromTag=Data.(matlab.lang.makeValidName(CellLabel)).ChromTag;
    PlotGain=RGGains(e);
    PlotSWeight=SWeights(e);
   
    if strcmp(PlotChromTag,'L-dominated')==1
        LDoms(l,:)=[PlotEccentricity,PlotCenterWeight,PlotSurroundWeight,PlotGain];
        l=l+1;
        color=[0.64 0.08 0.18];
    elseif strcmp(PlotChromTag,'M-dominated')==1
        MDoms(m,:)=[PlotEccentricity,PlotCenterWeight,PlotSurroundWeight,PlotGain];
        m=m+1;
        color=[0.47 0.67 0.19];
    elseif strcmp(PlotChromTag,'Achromatic')==1
        Achroms(a,:)=[PlotEccentricity,PlotCenterWeight,PlotSurroundWeight,PlotGain];
        a=a+1;
        color=[0.5 0.5 0.5];
    end
%     PlotLM=plot(PlotEccentricity,PlotGain,'ro');
%     set(PlotLM,...
%         'LineStyle','none',...
%         'Color',color,...
%         'MarkerFaceColor',color,...
%         'MarkerSize',12,...
%         'Marker','o')
    PlotS=plot(PlotEccentricity,PlotGain,'ro');
    hold on;
    set(PlotS,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',hsv2rgb([0.58 (PlotSWeight) .75]),...
        'MarkerSize',8,...
        'Marker','o')
end
figure;
hold on;
axis([0 10 0 15]);
axis normal
Unity=line([0 10],[1 1]);
set(Unity,...
    'LineStyle','--',...
    'Color','k')

PlotL=plot(LDoms(:,1),LDoms(:,4),'ro');
    set(PlotL,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[1-PlotSWeight 1-PlotSWeight 1],...
        'MarkerSize',8,...
        'Marker','o')
PlotM=plot(MDoms(:,1),MDoms(:,4),'go');
    set(PlotM,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[1-PlotSWeight 1-PlotSWeight 1],...
        'MarkerSize',8,...
        'Marker','o')
PlotA=plot(Achroms(:,1),Achroms(:,4),'go');
    set(PlotA,...
        'LineStyle','none',...
        'Color','k',...
        'MarkerFaceColor',[1-PlotSWeight 1-PlotSWeight 1],...
        'MarkerSize',8,...
        'Marker','o')

title('Center-Surround Weights for 5000 Model Cells');
xlabel('Eccentricity (mm)');
ylabel('(L-M)/(L+M)');

%% HISTOGRAMS

PureL=LDoms(find(LDoms(:,1)>.999),:)
PureM=MDoms(find(MDoms(:,1)==0),:)
MixedL=LDoms(find(LDoms(:,1)<.999),:)
MixedM=MDoms(find(MDoms(:,1)~=0),:)

figure;
subplot(2,3,1);
hist([PureL(:,1); PureM(:,1)],30);
subplot(2,3,2);
hist([MixedL(:,1); MixedM(:,1)],30);
subplot(2,3,3);
hist(Achroms(:,1),30);
subplot(2,3,4);
hist([PureL(:,2); PureM(:,2)],30);
subplot(2,3,5);
hist([MixedL(:,2); MixedM(:,2)],30);
subplot(2,3,6);
hist(Achroms(:,2),30);

figure;
subplot(1,2,1);
hist([MixedL(:,1); MixedM(:,1);Achroms(:,1)],30);
subplot(1,2,2);
hist([MixedL(:,2); MixedM(:,2);Achroms(:,2)],30)

figure;
subplot(1,2,1);
histogram([Achroms(:,1)],EDGES);
subplot(1,2,2);
histogram([Achroms(:,2);],EDGES)

%% ONE BIG PLOT WITH GRADIENT COLORING

figure;
hold on;
axis([0 1 0 1]);
axis square
Unity=line([0 1],[0 1]);
set(Unity,...
    'LineStyle','--',...
    'Color','k')
    
for e=2:length(bins)
    BinLabel=strcat('Bin',num2str(bins(e)));
    BinnedCells=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData;
    for c=1:size(BinnedCells,1)
        PlotCenterWeight=Histograms.(matlab.lang.makeValidName(BinLabel)).CenterWeights(c);
        PlotSurroundWeight=Histograms.(matlab.lang.makeValidName(BinLabel)).SurroundWeights(c);
        PlotRGResp=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData{c,7};
        PlotRGPhase=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData{c,9};
        if PlotRGPhase<5
            PlotL=plot(PlotCenterWeight,PlotSurroundWeight,'ro');
            hold on;
            set(PlotL,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',hsv2rgb([0.99 (PlotRGResp/100) .75]),... %for gain %'MarkerFaceColor',hsv2rgb([0.99 PlotRGResp/100 .75]),... %for resp
                'MarkerSize',8,...
                'Marker','o')
        elseif PlotRGPhase>175
            PlotM=plot(PlotCenterWeight,PlotSurroundWeight,'go');
            hold on;
            set(PlotM,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',hsv2rgb([0.236 (PlotRGResp/100) .75]),... %for gain %'MarkerFaceColor',hsv2rgb([0.236 PlotRGResp/100 .75]),... %for resp
                'MarkerSize',8,...
                'Marker','o')
        end
        %title(BinLabel);
        
    end
    clearvars BinnedCells
end

title('Center-Surround Weights for 5000 Model Cells');
xlabel('Center L/(L+M)');
ylabel('Surround L/(L+M)');

%%
figure;
hold on;
axis([0 1 0 1]);
axis square
Unity=line([0 1],[0 1]);
set(Unity,...
    'LineStyle','--',...
    'Color','k')

for e=1:length(Eccentricities)
    CellLabel=strcat(Cells{e});
    %CellLabel=strcat('Cell',num2str(CellNumber(e)));
    PlotEccentricity=Eccentricities(e);
    PlotCenterWeight=Data.(matlab.lang.makeValidName(CellLabel)).CenterWeights(1);
    PlotSurroundWeight=Data.(matlab.lang.makeValidName(CellLabel)).SurroundWeights(1);
    PlotChromTag=Data.(matlab.lang.makeValidName(CellLabel)).ChromTag;
    PlotRGResp=RGGains(e);
    PlotLumResp=LumGains(e);
    PlotRGPhase=RGPhases(e);
    if PlotRGPhase<5
            PlotL=plot(PlotCenterWeight,PlotSurroundWeight,'ro');
            hold on;
            set(PlotL,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',hsv2rgb([0.99 PlotRGResp/(PlotRGResp+PlotLumResp) .75]),... %for resp
                'MarkerSize',8,...
                'Marker','o')
        elseif PlotRGPhase>175
            PlotM=plot(PlotCenterWeight,PlotSurroundWeight,'go');
            hold on;
            set(PlotM,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',hsv2rgb([0.236 PlotRGResp/(PlotRGResp+PlotLumResp) .75]),... %for resp
                'MarkerSize',8,...
                'Marker','o')
    end
end

%% DIAMOND PLOT FOR EACH BIN
% figure;
% hold on;
for e=2:length(bins)
    figure;
    hold on;
    BinLabel=strcat('Bin',num2str(bins(e)));
    BinnedCells=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData;
    for c=1:length(BinnedCells)
        Center=[Histograms.(matlab.lang.makeValidName(BinLabel)).CenterWeights(c) 1-Histograms.(matlab.lang.makeValidName(BinLabel)).CenterWeights(c)];
        Surround=[Histograms.(matlab.lang.makeValidName(BinLabel)).SurroundWeights(c) 1-Histograms.(matlab.lang.makeValidName(BinLabel)).SurroundWeights(c)];
        CSRatio=Histograms.(matlab.lang.makeValidName(BinLabel)).CSRatio(c);
        DKL=(Center-CSRatio*Surround);
        normDKL=DKL/sum(abs(DKL))*Histograms.(matlab.lang.makeValidName(BinLabel)).Polarity(c);
        PlotBin=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData{c,4};
        if strcmp(PlotBin,'L-dominated')==1
            PlotL=plot(normDKL(1),normDKL(2),'ko');
            hold on;
%             set(PlotL,...
%                 'LineStyle','none',...
%                 'Color','k',...
%                 'MarkerFaceColor',[0.64 0.08 0.18],...
%                 'MarkerSize',5,...
%                 'Marker','o')
        elseif strcmp(PlotBin,'M-dominated')==1
            PlotM=plot(normDKL(1),normDKL(2),'ko');
            hold on;
%             set(PlotM,...
%                 'LineStyle','none',...
%                 'Color','k',...
%                 'MarkerFaceColor',[0.47 0.67 0.19],...
%                 'MarkerSize',5,...
%                 'Marker','o')
        elseif strcmp(PlotBin,'Achromatic')==1
            PlotA=plot(normDKL(1),normDKL(2),'ko');
            hold on;
%             set(PlotA,...
%                 'LineStyle','none',...
%                 'Color','k',...
%                 'MarkerFaceColor',[1 0.84 0],...
%                 'MarkerSize',5,...
%                 'Marker','o')
        end
        %title(BinLabel);
    end
    axis([-1 1 -1 1]);
    axis square
    Ord=line([0 0],[-1 1]);
    Abs=line([-1 1],[0 0]);
    Diag3=line([0 -1],[-1 0]);
    Diag1=line([0 1],[1 0]);
    Diag4=line([0 1],[-1 0]);
    Diag2=line([0 -1],[1 0]);
    set([Ord,Abs,Diag1,Diag2,Diag3,Diag4],...
        'LineStyle','-',...
        'Color','k')
    clearvars BinnedCells
end

%% 3D DIAMOND PLOT
% figure;
% hold on;
normDKL_all = [];
for e=2:length(bins)
    BinLabel=strcat('Bin',num2str(bins(e)));
    BinnedCells=Histograms.(matlab.lang.makeValidName(BinLabel)).AllData;
    for c=1:length(BinnedCells)
        Center=[Histograms.(matlab.lang.makeValidName(BinLabel)).CenterWeights(c) 1-Histograms.(matlab.lang.makeValidName(BinLabel)).CenterWeights(c)];
        Surround=[Histograms.(matlab.lang.makeValidName(BinLabel)).SurroundWeights(c) 1-Histograms.(matlab.lang.makeValidName(BinLabel)).SurroundWeights(c)];
        CSRatio=Histograms.(matlab.lang.makeValidName(BinLabel)).CSRatio(c);
        DKL=(Center-CSRatio*Surround);
        normDKL=DKL/sum(abs(DKL))*Histograms.(matlab.lang.makeValidName(BinLabel)).Polarity(c);
        normDKL_all.(matlab.lang.makeValidName(BinLabel))(c,:) = [normDKL(1)+normDKL(2),normDKL(2)-normDKL(1)];
    end
    clearvars BinnedCells
end
test = [Selectivity6_DataPoints(:,1)+Selectivity6_DataPoints(:,2),Selectivity6_DataPoints(:,2)-Selectivity6_DataPoints(:,1)]
%%
hold on;
subplot(2,3,4)
hist3(test,[20 20]);
set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
axis([-1 1 -1 1]);
axis square
Diag1=line([1 -1 0],[-1 1 0]);
Diag2=line([-1 1 0],[-1 1 0]);
    set([Diag1,Diag2],...
        'LineStyle','-',...
        'Color','k')
view(-45,60);
grid off
