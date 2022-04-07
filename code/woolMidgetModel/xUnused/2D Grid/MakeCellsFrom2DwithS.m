NumberOfEccentricities=[1 2 3 4 5 6 7 8 9 10]; %Sample cells at each eccentricity
NumberOfCells=linspace(1,100,100); %How many cells per sample?
    
for e=1:length(NumberOfEccentricities)
    EccentricityLabel=strcat('Ecc',num2str(NumberOfEccentricities(e)),'mm');
    for c=1:length(NumberOfCells)
        CellLabel=strcat('Cell',num2str(NumberOfCells(c)));
        Cell=DoGfrom2DwithS(NumberOfEccentricities(e));
        Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel))=Cell;
        disp(strcat(CellLabel,' @ ',EccentricityLabel,' complete.'));
    end
    disp(strcat(EccentricityLabel,' complete.'));
end

figure;
hold on;
% title(strcat(EccentricityLabel(4:end)));
% xlabel('Center L/(L+M)');
% ylabel('Surround L/(L+M)'); 
% axis([0 1 0 1]);
% axis square
% Unity=line([0 1],[0 1]);
% set(Unity,...
%     'LineStyle','--',...
%     'Color','k')
% EccMarkers=['o' 'd' '^' 'v' '<' '>' 's' 'p' 'h' '*'];
EccMarkers=['o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o' 'o'];
for e=1:length(NumberOfEccentricities)
    EccentricityLabel=strcat('Ecc',num2str(NumberOfEccentricities(e)),'mm');
    subplot(2,5,e);
    title(strcat(EccentricityLabel(4:end)));
    xlabel('Center L/(L+M)');
    ylabel('Surround L/(L+M)'); 
    axis([0 1 0 1]);
    axis square
    Unity=line([0 1],[0 1]);
    set(Unity,...
        'LineStyle','--',...
        'Color','k')
    for c=1:length(NumberOfCells)
        CellLabel=strcat('Cell',num2str(NumberOfCells(c)));
        PercentS=ceil((Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).CenterWeights(3)*100)/10);
        PlotCenterWeight=Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).CenterWeights(1);
        hold on;
        PlotSurroundWeight=Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).SurroundWeights(1);
        PlotChromTag=Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).ChromTag;
        if strcmp(PlotChromTag,'L-dominated')==1
            PlotL=plot(PlotCenterWeight,PlotSurroundWeight,'ro');
            set(PlotL,...
                'LineStyle','-',...
                'LineWidth',.1+PercentS/2,...
                'Color',[0.08 .17 .55],...
                'MarkerFaceColor',[0.64 0.08 0.18],...
                'MarkerSize',8,...
                'Marker',EccMarkers(e))
        elseif strcmp(PlotChromTag,'M-dominated')==1
            PlotM=plot(PlotCenterWeight,PlotSurroundWeight,'go');
            set(PlotM,...
                'LineStyle','-',...
                'LineWidth',.1+PercentS/2,...
                'Color',[0.08 .17 .55],...
                'MarkerFaceColor',[0.47 0.67 0.19],...
                'MarkerSize',8,...
                'Marker',EccMarkers(e))
        elseif strcmp(PlotChromTag,'Achromatic')==1
            PlotA=plot(PlotCenterWeight,PlotSurroundWeight,'ko');
            set(PlotA,...
                'LineStyle','-',...
                'LineWidth',.1+PercentS/2,...
                'Color',[0.08 .17 .55],...
                'MarkerFaceColor',[1 0.84 0],...
                'MarkerSize',8,...
                'Marker',EccMarkers(e))
        end
    end
end


% CellLAmpData=Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).ResponseFunctions.LIsoResponse.Amplitude;
% CellMAmpData=Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).ResponseFunctions.MIsoResponse.Amplitude;
% CellSAmpData=Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).ResponseFunctions.SIsoResponse.Amplitude;
% CenterSWeight=Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).CenterWeights(3);
% 
% PercentS=ceil((Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).CenterWeights(3)*100)/10);
% PercentS=ceil((CenterSWeight*100)/10);
% PercentS=ceil(((max(CellSAmpData)/(max(CellMAmpData)+max(CellLAmpData)+max(CellSAmpData)))*100)/10);
% PercentSLabel=[...
%     29 32 85;
%     32 65 127;
%     0 91 161;
%     6 115 186;
%     49 143 205;
%     87 170 222;
%     115 188 232;
%     140 207 242;
%     169 220 250;
%     210 240 255]/255;
% 
% PercentSLabel=[...
%     210 240 255;
%     169 220 250;
%     140 207 242;
%     115 188 232;
%     87 170 222;
%     49 143 205;
%     6 115 186;
%     0 91 161;
%     32 65 127;
%     29 32 85]/255;

