NumberOfEccentricities=[1 2 3 4 5 6 7 8 9 10]; %Sample cells at each eccentricity
NumberOfCells=linspace(1,100,100); %How many cells per sample?
    
for e=1:length(NumberOfEccentricities)
EccentricityLabel=strcat('Ecc',num2str(NumberOfEccentricities(e)),'mm');
    for c=1:length(NumberOfCells)
        CellLabel=strcat('Cell',num2str(NumberOfCells(c)));
        Cell=DoG(NumberOfEccentricities(e));
        Data.(matlab.lang.makeValidName(RetinaLabel)).(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).RetinaWeights.pL=pL;
        Data.(matlab.lang.makeValidName(RetinaLabel)).(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).RetinaWeights.pM=pM;
        Data.(matlab.lang.makeValidName(RetinaLabel)).(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).RetinaWeights.pS=pS;
        Data.(matlab.lang.makeValidName(RetinaLabel)).(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel))=Cell;
        disp(strcat(CellLabel,' @ ',EccentricityLabel,' complete.'));
    end
disp(strcat(EccentricityLabel,' complete.'));
end


figure;
hold on;
xlabel('Center L/(L+M)');
ylabel('Surround L/(L+M)'); 
axis([0 1 0 1]);
axis square
Unity=line([0 1],[0 1]);
set(Unity,...
    'LineStyle','--',...
    'Color','k')
EccMarkers=['o' 'd' '^' 'v' '<' '>' 's' 'p' 'h' '*'];
for r=1:NumberOfRetinas
RetinaLabel=strcat('LM',num2str(RetinaWeights(r)));
for e=1:length(NumberOfEccentricities)
    EccentricityLabel=strcat('Ecc',num2str(NumberOfEccentricities(e)),'mm');
    for c=1:length(NumberOfCells)
        CellLabel=strcat('Cell',num2str(NumberOfCells(c)));
        PlotCenterWeight=Data.(matlab.lang.makeValidName(RetinaLabel)).(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).CenterWeight(1);
        PlotSurroundWeight=Data.(matlab.lang.makeValidName(RetinaLabel)).(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).SurroundWeight(1);
        PlotChromTag=Data.(matlab.lang.makeValidName(RetinaLabel)).(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).ChromTag;
        if strcmp(PlotChromTag,'L-dominated')==1
            PlotL=plot(PlotCenterWeight,PlotSurroundWeight,'ro');
            set(PlotL,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',[0.64 0.08 0.18],...
                'MarkerSize',8,...
                'Marker',EccMarkers(e))
        elseif strcmp(PlotChromTag,'M-dominated')==1
            PlotM=plot(PlotCenterWeight,PlotSurroundWeight,'go');
            set(PlotM,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',[0.47 0.67 0.19],...
                'MarkerSize',8,...
                'Marker',EccMarkers(e))
        elseif strcmp(PlotChromTag,'Achromatic')==1
            PlotA=plot(PlotCenterWeight,PlotSurroundWeight,'ko');
            set(PlotA,...
                'LineStyle','none',...
                'Color','k',...
                'MarkerFaceColor',[1 0.84 0],...
                'MarkerSize',8,...
                'Marker',EccMarkers(e))
        end
    end
end
end

