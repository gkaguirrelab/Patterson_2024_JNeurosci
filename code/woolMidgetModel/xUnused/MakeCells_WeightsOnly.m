NumberOfEccentricities=[5]; %Sample cells at each eccentricity
NumberOfCells=linspace(1,100,100); %How many cells per sample?
    
for e=1:length(NumberOfEccentricities)
    EccentricityLabel=strcat('Ecc',num2str(NumberOfEccentricities(e)),'mm');
    for c=1:length(NumberOfCells)
        CellLabel=strcat('Cell',num2str(NumberOfCells(c)));
        Cell=DoG_WeightsOnly(NumberOfEccentricities(e));
        Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel))=Cell;
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
for e=1:length(NumberOfEccentricities)
    EccentricityLabel=strcat('Ecc',num2str(NumberOfEccentricities(e)),'mm');
    for c=1:length(NumberOfCells)
        CellLabel=strcat('Cell',num2str(NumberOfCells(c)));
        PlotCenterWeight=Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).CenterWeight(1);
        PlotSurroundWeight=Data.(matlab.lang.makeValidName(EccentricityLabel)).(matlab.lang.makeValidName(CellLabel)).SurroundWeight(1);
        PlotL=plot(PlotCenterWeight,PlotSurroundWeight,'ro');
        set(PlotL,...
            'LineStyle','none',...
            'Color','k',...
            'MarkerFaceColor',[0.5 0.5 0.5],...
            'MarkerSize',8,...
            'Marker',EccMarkers(e))
    end
end

