NumberOfCells=20;
MinEccentricity=0.25;
MaxEccentricity=1;
Eccentricities=(MaxEccentricity-MinEccentricity).*rand(NumberOfCells,1) + MinEccentricity;
CellNumber=linspace(1,length(Eccentricities),length(Eccentricities)); %How many cells per sample?
    
for e=1:length(Eccentricities)
    EccentricityLabel=strcat('Ecc',num2str(Eccentricities(e)),'mm');
    CellLabel=strcat('Cell',num2str(CellNumber(e)));
    Cell=DoGfrom2D_NoConv(Eccentricities(e));
    Cell.Eccentricity=Eccentricities(e);
    Data.(matlab.lang.makeValidName(EccentricityLabel))=Cell;
    disp(strcat(CellLabel,' @ ',EccentricityLabel,' complete.'));
end

figure;
xlabel('Center L/(L+M)');
ylabel('Eccentricity'); 
zlabel('Surround L/(L+M)'); 
axis([0 1 0.25 10 0 1]);
view([45 22]);
grid on;
hold on;
Unity=fill3([0,1,1,0],[0,0,10,10],[0,1,1,0],'k');
set(Unity,...
    'LineStyle','none',...
    'FaceAlpha',0.125)
for e=1:length(Eccentricities)
    EccentricityLabel=strcat('Ecc',num2str(Eccentricities(e)),'mm');
    CellLabel=strcat('Cell',num2str(CellNumber(e)));
    PlotEccentricity=Eccentricities(e);
    PlotCenterWeight=Data.(matlab.lang.makeValidName(EccentricityLabel)).CenterWeights(1);
    PlotSurroundWeight=Data.(matlab.lang.makeValidName(EccentricityLabel)).SurroundWeights(1);
    PlotChromTag=Data.(matlab.lang.makeValidName(EccentricityLabel)).ChromTag;
    if strcmp(PlotChromTag,'L-dominated')==1
        PlotL=plot3(PlotCenterWeight,PlotEccentricity,PlotSurroundWeight,'ro');
        set(PlotL,...
            'LineStyle','none',...
            'Color','k',...
            'MarkerFaceColor',[0.64 0.08 0.18],...
            'MarkerSize',8,...
            'Marker','o')
    elseif strcmp(PlotChromTag,'M-dominated')==1
        PlotM=plot3(PlotCenterWeight,PlotEccentricity,PlotSurroundWeight,'go');
        set(PlotM,...
            'LineStyle','none',...
            'Color','k',...
            'MarkerFaceColor',[0.47 0.67 0.19],...
            'MarkerSize',8,...
            'Marker','o')
    elseif strcmp(PlotChromTag,'Achromatic')==1
        PlotA=plot3(PlotCenterWeight,PlotEccentricity,PlotSurroundWeight,'ko');
        set(PlotA,...
            'LineStyle','none',...
            'Color','k',...
            'MarkerFaceColor',[1 0.84 0],...
            'MarkerSize',8,...
            'Marker','o')
    end
end

