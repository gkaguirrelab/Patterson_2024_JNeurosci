function PlotRFCharacteristics(Cell)
%% IMPORT DATA STRUCTURE

%Identifying information
% ImportCell=Data.Ecc1mm.Cell5;
% CellName='Data.Ecc1mm.Cell5';
%ImportCell=Cell;
CellName='Cell'
%CellName=strcat('Ecc',num2str(Cell.Eccentricity));

%Load variables
SortedTotalMosaic=ImportCell.Global.AllConeCoords;
RetinaWeights=ImportCell.RetinaWeights;
CenterWeights=ImportCell.CenterWeights;
SurroundWeights=ImportCell.SurroundWeights;
LCoords=ImportCell.Global.LConeCoords;
MCoords=ImportCell.Global.MConeCoords;
SCoords=ImportCell.Global.SConeCoords;

CenterX=ImportCell.Center.CenterConeDistances;
SurroundX=ImportCell.Surround.SurroundConeDistances;
ConesToCenter=ImportCell.Center.ConesToCenter;
ConesToSurround=ImportCell.Surround.ConesToSurround;
NormalizedCenterKernel=ImportCell.Center.CenterKernel;
NormalizedSurroundKernel=ImportCell.Surround.SurroundKernel;

fullfieldCenter=ImportCell.Center.CenterKernelInField;
fullfieldSurround=ImportCell.Surround.SurroundKernelInField;
DoG=ImportCell.Global.DoG.ReceptiveFieldKernel;
LDoG=ImportCell.Global.DoG.LDoG;
MDoG=ImportCell.Global.DoG.MDoG;
SDoG=ImportCell.Global.DoG.SDoG;

CenterLCoordWeights=ImportCell.Center.LCones.EachWeight;
CenterMCoordWeights=ImportCell.Center.MCones.EachWeight;
CenterSCoordWeights=ImportCell.Center.SCones.EachWeight;

SurroundLCoordWeights=ImportCell.Surround.LCones.EachWeight;
SurroundMCoordWeights=ImportCell.Surround.MCones.EachWeight;
SurroundSCoordWeights=ImportCell.Surround.SCones.EachWeight;

cnL=ImportCell.CenterConeCount(1);
cnM=ImportCell.CenterConeCount(2);
cnS=ImportCell.CenterConeCount(3);

snL=ImportCell.SurroundConeCount(1);
snM=ImportCell.SurroundConeCount(2);
snS=ImportCell.SurroundConeCount(3);

%CD to plot directory
PlotDirectory=strcat('/Users/lauren/Documents/SCHOOL/Dacey/RF Model/2D Grid/Analysis/Cell Plots/',CellName);
if (exist(PlotDirectory) == 0)
    mkdir(PlotDirectory);
    cd(PlotDirectory);
else
    cd(PlotDirectory);
end

%% PLOT THE DATA

% figure;
% %Plot the retinal universe and how the 1-mm2 patch compares
% subplot(1,2,1);
% hold on;
% plot(pixel_coords(:,1),pixel_coords(:,2),'b.')
% axis([-250 250 -250 250]);
% axis square
% set(gca, 'color', [0.25 0.25 0.25]);
% rectangle('Position',[-r,-r,2*r,2*r],'EdgeColor','r');

%Plot the 1-mm2 patch and its L,M,S cones
%subplot(1,2,2);
figure('Name',strcat(CellName,' Retinal Mosaic'));
hold on;
xlabel('x (mm)');
ylabel('y (mm)'); 
title(strcat(CellName,' ',num2str(RetinaWeights(1)),' L,',num2str(RetinaWeights(2)),' M,',num2str(RetinaWeights(3)),' S'));

MapL=plot(LCoords(:,1),LCoords(:,2),'r.');
MapM=plot(MCoords(:,1),MCoords(:,2),'g.');
MapS=plot(SCoords(:,1),SCoords(:,2),'b.');
axis square
set(gca, 'color', [0.25 0.25 0.25]);
set(MapL,...
    'DisplayName','L Cones',...
    'MarkerFaceColor',[.64 .08 .18],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',3,...
    'Marker','o')
set(MapM,...
    'DisplayName','M Cones',...
    'MarkerFaceColor',[.47 .67 .19],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',3,...
    'Marker','o')
set(MapS,...
    'DisplayName','S Cones',...
    'MarkerFaceColor',[0 .45 .74],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',3,...
    'Marker','o')

saveas(gcf,strcat('Cell Retinal Mosaic'));

%Plot the weighted center and surround inputs for this particular cell
figure('Name',strcat(CellName,' Response Components'));
subplot(1,2,1);
scatter(SortedTotalMosaic(1:ConesToCenter,1),SortedTotalMosaic(1:ConesToCenter,2),30,NormalizedCenterKernel,'fill')
axis([-max(CenterX)*1.1 max(CenterX)*1.1 -max(CenterX)*1.1 max(CenterX)*1.1]);
axis square
xlabel('x (mm)');
ylabel('y (mm)'); 
title(strcat('Center (',num2str(ConesToCenter),' cones)'));
subplot(1,2,2);
scatter(SortedTotalMosaic(1:ConesToSurround,1),SortedTotalMosaic(1:ConesToSurround,2),30,NormalizedSurroundKernel,'fill')
axis([-max(SurroundX)*1.1 max(SurroundX)*1.1 -max(SurroundX)*1.1 max(SurroundX)*1.1]);
axis square
xlabel('x (mm)');
ylabel('y (mm)'); 
title(strcat('Surround (',num2str(ConesToSurround),' cones)'));

saveas(gcf,strcat('Cell Response Components'));

%Plot the entire patch to see the receptive field in context (center,
%surround, DoG)
figure('Name',strcat(CellName,' Full Receptive Field'));
subplot(2,3,1);
scatter(SortedTotalMosaic(:,1),SortedTotalMosaic(:,2),30,fullfieldCenter,'fill')
axis([-max(SurroundX)*2 max(SurroundX)*2 -max(SurroundX)*2 max(SurroundX)*2]);
axis square
xlabel('x (mm)');
ylabel('y (mm)'); 
title(strcat('Center (',num2str(cnL),'L/',num2str(cnM),'M/',num2str(cnS),'S)'));
subplot(2,3,2);
scatter(SortedTotalMosaic(:,1),SortedTotalMosaic(:,2),30,fullfieldSurround,'fill')
axis([-max(SurroundX)*2 max(SurroundX)*2 -max(SurroundX)*2 max(SurroundX)*2]);
axis square
xlabel('x (mm)');
ylabel('y (mm)'); 
title(strcat('Surround (',num2str(snL),'L/',num2str(snM),'M/',num2str(snS),'S)'));
subplot(2,3,3)
scatter(SortedTotalMosaic(:,1),SortedTotalMosaic(:,2),30,DoG,'fill')
axis([-max(SurroundX)*2 max(SurroundX)*2 -max(SurroundX)*2 max(SurroundX)*2]);
axis square
xlabel('x (mm)');
ylabel('y (mm)'); 
title('DoG');

%Plot the same thing, but visualize the weights of the cone classes
subplot(2,3,4);
hold on;
CLPlot=plot3(LCoords(:,1),LCoords(:,2),CenterLCoordWeights,'r.');
CMPlot=plot3(MCoords(:,1),MCoords(:,2),CenterMCoordWeights,'g.');
CSPlot=plot3(SCoords(:,1),SCoords(:,2),CenterSCoordWeights,'b.');
view([57 -6]);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('Normalized amplitude');
title(strcat(num2str(CenterWeights(1)),' L/',num2str(CenterWeights(2)),' M/',num2str(CenterWeights(3)),' S'));
set(CLPlot,...
    'DisplayName','Surround L',...
    'MarkerFaceColor',[.64 .08 .18],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',2,...
    'Marker','o')
set(CMPlot,...
    'DisplayName','Surround M',...
    'MarkerFaceColor',[.47 .67 .19],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',2,...
    'Marker','o')
set(CSPlot,...
    'DisplayName','Surround S',...
    'MarkerFaceColor',[0 .45 .74],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',2,...
    'Marker','o')

subplot(2,3,5);
hold on;
SLPlot=plot3(LCoords(:,1),LCoords(:,2),SurroundLCoordWeights,'ro');
SMPlot=plot3(MCoords(:,1),MCoords(:,2),SurroundMCoordWeights,'go');
SSPlot=plot3(SCoords(:,1),SCoords(:,2),SurroundSCoordWeights,'bo');
view([57 -6]);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('Normalized amplitude');
title(strcat(num2str(SurroundWeights(1)),' L/',num2str(SurroundWeights(2)),' M/',num2str(SurroundWeights(3)),' S'));
set(SLPlot,...
    'DisplayName','Surround L',...
    'MarkerFaceColor',[.64 .08 .18],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',2,...
    'Marker','o')
set(SMPlot,...
    'DisplayName','Surround M',...
    'MarkerFaceColor',[.47 .67 .19],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',2,...
    'Marker','o')
set(SSPlot,...
    'DisplayName','Surround S',...
    'MarkerFaceColor',[0 .45 .74],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',2,...
    'Marker','o')

subplot(2,3,6);
hold on;
DoGLPlot=plot3(LCoords(:,1),LCoords(:,2),LDoG,'r.');
DoGMPlot=plot3(MCoords(:,1),MCoords(:,2),MDoG,'g.');
DoGSPlot=plot3(SCoords(:,1),SCoords(:,2),SDoG,'b.');
view([57 -6]);
xlabel('x (mm)');
ylabel('y (mm)');
zlabel('Normalized amplitude');
title('DoG');
set(DoGLPlot,...
    'DisplayName','Surround L',...
    'MarkerFaceColor',[.64 .08 .18],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',2,...
    'Marker','o')
set(DoGMPlot,...
    'DisplayName','Surround M',...
    'MarkerFaceColor',[.47 .67 .19],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',2,...
    'Marker','o')
set(DoGSPlot,...
    'DisplayName','Surround S',...
    'MarkerFaceColor',[0 .45 .74],...
    'LineStyle','none',...
    'Color','none',...
    'MarkerSize',2,...
    'Marker','o')

saveas(gcf,strcat('Cell Full Receptive Field'));
