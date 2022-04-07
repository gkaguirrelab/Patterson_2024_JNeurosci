
Cell=GenerateCell(10);
figure
k=Cell.Surround.LCones.Coords;
plot(k(:,1),k(:,2),'or','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.8196    0.1569    0.1569])
hold on
k=Cell.Surround.MCones.Coords;
plot(k(:,1),k(:,2),'og','MarkerSize',5,'MarkerEdgeColor','k','MarkerFaceColor',[0.4667    0.6745    0.1882])
xlim([-0.5 0.5]);
ylim([-0.5 0.5]);
axis equal
set(gca,'Color',[0.5 0.5 0.5])
centerCoords = mean(Cell.Center.CenterConeCoords);
viscircles(centerCoords(1:2),0.0625,'Color','y');
viscircles(centerCoords(1:2),0.375,'Color','y');
box off
h = gca; h.XAxis.Visible = 'off'; h.YAxis.Visible = 'off';



Cell=GenerateCell(2);
figure
k=Cell.Surround.LCones.Coords;
plot(k(:,1),k(:,2),'or','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.8196    0.1569    0.1569])
hold on
k=Cell.Surround.MCones.Coords;
plot(k(:,1),k(:,2),'og','MarkerSize',20,'MarkerEdgeColor','k','MarkerFaceColor',[0.4667    0.6745    0.1882])
xlim([-0.05 0.05]);
ylim([-0.05 0.05]);
axis equal
set(gca,'Color',[0.5 0.5 0.5])
centerCoords = mean(Cell.Center.CenterConeCoords);
viscircles(centerCoords(1:2),0.04,'Color','y');
viscircles(centerCoords(1:2),0.01,'Color','y');


Cell=GenerateCell(0.15);
figure
k=Cell.Surround.LCones.Coords;
plot(k(:,1),k(:,2),'or','MarkerSize',27,'MarkerEdgeColor','k','MarkerFaceColor',[0.8196    0.1569    0.1569])
hold on
k=Cell.Surround.MCones.Coords;
plot(k(:,1),k(:,2),'og','MarkerSize',27,'MarkerEdgeColor','k','MarkerFaceColor',[0.4667    0.6745    0.1882])
xlim([-0.02 0.02]);
ylim([-0.02 0.02]);
axis equal
set(gca,'Color',[0.5 0.5 0.5])
centerCoords = Cell.Center.CenterConeCoords;
viscircles(centerCoords(1:2),0.0025,'Color','y');
viscircles(centerCoords(1:2),0.0125,'Color','y');