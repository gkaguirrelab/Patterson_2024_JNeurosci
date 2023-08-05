% Colestone 1994, Figure 2

% Place to save figures
savePath = '~/Desktop/FlickerGrant/';

% Set up the figure
f1 = figure();
figuresize(200,200,'pt');
set(gcf,'color','w');
t = tiledlayout(1,2);
t.TileSpacing = 'tight';
t.Padding = 'none';


migraine = [
3.2753433756122483, 0.16684733770685123
5.325474261480487, 0.2872808061290815
8.298554750073702, 0.3055913306315554
10.331991915977309, 0.2687031821071737
12.450140831063102, 0.23013855787570636
15.24367503365383, 0.1564048830906824
20.159265192173127, 0.0809265369765898
25.498946874900103, 0.0020867149955778874
30.16182252799749, -0.05497313021456762
50.01758168379253, -0.13429884600220932];

control = [
3.282447086235493, 0.23379270662030316
5.249819743342936, 0.37431262720081837
8.399960219220507, 0.46123647187108185
10.439968317450617, 0.4862727895916432
12.390824847359013, 0.47114472744838265
15.277950437410983, 0.479416288098088
20.185548921479132, 0.32862440195636194
25.513331888912173, 0.13765108704531814
30.089009494109245, 0.05883683842255005
50.1117058495505, -0.04727270789897109];

migraine(:,2) = migraine(:,2) ./ max(migraine(:,2));
control(:,2) = control(:,2) ./ max(control(:,2));
W = ones(size(migraine(:,1)));
    p0A = [1.5 5 1.1 1.5];
    interpFreqs = logspace(log10(3),log10(50),101);

semilogx(migraine(:,1),migraine(:,2),'s',...
    'MarkerEdgeColor','none','MarkerFaceColor',[0.75 0.25 0.75],...
    'MarkerSize',6);

hold on
[~,~,~,yFitInterp] = fitWatsonModel(migraine(:,2)+0.5,W,migraine(:,1),p0A,interpFreqs);
semilogx(interpFreqs,yFitInterp-0.5,':k','LineWidth',2);

semilogx(control(:,1),control(:,2),'^',...
    'MarkerEdgeColor','none','MarkerFaceColor',[0.75 0.75 0.75],...
    'MarkerSize',6);
[~,~,~,yFitInterp] = fitWatsonModel(control(:,2)+0.5,W,control(:,1),p0A,interpFreqs);
semilogx(interpFreqs,yFitInterp-0.5,'-k','LineWidth',2);

    a = gca;
    a.TickDir = 'out';
    xlim([2 50]);
a.XTick = ([2 10 50]);
a.XTickLabel = {'2','10','50'};
xlabel('Frequency [Hz]')
ylabel('log relative sensitivity')
box off
title({'Detection threshold','(Coleston 1994)'})

plotNamesPDF = 'Figure X -- Coleston1994Detection.pdf';
saveas(f1,fullfile(savePath,plotNamesPDF));

