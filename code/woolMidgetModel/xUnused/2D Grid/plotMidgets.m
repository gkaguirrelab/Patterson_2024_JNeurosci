function plotMidgets(rig, filename1, filename2)

cd('/Users/lauren/Desktop/All fits')
try
    file1 = strcat(rig,'HG',filename1,'.LW.','spatial.sf.fit.mat');
    file2 = strcat(rig,'HG',filename2,'.LW.','spatial.sf.fit.mat');
    fitData1 = load(file1);
    fitData2 = load(file2);
catch
    try
        file1 = strcat(rig,'HG',filename1,'.LW.FREE.','spatial.sf.fit.mat');
        file2 = strcat(rig,'HG',filename2,'.LW.FREE.','spatial.sf.fit.mat');
        fitData1 = load(file1);
    fitData2 = load(file2);
    catch
        file1 = strcat(rig,'HG',filename1,'.spatial.sf.fit.mat');
        file2 = strcat(rig,'HG',filename2,'.spatial.sf.fit.mat');
        fitData1 = load(file1);
        fitData2 = load(file2);
    end
end
    
xMin = min([min(fitData1.handles.plotIndVar) min(fitData2.handles.plotIndVar) min(fitData1.handles.indVar) min(fitData2.handles.indVar)]);
ampMin = min([min(fitData1.handles.plotAmps) min(fitData2.handles.plotAmps) min(fitData1.handles.Amps) min(fitData2.handles.Amps)]);
phaseMin = min([min(fitData1.handles.plotPhases) min(fitData2.handles.plotPhases) min(fitData1.handles.Phases) min(fitData2.handles.Phases)]);
xMax = 1.1 * max([max(fitData1.handles.plotIndVar) max(fitData2.handles.plotIndVar) max(fitData1.handles.indVar) max(fitData2.handles.indVar)]);
ampMax = 1.1 * max([max(fitData1.handles.plotAmps) max(fitData2.handles.plotAmps) max(fitData1.handles.Amps) max(fitData2.handles.Amps)]);
phaseMax = 1.1 * max([max(fitData1.handles.plotPhases) max(fitData2.handles.plotPhases) max(fitData1.handles.Phases) max(fitData2.handles.Phases)]);

subplot(2,1,1);
fit1Amps = loglog(fitData1.handles.plotIndVar,fitData1.handles.plotAmps);
hold on;
plot1Amps = loglog(fitData1.handles.indVar,fitData1.handles.Amps);
fit2Amps = loglog(fitData2.handles.plotIndVar,fitData2.handles.plotAmps);
plot2Amps = loglog(fitData2.handles.indVar,fitData2.handles.Amps);
axis([xMin xMax 1 ampMax])

subplot(2,1,2);
fit1Phases = semilogx(fitData1.handles.plotIndVar,fitData1.handles.plotPhases);
hold on;
plot1Phases = semilogx(fitData1.handles.indVar,fitData1.handles.Phases);
fit2Phases = semilogx(fitData2.handles.plotIndVar,fitData2.handles.plotPhases);
plot2Phases = semilogx(fitData2.handles.indVar,fitData2.handles.Phases);
axis([xMin xMax phaseMin phaseMax]);

% set(fit1Amps,fit2Amps,fit1Phases,fit2Phases,...
%     'Color','k',...
%     'LineStyle','-')
% 
% set(plot1Amps,plot2Amps,plot1Phases,plot2Phases,...
%     'MarkerSize',6,...
%     'LineStyle','none')
