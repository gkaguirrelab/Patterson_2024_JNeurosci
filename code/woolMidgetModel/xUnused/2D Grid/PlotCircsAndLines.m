figure;

for s=1:length(SpatFreqs)
    subplot(2,6,s)
    SFLabel=strcat('SF',num2str(SpatFreqs(s)),'cpd');
    %Read the appropriate values from the RF FFTs
    PlotLOriAmp=semilogy(OriCircles.(matlab.lang.makeValidName(SFLabel)).Oris,TotalLData.(matlab.lang.makeValidName(SFLabel)).Amplitude,'ro');
    hold on;
    PlotMOriAmp=semilogy(OriCircles.(matlab.lang.makeValidName(SFLabel)).Oris,TotalMData.(matlab.lang.makeValidName(SFLabel)).Amplitude,'go');    
    axis([0 2*pi 0.1 60]);
    axis square;
    set(gca,'TickDir','in','TickLength', [.005 .005]);box off
    title(SFLabel);
    set(PlotLOriAmp,...
        'DisplayName','L',...
        'LineWidth',.5,...
        'LineStyle',':',...
        'Color','k',...
        'MarkerFaceColor',[204 0 0]/255)
    set(PlotMOriAmp,...
        'DisplayName','L',...
        'LineWidth',.5,...
        'LineStyle',':',...
        'Color','k',...
        'MarkerFaceColor',[119 172 48]/255)
end

figure;
for s=1:length(SpatFreqs)
    subplot(2,6,s)
    SFLabel=strcat('SF',num2str(SpatFreqs(s)),'cpd');
    %Read the appropriate values from the RF FFTs
    PlotLOriPhase=plot(OriCircles.(matlab.lang.makeValidName(SFLabel)).Oris,abs(TotalLData.(matlab.lang.makeValidName(SFLabel)).Phase),'ro');
    hold on;
    PlotMOriPhase=plot(OriCircles.(matlab.lang.makeValidName(SFLabel)).Oris,abs(TotalMData.(matlab.lang.makeValidName(SFLabel)).Phase),'go');    
    axis square;
    axis([0 2*pi 0 180]);
    set(gca,'TickDir','in','TickLength', [.005 .005]);box off
    title(SFLabel);
    set(PlotLOriPhase,...
        'DisplayName','L',...
        'LineWidth',.5,...
        'LineStyle',':',...
        'Color','k',...
        'MarkerFaceColor',[204 0 0]/255)
    set(PlotMOriPhase,...
        'DisplayName','L',...
        'LineWidth',.5,...
        'LineStyle',':',...
        'Color','k',...
        'MarkerFaceColor',[119 172 48]/255)
end