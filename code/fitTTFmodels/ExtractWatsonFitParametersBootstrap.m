%% ExtractTTFdataWatsonModelFit
%
% This function calculates Bootstrapped Watson Model fit parameters and
% plots them

function [p] = ExtractWatsonFitParametersBootstrap(w,yVal)

        
% Adjust the values for the zero frequency, obtain bootstrapped means and 95% CI
y0 = yVal{1}; % response at frequency = 0
yW = yVal(:,2:end); % responses to all other frequencies

bootVals = NaN*ones(length(w),1000);
for ff = 1:length(w)
    temp_data = yW{ff}-y0;
    bootVals(ff,:) = sort(bootstrp(1000,@mean,temp_data));
end

% Create a scaled-up, log-spaced, version of the frequency domain
wDelta = min(diff(log10(w))); 
upScale = 10;
wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);

options = optimoptions(@fmincon,... % The options for the search (mostly silence diagnostics)
    'Diagnostics','off',...
    'Display','off');

% Model guess
p0 = [1.5, 0.8, 0.015]; lb = [0 0 0]; ub = [8 1 0.025];

p = NaN*ones(3,1000); 

% Obtain Watson fit for bootstrapped values
for bb = 1:1000
    y = bootVals(:,bb)';
    
    % The objective function is the norm of the model fit error
    myObj = @(p) norm(y - watsonTTF(p,w));
        
    % Search
    p(:,bb) = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
    
    % plot the fit to check it
    figure(500)
    hold on
    yFit = watsonTTF(p(:,bb),wFit);
    plot(w,y,'ok')
    plot(wFit,yFit,'-b','LineWidth',1)
    ylabel('BOLD % change');
    xlabel('frequency [Hz]');
    semilogx([1 64],[0 0],':k','LineWidth',1)
    ylim([-2 8]);
    xlim([1 128])
    set(gca,'TickDir','out');
    box off
    
end

% plot bootstrapped parameters
p = sort(p,2);

figure
subplot(1,3,1)
hold on
errorbar(1,p(1,500),p(1,500)-p(1,25),p(1,975)-p(1,500),'ob','MarkerFaceColor','b','LineWidth',1.5)
set(gca,'TickDir','out');
box off
xticks(1)
xlim([0 2])
ylim([0 8])
xlabel('gain')
ylabel('fit parameter value')

subplot(1,3,2)
hold on
errorbar(1,p(2,500),p(2,500)-p(2,25),p(2,975)-p(2,500),'ob','MarkerFaceColor','b','LineWidth',1.5)
set(gca,'TickDir','out');
box off
xticks(1)
xlim([0 2])
ylim([0.5 1.2])
xlabel('surround gain')
ylabel('fit parameter value')

subplot(1,3,3)
hold on
errorbar(1,p(3,500),p(3,500)-p(3,25),p(3,975)-p(3,500),'ob','MarkerFaceColor','b','LineWidth',1.5)
set(gca,'TickDir','out');
box off
xticks(1)
xlim([0 2])
ylim([0 0.03])
xlabel('time constant')
ylabel('fit parameter value')

end
