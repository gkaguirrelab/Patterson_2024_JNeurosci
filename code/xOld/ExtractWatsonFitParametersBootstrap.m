%% ExtractTTFdataWatsonModelFit
%
% This function calculates Bootstrapped Watson Model fit parameters

function [p,Rsquared] = ExtractWatsonFitParametersBootstrap(w,yVal,varargin)

q = inputParser;
q.addParameter('p0',[1 0.9 0.01],@isnumeric);
q.addParameter('lb',[-10 -2 0],@isnumeric);
q.addParameter('ub',[10 2 0.05],@isnumeric);
q.addParameter('SurroundAmp',0.9,@isnumeric); % based on average Gs across subjects and directions V1 wide field
q.parse(varargin{:});

% Output:
% p = parameter values across 1000 bootstrapped samples
% Rsqaured = R-squared (1st row) and adjusted R-squared (2nd row) value for goodness of fit 

    
% Adjust the values for the zero frequency, obtain bootstrapped means and 95% CI
y0 = yVal{1}; % response at frequency = 0
yW = yVal(:,2:end); % responses to all other frequencies
p0 = q.Results.p0; lb = q.Results.lb; ub = q.Results.ub;

bootVals = NaN*ones(length(w),1000);
for ff = 1:length(w)
    temp_data = yW{ff}-y0;
    bootVals(ff,:) = sort(bootstrp(1000,@median,temp_data));
end

options = optimoptions(@fmincon,... % The options for the search (mostly silence diagnostics)
    'Diagnostics','off',...
    'Display','off');

p = NaN*ones(length(p0),1000);
Rsquared = NaN*ones(1,1000);

% Obtain Watson fit for bootstrapped values
for bb = 1:1000
    y = bootVals(:,bb)';
    
    % shift so that if min(y)=0
    y = y - min(y);

    myObj = @(p) norm(y - watsonTTF(p,w));
    
        
    % Search
    p(:,bb) = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
    
    yFit = watsonTTF(squeeze(p(:,bb)),w);
    
    R = corrcoef(y,yFit);
    Rsquared(1,bb) = R(1,2)^2;
    
%     % To plot individual bootstrapped values
%     wDelta = min(diff(log10(w))); % Create a scaled-up, log-spaced, version of the frequency domain
%     upScale = 10;
%     wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);
%     yFit2 = watsonTTF(squeeze(p(:,bb)),wFit);
%     figure(276)
%     clf
%     hold on
%     plot(w,y,'ok')
%     plot(wFit,yFit2,'-b','LineWidth',1)
%     str = {sprintf('[%2.1f, %2.2f, %2.3f]',p(1:3))};
%     str2 = {sprintf('Rsquared = %2.2f',Rsquared(1,bb))};
%     title(str)
%     xlabel(str2)
%     ylabel(num2str(bb))
%     ax = gca; ax.TickDir = 'out'; ax.Box = 'off'; ax.XScale = 'log';
%     pause
    
end
   
    
end
