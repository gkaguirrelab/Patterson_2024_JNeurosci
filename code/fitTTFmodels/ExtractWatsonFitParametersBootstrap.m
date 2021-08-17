%% ExtractTTFdataWatsonModelFit
%
% This function calculates Bootstrapped Watson Model fit parameters

function [p,Rsquared] = ExtractWatsonFitParametersBootstrap(w,yVal)

% Output:
% p = parameter values across 1000 bootstrapped samples
% Rsqaured = R-squared (1st row) and adjusted R-squared (2nd row) value for goodness of fit 
        
% Adjust the values for the zero frequency, obtain bootstrapped means and 95% CI
y0 = yVal{1}; % response at frequency = 0
yW = yVal(:,2:end); % responses to all other frequencies

bootVals = NaN*ones(length(w),1000);
for ff = 1:length(w)
    temp_data = yW{ff}-y0;
    bootVals(ff,:) = sort(bootstrp(1000,@median,temp_data));
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
Rsquared = NaN*ones(2,1000);

% Obtain Watson fit for bootstrapped values
for bb = 1:1000
    y = bootVals(:,bb)';
    
    % The objective function is the norm of the model fit error
    myObj = @(p) norm(y - watsonTTF(p,w));
        
    % Search
    p(:,bb) = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
    
    % Calculated R-squared and Adjusted R-squared for bootstrapped fits
    yFit = watsonTTF(p(:,bb),w);
    SST = sum(abs(diff([mean(y)*ones(1,6);y]))).^2;
    SSE = sum(abs(diff([y;yFit]))).^2;
    
    Rsquared(1,bb) = 1-(SSE./SST);
    
    n=6;
    param=length(p0);
    Rsquared(2,bb) = 1-(((n-1)/(n-param)).*(SSE./SST));
    
end

end
