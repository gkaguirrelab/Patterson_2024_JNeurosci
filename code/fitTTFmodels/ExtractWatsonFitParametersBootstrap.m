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

options = optimoptions(@fmincon,... % The options for the search (mostly silence diagnostics)
    'Diagnostics','off',...
    'Display','off');

% Model guess
p0 = [1.5 0.015]; lb = [-inf 0.005]; ub = [8 0.025];
% p0 = [1.5, 0.8, 0.015]; lb = [0 0 0]; ub = [8 1 0.025];
% p0 = [1.5, 0.8, 0.015, 0.018, 10]; lb = [0 0 0 0 0]; ub = [8 1 0.025 0.03 20];

p = NaN*ones(2,1000);
Rsquared = NaN*ones(1,1000);

% Obtain Watson fit for bootstrapped values
for bb = 1:1000
    y = bootVals(:,bb)';
    
    % The objective function is the norm of the model fit error
    myObj = @(p) norm(y - watsonTTF2param(p,w));
        
    % Search
    p(:,bb) = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
    
    % Calculated R-squared and Adjusted R-squared for bootstrapped fits
    yFit = watsonTTF2param(p(:,bb),w);
    
    R = corrcoef(y,yFit);
    Rsquared(1,bb) = R(1,2)^2;
    
end

end
