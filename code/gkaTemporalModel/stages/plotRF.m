function plotRF(cellEquation,figHandle,lineStyle,whichPanel,LineWidth)

% Handle arguments
if nargin == 1
    figHandle = figure();
    lineStyle = '-r';
    whichPanel = [1 2 3];
    LineWidth = 1;
end
if nargin == 2
    figure(figHandle);
    lineStyle = '-r';
    whichPanel = [1 2 3];
    LineWidth = 1;
end
if nargin == 3
    figure(figHandle);
    whichPanel = [1 2 3];
    LineWidth = 1;
end
if nargin == 4
    figure(figHandle);
    LineWidth = 1;
end

% Check if there is anything in the figure yet
newFigure = isempty(figHandle.Children);

% Define the support for the plots
myFreqs = logspace(log10(0.5),log10(100),101);
myTime = 0:0.001:0.1;

%% Panel 1 -- Gain by frequency
if iscell(cellEquation)
    gainVals = zeros(size(myFreqs)); angleVals = zeros(size(myFreqs));
    for ii=1:length(cellEquation)
        ttfComplex = double(subs(cellEquation{ii},myFreqs));
        gainVals = gainVals + abs(ttfComplex).*(1/length(cellEquation));
        angleVals = angleVals + unwrap(angle(ttfComplex)).*(1/length(cellEquation));
    end
else
    ttfComplex = double(subs(cellEquation,myFreqs));
    gainVals = abs(ttfComplex);
    angleVals = angle(ttfComplex);
end

if any(whichPanel==1)
    subplot(3,1,1)
    if ~newFigure
        hold on
    end
    semilogx(myFreqs,gainVals,lineStyle,'LineWidth',LineWidth);
 %   ylim([1e-2 1e1]);
    xlabel('frequency [Hz]'); ylabel('gain');
end


%% Panel 2 -- Phase by frequency
if any(whichPanel==2)
    subplot(3,1,2)
    if ~newFigure
        hold on
    end
    semilogx(myFreqs,unwrap(angleVals)*(180/pi),lineStyle,'LineWidth',LineWidth);
%    ylim([-1000 100]);
    xlabel('frequency [Hz]'); ylabel('phase [deg]');
end


%% Panel 3 -- Impulse response function
if any(whichPanel==3) && ~any(isnan(ttfComplex))
    subplot(3,1,3)

    % Use the inverse fourier transform to obtain the response in time
    % after converting units from frequency to radians/sec (w)
    syms f w x
    if iscell(cellEquation)
        irf = zeros(size(myTime));
        for ii=1:length(cellEquation)
            cellEquationTime = ifourier(subs(cellEquation{ii},f,w/(2*pi)));
            irf = irf + (double(subs(cellEquationTime,myTime)));
        end
    else
        cellEquationTime = ifourier(subs(cellEquation,f,w/(2*pi)),w,x);
        irf = double(subs(cellEquationTime,x,myTime));
    end

    % Scale the IRF to unit amplitude
    irf = irf ./ max(irf);
    if ~newFigure
        hold on
    end
    plot(myTime*1000,irf,lineStyle,'LineWidth',LineWidth);
    ylim([-1.1 1.1]);
    xlabel('time [msec]'); ylabel('relative response');
end

end
