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

% Use the inverse fourier transform to obtain the response in time after
% converting units from frequency to radians/sec (w)
syms w
cellEquationTime = ifourier(subs(cellEquation,w/(2*pi)));

% Define the support for the plots
myFreqs = logspace(log10(0.5),log10(100),100);
myTime = 0:0.001:0.2;


%% Panel 1 -- Gain by frequency
ttfComplex = eval(subs(cellEquation,myFreqs));
if any(whichPanel==1)
    subplot(3,1,1)
    if ~newFigure
        hold on
    end
    semilogx(myFreqs,abs(ttfComplex),lineStyle,'LineWidth',LineWidth);
    ylim([1e-2 1e1]);
    xlabel('frequency [Hz]'); ylabel('gain');
end


%% Panel 2 -- Phase by frequency
if any(whichPanel==2)
    subplot(3,1,2)
    if ~newFigure
        hold on
    end
    semilogx(myFreqs,unwrap(angle(ttfComplex))*(180/pi),lineStyle,'LineWidth',LineWidth);
    ylim([-1000 100]);
    xlabel('frequency [Hz]'); ylabel('phase [deg]');
end


%% Panel 3 -- Impulse response function
if any(whichPanel==3)
    subplot(3,1,3)
    irf = (eval(subs(cellEquationTime,myTime)));
    irf = irf ./ max(irf);
    if ~newFigure
        hold on
    end
    plot(myTime*1000,irf,lineStyle,'LineWidth',LineWidth);
    ylim([-1.1 1.1]);
    xlabel('time [msec]'); ylabel('relative response');
end

end
