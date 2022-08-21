function plotHandles = plotRF(cellEquation,figHandle,lineStyle)

if nargin == 1
    figHandle = figure();
    lineStyle = '-r';
end
if nargin == 2
    figure(figHandle);
    lineStyle = '-r';
end
if nargin == 3
    figure(figHandle);
end

% Check if there is anything in the figure yet
newFigure = isempty(figHandle.Children);

% Use the inverse fourier transform to obtain the response in time after
% converting units from frequency to radians/sec (w)
syms w
cellEquationTime = ifourier(subs(cellEquation,w/(2*pi)));

myFreqs = logspace(log10(0.5),log10(100),100);
myTime = 0:0.001:0.2;

ttfComplex = eval(subs(cellEquation,myFreqs));
subplot(3,1,1)
if ~newFigure
    hold on
end

loglog(myFreqs,abs(ttfComplex),lineStyle);
ylim([1e-2 1e1]);
%    axHandle.YTick = [1e-3 1e-2 1e-1 1e0];
xlabel('frequency [Hz]'); ylabel('gain');

% The phase change by frequency, completing the Bode plot
subplot(3,1,2)
if ~newFigure
    hold on
end
semilogx(myFreqs,unwrap(angle(ttfComplex))*(180/pi),lineStyle);
%    ylim([-1000 100]);
%    axHandle = gca;
xlabel('frequency [Hz]'); ylabel('phase [deg]');

% The impulse response function in time
subplot(3,1,3)
irf = (eval(subs(cellEquationTime,myTime)));
if ~newFigure
    hold on
end
plot(myTime*1000,irf,lineStyle);
xlabel('time [msec]'); ylabel('response');

end
