SF='SF21';
T=500;
figure;hold on
for f=1:1000;
fft.d=fft(Data.(matlab.lang.makeValidName(SF)).Bootstrap.Histograms(f,:));
m = length(Data.(matlab.lang.makeValidName(SF)).Bootstrap.Histograms(f,:));
x = 0:0.01:max(Data.(matlab.lang.makeValidName(SF)).RawData.TimeBins*(2*pi/T));
fft.a0 = fft.d(1)/m;
fft.a1 = 2*real(fft.d(2))/m;
fft.b1 = -2*imag(fft.d(2))/m;
fft.F1fit = fft.a0 + fft.a1*cos(x) + fft.b1*sin(x);
fit=plot(x,fft.F1fit);
set(fit,...
'Color',[.75 .75 .75]);
% max=max(fft.F1fit);
% xmax=x(find(fft.F1fit==max(fft.F1fit)));
% plot(xmax,max);
clear fft m
end

Data.(matlab.lang.makeValidName(SF)).fft.d=fft(Data.(matlab.lang.makeValidName(SF)).RawData.Histogram);
m = length(Data.(matlab.lang.makeValidName(SF)).RawData.Histogram);
x = 0:0.01:max(Data.(matlab.lang.makeValidName(SF)).RawData.TimeBins*(2*pi/T));
Data.(matlab.lang.makeValidName(SF)).fft.a0 = Data.(matlab.lang.makeValidName(SF)).fft.d(1)/m;
Data.(matlab.lang.makeValidName(SF)).fft.a1 = 2*real(Data.(matlab.lang.makeValidName(SF)).fft.d(2))/m;
Data.(matlab.lang.makeValidName(SF)).fft.b1 = -2*imag(Data.(matlab.lang.makeValidName(SF)).fft.d(2))/m;
Data.(matlab.lang.makeValidName(SF)).fft.F1fit = Data.(matlab.lang.makeValidName(SF)).fft.a0 + Data.(matlab.lang.makeValidName(SF)).fft.a1*cos(x) + Data.(matlab.lang.makeValidName(SF)).fft.b1*sin(x);
% max=max(Data.SF1.fft.F1fit);
% xmax=x(find(Data.SF1.fft.F1fit==max(Data.SF1.fft.F1fit)));
emp=plot(x,Data.(matlab.lang.makeValidName(SF)).fft.F1fit);
set(emp,...
'LineWidth',1.5,...
'LineStyle','-',...
'Color','red');

