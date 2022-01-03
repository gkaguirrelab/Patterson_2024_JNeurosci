x = 0:0.001:12-0.001;
f = [4,2,8];
hc = (-cos(x*2*pi*(1/3))+1)/2;
hc(and(x>1.5,x<10.5))=1;

d = [];
for ii = 1:length(f)
y = sin(x*2*pi*f(ii));
d = [d y.*hc];
end

% Add the attention event
x = 0:0.001:36-0.001;
d(and(x>30,x<30.250))=nan;

figure
plot(x,d)
box off
axis off

% Save the figure

resultsSaveDir = fullfile(localSaveDir,'Fig 1 - experimental design');
mkdir(resultsSaveDir);

fileName = fullfile(resultsSaveDir,'flickerBlocks.pdf');
print(fileName,'-dpdf');
