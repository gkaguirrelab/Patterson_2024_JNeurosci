
% Values obtained using the web plot digitized, applied to the 16 cell
% average from figure 9 of Lee 2007 JoV

midget = [0.6050875007272942, 207.30143890792442
    1.232846739442066, 163.0629189863448
    2.4088967285183034, 181.42087560797637
    4.806380863064389, 65.84651866079408
    9.792849742266032, 14.787139496981759
    19.539304046896113, 1.846754927214513
    25.118864315095795, 1];

parasol = [0.6165504968075404, 32.70727417075832
    1.232846739442066, 51.93751559698882
    2.465185075387883, 84.07873955784511
    4.929353553439602, 101.94553300259766
    9.856674331433043, 72.06768734806805
    20.14073128088476, 36.7158474278988
    53.366992312063125, 1];

finterp = logspace(log10(0.5),log10(100),100);

%% Midget
f = midget(:,1);
y = midget(:,2);
yfit = 10.^(spline(log10(f),log10(y),log10(finterp(finterp>1))));

% For the midget response, we do some hacking to force a plateau in the MTF
[a,b] = max(yfit);
yfit = [repmat(a,1,b+100-length(yfit)) yfit(b+1:end)];

% Force the MTF to have a floor above the CCF
a = find(yfit<1,1);
yfit(a:end)=1;

% Store the fit
midgetFit = yfit./max(yfit);

% plot
loglog(f,y./max(yfit),'or')
hold on
loglog(finterp,yfit./max(yfit),'-r')


%% Parasol
f = parasol(:,1);
y = parasol(:,2);
yfit = 10.^(spline(log10(f),log10(y),log10(finterp)));

% Force the MTF to have a floor above the CCF
a = find(yfit<1,1);
yfit(a:end)=1;

% Store the fit
parasolFit = yfit./max(yfit);

% plot
loglog(f,y./max(yfit),'ok')
hold on
loglog(finterp,yfit./max(yfit),'-k')
ylim([0 1.2])

%% Wool 2018
LminusMResponse = [0.539083557951483, 90.6832298136646
4.582210242587602, 69.25465838509317
12.331536388140162, 44.09937888198756
21.765498652291104, 28.881987577639748
38.67924528301887, 18.944099378881972
49.46091644204853, 16.149068322981364];

LplusMResponse = [0.6064690026954187, 65.21739130434781
4.851752021563343, 75.15527950310559
9.636118598382751, 82.91925465838509
14.285714285714288, 88.19875776397515
19.676549865229113, 92.2360248447205
29.31266846361186, 95.96273291925466
39.42048517520217, 97.51552795031056
49.595687331536396, 97.20496894409938];
% 
% semilogx(LminusMResponse(:,1),LminusMResponse(:,2),'-r');
% hold on
% semilogx(LplusMResponse(:,1),LplusMResponse(:,2),'-k');


% Make some mixes
for mRatio = 0:0.1:0.5
    myFit = parasolFit.*(1-mRatio) + midgetFit.*mRatio;
    finterp(find(myFit==max(myFit)))
end

