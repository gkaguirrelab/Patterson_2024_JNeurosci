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

f = parasol(:,1);
y = parasol(:,2);
finterp = logspace(log10(0.5),log10(100),100);
ymax = max(y);
y = y / ymax;

% Set up the p0 guess, and the bounds on the params
p0 = [1 8 1 1];
lb = [0 0 0 1];
ub = [100 8 100 100];

% The options for the search (mostly silence diagnostics)
options = optimoptions(@fmincon,...
    'Diagnostics','off',...
    'Display','off');

% The objective function is the norm of the model fit error
myObj = @(p) norm(y - betaTTF(p,f));

% The non-linear constraint
myNonlcon = @(p) betaTTF(p,f,y);

% Search
p = fmincon(myObj,p0,[],[],[],[],lb,ub,myNonlcon,options);

yFit = betaTTF(p,finterp);
y = y * ymax;
yFit = yFit * ymax;
loglog(f,y,'or')
hold on
loglog(finterp,yFit,'-r')



ylim([1 500])