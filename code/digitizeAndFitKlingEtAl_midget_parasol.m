
% Determine TSFs for midget and parasol cells based on impulse response function data from Kling et al (2020) Biorxiv

rON_parasol_x = [238.96,218.18,190.91,172.73,154.55,138.96,129.87,116.88,107.79,97.40,89.61,...
    93.51,84.42,79.22,75.32,71.43,68.83,66.23,63.64,59.74,55.84,51.95,49.35,...
    45.45,42.86,40.26,36.36,33.77,31.17,27.27,23.38,18.18,11.69,3.90];

rON_parasol_y = [0.00,0.00,-0.01,-0.02,-0.04,-0.07,-0.10,-0.17,-0.21,-0.28,-0.32,-0.31,...
    -0.30,-0.24,-0.15,-0.03,0.08,0.19,0.33,0.44,0.58,0.67,0.69,0.67,0.61,0.52,0.44,...
    0.33,0.26,0.18,0.10,0.05,0.01,0.01];

rOFF_parasol_x = [241.56,224.68,201.30,183.12,158.44,142.86,120.78,109.09,96.10,83.12,...
    76.62,68.83,63.64,53.25,49.35,46.75,41.56,37.66,31.17,27.27,22.08,16.88,...
    10.39,3.90];

rOFF_parasol_y = [0.02,0.02,0.03,0.05,0.08,0.11,0.16,0.19,0.21,0.19,0.14,0.05,-0.10,...
    -0.44,-0.58,-0.64,-0.66,-0.63,-0.52,-0.38,-0.22,-0.09,0.00,0.01];

rON_midget_x = [240.09,216.81,188.36,166.38,146.98,128.88,121.12,112.07,100.43,91.38,...
    84.91,78.45,73.28,62.93,57.76,53.88,50.00,47.41,42.24,38.36,33.19,22.84,...
    15.09,9.91,3.45];

rON_midget_y = [-0.02,-0.04,-0.05,-0.08,-0.12,-0.13,-0.13,-0.11,-0.06,0.00,0.09,0.18,...
    0.29,0.52,0.63,0.66,0.67,0.64,0.54,0.42,0.28,0.09,0.02,0.01,0.02];

rOFF_midget_x = [238.79,193.53,159.91,136.64,118.53,104.31,92.67,84.91,77.16,70.69,...
    62.93,56.47,52.59,47.41,40.95,34.48,29.31,21.55,16.38,8.62,0.86];

rOFF_midget_y = [0.03,0.03,0.02,-0.01,-0.03,-0.07,-0.13,-0.19,-0.28,-0.37,-0.48,-0.58,...
    -0.61,-0.63,-0.59,-0.49,-0.38,-0.19,-0.09,-0.03,-0.03];


% get interpolated fits

x_interp = fliplr([0:1:200]);
ON_parasol_y = spline(rON_parasol_x,rON_parasol_y,x_interp);
OFF_parasol_y = spline(rOFF_parasol_x,rOFF_parasol_y,x_interp);
ON_midget_y = spline(rON_midget_x,rON_midget_y,x_interp);
OFF_midget_y = spline(rOFF_midget_x,rOFF_midget_y,x_interp);

% plot interpolated values
subplot(1,4,1)
plot(x_interp.*-1,ON_parasol_y,'--k')
title('ON parasol')

subplot(1,4,2)
plot(x_interp.*-1,OFF_parasol_y,'--k')
title('OFF parasol')

subplot(1,4,3)
plot(x_interp.*-1,ON_midget_y,'--k')
title('ON midget')

subplot(1,4,4)
plot(x_interp.*-1,OFF_midget_y,'--k')
title('OFF midget')

% Send frequencies through parasol and midget filters

w = [0.5 1 2 4 6 8 12 16 24 32 48 64];
t = 0:0.001:2.4;
PP = 200; %size of parsing bins, must be 200 or less

for x=1:length(w)
    stimulus = sin(2*pi*w(x)*t);
    
    parse = 1:PP:1001;
    
    for y = 1:length(parse)-1
        ONp = ON_parasol_y(1:PP+1).*stimulus(parse(y):parse(y)+PP);
        OFFp = OFF_parasol_y(1:PP+1).*stimulus(parse(y):parse(y)+PP);
        temp = sum(ONp + OFFp);
        temp(temp<0) = 0;
        parasol(x,y) = temp;
        
        ONm = ON_midget_y(1:PP+1).*stimulus(parse(y):parse(y)+PP);
        OFFm = OFF_midget_y(1:PP+1).*stimulus(parse(y):parse(y)+PP);
        temp = sum(ONm + OFFm);
        temp(temp<0) = 0;
        midget(x,y) = temp;
    end
end

% get rectified sum of ON and OFF for parasol and midget cells

parasol = sum(parasol,2);
parasol = parasol./max(parasol);
midget = sum(midget,2);
midget = midget./max(midget);

% Determine and Plot TSFs
[wFit,yFit_parasol,yFit2_parasol,p1] = fitExp(w,parasol');
[~,yFit_midget,yFit2_midget,p2] = fitExp(w,midget');

figure
subplot(1,2,1)
hold on
plot(w,parasol,'.k','MarkerSize',12)
plot(wFit,yFit_parasol,'-k')
title('parasol')
ax = gca; ax.Box = 'off'; ax.XScale = 'log'; ax.TickDir = 'out';
xlabel(sprintf('%1.2g, %1.2g, %1.2g',p1))

subplot(1,2,2)
hold on
plot(w,midget,'.k','MarkerSize',12)
plot(wFit,yFit_midget,'-k')
title('midget')
ax = gca; ax.Box = 'off'; ax.XScale = 'log'; ax.TickDir = 'out';
xlabel(sprintf('%1.2g, %1.2g, %1.2g',p2))


%% Local functions
function [wFit,yFit,yFit2,p] = fitExp(w,Y)
            
            % TTF model guess
            p0 = [0.5 4 1];
            lb = [0 0 0]; 
            ub = [10 10 2];
            
            
            wDelta = min(diff(log10(w))); % Create a scaled-up, log-spaced, version of the frequency domain
            upScale = 500;
            wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);
    
            myObj=@(p)sqrt(sum((Y-watsonTemporalModelvep(w,p)).^2));
            
            p = fmincon(myObj,p0,[],[],[],[],lb,ub);
            
            % calculate model fit and undo scaling
            yFit = watsonTemporalModelvep(wFit,p);
            
            yFit(~isfinite(yFit))=nan;
            
            yFit2 = watsonTemporalModelvep(w,p);
end