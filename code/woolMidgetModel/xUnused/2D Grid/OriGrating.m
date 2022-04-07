theta = 90;                              % grating orientation
sigma = .05;                             % gaussian standard deviation in pixels
phase = 0;                            % phase (0 -> 1)
freq = 5;                    % cpd
phaseRad = (phase/360)*2*pi;
thetaRad=(theta/360)*2*pi;        % convert theta (orientation) to radians
X = linspace(-5,5,10001);                           % X is a vector from 1 to imageSize
[Xm, Ym] = meshgrid(X, X);             % 2D matrices
Xt = Xm*cos(thetaRad);                % compute proportion of Xm for given orientation
Yt = Ym*sin(thetaRad);                % compute proportion of Ym for given orientation
XYt=[Xt+Yt];                      % sum X and Y components
grating = sin(XYt*freq*2*pi+phaseRad);                   % make 2D sinewave
gauss=exp(-(((Xm.^2)+(Ym.^2))./(2*sigma^2))); % formula for 2D gaussian
gauss=exp(-((Xm-100).^2/2/sigma^2 + (Ym-100).^2/2/sigma^2));
gauss(gauss < .005) = 0;                 % trim around edges (for 8-bit colour displays)
gabor = grating .* gauss;                % use .* dot-product
gaborimage=imagesc( gabor, [-1 1] );                        % display
axis off; axis image;                    % use gray colormap
axis image; axis off; colormap gray(256);
set(gca,'pos', [0 0 1 1]);               % display nicely without borders
%set(gcf, 'menu', 'none', 'Color',[.5 .5 .5]); % without background

theta = 90;                              % grating orientation
phase = 0;                            % phase (0 -> 1)
freq = 1;                    % cpd
phaseRad = (phase/360)*2*pi;
thetaRad=(theta/360)*2*pi;        % convert theta (orientation) to radians
X = linspace(-5,5,10001);                           % X is a vector from 1 to imageSize
[Xm, Ym] = meshgrid(X, X);             % 2D matrices
Xt = Xm*cos(thetaRad);                % compute proportion of Xm for given orientation
Yt = Ym*sin(thetaRad);                % compute proportion of Ym for given orientation
XYt=[Xt+Yt];                      % sum X and Y components
grating = sin(XYt*freq*2*pi+phaseRad);                   % make 2D sinewave
axis off; axis image;                    % use gray colormap
axis image; axis off; colormap gray(256);
set(gca,'pos', [0 0 1 1]);               % display nicely without borders
%set(gcf, 'menu', 'none', 'Color',[.5 .5 .5]); % without background
imagesc(grating);
