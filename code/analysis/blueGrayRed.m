function [c] = blueGrayRed(m)

%   Creates a matrix useful for plotting pRF and ccRF correlation maps
%
%   Usage:
%   [mycolormap] = make_polar_colormap(show)
%
%   defaults:
%   show = 0; do not plot the resulting polar angle colormap
%
%   Written by Andrew S Bock Oct 2014

if nargin < 1, m = size(get(gcf,'colormap'),1); end


%% Create colormap
% Blue to Blue-Gray 
c(0*m/4+1:1*m/4,1) = linspace(0,0.375,m/4);
c(0*m/4+1:1*m/4,2) = linspace(0,0.375,m/4);
c(0*m/4+1:1*m/4,3) = linspace(1,0.875,m/4);
% Blue-Gray to Gray 
c(1*m/4+1:2*m/4,1) = linspace(0.375,0.75,m/4);
c(1*m/4+1:2*m/4,2) = linspace(0.375,0.75,m/4);
c(1*m/4+1:2*m/4,3) = linspace(0.875,0.75,m/4);
% Gray to Red-Gray 
c(2*m/4+1:3*m/4,1) = linspace(0.75,0.875,m/4);
c(2*m/4+1:3*m/4,2) = linspace(0.75,0.375,m/4);
c(2*m/4+1:3*m/4,3) = linspace(0.75,0.375,m/4);
% Red-Gray to Red
c(3*m/4+1:4*m/4,1) = linspace(0.875,1,m/4);
c(3*m/4+1:4*m/4,2) = linspace(0.375,0,m/4);
c(3*m/4+1:4*m/4,3) = linspace(0.375,0,m/4);

end