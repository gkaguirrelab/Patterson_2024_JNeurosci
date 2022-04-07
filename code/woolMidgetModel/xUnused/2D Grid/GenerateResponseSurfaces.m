function Cell=GenerateResponseSurfaces(Em)
%% IMPORT DATA
Cell=GenerateCell(Em);

cpAll=Cell.RetinaWeights;
AllCoords=Cell.Global.AllConeCoords(:,1:2);
LCoords=Cell.Global.LConeCoords;
MCoords=Cell.Global.MConeCoords;
SCoords=Cell.Global.SConeCoords;
LDoG=Cell.Global.DoG.LDoG;
MDoG=Cell.Global.DoG.MDoG;
SDoG=Cell.Global.DoG.SDoG;
TotalDoG=LDoG+MDoG;

%% TIDY UP THE VECTORS
%Clean out the NaNs from the coordinate arrays
LCoords(isnan(LCoords(:,1)),:) = [];
MCoords(isnan(MCoords(:,1)),:) = [];
SCoords(isnan(SCoords(:,1)),:) = [];

%Fix up the DoGs and get rid of extra zeros
[~,v]=max(cumsum(LDoG'~=0));
for h=1:v
    LDoG(LDoG==0)=[];
end
LDoG=padarray(LDoG,length(LCoords)-length(LDoG),'post');
[~,v]=max(cumsum(MDoG'~=0));
for h=1:v
    MDoG(MDoG==0)=[];
end
MDoG=padarray(MDoG,length(MCoords)-length(MDoG),'post');

if cpAll(3)~=0
    SCoords=(SCoords*1000)/200;
    SCoords(isnan(SCoords(:,1)),:) = [];
    for h=1:v
        SDoG(SDoG==0)=[];
    end
    SDoG=padarray(SDoG,length(SCoords)-length(SDoG),'post');
end
%% GENERATE RESPONSE SURFACES
%Assign a meshgrid to interpolate the data over
X=linspace(-.5,.5,1001); %Work within the 1mm2 patch that we've been using throughout
[xq,yq] = meshgrid(X,X);
yq=flipud(yq);

LResponseSurface=griddata(LCoords(:,1),LCoords(:,2),LDoG,xq,yq);
MResponseSurface=griddata(MCoords(:,1),MCoords(:,2),MDoG,xq,yq);
if cpAll(3)~=0
   SResponseSurface=griddata(SCoords(:,1),SCoords(:,2),MDoG,xq,yq);
end

%% SAVE THE DATA
Cell.Global.ResponseSurface.LResponse=LResponseSurface;
Cell.Global.ResponseSurface.MResponse=MResponseSurface;
if cpAll(3)~=0
    Cell.Global.ResponseSurface.SResponse=SResponseSurface;
end

%% PLOT THE SURFACES
figure;surf(X,X,LResponseSurface);
shading interp;
axis square;
title('L RESPONSE');
figure;surf(X,X,MResponseSurface);
shading interp;
axis square;
title('M RESPONSE');
