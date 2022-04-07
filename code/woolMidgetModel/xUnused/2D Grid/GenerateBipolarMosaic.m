function BipolarMosaic=GenerateBipolarMosaic(Em)
%% GENERATE THE BIPOLAR MOSAIC
%Load a universal hex grid
cd('/Users/lauren/Documents/SCHOOL/Dacey/RF Model/2D Grid');
load('HexGrid');

%Jitter the cartesian coordinates to make a unique mosaic
pixel_coords(:,1)=jitter(pixel_coords(:,1));
pixel_coords(:,2)=jitter(pixel_coords(:,2));

%Determine the bipolar density given the eccentricity
BipolarsPerMM=ceil(((52170*exp(-1.057.*Em))+(15550*exp(-0.1605.*Em)))/(2.67)); %From 'cones per midget gc data' Excel file
DFRadius=-0.000565*(Em^2)+(0.01863*Em)+0.00146; %In mm. Polynomial fit of data from Fig 9, Watanabe and Rodieck JCN 1989

%Determine the aperture width needed to encompass the correct number of
%bipolars
r=1;
BipolarCount=0;
while BipolarCount<BipolarsPerMM
    r=r+.1;
    counter=0;
    inc_pix_coords=[];
    for p=1:length(pixel_coords)
        if pixel_coords(p,1)>-r && pixel_coords(p,1)<r
            if pixel_coords(p,2)>-r && pixel_coords(p,2)<r
                inc_pix_coords(end+1,1)=pixel_coords(p,1);
                inc_pix_coords(end,2)=pixel_coords(p,2);
                counter=counter+1;
            end
        end
    end
    BipolarCount=counter;
end

%Adjust the aperture window to a 1-mm2 portion...now you have cones/mm2
RFactor=1/(2*r);
pixel_coordsADJ=pixel_coords*RFactor;
inc_pix_coordsADJ=inc_pix_coords*RFactor;

%Compute the distances of each bipolar from the center (0,0) and sort the
%coordinates from closest to farthest
distances=pdist2(inc_pix_coordsADJ,[0,0]);
SortedTotalBipolars=sortrows([inc_pix_coordsADJ,distances],3);

%Determine each bipolar's nearest neighbors and use their average as an estimate of
%intercone distance across the retina
NearestNeighbors=zeros(length(SortedTotalBipolars),1);
HowManyNeighbors=1;
parfor q=1:length(SortedTotalBipolars)
    dist=sort(pdist2(SortedTotalBipolars(q,:),SortedTotalBipolars));
    NearestNeighbors(q)=mean(dist(2:HowManyNeighbors+1));
end
InterBipolarDistance=mean(NearestNeighbors);

%% SAVE TO OUTPUT

BipolarMosaic.BipolarsPerMM=BipolarCount;
BipolarMosaic.AllCoords=SortedTotalBipolars;
BipolarMosaic.AllDistances=distances;
BipolarMosaic.InterBipolarDistance=InterBipolarDistance;