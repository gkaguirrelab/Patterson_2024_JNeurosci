function OriCircles=ComputeOriCircle(SpatFreqs,W,Freq)
% This function draws a circle in frequency space that represents a
% continuous vector of orientations (0-360 degrees) of a 2D grating at a particular
% spatial frequency. This can be used to retrieve amplitude and phase data from a
% model cell's FFT given a series of gratings with a single orientation.
%
% SpatFreqs: vector of spatial frequencies, in cpd
% W: size of a single dimension of the FFT domain matrix
% Freq: 1D frequency vector that defines frequency space

for t=1:length(SpatFreqs)
    SpatFreq=SpatFreqs(t);
    SFLabel=strcat('SF',num2str(SpatFreqs(t)),'cpd');
    [~, radInd]=min(abs(Freq-SpatFreq)); %Find the index value corresponding to SpatFreq
    [~, oriInd]=min(abs(Freq)); %Find the index value corresponding to SpatFreq=0 (i.e., origin)

    CircleSubstrate=zeros(W);
    Radius=abs(radInd-oriInd);

    %Compute the circle containing all orientations for the given SF
    CircleMatrix=MidpointCircle(CircleSubstrate,Radius,oriInd,oriInd,1);
    Circle=find(CircleMatrix);

    %Convert indices to subindices
    [iCircle,jCircle]=ind2sub(W,Circle);

    %True orientations  given the subindices (i.e.,Pythagorean Theorem)
    Orientations=atan2((iCircle-oriInd),(jCircle-oriInd));
    for o=1:length(Orientations)
        if Orientations(o)<0
        Orientations(o)=Orientations(o)+(2*pi);
        end
    end

    AllOriData=sortrows([iCircle, jCircle, Orientations],3);

    OriCircles.(matlab.lang.makeValidName(SFLabel)).Indices=sub2ind([W W],AllOriData(:,1),AllOriData(:,2));
    OriCircles.(matlab.lang.makeValidName(SFLabel)).SubIndices=[AllOriData(:,1),AllOriData(:,2)];
    OriCircles.(matlab.lang.makeValidName(SFLabel)).Oris=AllOriData(:,3);
end