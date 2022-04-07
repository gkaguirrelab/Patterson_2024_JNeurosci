function SFLines=ComputeSFLine(thetas,W,Freq)
% This function draws a line in frequency space that represents a
% continuous vector of spatial frequencies (0-64 cpd) of a 2D grating at a particular
% orientation. This can be used to retrieve amplitude and phase data from a
% model cell's FFT given a series of gratings with a single orientation.
%
% theta: orientation, in degrees
% W: size of a single dimension of the FFT domain matrix
% Freq: 1D frequency vector that defines frequency space
for z=1:12
    zInd(z,1)=2^(z-1);
end
for t=1:length(thetas)
    theta=thetas(t);
    OrientationLabel=strcat('Theta',num2str(thetas(t)),'deg');
    
    %Find the index value corresponding to SpatFreq=0 (i.e., origin). This is
    %P1 of the line
    [~, oriInd]=min(abs(Freq)); 

    %Choose an orientation for the continuous grating stimuli
    thetaRad=(theta/360)*2*pi; % convert to radians

    %Determine P2 of the line. This is much farther out than any sensible SF,
    %to be safe
    SFx=64*cos(thetaRad);
    SFy=64*sin(thetaRad);

    %Find the corresponding subindices for P2
    [~, SFi] = min(abs(Freq-SFx));
    [~, SFj] = min(abs(Freq-SFy));

    %Assemble the indices for P1 nd P2
    p1=[oriInd oriInd];
    p2=[SFi SFj];

    %Compute the line containing all SFs for the given theta
    Line=drawline(p1,p2,[W W]);
    
    %Convert indices to subindices
    [iLine,jLine]=ind2sub(W,Line);

    %True SFs  given the subindices (i.e.,Pythagorean Theorem)
    LineFreqs=sqrt(Freq(iLine).^2+Freq(jLine).^2);
    
    SFLines.(matlab.lang.makeValidName(OrientationLabel)).Indices=Line;
    SFLines.(matlab.lang.makeValidName(OrientationLabel)).zIndices=Line(zInd);
    SFLines.(matlab.lang.makeValidName(OrientationLabel)).SubIndices=[iLine',jLine'];
    SFLines.(matlab.lang.makeValidName(OrientationLabel)).zSubIndices=[iLine(zInd)',jLine(zInd)'];
    SFLines.(matlab.lang.makeValidName(OrientationLabel)).Freqs=LineFreqs;
    SFLines.(matlab.lang.makeValidName(OrientationLabel)).zFreqs=LineFreqs(zInd);
end