%% demoExtractTTFData
%
% This script demonstrates how to download the "results" files Flywheel and
% extract BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies


% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'L-M','S','LF'};
analysisIDs = { {'60e9ea50a85492ed8f96cabd','60e9ea334ef89230db2b7021','60e9ea6dbceb4c0bc9e0767e'} , ...
    {'60e9eaa1a74445f40c56b123', '60e9ea85bd00f64426dd9301','60e9eabf4ef89230db2b7027'} };

% Set up a temporary location to save files we will need
scratchSaveDir = tempdir();

% Create a flywheel object. You need to set you flywheelAPIKey in the
% "flywheelMRSupport" local hook.
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Download and unzip the retino maps
retinoMapID = '5dc88aaee74aa3005e169380';
retinoFileName = 'TOME_3021_cifti_maps.zip';
saveDir = fullfile(scratchSaveDir,'v0','output');
mkdir(saveDir);
tmpPath = fullfile(saveDir,retinoFileName);
fw.downloadOutputFromAnalysis(retinoMapID,retinoFileName,tmpPath);
command = ['unzip -q -n ' tmpPath ' -d ' saveDir];
system(command);

% Load the retino maps
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_varea.dtseries.nii'));
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_eccen.dtseries.nii'));
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_angle.dtseries.nii'));
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(saveDir,strrep(retinoFileName,'_cifti_maps.zip','_inferred_sigma.dtseries.nii'));
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Create a "subcortical" map
subcorticalMap = zeros(size(vArea));
subcorticalMap(1:26298)=1;

% This is the threshold for the goodness of fit to the fMRI time-series
% data
r2Thresh = 0.25;
attentionThresh = 0;

% This is the visual area and eccentricity range to grab. The visual areas
% are: V1 = 1, V2 = 2, V3 = 3, hV4/LO = [4 5], MT/MST = [8 9]
area = 1;

% Some values used for fitting peak frequency
freqs = [2 4 8 16 32 64];
nFreqs = length(freqs);
freqsIdx = 1:nFreqs;
deltaF10 = min(diff(log10(freqs)));
fitScaleUp = 10;
freqsFit = 10.^(log10(min(freqs))-deltaF10+deltaF10/fitScaleUp:deltaF10/fitScaleUp:log10(max(freqs))+deltaF10);
freqsFitIdx = 1/fitScaleUp:1/fitScaleUp:nFreqs+1;

% The fitting function
myFunc = @(f,A,B,C,D)  C.*ncbeta(f./D, 1e-6, A, B );

% x0 and bounds
x0 = [14.4936   25.5673    0.9136    6.5449];
lb = [0 1 0 0];
ub = [100 100 100 8];

% define some search options
options = optimoptions(@fmincon,...
    'Diagnostics','off',...
    'Display','off');


% Loop through the subjects and directions
for ss = 1:2
    
    % Set up the paths for this subject
    fileStem = [subjectNames{ss} '_agtcOL_'];
    resultsSaveDir = ['/Users/aguirre/Dropbox (Aguirre-Brainard Lab)/_Papers/Patterson_2021_EccentricityFlicker/matlabFigures/' subjectNames{ss}];
    
    allResults = cell(3,1);
    
    for dd = 1:3
        
        % Download the results file for this subject / direction
        fileStem = [subjectNames{ss} '_agtcOL_'];
        fileName = [fileStem 'results.mat'];
        tmpPath = fullfile(saveDir,[analysisLabels{dd} '_' fileName]);
        fw.downloadOutputFromAnalysis(analysisIDs{ss}{dd},fileName,tmpPath);
        
        % Load the result file into memory and delete the downloaded file
        clear results
        load(tmpPath,'results')
        
        % Store the results file in a cell arrray
        allResults{dd}=results;
        
    end
    
    % Set the eccentricity divisions
    eccenDivs = [0 90./(2.^(5:-1:0))];
    
    figHandle0 = figure();
    
    dkl = [];
    rgb=[];
    for ee = 1:length(eccenDivs)-1
        % All vertices in V1 that are within the eccentricity range, and
        % for which there was a response to the attention task above
        % threshold.
        goodIdxA = ( (vArea>=1) .* (vArea<=1) .* (eccenMap > eccenDivs(ee)) .* (eccenMap < eccenDivs(ee+1)) ) .* ...
            ( ((allResults{1}.attention-allResults{1}.f0Hz)>attentionThresh) .* ((allResults{2}.attention-allResults{2}.f0Hz)>attentionThresh) .* ((allResults{3}.attention-allResults{3}.f0Hz)>attentionThresh) );
        
        % Also require that at least one of the three modulation directions
        % was associated with an fMRI time-series fit above the R2
        % threshold.
        goodIdxB = (allResults{1}.R2 > r2Thresh) + (allResults{2}.R2 > r2Thresh) + (allResults{3}.R2 > r2Thresh);
        goodIdxB(goodIdxB>0)=1;
        goodIdx = logical(goodIdxA .* goodIdxB);
        
        y = [];
        for dd = 1:3
            
            % Get the beta values for these indices
            yVals = allResults{dd}.params(goodIdx,1:nFreqs+1);
            
            % Get the mean betas across vertices / voxels
            yVals = nanmean(yVals);
            
            % Obtain the set of beta valyes, relative to the baseline condition
            yVals = yVals(2:end)-yVals(1);
            
            % Define the objective and non-linear constraint
            myObj = @(p) norm(yVals - myFunc(freqsIdx,p(1),p(2),p(3),p(4)));
            myNonlcon = @(p) betaNonlcon(p,yVals,freqsFitIdx);
            
            % Fit
            p = fmincon(myObj,x0,[],[],[],[],lb,ub,myNonlcon,options);
            R2 = corr(yVals',myFunc(freqsIdx,p(1),p(2),p(3),p(4))').^2;
            myFit = myFunc(freqsFitIdx,p(1),p(2),p(3),p(4));
            myFit(~isfinite(myFit))=nan;
            
            % Store the interpolated maximum response
            y(dd)=max(myFit);
            
            % Add this fit to the plot
            subplot(3,length(eccenDivs)-1,ee+((dd-1)*(length(eccenDivs)-1)))
            semilogx(freqs,yVals,'ok')
            hold on
            semilogx(freqsFit,myFit,'-r')
            semilogx([1 64],[0 0],':k','LineWidth',1)
            ylim([-1 6]);
            axis off
            if dd==1
                ax = gca;
                ax.Title.Visible = 'on';
                str = {sprintf('%2.1f-%2.1f°',eccenDivs(ee:ee+1)),...
                    sprintf('n=%d',sum(goodIdx))};
                title(str);
            end
            if ee==1
                ax = gca;
                ax.YLabel.Visible = 'on';
                ylabel(analysisLabels{dd});
                semilogx([1 1],[0 5],'-k','LineWidth',1)
            end
            
        end
        
        % Normalize the set of responses across the three channels to have
        % unit magnitude
        dkl(ee,:) = y/norm(y);
        
    end
    
    % Save the TTF plot
    outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_TTFbyEccenbyChannel.pdf']);
    print(figHandle0,outFile,'-dpdf','-bestfit');
    
    % Plot the DKL surface
    figHandle1 = figure();
    
    % Render the DKL surface
    [X,Y,Z] = sphere(40);
    XVec = X(:); YVec = Y(:); ZVec = Z(:);
    dklSurf = [XVec, YVec, ZVec];
    rgbSurf = dklToRGB(dklSurf);
    posQuadrant = logical((XVec>=-1e-6) .* (YVec>=-1e-6) .* (ZVec>=-1e-6));
    rgbSurf(~posQuadrant,:)=nan;
    rgbSurf = reshape(rgbSurf,41,41,3);
    X(~posQuadrant)=nan;
    Y(~posQuadrant)=nan;
    Z(~posQuadrant)=nan;
    surf(X,Y,Z,rgbSurf,'EdgeColor','none','FaceColor','interp');
    hold on
    
    plotOnSphere(dkl(:,1),dkl(:,2),dkl(:,3),'Color','k','LineWidth',2);
    
    % Add circles around each location
    rVals = linspace(0.9999,0.999,length(eccenDivs)-1);
    for ii = 1:size(dkl,1)
        addCircleAtCoordinate(dkl(ii,:),rVals(ii));
    end
    
    axis equal
    xlim([0 1])
    ylim([0 1])
    zlim([0 1])
    zlabel('luminance');
    xlabel('red-green');
    ylabel('yellow-violet');
    view(135,45);
    
    % Save the DKL surface plot
    outFile = fullfile(resultsSaveDir, [subjectNames{ss} '_RelativeChannelResponses.pdf']);
    print(figHandle1,outFile,'-dpdf','-bestfit');
    
end




function [c, ceq] = betaNonlcon(p,y,f)

% Evaluate the function
yFit = p(3).*ncbeta(f./p(4), 1e-6, p(1), p(2) );
yFit = yFit(isfinite(yFit));

% The rate of change should not exceed 1 unit
d = diff(yFit);
d = d(isfinite(d));
c = max(abs(d))-1;

% The interpolated peak should not be more than 10% the size of the
% largest y value
peak = (max(yFit)-max(y))./max(y) - 0.1;
if peak > 0
    ceq = peak;
else
    ceq = 0;
end

end


function [ pz ] = ncbeta (x, a, b, lam )
%This function computes the probability density function for the
%noncentral beta distribution using a transformation of variables to put
%the desired density function in terms of a noncentral F, which is included
%in Matlabs statsitics toolbox already.
%
% USAGES
% [ pz ] = ncbeta (x, a, b, lambda )
%
% INPUT
% x:    Vector of possible argument values for the pdf.
% a:    The first degree of freedom/shape parameter.
% b:    The second degree of freedom/shape parameter.
% lam:  The noncentrality parameter.
%
% OUTPUT
% pz:   The probability density function for x, given the input parameters.
%
%-----------------------------------------------------------------------
% Latest Edit: 18.Feb.2014
% Joshua D Carmichael
% josh.carmichael@gmail.com
%-----------------------------------------------------------------------
const   = a./b;
pz      = @(r) ncfpdf( r./(const*(1-r)), a, b, lam).* 1./(const*(1-r).^2);
pz      = pz(x);

end


function h = addCircleAtCoordinate(X,r)

theta=-pi:0.1:pi;

    x = r*ones(size(theta));
    y = ((1-r^2)^0.5)*cos(theta);
    z = ((1-r^2)^0.5)*sin(theta);
    x = [x x(1)];
    y = [y y(1)];
    z = [z z(1)];
    
    h = plot3(x,y,z,'-k','LineWidth',2);
    
    [az,el] = cart2sph(X(1),X(2),X(3));

    rotate(h,[0 0 1],rad2deg(az),[0 0 0]);
    R = [cos(az) -sin(az); sin(az) cos(az)];
    direction =  [(R*[0;1])' 0];
    rotate(h,direction,-rad2deg(el),[0 0 0]);
end
    
