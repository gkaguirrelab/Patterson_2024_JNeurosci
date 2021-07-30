%% channelWeightPlotByEccen
%
% This routine loads the forwardModel results.mat files, and then creates
% plots of the TTF data and fitting across eccentricity bands

% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% Define where we want to save these figures
resultsSaveDir = fullfile(localSaveDir,'Fig X - channel weight by eccen');
mkdir(resultsSaveDir);

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'L-M','S','LF'};

% Load the retino maps
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
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
deltaF10 = min(diff(log10(freqs)));
fitScaleUp = 10;
freqsFit = 10.^(log10(min(freqs))-deltaF10+deltaF10/fitScaleUp:deltaF10/fitScaleUp:log10(max(freqs))+deltaF10);

% Which TTF model to use
ttfModel = 'watson';

% x0 and bounds
switch ttfModel
    case 'watson'
        p0 = [1.72, 0.87, 0.013];
        lb = [0 0 0];
        ub = [5 1 0.025];
    case 'beta'
        p0 = [1.8518    5.2050    6.7683   11.7180];
        lb = [0 0 0 1];
        ub = [100 8 100 100];
end

% define some search options
options = optimoptions(@fmincon,...
    'Diagnostics','off',...
    'Display','off');

% Set up a variable to hold the model fit results
pVals = [];

% Loop through the subjects and directions
for ss = 1:2
        
    allResults = cell(3,1);
    
    % Load the results files for this subject
    for dd = 1:3
        
        % Load the results file for this subject / direction
        filePath = fullfile(localSaveDir,'resultsFiles',[analysisLabels{dd} '_' subjectNames{ss} '_agtcOL_results.mat']);
        load(filePath,'results')
        
        % Store the results file in a cell arrray
        allResults{dd}=results;
        
    end
    
    % Set the eccentricity divisions
    eccenDivs = [0 90./(2.^(5:-1:0))];
    
    % Create a figure
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
            switch ttfModel
                case 'watson'
                    myObj = @(p) norm(yVals - watsonTTF(p,freqs));
                    myNonlcon = [];
                case 'beta'
                    myObj = @(p) norm(yVals - betaTTF(p,freqs));
                    myNonlcon = @(p) betaTTF(p,freqsFit,yVals);
            end
            
            % Search
            p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
            
            % Save the params
            pVals(ss,dd,ee,:)=p;
            
            % Get the fitted response
            switch ttfModel
                case 'watson'
                    myFit = watsonTTF(p,freqsFit);
                case 'beta'
                    myFit = betaTTF(p,freqsFit);
            end
            myFit(~isfinite(myFit))=nan;
            
            % Store the interpolated maximum response
            y(dd)=max(myFit);
            
            % Add this fit to the plot
            subplot(3,length(eccenDivs)-1,ee+((dd-1)*(length(eccenDivs)-1)))
            semilogx(freqs,yVals,'ok')
            hold on
            semilogx(freqsFit,myFit,'-r')
            semilogx([1 64],[0 0],':k','LineWidth',1)
            ylim([-1 7]);
            axis off
            
            % Report the parameters
           str = {sprintf('[%2.1f, %2.2f]',p(1:2)),...
               sprintf('%2.1f ms',p(3)*1000)};
           text(1,6,str);
            
            % Add some chart stuff
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
                semilogx([1 1],[0 7],'-k','LineWidth',1)
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
    
