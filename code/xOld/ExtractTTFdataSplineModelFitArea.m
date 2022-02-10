%% ExtractTTFdataWatsonModelFit
%
% This script downloads the "results" files Flywheel and
% extracts BOLD fMRI response amplitudes for each of the stimulus temporal
% frequencies, plots them, and calculates 95% CI and Watson model fits


% Get the localSaveDir pref
localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
directions = {'LminusM','S','LMS'};
freqs = [0,2,4,8,16,32,64];
nFreqs = length(freqs);
Color = {'r','b','k'};

% TTF model guess
p0 = [1.5 0.9 0.015]; 
lb = [-2 0.5 0.005]; 
ub = [6 2 0.05];

% Frequency components for model fitting
deltaF10 = min(diff(log10(freqs(2:end))));
fitScaleUp = 10;
wFit = 10.^(log10(min(freqs(2:end)))-deltaF10+deltaF10/fitScaleUp:deltaF10/fitScaleUp:log10(max(freqs(2:end)))+deltaF10);


% Create a flywheel object. You need to set you flywheelAPIKey in the
% "flywheelMRSupport" local hook.
fw = flywheel.Flywheel(getpref('flywheelMRSupport','flywheelAPIKey'));

% Load the retino maps
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_varea.dtseries.nii');
vArea = cifti_read(tmpPath); vArea = vArea.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_eccen.dtseries.nii');
eccenMap = cifti_read(tmpPath); eccenMap = eccenMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_angle.dtseries.nii');
polarMap = cifti_read(tmpPath); polarMap = polarMap.cdata;
tmpPath = fullfile(localSaveDir,'retinoFiles','TOME_3021_inferred_sigma.dtseries.nii');
sigmaMap = cifti_read(tmpPath); sigmaMap = sigmaMap.cdata;

% Load the subcortical ROIs
projectID = '5ca7803af546b60029ef118e';
subCorticalROIsFullNames = {'LGN_bilateral.dtseries.nii','thalamus_bilateral.dtseries.nii','midbrain_bilateral.dtseries.nii'};
subCorticalROIsLabels = {'LGN','thalamus','midbrain'};
for rr = 1:length(subCorticalROIsFullNames)
    tmpPath = fullfile(localSaveDir,'retinoFiles',subCorticalROIsFullNames{rr});    
    fw.downloadFileFromProject(projectID,subCorticalROIsFullNames{rr},tmpPath);
    tmpRegion = cifti_read(tmpPath); tmpRegion = tmpRegion.cdata;
    str = [subCorticalROIsLabels{rr} 'ROI = tmpRegion;'];
    eval(str);
end

% This is the threshold for the goodness of fit to the fMRI time-series
% data
r2Thresh = 0.1;

eccenRange = [0 90];

% Create a figure
figure;

% Create a data variable to hold the results. This will be a 2 x 3 x 6
% (subjects x directons x frequencies) cell array for the selected subject. Each cell
% will have the 12 measurements
data = cell(2,3,6);

areaLabels = {'LGN';'V1'; 'V23'; 'hV4'; 'MT'};

R_squared =  NaN*ones(2,3,length(areaLabels),1000);
        
% Loop through the directions and visual areas
for ss = 1:2
    for dd = 1:3
        pBoot = NaN*ones(length(areaLabels),length(p0),1000);
        Rsquared = NaN*ones(length(areaLabels),1000);
        for aa = 1:length(areaLabels)
            areaLabel = areaLabels{aa};
            % Load the results file for this subject
            filePath = fullfile(localSaveDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
            load(filePath,'results')
           
               % Grab the stimLabels
                stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};
        
        % Find the vertices that we wish to analyze
        switch areaLabel
            case 'V1'
                areaIdx = (vArea==1) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));
            case 'V23'
                areaIdx = (vArea>=2) .* (vArea<=3) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));
            case 'hV4'
                areaIdx = (vArea>=4) .* (vArea<=5) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));
            case 'MT'
                areaIdx = (vArea>=8) .* (vArea<=9) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));
            case 'LGN'
                areaIdx = LGNROI;
            case 'midbrain'
                areaIdx = midbrainROI;
            case 'thalamus'
                areaIdx = thalamusROI;
        end
        
        goodIdx = logical( (results.R2 > r2Thresh) .* areaIdx );
        
        % Loop through the frequencies and obtain the set of values
        vals = cell(1,nFreqs);
        for ff = 1:nFreqs
            subString = sprintf(['f%dHz_' directions{dd}],freqs(ff));
            idx = find(contains(stimLabels,subString));
            vals{ff} = mean(results.params(goodIdx,idx),'omitnan');
        end

            % Prepare to plot into this subplot
            figure(1)
            subplot(2,3,dd+(ss-1)*3);

            % Adjust the values for the zero frequency, obtain bootstrapped
            % medians and 95% CI
            mBoot = NaN*ones(nFreqs-1,1);
            lBoot = NaN*ones(nFreqs-1,1);
            uBoot = NaN*ones(nFreqs-1,1);
            for ff = 2:nFreqs
                data{ss,dd,aa,ff-1} = vals{ff}-vals{1};
                temp_data = vals{ff}-vals{1};
                bootVals = sort(bootstrp(1000,@mean,temp_data));
                mBoot(ff-1,:) = bootVals(500);
                lBoot(ff-1,:) = bootVals(25);
                uBoot(ff-1,:) = bootVals(975);

                semilogx(zeros(1,length(data{ss,dd,aa,ff-1}))+freqs(ff),data{ss,dd,aa,ff-1},'.','Color',[0.9 0.9 0.9]);
                hold on
            end

            % Plot bootstrapped mean across the frequencies
            errorbar(freqs(2:end),mBoot,abs(diff([mBoot lBoot],1,2)),abs(diff([uBoot mBoot],1,2)),'ob','MarkerSize',8,'MarkerFaceColor','b','LineWidth',2)

            % Obtain Spline fit for bootstrapped mean and plot it 
            Y = mBoot'; % select bootstrapped mean to be fit to the model
            w = freqs(2:end);
            
            y = Y - min(Y);

            wDelta = min(diff(log10(w))); % Create a scaled-up, log-spaced, version of the frequency domain
            upScale = 10;
            wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);

            yFit = spline(w,y,wFit);
            yFit(~isfinite(yFit))=nan;
            
            yFit = yFit+min(Y);
                
            plot(wFit,yFit,'-b','LineWidth',1)

            % Clean up the plot
            title([shortNames{ss} ' ' directions{dd}])
            ylabel('BOLD % change');
            xlabel('frequency [Hz]');
            semilogx([1 64],[0 0],':k','LineWidth',1)
            ylim([-2 8]);
            xlim([1 128])
            set(gca,'TickDir','out');
            box off

            % Get bootstrapped fit parameters
            y0 = vals{1}; % response at frequency = 0
            yW = vals(:,2:end); % responses to all other frequencies
            bootVals = NaN*ones(length(w),1000);
            for ff = 1:length(w)
                temp_data = yW{ff}-y0;
                bootVals(ff,:) = sort(bootstrp(1000,@median,temp_data));
            end

            Rsquared = NaN*ones(1,1000);

            % Obtain Spline fit for bootstrapped means
            for bb = 1:1000
                y = bootVals(:,bb)';

                % shift so that if min(y)=0
                y = y - min(y);


               % Spline
                yFit_boot = spline(w,y,w);

                R = corrcoef(y,yFit_boot);
                Rsquared(1,bb) = R(1,2)^2;
            end

            title([shortNames{ss} ', ' directions{dd} ', visual area ' areaLabels{aa}])
            
            R_squared(ss,dd,aa,:) = sort(Rsquared,2);
            clear Rsquared
        end
    end
end

save param_by_areaSpline R_squared areaLabels data