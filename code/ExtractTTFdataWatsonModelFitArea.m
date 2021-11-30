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

p_Boot = NaN*ones(2,3,length(areaLabels),length(p0),1000);
R_squared =  NaN*ones(2,3,length(areaLabels),1000);
p_mean = NaN*ones(2,3,length(areaLabels),length(p0));
        
% Loop through the directions and visual areas
for ss = 1:2
    for dd = 1:3
        pBoot = NaN*ones(length(areaLabels),length(p0),1000);
        Rsquared = NaN*ones(length(areaLabels),1000);
        for aa = 1:length(areaLabels)
            areaLabel = areaLabels{aa};
            % Load the results file for this subject
            filePath = fullfile(localSaveDir,'resultsFiles',[subjectNames{ss} '_mtSinai_results.mat']);
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

            % Obtain Watson fit for bootstrapped mean and plot it 
            Y = mBoot'; % select bootstrapped mean to be fit to the model
            w = freqs(2:end);
            
            y = Y - min(Y);

            wDelta = min(diff(log10(w))); % Create a scaled-up, log-spaced, version of the frequency domain
            upScale = 10;
            wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);

            options = optimoptions(@fmincon,... % The options for the search (mostly silence diagnostics)
                'Diagnostics','off',...
                'Display','off');

            % fits model, and provides high
            % resolution fit

            myObj = @(p) norm(y - watsonTTF(p,w));
            p0(1,1) = max(y);
            p = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
            yFit = watsonTTF(p,wFit);
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
            [pBoot(aa,:,:),Rsquared(aa,:)] = ExtractWatsonFitParametersBootstrap(w,vals,'p0',p0,'lb',lb,'ub',ub);
            title([shortNames{ss} ', ' directions{dd} ', visual area ' areaLabels{aa}])
            
            pBoot(aa,:,:) = sort(pBoot(aa,:,:),3);
            Rsquared(aa,:) = sort(Rsquared(aa,:),2);
            p_mean(ss,dd,aa,:) = p;
            clear p
        end
        % plot bootstrapped parameters across eccentricity
        
                figure(10)
        switch ss
            case 1
                loc = 1;
            case 2
                loc = 5;
        end
        subplot(2,4,loc)
        hold on
        errorbar(1:length(areaLabels),squeeze(pBoot(:,1,500)),abs(diff(squeeze(pBoot(:,1,[25 500])),[],2)),abs(diff(squeeze(pBoot(:,1,[500 975])),[],2)),'o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        set(gca,'TickDir','out');
        box off
        xticks(1:length(areaLabels))
        xticklabels(areaLabels)
        xlim([0 7])
        ylim([-2 6])
        xlabel('gain')
        ylabel('fit parameter value')
        
        if length(p0)==3
            subplot(2,4,loc+1)
            hold on
            errorbar(1:length(areaLabels),squeeze(pBoot(:,2,500)),abs(diff(squeeze(pBoot(:,2,[25 500])),[],2)),abs(diff(squeeze(pBoot(:,2,[500 975])),[],2)),'o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
            set(gca,'TickDir','out');
            box off
            xticks(1:length(areaLabels))
            xticklabels(areaLabels)
            xlim([0 7])
            ylim([0 2])
            xlabel('surround gain')
        end
        
        subplot(2,4,loc+2)
        hold on
        errorbar(1:length(areaLabels),squeeze(pBoot(:,end,500)),abs(diff(squeeze(pBoot(:,end,[25 500])),[],2)),abs(diff(squeeze(pBoot(:,end,[500 975])),[],2)),'o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        set(gca,'TickDir','out');
        box off
        xticks(1:length(areaLabels))
        xticklabels(areaLabels)
        xlim([0 7])
        ylim([0 0.04])
        xlabel('time constant')
        
        subplot(2,4,loc+3)
        hold on
        errorbar(1:length(areaLabels),Rsquared(:,500),abs(diff(Rsquared(:,[25 500]),[],2)),abs(diff(Rsquared(:,[500 975]),[],2)),'o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        set(gca,'TickDir','out');
        box off
        xticks(1:length(areaLabels))
        xticklabels(areaLabels)
        xlim([0 7])
        ylim([0 1])
        xlabel('R squared')
        
        p_Boot(ss,dd,:,:,:) = pBoot;
        R_squared(ss,dd,:,:) = Rsquared;
        
        clear pBoot Rsquared
    end
end

save param_by_area p_mean p_Boot R_squared areaLabels data