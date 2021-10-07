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
analysisIDs = {'6117d4db18adcc19d6e0f820','611d158fa296f805e7a2da75'};
shortNames = {'gka','asb'};
directions = {'LminusM','S','LMS'};
freqs = [0,2,4,8,16,32,64];
nFreqs = length(freqs);
Color = {'r','b','k'};

% Frequency components for model fitting
deltaF10 = min(diff(log10(freqs(2:end))));
fitScaleUp = 10;
freqsFit = 10.^(log10(min(freqs(2:end)))-deltaF10+deltaF10/fitScaleUp:deltaF10/fitScaleUp:log10(max(freqs(2:end)))+deltaF10);

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

% This is the threshold for the goodness of fit to the fMRI time-series
% data
r2Thresh = 0.1;


% Create a figure
figure;

eccenDivs = [0 90./(2.^(5:-1:0))]; % eccentricity bins

data = cell(2,3,length(eccenDivs)-1,6);

% Loop through the directions and eccentricities

p = NaN*ones(2,3,2);
        
for ss = 1:2
    for dd = 1:3
        
        pBoot = NaN*ones(length(eccenDivs)-1,2,1000);
        Rsquared = NaN*ones(length(eccenDivs)-1,1000);
        
        for ee = 1:length(eccenDivs)-1
            eccenRange = [eccenDivs(ee) eccenDivs(ee+1)];
        
            % Load the results file for this subject
            filePath = fullfile(localSaveDir,'resultsFiles',[subjectNames{ss} '_mtSinai_results.mat']);
            load(filePath,'results')

            % Grab the stimLabels
            stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

            % Find the vertices that we wish to analyze (V1 for
            % eccentricity)
            areaIdx = (vArea == 1) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));
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
                data{ss,dd,ee,ff-1} = vals{ff}-vals{1};
                temp_data = vals{ff}-vals{1};
                bootVals = sort(bootstrp(1000,@mean,temp_data));
                mBoot(ff-1,:) = bootVals(500);
                lBoot(ff-1,:) = bootVals(25);
                uBoot(ff-1,:) = bootVals(975);

                semilogx(zeros(1,length(data{ss,dd,ee,ff-1}))+freqs(ff),data{ss,dd,ee,ff-1},'.','Color',[0.9 0.9 0.9]);
                hold on
            end

            % Plot bootstrapped mean across the frequencies
            errorbar(freqs(2:end),mBoot,abs(diff([mBoot lBoot],1,2)),abs(diff([uBoot mBoot],1,2)),'ob','MarkerSize',8,'MarkerFaceColor','b','LineWidth',2)

            % Obtain Watson fit for bootstrapped mean and plot it 
            y = mBoot'; % select bootstrapped mean to be fit to the model
            w = freqs(2:end);

            wDelta = min(diff(log10(w))); % Create a scaled-up, log-spaced, version of the frequency domain
            upScale = 10;
            wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);

            p0 = [1.5 0.015]; lb = [-inf 0.005]; ub = [8 0.025];
            
            options = optimoptions(@fmincon,... % The options for the search (mostly silence diagnostics)
                'Diagnostics','off',...
                'Display','off');

            % The objective function is the norm of the model fit error
            myObj = @(p) norm(y - watsonTTF2param(p,w));

            % Search
            p(ss,dd,:) = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

            % Obtain the high-resolution model fit and plot it
            yFit = watsonTTF2param(squeeze(p(ss,dd,:)),wFit);
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
            [pBoot(ee,:,:),Rsquared(ee,:)] = ExtractWatsonFitParametersBootstrap(w,vals);
            
            pBoot(ee,:,:) = sort(pBoot(ee,:,:),3);
            Rsquared(ee,:) = sort(Rsquared(ee,:),2);
            title([shortNames{ss} ', ' directions{dd} ', eccentricity' eccenDivs(ee) ' to ' eccenDivs(ee+1)])
        end
        % plot bootstrapped parameters across eccentricity

        figure(10)
        switch ss
            case 1
                loc = 1;
            case 2
                loc = 4;
        end
        
        subplot(2,3,loc)
        hold on
        errorbar(1:length(eccenDivs)-1,squeeze(pBoot(:,1,500)),abs(diff(squeeze(pBoot(:,1,[25 500])),[],2)),abs(diff(squeeze(pBoot(:,1,[500 975])),[],2)),'-o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        set(gca,'TickDir','out');
        box off
        xticks(1:length(eccenDivs))
        xticklabels(eccenDivs(1:end-1))
        xlim([0 7])
        ylim([-2 8])
        xlabel('gain')
        ylabel('fit parameter value')
        
        subplot(2,3,loc+1)
        hold on
        errorbar(1:length(eccenDivs)-1,squeeze(pBoot(:,2,500)),abs(diff(squeeze(pBoot(:,2,[25 500])),[],2)),abs(diff(squeeze(pBoot(:,2,[500 975])),[],2)),'-o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        set(gca,'TickDir','out');
        box off
        xticks(1:length(eccenDivs))
        xticklabels(eccenDivs(1:end-1))
        xlim([0 7])
        ylim([0 0.05])
        xlabel('time constant')
        
        subplot(2,3,loc+2)
        hold on
        errorbar(1:length(eccenDivs)-1,Rsquared(:,500),abs(diff(Rsquared(:,[25 500]),[],2)),abs(diff(Rsquared(:,[500 975]),[],2)),'-o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        set(gca,'TickDir','out');
        box off
        xticks(1:length(eccenDivs))
        xticklabels(eccenDivs(1:end-1))
        xlim([0 7])
        ylim([0 1])
        xlabel('R squared')
        
        p_Boot(ss,dd,:,:,:) = pBoot;
        R_squared(ss,dd,:,:) = Rsquared;
    end
end

save param_by_V1eccen p_Boot R_squared eccenDivs data