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
areaLabel = {'LGN','V1','V2','V3','hV4/LO','MT/MST'};

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

% This is the visual area and eccentricity range to grab. The visual areas
% are: V1 = 1, V2 = 2, V3 = 3, hV4/LO = [4 5], MT/MST = [8 9]
area = 1;
% eccenRange = [0 90];

% Create a figure
figure;

% Create a data variable to hold the results. This will be a 2 x 3 x 6
% (subjects x directons x frequencies) cell array for the selected subject. Each cell
% will have the 12 measurements
data = cell(2,3,6);

eccenRange = [0 90];

areas = {0; 1; 2; 3; [4 5]; [8 9]};

% Loop through the directions and visual areas
p = NaN*ones(2,3,3);
for ss = 1:2
    for dd = 1:3
        pBoot = NaN*ones(length(areas),3,1000);
        for aa = 1:length(areas)
            area = areas{aa};
            % Load the results file for this subject
            filePath = fullfile(localSaveDir,'resultsFiles',[subjectNames{ss} '_mtSinai_results.mat']);
            load(filePath,'results')
            stimLabels = results.model.opts{6};

            % Find the vertices that we wish to analyze
            switch area(1)
                case 0
                    goodIdx = logical( (results.R2 > r2Thresh) .* (subcorticalMap==1)  );
                case {1, 2, 3}
                    goodIdx = logical( (results.R2 > r2Thresh) .* (vArea==area) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2)) );
                otherwise
                    goodIdx = logical( (results.R2 > r2Thresh) .* (vArea==area(1) | vArea==area(2)) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2)) );
            end

            % Loop through the frequencies and obtain the set of values
            vals = cell(1,nFreqs);
            for ff = 1:nFreqs
                subString = sprintf(['f%dHz_' directions{dd}],freqs(ff));
                idx = find(contains(stimLabels,subString));
                vals{ff} = mean(results.params(goodIdx,idx));
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
                data{ss,dd,ff-1} = vals{ff}-vals{1};
                temp_data = vals{ff}-vals{1};
                bootVals = sort(bootstrp(1000,@median,temp_data));
                mBoot(ff-1,:) = bootVals(500);
                lBoot(ff-1,:) = bootVals(25);
                uBoot(ff-1,:) = bootVals(975);

                semilogx(zeros(1,length(data{ss,dd,ff-1}))+freqs(ff),data{ss,dd,ff-1},'.','Color',[0.9 0.9 0.9]);
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


            p0 = [1.5, 0.8, 0.015]; lb = [0 0 0]; ub = [8 1 0.025];  % Set up the p0 guess, and the bounds on the params based on post-receptoral channel

            options = optimoptions(@fmincon,... % The options for the search (mostly silence diagnostics)
                'Diagnostics','off',...
                'Display','off');

            % The objective function is the norm of the model fit error
            myObj = @(p) norm(y - watsonTTF(p,w));

            % Search
            p(ss,dd,:) = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);

            % Obtain the high-resolution model fit and plot it
            yFit = watsonTTF(squeeze(p(ss,dd,:)),wFit);
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
            [pBoot(aa,:,:)] = ExtractWatsonFitParametersBootstrap(w,vals);
            title([shortNames{ss} ', ' directions{dd} ', visual area' areaLabel{aa}])
        end
        % plot bootstrapped parameters across eccentricity

        figure
        subplot(1,3,1)
        hold on
        errorbar(1:length(areas),squeeze(pBoot(:,1,500)),squeeze(pBoot(:,1,500))-squeeze(pBoot(:,1,25)),squeeze(pBoot(:,1,975))-squeeze(pBoot(:,1,500)),'-o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        set(gca,'TickDir','out');
        box off
        xticks(1:length(areas))
        xticklabels(areaLabel)
        xlim([0 7])
        ylim([0 8])
        xlabel('gain')
        ylabel('fit parameter value')

        subplot(1,3,2)
        hold on
        errorbar(1:length(areas),squeeze(pBoot(:,2,500)),squeeze(pBoot(:,2,500))-squeeze(pBoot(:,2,25)),squeeze(pBoot(:,2,975))-squeeze(pBoot(:,2,500)),'-o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        set(gca,'TickDir','out');
        box off
        xticks(1:length(areas))
        xticklabels(areaLabel)
        xlim([0 7])
        ylim([0.5 1.2])
        xlabel('surround gain')
        ylabel('fit parameter value')

        subplot(1,3,3)
        hold on
        errorbar(1:length(areas),squeeze(pBoot(:,3,500)),squeeze(pBoot(:,3,500))-squeeze(pBoot(:,3,25)),squeeze(pBoot(:,3,975))-squeeze(pBoot(:,3,500)),'-o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        set(gca,'TickDir','out');
        box off
        xticks(1:length(areas))
        xticklabels(areaLabel)
        xlim([0 7])
        ylim([0 0.03])
        xlabel('time constant')
        ylabel('fit parameter value')
        
        clear p pBoot
    end
end
