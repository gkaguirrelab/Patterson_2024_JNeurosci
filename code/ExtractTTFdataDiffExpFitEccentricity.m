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
w = freqs(2:end);
Color = {'r','b','k'};

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

% This is the threshold for the goodness of fit to the fMRI time-series
% data
r2Thresh = 0.1;

% Create a figure
figure;

eccenDivs = [0 90./(2.^(5:-1:0))]; % eccentricity bins

data = cell(2,3,length(eccenDivs)-1,nFreqs-1);

% Loop through the directions and eccentricities
        
for ss = 1:2
    for dd = 1:3
        for ee = 1:length(eccenDivs)-1
            eccenRange = [eccenDivs(ee) eccenDivs(ee+1)];
        
            % Load the results file for this subject
            filePath = fullfile(localSaveDir,[subjectNames{ss} '_resultsFiles'],[subjectNames{ss} '_mtSinai_results.mat']);
            load(filePath,'results')
            

            % Grab the stimLabels
            stimLabels = results.model.opts{find(strcmp(results.model.opts,'stimLabels'))+1};

            % Find the vertices that we wish to analyze (V1 for
            % eccentricity)
            areaIdx = (vArea==1) .* (eccenMap > eccenRange(1)) .* (eccenMap < eccenRange(2));            
            goodIdx = logical( (results.R2 > r2Thresh) .* areaIdx );
            
        % Loop through the frequencies and obtain the set of values
        vals = cell(1,nFreqs);
        for ff = 1:nFreqs
            subString = sprintf(['f%dHz_' directions{dd}],freqs(ff));
            idx = find(contains(stimLabels,subString));
            vals{ff} = mean(results.params(goodIdx,idx),'omitnan');
        end
        

            % Prepare to plot into this subplot
            figure(ss)
            subplot(3,length(eccenDivs)-1,ee+(dd-1)*(length(eccenDivs)-1));

            % Adjust the values for the zero frequency, obtain bootstrapped
            % medians and 95% CI
            mBoot = NaN*ones(nFreqs-1,1);
            lBoot = NaN*ones(nFreqs-1,1);
            uBoot = NaN*ones(nFreqs-1,1);
            for ff = 2:nFreqs
                data{ss,dd,ee,ff-1} = vals{ff}-vals{1};
                temp_data = vals{ff}-vals{1};
                bootVals = sort(bootstrp(1000,@mean,temp_data));
                BootVals(ff-1,:) = bootVals;
                mBoot(ff-1,:) = bootVals(500);
                lBoot(ff-1,:) = bootVals(25);
                uBoot(ff-1,:) = bootVals(975);

                semilogx(zeros(1,length(data{ss,dd,ee,ff-1}))+freqs(ff),data{ss,dd,ee,ff-1},'.','Color',[0.9 0.9 0.9]);
                hold on
            end
            Boot_Vals(ss,dd,ee,:,:) = BootVals;

            %fit bootstrapped values with R2 values
            for xx = 1:1000
                yy = BootVals(:,xx)';
                [wFit,yFit,yFit2,pBoot(ss,dd,ee,xx,:)] = fitExp(w,yy);
                maxBoot(ss,dd,ee,xx,:) = max(yFit);
               if max(yFit)<0
                   temp = wFit(find(yFit-min(yFit)>max(yFit-min(yFit))*0.95));
                else
                   temp = wFit(find(yFit>max(yFit)*0.95));
               end
               
               if length(temp)>1
                    temp = temp(round(length(temp)/2));
               end
                peakFreqBoot(ss,dd,ee,xx,:) = temp;
                R = corrcoef(yy,yFit2);
                r2Boot(ss,dd,ee,xx) = R(1,2)^2;
                yFitBoot(ss,dd,ee,xx,:) = yFit;
                
%                 figure(46)
%                 hold on
%                 plot(w,yy,'ok')
%                 plot(wFit,yFit)
%                 title(sprintf('max response = %.3g, peak frequency = %.3g',[max(yFit) temp(1)]))
%                 xlabel(sprintf('p1 = %.3g, p2 = %.3g, p3 = %.3g',squeeze(squeeze(squeeze(pBoot(ss,dd,ee,xx,:))))))
%                 ax=gca; ax.TickDir='out'; ax.Box ='off'; ax.XScale = 'log';
%                 pause
%                 clf
            end
            
            % Plot bootstrapped mean across the frequencies
            errorbar(freqs(2:end),mBoot,abs(diff([mBoot lBoot],1,2)),abs(diff([uBoot mBoot],1,2)),'ob','MarkerSize',8,'MarkerFaceColor','b','LineWidth',2)

            % Obtain difference of exponentials fit for bootstrapped mean and plot it 
            Y = mBoot'; % select bootstrapped mean to be fit to the model
    
            [wFit,yFit,yFit2,p] = fitExp(w,Y);

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
            
            P(ss,dd,ee,:) = p;
            clear p
        end
    end
end

save_file = [localSaveDir '/expFitsEcc.mat'];
save(save_file,'data','eccenDivs','Boot_Vals','peakFreqBoot','maxBoot','pBoot','r2Boot','yFitBoot');

%% Local functions
function [wFit,yFit,yFit2,p] = fitExp(w,Y)
            
            % TTF model guess
            p0 = [0.5 4 1];
            lb = [0 0 0]; 
            ub = [10 10 2];
            
            % set minimum to 0
            scaledY = Y-min(Y);
            
            % Find the maximum interpolated VEP response
            wDelta = min(diff(log10(w))); % Create a scaled-up, log-spaced, version of the frequency domain
            upScale = 10;
            wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);
    
            % Scale the Y vector so that the max is 1
            scaledY=scaledY./max(scaledY);
            myObj=@(p)sqrt(sum((scaledY-watsonTemporalModelvep(w,p)).^2));
            
            p = fmincon(myObj,p0,[],[],[],[],lb,ub);
            
            % calculate model fit and undo scaling
            yFit = (watsonTemporalModelvep(wFit,p).*(max(Y)-min(Y)))+min(Y);
            
            yFit(~isfinite(yFit))=nan;
            
            yFit2 = (watsonTemporalModelvep(w,p).*(max(Y)-min(Y)))+min(Y);
end