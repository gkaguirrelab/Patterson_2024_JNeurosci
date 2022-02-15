% Plots TTF figures for MRI flicker paper across area
% the variable data is baseline subtracted

localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');
open_file = [localSaveDir '/expFitsArea.mat'];
load(open_file)

freqs = [2,4,8,16,32,64];
w = freqs;
nFreqs = length(freqs);

maxBoot = sort(maxBoot,4);
peakFreqBoot = sort(peakFreqBoot,4);

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'L-M','S','LF'};
plotColors = {'r','b','k'};
plotColors2 = [1 0.9 0.9;0.9 0.9 1;0.9 0.9 0.9];

% Loop through the subjects and directions
for ss = 1:2
    
    % Create a figure
    figHandle0 = figure();
        
        y = [];
        for aa = 1:length(areaLabels)
            for dd = 1:3
 
                for ff = 1:nFreqs
                    vals(ss,dd,aa,ff) = nanmean(data{ss,dd,aa,ff});
                    bootVals(ss,dd,aa,ff,:) = sort(bootstrp(1000,@nanmean,data{ss,dd,aa,ff}));
                end
            
                yVals = squeeze(squeeze(vals(ss,dd,aa,:)));
                yBoot = squeeze(squeeze(squeeze(bootVals(ss,dd,aa,:,:))));
                
                [wFit,yFit,yFit2,p] = fitExp(w,yBoot(:,500)');
                
                % Add this fit to the plot
                subplot(3,length(areaLabels),aa+((dd-1)*(length(areaLabels))))
                hold on
                XX = [freqs fliplr(freqs)];
                YY = [yBoot(:,25)' fliplr(yBoot(:,975)')];
                fill(XX,YY,plotColors2(dd,:),'LineStyle','none')
                plot(freqs,yBoot(:,500),['.' plotColors{dd}])
                plot(wFit,yFit,plotColors{dd})
                plot([1 64],[0 0],':k','LineWidth',1)
                ylim([-1 7]);
                xlim([0.5 70])
                set(gca,'Box','off')
                set(gca,'XScale','log')
                set(gca,'TickDir','out')
                axis off

                % Report the parameters
                temp = wFit(find(yFit==max(yFit))); 
                temp = temp(end);
               str = {sprintf('[%2.1f, %2.0f]',[max(yFit) temp])};
               text(1,6,str);

                % Add some chart stuff
                if dd==1
                    ax = gca;
                    ax.Title.Visible = 'on';
                    title(areaLabels(aa));
                end
                if aa==1
                    ax = gca;
                    ax.YLabel.Visible = 'on';
                    ylabel(analysisLabels{dd});
                    semilogx([1 1],[0 7],'-k','LineWidth',1)
                end
                
            end

        end
end

% Create an average peak frequency and maximum response plot
Xtic = 1:length(areaLabels);
Xtic_shift = [-0.2 0 0.2];
paramNames = {'maximum response','peak frequency'};


for pp = 1:2
    figure
    for aa=1:length(areaLabels)
        for ss=1:2
            for dd=1:3
                subplot(2,1,ss)
                hold on
                xx = aa+Xtic_shift(dd);
                switch pp
                    case 1
                        yy = squeeze(squeeze(maxBoot(ss,dd,aa,:)));
                        set(gca,'YLim',[0.1 10])
                        set(gca,'YScale','log')
                        set(gca,'YTick',[0.2:0.2:1 2:2:10])


                    case 2
                        yy = squeeze(squeeze(peakFreqBoot(ss,dd,aa,:)));
                        set(gca,'YLim',[12 32])
                        set(gca,'YScale','linear')
                        set(gca,'YTick',15:5:30)
                end


                errorbar(xx,yy(500),abs(diff(yy([25 500]))),abs(diff(yy([500 975]))),'o','Color',plotColors{dd},'MarkerFaceColor',plotColors{dd},'Linewidth',1.5)

                xlabel('area');
                ylabel(paramNames{pp});
                title(shortNames{ss})
                set(gca,'Box','off')
                set(gca,'TickDir','out')
                set(gca,'XLim',[0 length(areaLabels)+1])
                set(gca,'XTick',Xtic)
                set(gca,'XTickLabel',areaLabels);
            end
        end
        
    end
end


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
            upScale = 500;
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