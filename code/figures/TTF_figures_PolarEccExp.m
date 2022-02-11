% Plots TTF figures for MRI flicker paper
% the variable data is baseline subtracted

localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');
open_file = [localSaveDir '/expFitsPolarEcc.mat'];
load(open_file)

freqs = [2,4,8,16,32,64];
w = freqs;
nFreqs = length(freqs);

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'L-M','S','LF'};

maxBoot = sort(maxBoot,4);
peakFreqBoot = sort(peakFreqBoot,4);

% Loop through the subjects and directions
for ss = 1:2
        for ee = 1:length(eccenDivs)
            figHandle0 = figure();
            for dd = 1:3
                    for pp = 1:size(polarDivs,1)


                        for ff = 1:nFreqs
                            vals(ss,dd,ee,ff) = nanmean(data{ss,dd,ee,pp,ff});
                            bootVals(ss,dd,ee,ff,:) = sort(bootstrp(1000,@nanmean,data{ss,dd,ee,pp,ff}));
                        end


                        yVals = squeeze(squeeze(squeeze(vals(ss,dd,ee,pp,:))));
                        yBoot = squeeze(squeeze(squeeze(squeeze(bootVals(ss,dd,ee,pp,:,:)))));

                        [wFit,yFit,yFit2,p] = fitExp(w,yBoot(:,500)');

                        % Add this fit to the plot
                        subplot(3,length(polarDivs),pp+((dd-1)*(length(polarDivs))))
                        hold on
                        XX = [freqs fliplr(freqs)];
                        YY = [yBoot(:,25)' fliplr(yBoot(:,975)')];
                        fill(XX,YY,[0.9 0.9 0.9],'LineStyle','none')
                        plot(freqs,yBoot(:,500),'.k')
                        plot(wFit,yFit,'-r')
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
                            str = {sprintf('%2.1f-%2.1f°',polarDivs(pp,:))};
                            title(str);
                        end
                        if ee==1
                            ax = gca;
                            ax.YLabel.Visible = 'on';
                            ylabel(analysisLabels{dd});
                            semilogx([1 1],[0 7],'-k','LineWidth',1)
                        end

                        polar_tick_labels(pp) = {sprintf('%2.1f-%2.1f°',polarDivs(pp,:))};
                end

                
            end
        
        end
end

% Create an average params plot
colors = {'r','b','k'};
Y_lim = [0 5.5; 0 1; 0 0.03];
paramNames = {'maximum response','peak frequency'};
P1 = 1:2;
P2 = 3:4;
xx = 1:length(eccenDivs)-1;
ft = fittype('smoothingspline');

figHandle = figure();
for ee = 1:length(eccenDivs)
    for pp = 1:2
        for ss=1:2
            for dd = 1:3

                switch pp
                    case 1
                        subplot(2,2,P1(ss))
                        hold on
                        set(gca,'TickDir','out')
                        set(gca,'Box','off')
                        for ff = 1:8
                            sortVals = sort(squeeze(squeeze(maxBoot(ss,dd,ee,ff,:))));
                            meanVal(ff) = mean(sortVals);
                            lowVal(ff) = sortVals(round(length(sortVals)*0.025));
                            hiVal(ff) = sortVals(round(length(sortVals)*0.975));
                            rangeVal(ff) = hiVal(ff) - lowVal(ff);
                        end
                        semilogy(meanVal,['o' colors{dd}])
                        hold on
                        for ff = 1:8
                            semilogy([ff, ff],[lowVal(ff), hiVal(ff)],['-' colors{dd}],'LineWidth',2)
                        end
                        f = fit([1:8]',log10(meanVal'),ft,'Weights',1./rangeVal);
                        semilogy(1:0.1:8,10.^(f(1:0.1:8)),['-' colors{dd}])
                        set(gca,'YLim',[0.1 10])
                        set(gca,'Yscale','log')

                    case 2                  
                        subplot(2,2,P2(ss))
                        hold on
                        for ff = 1:8
                            sortVals = sort(squeeze(squeeze(peakFreqBoot(ss,dd,ee,ff,:))));
                            meanVal(ff) = mean(sortVals);
                            lowVal(ff) = sortVals(round(length(sortVals)*0.025));
                            hiVal(ff) = sortVals(round(length(sortVals)*0.975));
                            rangeVal(ff) = hiVal(ff) - lowVal(ff);
                        end
                        rangeVal(rangeVal==0) = 1;       
                        plot(meanVal,['o' colors{dd}])
                        for ff = 1:8
                            plot([ff, ff],[lowVal(ff), hiVal(ff)],['-' colors{dd}],'LineWidth',2)
                        end
                        f = fit([1:8]',meanVal',ft,'Weights',1./rangeVal);
                        plot(1:0.1:8,f(1:0.1:8),['-' colors{dd}])
                        set(gca,'YLim',[2 25])
                        set(gca,'TickDir','out')
                end
            end 
        end
    end
end




%% Local functions
function [wFit,yFit,yFit2,p] = fitExp(w,Y)
            
            % TTF model guess
            p0 = [0.5 4 1 0];
            lb = [0 0 0 0]; 
            ub = [10 10 2 0];
            
            % set minimum to 0
            scaledY = Y-min(Y);
            
            % Find the maximum interpolated VEP response
            wDelta = min(diff(log10(w))); % Create a scaled-up, log-spaced, version of the frequency domain
            upScale = 500;
            wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);
    
            % Scale the Y vector so that the max is 1
            scaledY=scaledY./max(scaledY);
            myObj=@(p)sqrt(sum((scaledY-watsonTemporalModelTimeScale(w,p)).^2));
            
            p = fmincon(myObj,p0,[],[],[],[],lb,ub);
            
            % calculate model fit and undo scaling
            yFit = (watsonTemporalModelTimeScale(wFit,p).*(max(Y)-min(Y)))+min(Y);
            
            yFit(~isfinite(yFit))=nan;
            
            yFit2 = (watsonTemporalModelTimeScale(w,p).*(max(Y)-min(Y)))+min(Y);
end
