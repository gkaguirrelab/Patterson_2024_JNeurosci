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

maxBoot = sort(maxBoot,5);
peakFreqBoot = sort(peakFreqBoot,5);

%% Loop through the subjects and directions to create TSFs with fits
% for ss = 1:2
%         for ee = 1:length(eccenDivs)
%             figHandle0 = figure();
%             for dd = 1:3
%                     for pp = 1:size(polarDivs,1)
% 
% 
%                         yBoot = squeeze(squeeze(squeeze(squeeze(Boot_Vals(ss,dd,ee,pp,:,:)))));
% 
%                         [wFit,yFit,yFit2,p] = fitExp(w,yBoot(:,500)');
% 
%                         % Add this fit to the plot
%                         subplot(3,length(polarDivs),pp+((dd-1)*(length(polarDivs))))
%                         hold on
%                         XX = [freqs fliplr(freqs)];
%                         YY = [yBoot(:,25)' fliplr(yBoot(:,975)')];
%                         fill(XX,YY,[0.9 0.9 0.9],'LineStyle','none')
%                         plot(freqs,yBoot(:,500),'.k')
%                         plot(wFit,yFit,'-r')
%                         plot([1 64],[0 0],':k','LineWidth',1)
%                         ylim([-1 7]);
%                         xlim([0.5 70])
%                         set(gca,'Box','off')
%                         set(gca,'XScale','log')
%                         set(gca,'TickDir','out')
%                         axis off
% 
%                         % Report the parameters
%                         temp = wFit(find(yFit==max(yFit))); 
%                         temp = temp(end);
%                        str = {sprintf('[%2.1f, %2.0f]',[max(yFit) temp])};
%                        text(1,6,str);
% 
%                         % Add some chart stuff
%                         if dd==1
%                             ax = gca;
%                             ax.Title.Visible = 'on';
%                             str = {sprintf('%2.1f-%2.1f°',polarDivs(pp,:))};
%                             title(str);
%                         end
%                         if ee==1
%                             ax = gca;
%                             ax.YLabel.Visible = 'on';
%                             ylabel(analysisLabels{dd});
%                             semilogx([1 1],[0 7],'-k','LineWidth',1)
%                         end
%                     polar_tick_labels(pp) = {sprintf('%2.1f-%2.1f°',polarDivs(pp,:))};
%                     end
%                 
%             end
%         
%         end
% end

%% Create an average params plot
colors = {'r','b','k'};
linetype = {'-','--'};
Y_lim = [0 5.5; 0 1; 0 0.03];
paramNames = {'maximum response','peak frequency'};
P1 = 1:2;
P2 = 3:4;
xx = 1:length(eccenDivs)-1;
ft = fittype('smoothingspline');

polar_tick_labels2 = {'0','45','90','135','180','-135','-90','-45'};

polarTics = [0 45 90 135 180 225 270 315];

polarTics2 = 0:1:359;

figHandle = figure();
for ee = 1:length(eccenDivs)
    for pp = 1:2
        for ss=1:2
            for dd = 3

                switch pp
                    case 1
                        subplot(2,2,P1(ss))
                        for ff = 1:length(polarTics)
                            sortVals = sort(squeeze(squeeze(maxBoot(ss,dd,ee,ff,:))));
                            meanVal(ff) = mean(sortVals);
                            lowVal(ff) = sortVals(round(length(sortVals)*0.025));
                            hiVal(ff) = sortVals(round(length(sortVals)*0.975));
                            rangeVal(ff) = hiVal(ff) - lowVal(ff);
                        end
                        polarplot(deg2rad(polarTics),meanVal,[linetype{ee} 'o' colors{dd}])
                        hold on
                        set(gca,'RLim',[0.1 10])
                        set(gca,'ThetaTick',polarTics)
                        set(gca,'ThetaTickLabels',polar_tick_labels2)
                        set(gca,'Box','off')
                        set(gca,'ThetaZeroLocation','top')

                    case 2                  
                        subplot(2,2,P2(ss))
                        
                        for ff = 1:length(polarTics)
                            sortVals = sort(squeeze(squeeze(peakFreqBoot(ss,dd,ee,ff,:))));
                            meanVal(ff) = mean(sortVals);
                            lowVal(ff) = sortVals(round(length(sortVals)*0.025));
                            hiVal(ff) = sortVals(round(length(sortVals)*0.975));
                            rangeVal(ff) = hiVal(ff) - lowVal(ff);
                        end
                        rangeVal(rangeVal==0) = 1;       
                        polarplot(deg2rad(polarTics),meanVal,[linetype{ee} 'o' colors{dd}])
                        hold on
                        set(gca,'RLim',[2 25])
                        set(gca,'ThetaTick',polarTics)
                        set(gca,'ThetaTickLabels',polar_tick_labels2)
                        set(gca,'Box','off')
                        set(gca,'ThetaZeroLocation','top')
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
