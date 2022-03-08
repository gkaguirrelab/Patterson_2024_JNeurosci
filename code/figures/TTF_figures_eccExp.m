% Plots TTF figures for MRI flicker paper
% the variable data is baseline subtracted

localSaveDir = getpref('mriSinaiAnalysis','localSaveDir');
open_file = [localSaveDir '/expFitsEcc.mat'];
load(open_file)

freqs = [2,4,8,16,32,64];
w = freqs;
nFreqs = length(freqs);

wDelta = min(diff(log10(w))); % Create a scaled-up, log-spaced, version of the frequency domain
upScale = 500;
wFit = 10.^(log10(min(w))-wDelta+wDelta/upScale:wDelta/upScale:log10(max(w))+wDelta);

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'L-M','S','LF'};
plotColors = {'r','b','k'};
plotColors2 = [1 0.9 0.9;0.9 0.9 1;0.9 0.9 0.9];

maxBoot = sort(maxBoot,4);
peakFreqBoot = sort(peakFreqBoot,4);

% Loop through the subjects and directions
for ss = 1:2
    
    % Create a figure
    figHandle0 = figure();
    
    dkl = [];
    rgb=[];
        
        y = [];
        for ee = 1:length(eccenDivs)-1
            for dd = 1:3

            
                for ff = 1:nFreqs
                    vals(ss,dd,ee,ff) = nanmean(data{ss,dd,ee,ff});
                    bootVals(ss,dd,ee,ff,:) = sort(bootstrp(1000,@nanmean,data{ss,dd,ee,ff}));
                end
                
                
                yVals = squeeze(squeeze(vals(ss,dd,ee,:)));
                yBoot = squeeze(squeeze(squeeze(bootVals(ss,dd,ee,:,:))));
                
                [wFit,yFit,yFit2,p] = fitExp(w,yBoot(:,500)');
              
                % Add this fit to the plot
                subplot(3,length(eccenDivs)-1,ee+((dd-1)*(length(eccenDivs)-1)))
                hold on
                XX = [freqs fliplr(freqs)];
                YY = [yBoot(:,25)' fliplr(yBoot(:,975)')];
                fill(XX,YY,plotColors2(dd,:),'LineStyle','none')
                plot(freqs,yBoot(:,500),['.' plotColors{dd}])
                plot(wFit,yFit,['-' plotColors{dd}])
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
                    str = {sprintf('%2.1f-%2.1f°',eccenDivs(ee:ee+1))};
                    title(str);
                end
                if ee==1
                    ax = gca;
                    ax.YLabel.Visible = 'on';
                    ylabel(analysisLabels{dd});
                    semilogx([1 1],[0 7],'-k','LineWidth',1)
                end
                
                ecc_tick_labels(ee) = {sprintf('%2.1f-%2.1f°',eccenDivs(ee:ee+1))};

                
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
for pp = 1:2
    for ss=1:2
        for dd = 1:3
           
            switch pp
                case 1
                    subplot(2,2,P1(ss))
                    hold on
                    set(gca,'TickDir','out')
                    set(gca,'Box','off')
                    for ff = 1:6
                        sortVals = sort(squeeze(maxBoot(ss,dd,ff,:)));
                        meanVal(ff) = mean(sortVals);
                        lowVal(ff) = sortVals(round(length(sortVals)*0.025));
                        hiVal(ff) = sortVals(round(length(sortVals)*0.975));
                        rangeVal(ff) = hiVal(ff) - lowVal(ff);
                    end
                    semilogy(meanVal,['o' colors{dd}])
                    hold on
                    for ff = 1:6
                    semilogy([ff, ff],[lowVal(ff), hiVal(ff)],['-' colors{dd}],'LineWidth',2)
                    end
                    f = fit([1:6]',log10(meanVal'),ft,'Weights',1./rangeVal);
                    semilogy(1:0.1:6,10.^(f(1:0.1:6)),['-' colors{dd}])
                    set(gca,'YLim',[0.1 10])
                    set(gca,'Yscale','log')

                case 2                  
                    subplot(2,2,P2(ss))
                    hold on
                    for ff = 1:6
                        sortVals = sort(squeeze(peakFreqBoot(ss,dd,ff,:)));
                        meanVal(ff) = mean(sortVals);
                        lowVal(ff) = sortVals(round(length(sortVals)*0.025));
                        hiVal(ff) = sortVals(round(length(sortVals)*0.975));
                        rangeVal(ff) = hiVal(ff) - lowVal(ff);
                    end
                    rangeVal(rangeVal==0) = 1;       
                    plot(meanVal,['o' colors{dd}])
                    for ff = 1:6
                        plot([ff, ff],[lowVal(ff), hiVal(ff)],['-' colors{dd}],'LineWidth',2)
                    end
                    f = fit([1:6]',meanVal',ft,'Weights',1./rangeVal);
                    plot(1:0.1:6,f(1:0.1:6),['-' colors{dd}])
                    set(gca,'YLim',[2 25])
                    set(gca,'TickDir','out')
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