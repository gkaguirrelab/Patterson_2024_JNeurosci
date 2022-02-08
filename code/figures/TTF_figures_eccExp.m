% Plots TTF figures for MRI flicker paper
% the variable data is baseline subtracted

load expFitsEcc
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

for ss = 1:2
    for dd = 1:3
        for ee = 1:length(eccenDivs)-1
            for xx = 1:1000
                Y = squeeze(squeeze(squeeze(Boot_Vals(ss,dd,ee,:,xx))))';
                P = pBoot(ss,dd,ee,xx,:);
                yFit = (watsonTemporalModelvep(wFit,P).*(max(Y)-min(Y)))+min(Y);
                maxBoot(ss,dd,ee,xx,:) = max(yFit);

                temp = wFit(find(yFit==max(yFit)));
                peakFreqBoot(ss,dd,ee,xx,:) = temp(end);
            end
        end
    end
end

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
                temp = wFit(find(yFit>max(yFit)*0.95)); 
                temp = temp(round(length(temp)/2));
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
plotColors = {'r','b','k'};
plotColors2 = [1 0.9 0.9 ; 0.9 0.9 1; 0.9 0.9 0.9];
Y_lim = [0 5.5; 0 1; 0 0.03];
paramNames = {'maximum response','peak frequency'};
P1 = 1:2;
P2 = 3:4;
xx = 1:length(eccenDivs)-1;

figHandle = figure();
for pp = 1:2
    for ss=1:2
        for dd = 1:3
           
            

            switch pp
                case 1
                    yy = squeeze(squeeze(maxBoot(ss,dd,:,:)));
                    subplot(2,2,P1(ss))
                    hold on
                    set(gca,'YLim',[0.1 10])
                    set(gca,'YScale','log')
                    set(gca,'YTick',[0.2:0.2:1 2:2:10])
                    
                case 2
                    yy = squeeze(squeeze(peakFreqBoot(ss,dd,:,:)));                   
                    subplot(2,2,P2(ss))
                    hold on
                    set(gca,'YLim',[2 25])
                    set(gca,'YScale','linear')
                    set(gca,'YTick',5:5:25)
            end
            meanVal = yy(:,500);
            rangeVal = abs(diff(yy(:,[25 975]),[],2));   
            f = fit(eccenDivs(2:end)',meanVal,'smoothingspline','Weights',1./rangeVal);
            xx = eccenDivs(2:end);
            yy(yy<0) = 0.001;
            xxFit = linspace(1,91,60);
            yyFit = f(xxFit);
            
            errorbar(xx,yy(:,500),abs(diff(yy(:,[25 500]),[],2)),abs(diff(yy(:,[500 975]),[],2)),'o','Color',plotColors{dd},'MarkerFaceColor',plotColors{dd},'Linewidth',1.5)
            plot(xxFit,yyFit,'-','Color',plotColors{dd});
            title(shortNames{ss});
            xlabel('Eccentricity');
            ylabel(paramNames{pp});
            set(gca,'Box','off')
            set(gca,'TickDir','out')
            set(gca,'XLim',[2.5 95])
            set(gca,'XScale','log')
            set(gca,'XTick',xx)
            set(gca,'XTickLabel',ecc_tick_labels);
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
