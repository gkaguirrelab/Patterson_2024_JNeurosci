% Plots TTF figures for MRI flicker paper across area
% the variable data is baseline subtracted

load param_by_area
freqs = [2,4,8,16,32,64];
w = freqs;
nFreqs = length(freqs);

% TTF model guess
p0 = [1.5 0.9 0.015]; lb = [-10 -2 0.001]; ub = [10 2 0.1];

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'L-M','S','LF'};

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
                p = squeeze(squeeze(squeeze(p_Boot(ss,dd,aa,:,500))));
                
            
                 options = optimoptions(@fmincon,... % The options for the search (mostly silence diagnostics)
                        'Diagnostics','off',...
                        'Display','off');
                
                % Frequency components for model fitting
                deltaF10 = min(diff(log10(freqs(2:end))));
                fitScaleUp = 10;
                freqsFit = 10.^(log10(min(freqs(2:end)))-deltaF10+deltaF10/fitScaleUp:deltaF10/fitScaleUp:log10(max(freqs(2:end)))+deltaF10);
                
                % runs model, and gets fitted response
                y = yBoot(:,500)'-min(yBoot(:,500));
                myObj = @(P) norm(y - watsonTTF(P,w));
                p0(1,1) = max(yBoot(:,500));
                P = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
                myFit = watsonTTF(P,freqsFit);
                myFit(~isfinite(myFit))=nan;
                paramNames = {'amplitude [%change]','surround amplitude [%change]','time constant [secs]'};
                Y_lim = [-1 5.5;0.6 1.2;0 0.05];

                myFit = myFit+(min(yBoot(:,500)));
                

                % Add this fit to the plot
                subplot(3,length(areaLabels),aa+((dd-1)*(length(areaLabels))))
                hold on
                XX = [freqs fliplr(freqs)];
                YY = [yBoot(:,25)' fliplr(yBoot(:,975)')];
                fill(XX,YY,[0.9 0.9 0.9],'LineStyle','none')
                plot(freqs,yBoot(:,500),'.k')
                plot(freqsFit,myFit,'-r')
                plot([1 64],[0 0],':k','LineWidth',1)
                ylim([-1 7]);
                xlim([0.5 70])
                set(gca,'Box','off')
                set(gca,'XScale','log')
                set(gca,'TickDir','out')
                axis off

                % Report the parameters
               str = {sprintf('[%2.1f, %2.2f, %2.3f]',p(1:end))};
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

% Create an average params plot
plotColors = {'r','b','k'};
plotColors2 = [1 0.9 0.9 ; 0.9 0.9 1; 0.9 0.9 0.9];
area_tic = 1:length(areaLabels);
for pp = 1:length(p0)
    figHandle = figure();
    for ss=1:2
        subplot(2,1,ss)
        hold on
        for dd=1:3
            p_val =  squeeze(squeeze(squeeze(p_Boot(ss,dd,:,pp,:))));
            Rsq = squeeze(squeeze(squeeze(R_squared(ss,dd,:,500))));
            Rsq_lo = squeeze(squeeze(squeeze(R_squared(ss,dd,:,25))));
            yy = p_val((Rsq>0.5),:);
            xx = area_tic;
            xx = xx(Rsq>0.5);
           
            errorbar(xx,yy(:,500),abs(diff(yy(:,[25 500]),[],2)),abs(diff(yy(:,[500 975]),[],2)),'o','Color',plotColors{dd},'MarkerFaceColor',plotColors{dd},'Linewidth',1.5)
        end
        title(subjectNames{ss});
        xlabel('Visual area');
        ylabel(paramNames{pp});
        set(gca,'Box','off')
        set(gca,'TickDir','out')
        set(gca,'XLim',[0 6])
        set(gca,'YLim',Y_lim(pp,:))
        set(gca,'XTick',area_tic)
        set(gca,'XTickLabel',areaLabels);
        
    end
end


