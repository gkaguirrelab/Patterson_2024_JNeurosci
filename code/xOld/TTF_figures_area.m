% Plots TTF figures for MRI flicker paper across area
% the variable data is baseline subtracted

load param_by_area
freqs = [2,4,8,16,32,64];
w = freqs;
nFreqs = length(freqs);

% For only LGN and V1
areaLabels = areaLabels(1:2);
p_Boot = p_Boot(:,:,1:2,:,:);
p_mean = p_mean(:,:,1:2,:);
R_squared = R_squared(:,:,1:2,:);

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
subject_tic = 1:length(shortNames);
subject_tic_shift = [-0.2 0.2];
paramNames = {'amplitude','surround amplitude','time constant [secs]'};
Y_lim = [0 1.5;0.6 1.2;0 0.05];
figure
for pp = 1:length(p0)
    for aa=1:length(areaLabels)
        for dd=1:3
            subplot(length(p0),3,pp+((dd-1)*(length(p0))))
            hold on
            Rsq = squeeze(squeeze(squeeze(R_squared(:,dd,aa,500))));
            Rsq_lo = squeeze(squeeze(squeeze(R_squared(:,dd,aa,25))));
            yy = squeeze(squeeze(squeeze(p_Boot(:,dd,aa,pp,:))));
            xx = subject_tic+subject_tic_shift(aa);
            
            % Normalize amplitude by sum of squares
            if pp == 1
                lm = squeeze(squeeze(squeeze(p_Boot(:,1,aa,pp,500))));
                s = squeeze(squeeze(squeeze(p_Boot(:,2,aa,pp,500))));
                lms = squeeze(squeeze(squeeze(p_Boot(:,3,aa,pp,500))));
                norm = sqrt(lm.^2+s.^2+lms.^2);
                for YY = 1:2
                    yy(YY,:) = yy(YY,:)./norm(YY);
                end
            end

            errorbar(xx,yy(:,500),abs(diff(yy(:,[25 500]),[],2)),abs(diff(yy(:,[500 975]),[],2)),'o','Color',plotColors{dd},'MarkerFaceColor',plotColors{dd},'Linewidth',1.5)
            title(analysisLabels{dd});
            xlabel('Subject');
            ylabel(paramNames{pp});
            set(gca,'Box','off')
            set(gca,'TickDir','out')
            set(gca,'XLim',[0 length(shortNames)+1])
            set(gca,'YLim',Y_lim(pp,:))
            set(gca,'XTick',subject_tic)
            set(gca,'XTickLabel',shortNames);
        end
        
    end
end


