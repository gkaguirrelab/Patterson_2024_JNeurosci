% Determines max of TTF fits for LGN and V1 from Watson fit and Spline for MRI flicker paper

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
                freqsFit = freqsFit(1:50);
                
                % runs model, and gets fitted response
                y = yBoot(:,500)'-min(yBoot(:,500));
                
                myObj = @(P) norm(y - watsonTTF(P,w));
                p0(1,1) = max(yBoot(:,500));
                P = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
                watFit = watsonTTF(P,freqsFit);
                watFit(~isfinite(watFit))=nan;
                max_wat = max(watFit);
                max_wat_freq = freqsFit(watFit==max_wat);

                watFit = watFit+(min(yBoot(:,500)));
                
                splineFit = pchip(w,y,freqsFit)+(min(yBoot(:,500)));
                max_spline = max(splineFit);
                max_spline_freq = freqsFit(splineFit==max_spline);

                % Add this fit to the plot
                subplot(3,length(areaLabels),aa+((dd-1)*(length(areaLabels))))
                hold on
                XX = [freqs fliplr(freqs)];
                YY = [yBoot(:,25)' fliplr(yBoot(:,975)')];
                fill(XX,YY,[0.9 0.9 0.9],'LineStyle','none')
                plot(freqs,yBoot(:,500),'.k')
                plot(freqsFit,watFit,'-r')
                plot(freqsFit,splineFit,'-b')
                plot([1 64],[0 0],':k','LineWidth',1)
                ylim([-1 7]);
                xlim([0.5 70])
                set(gca,'Box','off')
                set(gca,'XScale','log')
                set(gca,'TickDir','out')
                axis off

                % Report the parameters
               strSpline = {sprintf('[Spline: %2.2f, %2.1f]',[max_spline max_spline_freq])};
               text(1,6,strSpline);
               
               strWat = {sprintf('[Watson: %2.2f, %2.1f]',[max_wat max_wat_freq])};
               text(1,7,strWat);

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
