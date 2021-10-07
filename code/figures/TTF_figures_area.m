% Plots TTF figures for MRI flicker paper across area
% the variable data is baseline subtracted

load param_by_area
freqs = [2,4,8,16,32,64];
nFreqs = length(freqs);

% These variables define the subject names, stimulus directions, and the
% analysis IDs
subjectNames = {'HEROgka1','HEROasb1'};
shortNames = {'gka','asb'};
analysisLabels = {'L-M','S','LF'};

% Loop through the subjects and directions
for ss = 1:2
    
    % Create a figure
    figHandle0 = figure();
    
    dkl = [];
    rgb=[];
        
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

                % Frequency components for model fitting
                deltaF10 = min(diff(log10(freqs(2:end))));
                fitScaleUp = 10;
                freqsFit = 10.^(log10(min(freqs(2:end)))-deltaF10+deltaF10/fitScaleUp:deltaF10/fitScaleUp:log10(max(freqs(2:end)))+deltaF10);


                % Get the fitted response
                myFit = watsonTTF2param(p,freqsFit);
                myFit(~isfinite(myFit))=nan;
                
                y(dd) = p(1);
                y2(dd) = p(2);

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
               str = {sprintf('[%2.1f, %2.3f]',p(1:2))};
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

                % Normalize the set of responses across the three channels to have
                
            end

        end
end

% Create an average params plot
plotColors = {'r','b','k'};
plotColors2 = [1 0.9 0.9 ; 0.9 0.9 1; 0.9 0.9 0.9];
Y_lim = [-3 5 ; 0 0.03];
area_tic = 1:length(areaLabels);
paramNames = {'amplitude [%change]','time constant [secs]'};
for pp = 1:2
    figHandle = figure();
    for ss=1:2
        subplot(2,1,ss)
        hold on
        for dd=1:3
            p_val =  squeeze(squeeze(squeeze(p_Boot(ss,dd,:,pp,:))));
            Rsq = squeeze(squeeze(squeeze(R_squared(ss,dd,:,500))));
            Rsq_lo = squeeze(squeeze(squeeze(R_squared(ss,dd,:,25))));
            yy = p_val((Rsq>0.5 & Rsq_lo>0.2),:);
            xx = area_tic;
            xx = xx(Rsq>0.5 & Rsq_lo>0.2);
           
            errorbar(xx,yy(:,500),abs(diff(yy(:,[25 500]),[],2)),abs(diff(yy(:,[500 975]),[],2)),'o','Color',plotColors{dd},'MarkerFaceColor',plotColors{dd},'Linewidth',1.5)
        end
        title(subjectNames{ss});
        xlabel('Visual area');
        ylabel(paramNames{pp});
        set(gca,'Box','off')
        set(gca,'TickDir','out')
        set(gca,'XLim',[0 5])
        set(gca,'YLim',Y_lim(pp,:))
        set(gca,'XTick',area_tic)
        set(gca,'XTickLabel',areaLabels);
        
    end
end



%% statistical testing

subject = repmat([1;2],12,1);
channel = repmat([1 1 2 2 3 3]',4,1);
area = [ones(6,1);2*ones(6,1);3*ones(6,1);4*ones(6,1)];

channel2 = cell(24,1);
channel2(channel==1) = {'LM'}; channel2(channel==2) = {'S'}; channel2(channel==3) = {'LMS'};
channel2 = categorical(channel2);

subject2 = cell(24,1);
subject2(subject==1) = {'GKA'}; subject2(subject==2) = {'ASB'};
subject2 = categorical(subject2);

area2 = cell(24,1); area2 (area==1) = areaLabels(1); area2 (area==2) = areaLabels(2); area2 (area==3) = areaLabels(3); area2 (area==4) = areaLabels(4);

gain = reshape(squeeze(squeeze(p_Boot(:,:,:,1,500))),2,12,1);
gain = reshape(gain,24,1);
timeConstant = reshape(squeeze(squeeze(p_Boot(:,:,:,2,500))),2,12,1);
timeConstant = reshape(timeConstant,24,1);
R_sqr = reshape(squeeze(R_squared(:,:,:,500)),2,12,1);
R_sqr = reshape(R_sqr,24,1);

R_sqrL = reshape(squeeze(R_squared(:,:,:,25)),2,12,1);
R_sqrL = reshape(R_sqrL,24,1);

Tbl = table(subject2,channel2,area2,gain,timeConstant,R_sqr,R_sqrL,'VariableNames',{'Subject','Channel','Area','Gain','TimeConstant','R_squared','R_squared25'});
Tbl = Tbl(Tbl.R_squared>0.5 & Tbl.R_squared25>0.2,:);

[p,tbl,statsT] = anovan(Tbl.TimeConstant,{Tbl.Channel,Tbl.Area},'model','interaction','VarNames',{'Channel','Area'});
multcompare(statsT,'Dimension',[1 2]);

[p,tbl,statsG] = anovan(Tbl.Gain,{Tbl.Channel,Tbl.Area},'model','interaction','VarNames',{'Channel','Area'});
multcompare(statsG,'Dimension',[1 2]);

