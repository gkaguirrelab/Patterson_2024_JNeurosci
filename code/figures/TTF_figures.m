% Plots TTF figures for MRI flicker paper
% the variable data is baseline subtracted

load param_by_V1eccen
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
        for ee = 1:length(eccenDivs)-1
            for dd = 1:3

            
                for ff = 1:nFreqs
                    vals(ss,dd,ee,ff) = nanmean(data{ss,dd,ee,ff});
                    bootVals(ss,dd,ee,ff,:) = sort(bootstrp(1000,@nanmean,data{ss,dd,ee,ff}));
                end
            
                yVals = squeeze(squeeze(vals(ss,dd,ee,:)));
                yBoot = squeeze(squeeze(squeeze(bootVals(ss,dd,ee,:,:))));
                p = squeeze(squeeze(squeeze(p_Boot(ss,dd,ee,:,500))));

                % Frequency components for model fitting
                deltaF10 = min(diff(log10(freqs(2:end))));
                fitScaleUp = 10;
                freqsFit = 10.^(log10(min(freqs(2:end)))-deltaF10+deltaF10/fitScaleUp:deltaF10/fitScaleUp:log10(max(freqs(2:end)))+deltaF10);


                % Get the fitted response
                myFit = watsonTTF2param(p,freqsFit);
                myFit(~isfinite(myFit))=nan;

                % Store the interpolated maximum response
%                 y(dd) = max(myFit);
                
                y(dd) = p(1);
                y2(dd) = p(2);

                % Add this fit to the plot
                subplot(3,length(eccenDivs)-1,ee+((dd-1)*(length(eccenDivs)-1)))
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
                    str = {sprintf('%2.1f-%2.1f°',eccenDivs(ee:ee+1))};
                    title(str);
                end
                if ee==1
                    ax = gca;
                    ax.YLabel.Visible = 'on';
                    ylabel(analysisLabels{dd});
                    semilogx([1 1],[0 7],'-k','LineWidth',1)
                end

                % Normalize the set of responses across the three channels to have
                
            end
            
            % unit magnitude
            dkl(ee,:) = y/norm(y);
            dkl2(ee,:) = y2/norm(y2);
        
        end
            % Plot the DKL surface
            figHandle1 = figure();

            % Render the DKL surface
            [X,Y,Z] = sphere(40);
            XVec = X(:); YVec = Y(:); ZVec = Z(:);
            dklSurf = [XVec, YVec, ZVec];
            rgbSurf = dklToRGB(dklSurf);
            posQuadrant = logical((XVec>=-1e-6) .* (YVec>=-1e-6) .* (ZVec>=-1e-6));
            rgbSurf(~posQuadrant,:)=nan;
            rgbSurf = reshape(rgbSurf,41,41,3);
            X(~posQuadrant)=nan;
            Y(~posQuadrant)=nan;
            Z(~posQuadrant)=nan;
            surf(X,Y,Z,rgbSurf,'EdgeColor','none','FaceColor','interp');
            hold on

            plotOnSphere(dkl(:,1),dkl(:,2),dkl(:,3),'Color','k','LineWidth',2);

            % Add circles around each location
            rVals = linspace(0.9999,0.999,length(eccenDivs)-1);
            for ii = 1:size(dkl,1)
                addCircleAtCoordinate(dkl(ii,:),rVals(ii));
            end

            axis equal
            xlim([0 1])
            ylim([0 1])
            zlim([0 1])
            zlabel('luminance');
            xlabel('red-green');
            ylabel('yellow-violet');
            view(135,45);
            
            
            
            
             % Plot the DKL surface
            figHandle1 = figure();

            % Render the DKL surface
            [X,Y,Z] = sphere(40);
            XVec = X(:); YVec = Y(:); ZVec = Z(:);
            dklSurf = [XVec, YVec, ZVec];
            rgbSurf = dklToRGB(dklSurf);
            posQuadrant = logical((XVec>=-1e-6) .* (YVec>=-1e-6) .* (ZVec>=-1e-6));
            rgbSurf(~posQuadrant,:)=nan;
            rgbSurf = reshape(rgbSurf,41,41,3);
            X(~posQuadrant)=nan;
            Y(~posQuadrant)=nan;
            Z(~posQuadrant)=nan;
            surf(X.*0.5,Y.*0.5,Z.*0.5,rgbSurf,'EdgeColor','none','FaceColor','interp');
            hold on

            plotOnSphere(dkl2(:,1),dkl2(:,2),dkl2(:,3),'Color','k','LineWidth',2);

            % Add circles around each location
            rVals = linspace(0.9999,0.999,length(eccenDivs)-1);
            for ii = 1:size(dkl2,1)
                addCircleAtCoordinate(dkl2(ii,:),rVals(ii));
            end

            axis equal
            xlim([0 1])
            ylim([0 1])
            zlim([0 1])
            zlabel('luminance');
            xlabel('red-green');
            ylabel('yellow-violet');
            view(135,45);

    
end

% Create an average params plot
plotColors = {'r','b','k'};
plotColors2 = [1 0.9 0.9 ; 0.9 0.9 1; 0.9 0.9 0.9];
Y_lim = [0 7 ; 0 0.03];
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
            xx = eccenDivs(2:end);
            xx = xx(Rsq>0.5 & Rsq_lo>0.2);
            
            XX = [xx fliplr(xx)];
            YY = [yy(:,25)' fliplr(yy(:,975)')];
            fill(XX,YY,plotColors2(dd,:),'LineStyle','none')
            plot(xx,yy(:,500),'-','Color',plotColors{dd},'Linewidth',1.5)
        end
        title(subjectNames{ss});
        xlabel('Eccentricity band');
        ylabel(paramNames{pp});
        set(gca,'Box','off')
        set(gca,'XScale','log')
        set(gca,'TickDir','out')
        set(gca,'XLim',[2.7 95])
        set(gca,'YLim',Y_lim(pp,:))
        set(gca,'XTick',xx)
        
    end
end


%% LOCAL FUNCTIONS

function h = addCircleAtCoordinate(X,r)

theta=-pi:0.1:pi;

    x = r*ones(size(theta));
    y = ((1-r^2)^0.5)*cos(theta);
    z = ((1-r^2)^0.5)*sin(theta);
    x = [x x(1)];
    y = [y y(1)];
    z = [z z(1)];
    
    h = plot3(x,y,z,'-k','LineWidth',2);
    
    [az,el] = cart2sph(X(1),X(2),X(3));

    rotate(h,[0 0 1],rad2deg(az),[0 0 0]);
    R = [cos(az) -sin(az); sin(az) cos(az)];
    direction =  [(R*[0;1])' 0];
    rotate(h,direction,-rad2deg(el),[0 0 0]);
end