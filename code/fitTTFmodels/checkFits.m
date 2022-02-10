% Plots TTF figures for MRI flicker paper
% the variable data is baseline subtracted

load expFitsEcc Boot_Vals eccenDivs
w = [2,4,8,16,32,64];

for ss = 1
    for dd = 3
        for ee = 5
            for xx = 1:1000
                yy = squeeze(squeeze(squeeze(Boot_Vals(ss,dd,ee,:,xx))))';
                [wFit,yFit,yFit2,P] = fitExp(w,yy);

                peak = wFit(find(yFit==max(yFit)));
                R = corrcoef(yy,yFit2);
                
                figure(46)
                hold on
                plot(w,yy,'ok')
                plot(wFit,yFit)
                title(sprintf('max response = %.3g, peak frequency = %.3g',[max(yFit) peak(end)]))
                xlabel(sprintf('boot %.3g: p1 = %.3g, p2 = %.3g, p3 = %.3g, p4 = %.3g, r2 = %.3g',[xx P R(1,2)^2]))
                ax=gca; ax.TickDir='out'; ax.Box ='off'; ax.XScale = 'log';
                pause
                clf
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