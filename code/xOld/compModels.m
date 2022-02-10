% Compare model parameters and R squared values of Watson models with
% different parameters

load('model2param.mat')
pBoot2 = squeeze(p_Boot); Rsquared2 = squeeze(R_squared);

load('model3param.mat')
pBoot3 = squeeze(p_Boot); Rsquared3 = squeeze(R_squared);

clear p_Boot R_squared

Color = {'r','b','k'};

figure(100)
for ss = 1:2
    switch ss
        case 1
            loc = 0;
        case 2
            loc = 3;
    end

    for dd = 1:3
        Rsquared = [squeeze(squeeze(Rsquared2(ss,dd,:,:))) squeeze(squeeze(Rsquared3(ss,dd,:,:)))];
        param_group = [zeros(2,1000) ones(2,1000)];

        [h,p] = kstest2(squeeze(squeeze(Rsquared2(ss,dd,2,:))),squeeze(squeeze(Rsquared3(ss,dd,2,:))));
        
        figure(100)
        subplot(2,3,loc+dd)
        hold on
        errorbar(1,Rsquared(1,500),abs(diff(Rsquared(1,[25 500]))),abs(diff(Rsquared(1,[500 975]))),'-o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        errorbar(2,Rsquared(1,1500),abs(diff(Rsquared(1,[1025 1500]))),abs(diff(Rsquared(1,[1500 1975]))),'-o','MarkerFaceColor',Color{dd},'Color',Color{dd},'LineWidth',1.5)
        set(gca,'TickDir','out');
        box off
        xticks(1:2)
        xticklabels({'2 params','3 params'})
        xlim([0 3])
        ylim([0 1])
        ylabel('Adj R squared')
%         title(['p = ' num2str(p)])
    end
end