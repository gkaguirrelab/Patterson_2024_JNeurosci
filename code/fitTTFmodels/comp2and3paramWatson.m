% Compare R squared values for 3 parameter and 2 parameter model across
% eccentricity

load param_by_V1eccen2param R_squared eccenDivs
R_squared2param = R_squared;
clear R_squared

load param_by_V1eccen3param R_squared
R_squared3param = R_squared;
clear R_squared

shortNames = {'gka','asb'};
analysisLabels = {'L-M','S','LMS'};

% Creates table with Subject, Channel, Eccentricity, 2 parameter Rsquared
% and 3 parameter Rsquared with 95% CI
Tbl_Rsquared = table('Size',[36,5],'VariableTypes',{'cellstr','cellstr','cellstr','cellstr','cellstr'},...
    'VariableNames',{'Subject','Channel','Eccentricity','param2Rsq','param3Rsq'});


counter = 1;
for ss = 1:2
    for d = 1:3
        for ee = 1:length(eccenDivs)-1
            Tbl_Rsquared.Subject(counter) = shortNames(ss);
            Tbl_Rsquared.Channel(counter) = analysisLabels(d);
            Tbl_Rsquared.Eccentricity(counter) = {[num2str(eccenDivs(ee)) ' - ' num2str(eccenDivs(ee+1))]};
            
            Tbl_Rsquared.param2Rsq(counter) = {sprintf('%1.2f [%1.2f, %1.2f]',squeeze(squeeze(R_squared2param(ss,d,ee,[500 25 975]))))};
            Tbl_Rsquared.param3Rsq(counter) = {sprintf('%1.2f [%1.2f, %1.2f]',squeeze(squeeze(R_squared3param(ss,d,ee,[500 25 975]))))};
            counter = counter+1;
        end
    end
end

figure
TString = evalc('disp(Tbl_Rsquared)');
TString = strrep(TString,'<strong>','\bf');
TString = strrep(TString,'</strong>','\rm');
TString = strrep(TString,'_','\_');
FixedWidth = get(0,'FixedWidthFontName');
annotation(gcf,'Textbox','String',TString,'Interpreter','Tex','FontName',FixedWidth,'Units','Normalized','Position',[0 0 1 1]);