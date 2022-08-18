
% Figure 1 "P cell" (midget), 7 degree eccentricity, 4600 td
%   Purpura K, Tranchina D, Kaplan E, Shapley RM. Light adaptation in the
%   primate retina: analysis of changes in gain and dynamics of monkey
%   retinal ganglion cells. Visual neuroscience. 1990 Jan;4(1):75-93.
midgetVals = [0.5255334185061452, 0.009137108418713314
1.041915028370383, 0.011912945088493911
2.067069488534802, 0.013269502988261829
4.138242167052262, 0.017529955431378554
8.362496490940943, 0.025722198606226675
16.43887918132244, 0.024799833930571988
25.351402035687936, 0.014508478481605257
32.45627863356696, 0.00859316131505289
47.850188455324336, 0.0019800923497398522];


% Figure 2 "M cell" (parasol), 4 degree eccentricity, 3600 td
%   Purpura K, Tranchina D, Kaplan E, Shapley RM. Light adaptation in the
%   primate retina: analysis of changes in gain and dynamics of monkey
%   retinal ganglion cells. Visual neuroscience. 1990 Jan;4(1):75-93.
parasolVals = [0.5300853821416367, 0.02942485367678097
1.0418195539046595, 0.041693355206741535
2.1077298093501895, 0.06842493667536828
4.225098476153642, 0.09315238576582831
8.392217275806795, 0.10380266849428775
10.411189684250672, 0.11101824423044634
12.551522527317879, 0.09213161239459439
16.529813989472647, 0.07445791988530988
25.514860302470176, 0.033906683040577416];


%xVals = [2,4,8,16,32,64];
xVals = midgetVals(:,1);
xValsFine = logspace(log10(min(xVals)),log10(max(xVals)),100);

%yVals = GKA_V1m(1,:);
yVals = midgetVals(:,2);

myAmplitudeFunc = @(f,k,fc,fcl) ...
    (1./sqrt(f.^2+fc^2)).^4 .* ...  % Four early stages
    (sqrt(f.^2 + ((1-k)*fc)^2)./sqrt(f.^2+fc^2)).^2 .* ...  % Inhibitory feedback
    (1./sqrt(f.^2+fcl^2)).^2;   % Two late stages with fixed corner frequency

mYVals = @(x,p) p(1) .* (myAmplitudeFunc(x,p(2),p(3),30)./max(myAmplitudeFunc(xValsFine,p(2),p(3),30)));

myObj = @(p) norm(yVals - mYVals(xVals,p));

p = fmincon(myObj,[1,0.8,20,2],[],[],[],[],[0 0.6 10],[4,0.9,60]);

semilogx(xValsFine,mYVals(xValsFine,p),'-k');
hold on
semilogx(xVals,yVals,'ok')