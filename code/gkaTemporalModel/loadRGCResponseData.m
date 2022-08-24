function [midgetData, parasolData] = loadRGCResponseData()

%
%   Solomon SG, Lee BB, White AJ, Rüttiger L, Martin PR. Chromatic
%   organization of ganglion cell receptive fields in the peripheral
%   retina. Journal of Neuroscience. 2005 May 4;25(18):4527-39.
%
% Figures 9A and 9B. Responses were measured between 0 and 47° eccentricity
%
% The specific frequencies studied are not specified in the paper, but I
% deduced these from the tendencies of the values on the plots and some
% hints in the methods. The measurements at 0-10 degrees eccentricity did
% not include 15 Hz.
%
% The data point for midget 20-30 degrees, chromatic gain at 78 Hz was
% below the range of the figure, but the presence of a corresponding phase
% value at this frequency suggests that it exists. I assigned a value by
% extrapolating the fit.
%
% Responses were to a spatially uniform, 4.7° diameter modulations at close
% to 2000 Trolands.
%

midgetData.e0.luminance.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 20, 30, 40, 50, 60 78];
midgetData.e0.luminance.g = [0.240956222	0.344570704	0.476634602	0.92272474	1.234188352	1.102638087	1.145143666	1.036248506	0.512739655	0.419529869	0.200628096];
midgetData.e0.luminance.p = [24.79881754	23.96124158	8.769912958	-17.8190179	-70.12645755	-171.0133027	-248.7436361	-346.1980621	-449.4334045	-495.3851207	-618.5580555];
midgetData.e0.chromatic.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 20, 30, 40, 50, 60];
midgetData.e0.chromatic.g = [0.883528514	1.041249317	1.335647522	1.994882817	2.737931338	2.222914502	1.145143666	0.642829544	0.325907056	0.192572001];
midgetData.e0.chromatic.p = [-0.870422073	-7.488914436	-11.2169486	-37.82230251	-98.71900148	-219.559862	-328.7403515	-429.0688126	-532.2713089	-566.8089998];

midgetData.e20.luminance.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 15, 20, 30, 40, 50, 60 78];
midgetData.e20.luminance.g = [0.16537005	0.156372507	0.311403691	0.599490731	1.5	1.9	2.112368642	2.253121364	2.406294456	1.4449392	0.685169767	0.202168197];
midgetData.e20.luminance.p = [58.33333333	83.33333333	86.11111111	61.11111111	2.777777778	-61.11111111	-86.11111111	-194.4444444	-272.2222222	-405.5555556	-497.2222222	-616.6666667];
midgetData.e20.chromatic.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 15, 20, 30, 40, 50, 60 78];
midgetData.e20.chromatic.g = [0.674126896	0.909586301	1.24836258	1.47089674	1.539420942	1.912429335	1.908340174	1.426084262	1.201997763	0.599140868	0.104648486	0.06];
midgetData.e20.chromatic.p = [26.81818182	10.81168831	-5.12987013	-35.35714286	-85.68181818	-132.4350649	-193.1168831	-311.3636364	-403.3766234	-552.6298701	-650.2597403	-733.7012987];

midgetData.e30.luminance.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 15, 20, 30, 40, 50, 60 78];
midgetData.e30.luminance.g = [0.230838514	0.236127818	0.427498386	0.722467765	1.882228718	2.284530668	2.916885902	5.182037328	6.969143018	8.450342831	4.946604595	1.879437669];
midgetData.e30.luminance.p = [40.55080722	35.99240266	48.40139285	41.05729661	8.452041785	-17.98037354	-69.29408041	-149.2877493	-228.7749288	-373.156062	-455.2389997	-568.5660019];
midgetData.e30.chromatic.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 15, 20, 30, 40, 50, 60 78];
midgetData.e30.chromatic.g = [0.415777936	0.644266457	0.900089539	1.057382233	1.119496026	0.808479799	1.032163881	0.886480472	1.074358692	0.761510012	0.494562596	0.344497299];
midgetData.e30.chromatic.p = [23.58341247	-0.696422919	8.990186768	-29.34472934	-78.91737892	-105.3181387	-128.4900285	-219.784742	-285.1535296	-426.7173156	-539.8227287	-655.8721114];

%   Solomon SG, Martin PR, White AJ, Rüttiger L, Lee BB. Modulation
%   sensitivity of ganglion cells in peripheral retina of macaque. Vision
%   Research. 2002 Dec 1;42(27):2893-8..
%
% Figures 1A and 1B
%
% Mean response of cells from 3-10° eccentricity at 2000 Trolands

parasolData.e0.luminance.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 20, 30, 40, 50, 60 78];
parasolData.e0.luminance.g = [0.262436577	0.591149071	1.644125739	3.702102437	5.900662716	7.471381688	4.057508719	2.615399798	1.591726565	0.672262083	0.495504391];

parasolData.e20.luminance.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 15, 20, 30, 40, 50, 60 78];
parasolData.e20.luminance.g = [0.229534081	0.744221684	1.775076399	5.228815432	8.177641335	8.049716059	7.052555544	5.956457628	5.022726615	3.495853844	2.340097948	0.898244473];

parasolData.e30.luminance.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 15, 20, 30, 40, 50, 60 78];
parasolData.e30.luminance.g = [0.27808972	0.638310685	2.413000771	5.755355979	8.33404175	7.454962426	7.1900742	7.6431795	9.831543136	9.482210723	8.303539014	3.645080274];

%{
    figure
    loglog(midgetData.e0.chromatic.f,midgetData.e0.chromatic.g,'-r','LineWidth',4);
    hold on
    loglog(midgetData.e20.chromatic.f,midgetData.e20.chromatic.g,'-r','LineWidth',2);
    loglog(midgetData.e30.chromatic.f,midgetData.e30.chromatic.g,'-r','LineWidth',1);
    
    loglog(midgetData.e0.luminance.f,midgetData.e0.luminance.g,'-k','LineWidth',4);
    loglog(midgetData.e20.luminance.f,midgetData.e20.luminance.g,'-k','LineWidth',2);
    loglog(midgetData.e30.luminance.f,midgetData.e30.luminance.g,'-k','LineWidth',1);
    
    figure
    loglog(parasolData.e0.luminance.f,parasolData.e0.luminance.g,'-k','LineWidth',4);
    hold on
    loglog(parasolData.e20.luminance.f,parasolData.e20.luminance.g,'-k','LineWidth',2);
    loglog(parasolData.e30.luminance.f,parasolData.e30.luminance.g,'-k','LineWidth',1);
%}

end
