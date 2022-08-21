function [midgetData, parasolData] = rgcFlickerResponseData()

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
% Responses are to spatially uniform, wide field modulations at close to
% 2000 Trolands.
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
midgetData.e20.chromatic.p = [55.55555556	83.33333333	86.11111111	63.88888889	2.777777778	-58.33333333	-83.33333333	-191.6666667	-275	-402.7777778	-497.2222222	-613.8888889];

midgetData.e30.luminance.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 15, 20, 30, 40, 50, 60 78];
midgetData.e30.luminance.g = [0.230838514	0.236127818	0.427498386	0.722467765	1.882228718	2.284530668	2.916885902	5.182037328	6.969143018	8.450342831	4.946604595	1.879437669];
midgetData.e30.luminance.p = [40.55080722	35.99240266	48.40139285	41.05729661	8.452041785	-17.98037354	-69.29408041	-149.2877493	-228.7749288	-373.156062	-455.2389997	-568.5660019];
midgetData.e30.chromatic.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 15, 20, 30, 40, 50, 60 78];
midgetData.e30.chromatic.g = [0.415777936	0.644266457	0.900089539	1.057382233	1.119496026	0.808479799	1.032163881	0.886480472	1.074358692	0.761510012	0.494562596	0.344497299];
midgetData.e30.chromatic.p = [23.58341247	-0.696422919	8.990186768	-29.34472934	-78.91737892	-105.3181387	-128.4900285	-219.784742	-285.1535296	-426.7173156	-539.8227287	-655.8721114];

%   Smith VC, Pokorny J, Lee BB, Dacey DM. Sequential processing in vision:
%   The interaction of sensitivity regulation and temporal dynamics. Vision
%   Research. 2008 Nov 1;48(26):2649-56.
%
% Figures 3B and 4
%
% Mean response of cells from 3-10° eccentricity at 2000 Trolands

midgetData.e3.chromatic.f = [0.609375, 0.8618, 1.2188, 1.625, 2.4375, 3.25, 4.8750, 6.57, 9.7500, 13.28, 20, 26, 39, 53];
midgetData.e3.chromatic.g = 10.^[-0.94	-0.93	-0.86	-0.89	-0.8	-0.74	-0.66	-0.63	-0.74	-0.52	-0.68	-0.74	-1.24	-1.59];
midgetData.e3.chromatic.p = rad2deg([0.163729128	0.243970315	0.148886827	0.050092764	0.027365492	-0.070037106	-0.343692022	-0.798701299	-1.572820037	-2.45593692	-3.979128015	-5.684601113	-9.1716141	-11.91280148]);
midgetData.e3.luminance.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 20, 39, 53];
midgetData.e3.luminance.g = 10.^[-1.74	-1.73	-1.55	-1.33	-1	-0.84	-1.02	-1.4];
midgetData.e3.luminance.p = rad2deg([0.44851577	0.327458256	0.313543599	0.156307978	-1.250463822	-3.586734694	-6.243042672	-8.269944341]);

parasolData.e3.luminance.f = [0.609375, 0.8618, 1.2188, 1.625, 2.4375, 3.25, 4.8750, 6.57, 9.7500, 13.28, 20, 26, 39, 53, 78];
parasolData.e3.luminance.g = 10.^[-1.526304534	-1.393926433	-1.189478186	-1.033219276	-0.781151982	-0.720202452	-0.468135158	-0.430852581	-0.357356715	-0.332121471	-0.198959224	-0.150057029	-0.314442543	-0.527373824	-1.537140006];
parasolData.e3.luminance.p = rad2deg([1.913043478	1.8	1.782608696	1.7	1.6	1.47826087	0.782608696	0.130434783	-0.47826087	-1.086956522	-2.260869565	-3.434782609	-5.956521739	-8.47826087	-12.73913043]);

% figure
% loglog(midget.e0.chromatic.f,midget.e0.chromatic.g,'-r','LineWidth',4);
% hold on
% loglog(midget.e20.chromatic.f,midget.e20.chromatic.g,'-r','LineWidth',2);
% loglog(midget.e30.chromatic.f,midget.e30.chromatic.g,'-r','LineWidth',1);
% 
% loglog(midget.e0.luminance.f,midget.e0.luminance.g,'-k','LineWidth',4);
% loglog(midget.e20.luminance.f,midget.e20.luminance.g,'-k','LineWidth',2);
% loglog(midget.e30.luminance.f,midget.e30.luminance.g,'-k','LineWidth',1);
% 
% figure
% loglog(parasol.e3.luminance.f,parasol.e3.luminance.g,'-b','LineWidth',3);
% hold on
% loglog(midget.e0.chromatic.f,midget.e0.chromatic.g,'-r','LineWidth',4);
% loglog(midget.e0.luminance.f,midget.e0.luminance.g,'-k','LineWidth',4);

end
