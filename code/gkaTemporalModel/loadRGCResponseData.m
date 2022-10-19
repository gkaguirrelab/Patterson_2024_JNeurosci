function rgcData = loadRGCResponseData()

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

% Note that there is no measurement at 15 Hz in this first eccentricity
% range.
rgcData.midget.e0.LMS.f = [0.625, 1.25, 2.5, 5, 10, 20, 30, 40, 50, 60, 80];
rgcData.midget.e0.LMS.g = [0.240956222	0.344570704	0.476634602	0.92272474	1.234188352	1.102638087	1.145143666	1.036248506	0.512739655	0.419529869	0.200628096];
rgcData.midget.e0.LMS.p = [24.79881754	23.96124158	8.769912958	-17.8190179	-70.12645755	-171.0133027	-248.7436361	-346.1980621	-449.4334045	-495.3851207	-618.5580555];
rgcData.midget.e0.LminusM.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 20, 30, 40, 50, 60];
rgcData.midget.e0.LminusM.g = [0.883528514	1.041249317	1.335647522	1.994882817	2.737931338	2.222914502	1.145143666	0.642829544	0.325907056	0.192572001];
rgcData.midget.e0.LminusM.p = [-0.870422073	-7.488914436	-11.2169486	-37.82230251	-98.71900148	-219.559862	-328.7403515	-429.0688126	-532.2713089	-566.8089998];

rgcData.midget.e20.LMS.f = [0.625, 1.25, 2.5, 5, 10, 15, 20, 30, 40, 50, 60, 80];
rgcData.midget.e20.LMS.g = [0.16537005	0.156372507	0.311403691	0.599490731	1.5	1.9	2.112368642	2.253121364	2.406294456	1.4449392	0.685169767	0.202168197];
rgcData.midget.e20.LMS.p = [58.33333333	83.33333333	86.11111111	61.11111111	2.777777778	-61.11111111	-86.11111111	-194.4444444	-272.2222222	-405.5555556	-497.2222222	-616.6666667];
rgcData.midget.e20.LminusM.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 15, 20, 30, 40, 50, 60 78];
rgcData.midget.e20.LminusM.g = [0.674126896	0.909586301	1.24836258	1.47089674	1.539420942	1.912429335	1.908340174	1.426084262	1.201997763	0.599140868	0.104648486	0.06];
rgcData.midget.e20.LminusM.p = [26.81818182	10.81168831	-5.12987013	-35.35714286	-85.68181818	-132.4350649	-193.1168831	-311.3636364	-403.3766234	-552.6298701	-650.2597403	-733.7012987];

rgcData.midget.e30.LMS.f = [0.625, 1.25, 2.5, 5, 10, 15, 20, 30, 40, 50, 60, 80];
rgcData.midget.e30.LMS.g = [0.230838514	0.236127818	0.427498386	0.722467765	1.882228718	2.284530668	2.916885902	5.182037328	6.969143018	8.450342831	4.946604595	1.879437669];
rgcData.midget.e30.LMS.p = [40.55080722	35.99240266	48.40139285	41.05729661	8.452041785	-17.98037354	-69.29408041	-149.2877493	-228.7749288	-373.156062	-455.2389997	-568.5660019];
rgcData.midget.e30.LminusM.f = [0.609375, 1.2188, 2.4375, 4.8750, 9.7500, 15, 20, 30, 40, 50, 60 78];
rgcData.midget.e30.LminusM.g = [0.415777936	0.644266457	0.900089539	1.057382233	1.119496026	0.808479799	1.032163881	0.886480472	1.074358692	0.761510012	0.494562596	0.344497299];
rgcData.midget.e30.LminusM.p = [23.58341247	-0.696422919	8.990186768	-29.34472934	-78.91737892	-105.3181387	-128.4900285	-219.784742	-285.1535296	-426.7173156	-539.8227287	-655.8721114];

%   Solomon SG, Martin PR, White AJ, Rüttiger L, Lee BB. Modulation
%   sensitivity of ganglion cells in peripheral retina of macaque. Vision
%   Research. 2002 Dec 1;42(27):2893-8.
%
% Figures 1A and 1B
%

rgcData.parasol.e0.LMS.f = [0.625, 1.25, 2.5, 5, 10, 20, 30, 40, 50, 60, 80];
rgcData.parasol.e0.LMS.g = [0.262436577	0.591149071	1.644125739	3.702102437	5.900662716	7.471381688	4.057508719	2.615399798	1.591726565	0.672262083	0.495504391];

rgcData.parasol.e20.LMS.f = [0.625, 1.25, 2.5, 5, 10, 15, 20, 30, 40, 50, 60, 80];
rgcData.parasol.e20.LMS.g = [0.229534081	0.744221684	1.775076399	5.228815432	8.177641335	8.049716059	7.052555544	5.956457628	5.022726615	3.495853844	2.340097948	0.898244473];

rgcData.parasol.e30.LMS.f = [0.625, 1.25, 2.5, 5, 10, 15, 20, 30, 40, 50, 60, 80];
rgcData.parasol.e30.LMS.g = [0.27808972	0.638310685	2.413000771	5.755355979	8.33404175	7.454962426	7.1900742	7.6431795	9.831543136	9.482210723	8.303539014	3.645080274];

% Yeh T, Lee BB, Kremers J. Temporal response of ganglion cells of the macaque retina to cone-specific modulation. JOSA A. 1995 Mar 1;12(3):456-64.
%
% Figures 3E and 4C (just the S-cone modulation)
%
% Responses were to a spatially uniform, 4.7° diameter modulations at close
% to 400 Trolands.
%
% Mean response from 3-10° ecentricity for "S responsive" RGCs.

rgcData.bistratified.e0.S.f = [0.625, 1.25, 2.5, 5, 10, 20, 40];
rgcData.bistratified.e0.S.g = [0.826341879	1.085182942	1.444652234	2.058777757	2.933969679	2.327361819	0.655492984];
rgcData.bistratified.e0.S.p = [15.61201158	13.08774961	0.695086632	-32.56846705	-108.625578	-258.2246837	-513.1842894];

rgcData.bistratified.e20.S.f = [0.625, 1.25, 2.5, 5, 10, 20, 40].*1.05;
rgcData.bistratified.e20.S.g = [0.826341879	1.085182942	1.444652234	2.058777757	2.933969679	2.327361819	0.655492984];
rgcData.bistratified.e20.S.p = [15.61201158	13.08774961	0.695086632	-32.56846705	-108.625578	-258.2246837	-513.1842894];

rgcData.bistratified.e30.S.f = [0.625, 1.25, 2.5, 5, 10, 20, 40].*1.10;
rgcData.bistratified.e30.S.g = [0.826341879	1.085182942	1.444652234	2.058777757	2.933969679	2.327361819	0.655492984];
rgcData.bistratified.e30.S.p = [15.61201158	13.08774961	0.695086632	-32.56846705	-108.625578	-258.2246837	-513.1842894];

%{
    figure
    loglog(rgcData.midget.e0.LminusM.f,rgcData.midget.e0.LminusM.g,'-r','LineWidth',4);
    hold on
    loglog(rgcData.midget.e20.LminusM.f,rgcData.midget.e20.LminusM.g,'-r','LineWidth',2);
    loglog(rgcData.midget.e30.LminusM.f,rgcData.midget.e30.LminusM.g,'-r','LineWidth',1);
    
    loglog(rgcData.midget.e0.LMS.f,rgcData.midget.e0.LMS.g,'-k','LineWidth',4);
    loglog(rgcData.midget.e20.LMS.f,rgcData.midget.e20.LMS.g,'-k','LineWidth',2);
    loglog(rgcData.midget.e30.LMS.f,rgcData.midget.e30.LMS.g,'-k','LineWidth',1);
    
    figure
    loglog(rgcData.parasol.e0.LMS.f,rgcData.parasol.e0.LMS.g,'-k','LineWidth',4);
    hold on
    loglog(rgcData.parasol.e20.LMS.f,rgcData.parasol.e20.LMS.g,'-k','LineWidth',2);
    loglog(rgcData.parasol.e30.LMS.f,rgcData.parasol.e30.LMS.g,'-k','LineWidth',1);

    figure
    loglog(rgcData.bistratified.e01.S.f,rgcData.bistratified.e01.S.g,'-b','LineWidth',4);


%}



end
