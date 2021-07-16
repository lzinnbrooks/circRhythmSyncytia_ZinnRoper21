# circRhythmSyncytia_ZinnRoper21
Public repository of Matlab scripts and figures from "Circadian rhythm shows potential...in multinucleate cells" - in PLoS Comp. Biol. 2021

Matlab scripts:

-main_GillespieSim.m is the most general code, and the most well-commented. It runs stochastic simulations of our syncytial circadian rhythm model. There are many additional options that you may toggle on and off, such as uniform/non-uniform compartments, bursty/non-bursty translation, and whether of not to display plots. This code was used (with minor modifications) to generate figures 2-5.

-powerSpectrum.m is a function called by main_GillespieSim.m, twoComp_divOfLabor.m, and eightComp_limitCycleConsistency.m. It computes the power spectrum, given a time vector and a signal.

-twoComp_divOfLabor.m is used to study division of labor (i.e., number of mRNAs transcribed in each compartment over each circadian cycle) in a two-compartment model syncytium. It was used to help generate figures 6-9.

-colorMap_mRNA.m is fed in a data file (scatterplot data from twoComp_divOfLabor.m) and uses this file to generate a color map for mRNA transcription, in which each point in the scatterplot is shaded depending on its frequency. This code was used to generate figure 6(a) and figure 8.

-eightComp_limitCycleConsistency.m is used to study limit cycle consistency in an eight-compartment model syncytium, with and without protein sharing. It computes the times that the first five peaks in total nuclear protein expression occurred. This code was run repeatedly (100 times each) in order to generate figures 10(a) and (b).

Figures:

Tif files are included of each figure. See our manuscript for parameter values used to generate figures.
