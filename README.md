# mriSinaiAnalysis
Analysis code for the paper:

CP Gentile, M Spitchan, HO Taskin, A Bock, GK Aguirre. (in revision) Temporal sensitivity for achromatic and chromatic flicker across the visual cortex. J Neuroscience.

Overview of analysis pipeline

Raw fMRI data were stored on Flywheel

Pre-processing using fmriprep and ICA aroma pursued using Flywheel gears

The time-series data were analyzed using forwardModel. To do so, the non-linear fitting routine requires a "stimulus" matrix, and details regarding the order and timing of the stimulus events. These elements were created using the functions:

- makeStimStruct (downloads the raw "results" files and subject responses created by the computer that presented the stimuli at the time of scanning)
- loadStimStructCellArray (loads the "stimStruct" that summarizes these raw files)
- getAttentionEvents (finds where in the stimStruct attention events took place)
- stimConstructor_gka_asb (assembles the "stimulus" matrix from the "stimStruct")
- stimConstructor_cgp (subject cgp was studied with a different stimulus generation system; the stimulus matrix for this subject was generated in one step from the raw, "results" files generated by the stimulus presentation computer at the time of scanning).

The result of this stage of analysis were "stimulus.mat" files that were uploaded to Flywheel and used as inputs, along with the pre-processed fMRI data, to the forwardModel. The forwardModel gear runs on Flywheel were defined using the CSV files located within the "submitGears" directory.

The output of the forwardModel were the files "HEROxxx1_mtSinai_results.mat", where xxx is the subject ID. These files are stored in the data directory of this repo. The "mtSinai_results" files contains, for each grayordinate in the HCP CIFTI space, the parameters of the fit of the forwardModel. These parameters are effectively a beta weight for the amplitude of response to each trial type from each acquisition, with an additional parameter for the response to the attention event in each acquisition, and three parameters for the shape of the HRF found to best fit the responses for that grayordinate. The results file also contains an R2 value of the model fit to the time-series data for that coordinate, excluding the effect of the attention events.

These model parameters then served as inputs to subsequent analyses. The primary analysis was then pursued using the routine:

- fitWatsonModelVertex

located in the "temporalModel" directory. This function fits the Watson double-exponential TTF to the response parameters at each vertex. The result of this fitting operation was stored in the data results directory for each subject with the name "HEROxxx1_WatsonFit_results.mat".

The figures for the paper were then generated using the routines located in the figures directory. These figures make use of the mtSinai and Watson fit results, as well as retinotopic regions of interest defined using cortical topology (these retino files and regions of interest masks are stored in the data directory).
