To submit these jobs to Flywheel, use the function submitGears that is included in the flywheelMRSupport toolbox. The syntax is:

	submitGears('myGearSubmission.csv')

The order of execution for these should be:

mtSinai_forwardModel_defineHRFParams.csv
	This runs the eventGain model on the average V1 time-series for each of the three stimulus types LF, L-M, S). The model fit includes finding the 3 parameters of a FLOBS HRF. I take the results of these three analyses and obtain the average HRF parameters across the three stimulus types in the average V1 region. These parameters could be used in subsequent analyses if we wish to lock the HRF parameters.

mtSinai_forwardModel_eventGain.csv

