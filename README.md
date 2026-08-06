This is the Git Repo for the BNL EICROC analysis efforts.

Scripts:

Full analysis of "real" data:

1) analyzeCSVData_radSource.cxx is a ROOT macro used for analysis of EICROC digital data collected using a Sr-90 source, but can also be used for other analyses, one just needs to make some modifications after the data are read-in.

ASIC calibration:
1) analyzeCSVData_ADC.cxx --> used to extract the ADC information, find the maximum of the ADC vs. time bin waveforms, and plot the ADC_max vs. your input charge values.
2) analyzeCSVData_ADCCorr.cxx --> used to analyze Vref_corr scan data to calculate Vref_corr offset values to account for high ADC baseline (and avoid saturation of the ADC for high input charge).
3) analyzeCSVData_Vth_global_scan.cxx --> used to analyze data from a Vth Global scan, and produce the s-curves from the data. Used both as input for "pivot" Vth global in a Vth_corr scan (to assign discriminator Vth trim values per pixel), and to then check that calibration after it is complete.
4) analyzeCSVData_thresholdCorr.cxx --> used to analyze Vth_corr scan data to calculate the per-pixel offsets to adjust for variable gain (seen as unaligned s-curves)


To run, you need ROOT, and simply type "root -b -q analyzeCSVData_radSource.cxx" to run in batch mode (no drawing of histograms, runs faster), or you can use "-l" in place of "-b -q" to plot histograms in real-time for QA.

There are several flags to be potentially toggled in the code, and the input file names are set in the file loop. This can fixed later.



