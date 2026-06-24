This is the Git Repo for the BNL EICROC analysis efforts.

Scripts:

analyzeCSVData_radSource.cxx is a ROOT macro used for analysis of EICROC digital data collected using a Sr-90 source, but can also be used for other analyses, one just needs to make some modifications after the data are read-in.

To run, you need ROOT, and simply type "root -b -q analyzeCSVData_radSource.cxx" to run in batch mode (no drawing of histograms, runs faster), or you can use "-l" in place of "-b -q" to plot histograms in real-time for QA.

There are several flags to be potentially toggled in the code, and the input file names are set in the file loop. This can fixed later.
