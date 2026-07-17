#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TPaveText.h>
#include <sstream>


#include <iostream>
#include <fstream>

//Lightweight EICROC digital output analyzer
//
//
// Author: Alex Jentsch


// Event number (1024 of them, pixels (0,0) ... (31,31)
//          
//                                                                                                           TDC ADC  HB
//57876376839;   0;  68;   0;  0;  79;  0;  0;  74;  0;  0;  63;  0;  0;  57;  0;  0;  59;  0;  0;  63;  0;  0;  73;  0

//Data format is read form the end to the beginning (it's stored backwards).
//Format (right to left) is hit bit, ADC, and then TDC.
//each event is 25ns, and there are 8 time bins (the triplets stored back to front).
//There are 32 "copies" of the event number, where each one is the various pixels.

/****************

TDC = 10 bit (0 - 1023), 1 unit is 25ns / 1024
ADC = 8 bit (0 - 255), 40 MHz, Maximum ADC for each time sample

1 event is 8 time samples --> 25ns --> 10 bit

****************/


/****************

Scope of code:

Inputs: 
1) Sr-90 source data file (or truly any EICROC digital output file, just modity) 
2) pedestal file

The code will first try and find the pedestal file. 
If one is not found, it will let you know, and then skip the entire block related to the pedestal file

After that, the code will begin processing the Sr-90 (or other) data file, by reading the data lines 
for each pixel, and for each unique event ID into a buffer array, and then storing information in 
ROOT histograms for analysis. 

This part can easily be modified to do things differently.

*********************/


#include <string>

using namespace std;

int getPixelIndex(int i, int j);
int getPixelCanvasIndex(int i);

// string for input datafile
TString inputFileString = "/Users/rkfuentes/Documents/phd/research/BNL_summer_2026/eicroc1/converted_output.csv";

Color_t markerColor[16] = {kGray, kGray+1, kGray+2, kGray+3, kGreen, kGreen+1, kGreen+2, kGreen+3, kBlue, kBlue+1, kBlue+2, kBlue+3, kRed, kRed+1, kRed+2, kRed+3};

int markerStyle[26] = {20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 29, 29, 29, 29};

int numGoodEvents[32][32];
int numEvents = 1999;
int pixels_of_interest[16] = {0, 1, 2, 3, 32, 33, 34, 35, 64, 65, 66, 67, 96, 97, 98, 99}; // {0, 7, 15, 31, 32, 39, 47, 63, 64, 71, 78, 85, 86, 93, 100, 107};
int pixel_colors[16] = {kBlack, kRed, kBlue, kGreen+2, kMagenta, kCyan+2, kOrange+1, kSpring-1, 
                        kViolet, kPink+9, kTeal-1, kAzure+1, kYellow+2, kGray+2, kRed-7, kBlue-7};

void eicroc1_plotter(TString inputFileName = "", TString inputFileName_pedestal = ""){
	// set up input file stream
	ifstream inputCSVFile;
	ifstream inputCSVFile_pedestal;
	
	bool subtractPedestals = false;
	bool usePedestalFile = false;
	
	bool subtractADC_baseline = false; //this is to set the middle of the ADC distributions to 0 to try and mimic Arzoo's analysis
	
	double thresholdStepSize   = 10;
	const int numOfAnalyzedEvents = 100; //This sets the size of buffer array for the event by event analysis - meaning, number of event user deems "good"
	
	int trigger_pixel = 5; //This sets the trigger pixel to focus on, but this logic can be disabled in the main code block, if desired
	double trigger_bin_ADC_sum = 0.0;
	
	TGraph * s_curve[1024]; // storage graphs for s_curves -- not used in this code
	for(int pixel = 0; pixel < 1024; pixel++){s_curve[pixel] = new TGraph();}
	
	std::vector<double> threshold_values;
	std::vector<double> efficiency;
	
	int pad = 1; //Sets a pad number that you can use later for drawing histograms on a divide canvas
	
	//various histograms for storing information for a 4x4 pixel array - can be modified to handle 32x32 pixel array
	// 32x32 for eicroc1
	TH1D * adc_distributions[32][32];
	TH1D * tdc_distributions[32][32];

	double adc_event_buffer[numOfAnalyzedEvents][1024][8] = {0.0};
	
	TH1D * adc_distributions_PEDESTAL[32][32];

	// number of thresholds, number of time bins, number of pixels
	TH1D * adc_RAW_distributions[numOfAnalyzedEvents][8][1024];
	TH1D * adc_RAW_distributions_pedestal[8][1024];

	// noise distr histogram
	TH1D * h_noise_distributions[1024];
	TH1D * h_noise_per_timebin[2][8]; // for 2 pixels only
	
	// array which identifies pixels where hitbit is 1
	TH1D * hitPixel = new TH1D("hit_pixel", "hit_pixel", 1024, 0.0, 1024.0);
	// array of ADC means across ALL pixels including hitbit 0
	TH2D * ADC_mean_map = new TH2D("ADC_mean_map", "ADC_mean_map", 32, 0, 32, 32, 0, 32);
	
	// number of thresholds, AVERAGED OVER TIME BINS, number of pixels
	TGraphErrors * adc_mean_distributions[1024];
	TGraph * adc_mean_distributions_PEDESTAL[1024];

	TGraph * tdc_mean_distributions[numOfAnalyzedEvents][1024];
	TGraph * tdc_mean_distributions_PEDESTAL[numOfAnalyzedEvents][1024];

	// 16 scatter plots for pedestal subtracted means per timebin
	TMultiGraph *mg_drift = new TMultiGraph();
	mg_drift->SetTitle("SUM(ADCn)/Nevents - ADC(mean);Time Bin (25ns);Average Delta [DACu]");

	TGraph *g_pixel_drift[16];
	for (int i = 0; i < 16; i++) {
		g_pixel_drift[i] = new TGraph();
	}
	
	TGraphErrors * pedestal_waveforms[1024];

	// loops to generate all required 1d histograms for pedestal and adc max distr threshold
	for(int tBin = 0; tBin < 8; tBin++){
		for(int pixel = 0; pixel < 1024; pixel++){
			
			TString title;
			title.Form("adc_max_distribution_pixel_%d_timeBin_%d_PEDESTAL", pixel, tBin);

			adc_RAW_distributions_pedestal[tBin][pixel] = new TH1D(title, "; ADC value [DACu]; counts", 256, 0, 255);
			adc_RAW_distributions_pedestal[tBin][pixel]->SetTitle(title);
				
		}
	}
	
	for(int ev = 0; ev < numOfAnalyzedEvents; ev++){
		for(int tBin = 0; tBin < 8; tBin++){
			for(int pixel = 0; pixel < 1024; pixel++){

				TString title;
				title.Form("adc_max_distribution_thresh_%d_pixel_%d_timeBin_%d", ev, pixel, tBin);

				adc_RAW_distributions[ev][tBin][pixel] = new TH1D(title, "; ADC value [DACu]; counts", 256, 0, 255);
				adc_RAW_distributions[ev][tBin][pixel]->SetTitle(title);
				
	
			}
		}
	}	

	// initialize arrays for ADC mean and noise distr for each pixel across all events
	for(int pixel = 0; pixel < 1024; pixel++){

		TString title;
		title.Form("adc_MEAN_distribution_thresh_pixel_%d", pixel);

		adc_mean_distributions[pixel] = new TGraphErrors();
		
		adc_mean_distributions[pixel]->GetXaxis()->SetTitle("time bin [time bin = 25ns"); 
		adc_mean_distributions[pixel]->GetYaxis()->SetTitle("ADC value [DACu]");
		adc_mean_distributions[pixel]->SetTitle(title);

		TString noise_title;
		noise_title.Form("h_noise_distribution_pixel_%d",pixel);

		h_noise_distributions[pixel] = new TH1D(noise_title, "Raw ADC Noise Distribution;ADC Channel;Counts", 256, 0, 256);

		if (pixel <= 1) {
			for (int tbin = 0; tbin < 8; tbin++){
				TString noise_tbin_title;
				noise_tbin_title.Form("h_noise_distribution_pixel_%d_tbin_%d",pixel, tbin);
				h_noise_per_timebin[pixel][tbin] = new TH1D(noise_tbin_title, "Raw ADC Noise;ADC Channel;Counts", 25, 0, 256);
			}
		}

	}

	for(int i = 0; i < 32; i++){
		for(int j = 0; j < 32; j++){
		
			adc_distributions[i][j] = new TH1D(Form("adc_max_distribution_pixel_%d_%d", i, j), "counts; ADC value [DACu]", 256, 0, 255);
			tdc_distributions[i][j] = new TH1D(Form("tdc_distribution_pixel_%d_%d", i, j), "counts; TDC value [DACu]", 1024, 0, 1024);
			
			adc_distributions_PEDESTAL[i][j] = new TH1D(Form("adc_max_distribution_pixel_%d_%d_PEDESTAL", i, j), "counts; ADC value [DACu]", 256, 0, 255);
		}
	}
	
	double numFiles = 0;
	
	int charge_DACu = 45;
	
	int numTriggeredEvents = 0;
	
	//older data structure to hold pre-calculated pedstals -- not used right now, but left here for comparison testing, if needed.
	//this will be replaced by a 32x32 pedestal map, once it is determined 
	const double pedestal_map[16][8] = {
	    {31.517, 32.544, 30.381, 29.391, 30.371, 30.646, 30.762, 32.269},  // pixel 0
	    {34.122, 35.456, 33.299, 32.568, 32.942, 32.759, 33.156, 35.034},  // pixel 1
	    {34.728, 35.983, 33.714, 32.918, 33.503, 33.347, 33.656, 35.609},  // pixel 2
	    {32.065, 33.088, 30.820, 29.915, 30.707, 31.323, 31.272, 32.803},  // pixel 3
	    {30.578, 31.609, 29.554, 28.476, 29.490, 29.867, 29.759, 31.293},  // pixel 4
	    {28.561, 29.459, 27.265, 26.381, 27.622, 28.367, 27.748, 29.224},  // pixel 5
	    {33.000, 34.156, 31.898, 30.976, 31.765, 31.793, 32.095, 33.833},  // pixel 6
	    {33.184, 34.194, 31.905, 30.952, 32.003, 32.677, 32.241, 33.912},  // pixel 7
	    {31.711, 31.054, 29.588, 31.969, 41.320, 42.265, 33.629, 35.279},  // pixel 8
	    {33.878, 34.687, 32.401, 31.656, 33.313, 34.646, 33.017, 34.789},  // pixel 9
	    {34.119, 35.031, 32.728, 31.976, 33.146, 34.327, 33.024, 34.898},  // pixel 10
	    {33.993, 34.711, 32.286, 31.310, 33.279, 35.629, 33.299, 34.932},  // pixel 11
	    {32.537, 33.344, 31.020, 30.092, 31.293, 31.918, 31.690, 33.320},  // pixel 12
	    {32.990, 33.473, 31.031, 30.381, 32.503, 34.827, 32.143, 33.932},  // pixel 13
	    {33.350, 31.949, 29.912, 33.054, 44.241, 46.170, 35.704, 37.537},  // pixel 14
	    {35.044, 34.439, 32.612, 34.776, 41.293, 42.065, 36.241, 38.031},  // pixel 15
	};
	
	
	/********************
	
	this for loop is only for if you are analyzing many files with the same name, but a different index.
	CURRENTLY only one file is being processed at a time
	 
	to loop over multiple filenames, they must all have the format:
	inputFileName.Form("[path_to_data_files]/[my_file_name]_%.0f_1.csv", fileIdx);
	
	then you can loop over files from the same run if you broke them up into smaller pieces, or did a scan, etc.
	
	***********************/
	
	for(int fileIdx = 0; fileIdx < 1; fileIdx++){ //BEGIN for loop over file indices
		
		//Input files go here --> one is for a pedestal file, the other is for the main data file
		// pedestal file is not currently being used
		inputFileName.Form(inputFileString);
		
		cout << inputFileName << endl;
		cout << inputFileName_pedestal << endl;

		// read in input data file
		inputCSVFile.open(inputFileName.Data());
		inputCSVFile_pedestal.open(inputFileName_pedestal.Data());

		if(!inputCSVFile){ 
			cout << "File does not exist, skipping...." << endl; 
			return;
		}
		
		if(!inputCSVFile_pedestal){ 
			cout << "PEDESTAL File does not exist, using internal pedestal map..." << endl; 
			usePedestalFile = false;
		} else {
			usePedestalFile = true;
		}

		numFiles++;

		// initialize arrays for events, ADC values, hitbits, etc
		string eventNumber[32][32];
		double numGoodEvents[32][32];
		double ADC_values[8];
		int TDC_values[8];
		int hitBits[8];
		
		double ADC_values_event[1024][8];
		int TDC_values_event[1024][8];
		int hitBits_event[1024][8];

		double adc_event_buffer[numOfAnalyzedEvents][1024][8] = {0.0};

		double ADC_sum_all_timebins = 0.0;

		double ADC_max_pixel[32][32];
		double ADC_max_whole_ASIC = 0;
		int pixel_with_max_ADC_ASIC = 0;

		string eventNumber_str;
		string pixelNumber_str;
		string ADC_values_str[8];
		string TDC_values_str[8];
		string hitBits_str[8];

		double tmp_ADC_MEAN_pedestal[8];
		double tmp_ADC_MEAN_err_pedestal[8];
		
		// pedestal currently not being used
		double tmp_ADC_MEAN_pedestal_BINS[8]      = {0, 1, 2, 3, 4, 5, 6, 7};
		double tmp_ADC_MEAN_err_pedestal_x_err[8] = {0, 0, 0, 0, 0, 0, 0, 0};
		
		double avg_ADC_offset_pedestal_pixel[1024];
		double max_ADC_pedestal_pixel[1024];
		
		// array to indicate which pixels have hitbit 1
		int pixel_has_hit[1024];


		// initialize to start and pixel 0,0 timebin 7, event 0
		int timeBin     = 7;
		int pixelRow    = 0;
		int pixelColumn = 0;
		
		int ADC_max_time_bin = 0;

		string tmp;
		
		// for current example, 10 total events in input file
		int numEvents = 10;
		
		for(int i = 0; i < 32; i++){
			for(int j = 0; j < 32; j++){
				numGoodEvents[i][j]++;
			}
		}
		
		int numHitBitSet = 0;
		bool hitBitInTriggerPixel = false; // currently pixel of interest is pixel 5, checks if there is a hitbit in that pixel
		bool hitBitInCentralTimeBins = false; // not currently used
		bool hitBitPixelFifteen = false; // not currently used
		
		bool pedestalNotSynched = true;

		int hitBitPixel[32][32];

		///////////////////////////////////////////////////////////////////
		//////// BEGINNING OF PEDESTAL CODE BLOCK /////////////////////////
		///////////////////////////////////////////////////////////////////

		// not currently used!
		if (usePedestalFile) {
			// read in lines from inputcsv file
			string rawLine;
			while (getline(inputCSVFile, rawLine)) { // read entire line
				if (rawLine.empty()) continue;

				stringstream ss(rawLine);
				string token;
				
				if (!getline(ss, eventNumber_str, ';')) continue;
				
				// parse the 8 triplets (TDC; ADC; HitBit) from back to front
				bool lineParseSuccess = true;
				for (int tBin = 7; tBin >= 0; tBin--) {
					if (getline(ss, TDC_values_str[tBin], ';') &&
						getline(ss, ADC_values_str[tBin], ';') &&
						getline(ss, hitBits_str[tBin], ';')) {
						// Token read successfully
					} else {
						lineParseSuccess = false;
						break; // Missing columns, skip corrupt row cleanly
					}
				}

				if (!lineParseSuccess) {
					cout << "WARNING: Skipped a malformed line: " << rawLine << endl;
					continue; 
				}

				// continue through 32x32 pixel array, if exists
					pixelRow++;
					if(pixelRow == 32 && pixelColumn < 31) { pixelColumn++; pixelRow = 0;}
					if(pixelRow == 32 && pixelColumn == 3){ pixelColumn++;}
					if(pixelRow == 32 && pixelColumn == 32){ 
						
						if(numHitBitSet >= 1 && numTriggeredEvents < numOfAnalyzedEvents){
							
							numTriggeredEvents++;
							
							for(int i = 0; i < 32; i++){
								for(int j = 0; j < 32; j++){
									
									int pixel = getPixelIndex(i, j);
																	
									for(int tBin = 0; tBin < 8; tBin++){
																	
										if(pixel != -1){
											adc_RAW_distributions_pedestal[tBin][pixel]->Fill(ADC_values_event[pixel][tBin]);
										}
										
										// reset internal counter and flags
										pixelRow = 0; pixelColumn = 0; numHitBitSet = 0; 
										hitBitInTriggerPixel = false; // Now it's safe to reset
										
										numEvents++;

									}
								}
							}
						}
						
						pixelRow = 0; pixelColumn = 0; numHitBitSet = 0;
						
						numEvents++;

					}
				} //END while loop over PEDESTAL file

		inputCSVFile_pedestal.close();

		cout << "number of pedestal events --> " << numTriggeredEvents << endl;
		cout << "\n\n";

		} else {
		cout << "Skipping pedestal file reading" << endl;
		}
		
	if(usePedestalFile){ //BEGIN conditional to print pedestal plots --> only used if you had a pedestal file, and set the flag

			TCanvas * adcRAWCan[8]; 
		
			for(int tBin = 0; tBin < 8; tBin++){
				adcRAWCan[tBin] = new TCanvas(Form("can_tBin_%d", tBin), Form("can_tBin_%d", tBin), 2000, 1600);
				adcRAWCan[tBin]->Divide(4,4);
			
			}
		
			pad = 1;
		
			adcRAWCan[0]->cd(1);

			for(int pixel_id = 0; pixel_id < 16; pixel_id++){
				int pixel = pixels_of_interest[pixel_id];
			
				int idx = getPixelCanvasIndex(pad);
			
				avg_ADC_offset_pedestal_pixel[pixel] = 0.0;
				max_ADC_pedestal_pixel[pixel] = 0.0;
					
				for(int tBin = 0; tBin < 8; tBin++){
		
					adcRAWCan[tBin]->cd(idx);
		
					adc_RAW_distributions_pedestal[tBin][pixel]->SetLineColor(markerColor[tBin]);
				
					adc_RAW_distributions_pedestal[tBin][pixel]->Draw();
				
					if(tBin == 0) {cout << "Pixel " << pixel << " ADC pedestals -->\t" << adc_RAW_distributions_pedestal[tBin][pixel]->GetMean() << "\t";}
					if(tBin > 0 && tBin < 7) {cout << adc_RAW_distributions_pedestal[tBin][pixel]->GetMean() << "\t";}
					if(tBin == 7){cout << adc_RAW_distributions_pedestal[tBin][pixel]->GetMean() << endl;}
		
					tmp_ADC_MEAN_pedestal[tBin]     = adc_RAW_distributions_pedestal[tBin][pixel]->GetMean();
					tmp_ADC_MEAN_err_pedestal[tBin] = adc_RAW_distributions_pedestal[tBin][pixel]->GetRMS();
		
					avg_ADC_offset_pedestal_pixel[pixel]+= adc_RAW_distributions_pedestal[tBin][pixel]->GetMean();
		
					if(adc_RAW_distributions_pedestal[tBin][pixel]->GetMean() > max_ADC_pedestal_pixel[pixel]) { max_ADC_pedestal_pixel[pixel] = adc_RAW_distributions_pedestal[tBin][pixel]->GetMean();}
		
				}
			
				avg_ADC_offset_pedestal_pixel[pixel] = avg_ADC_offset_pedestal_pixel[pixel]/8.0;
			
				pedestal_waveforms[pixel] = new TGraphErrors(8, tmp_ADC_MEAN_pedestal_BINS, tmp_ADC_MEAN_pedestal, tmp_ADC_MEAN_err_pedestal_x_err, tmp_ADC_MEAN_err_pedestal);
				pad++;
			}
		
		
			for(int tBin = 0; tBin < 8; tBin++){
				adcRAWCan[tBin]->SaveAs(Form("pedestals_for_timeBin_%d.png", tBin));
			}
		
			//DRAW PEDESTAL FITS HERE
		
			TCanvas * pedestal_fit_canvas = new TCanvas("ped_fit_can", "ped_fit_can", 2000, 1600);
			pedestal_fit_canvas->Divide(4,4);
		
			pad = 1;
			TF1 * pedFits[1024];
		
			for(int pixel = 0; pixel < 1024; pixel++){
			
			
				int idx = getPixelCanvasIndex(pad);
			
				pedestal_fit_canvas->cd(idx);
			
				pedFits[pixel] = new TF1(Form("pixel_%d_ped_fit", pixel), "[0] + [1]*sin([2]*x + [3])", 0, 7);
			
				double amplitude =  max_ADC_pedestal_pixel[pixel] - avg_ADC_offset_pedestal_pixel[pixel];
			
				cout << "pixel " << pixel << " avg ADC = " << avg_ADC_offset_pedestal_pixel[pixel] << " and max ADC = " << max_ADC_pedestal_pixel[pixel] << " possible amp = " << amplitude << endl;
			
				pedFits[pixel]->SetParLimits(1, 25, 35);
				pedFits[pixel]->SetParLimits(0, 35, 45);
			
				pedFits[pixel]->SetParLimits(2, 1.0, 2.0);
				pedFits[pixel]->SetParLimits(3, -TMath::PiOver2(), TMath::PiOver2());
			
			
				pedestal_waveforms[pixel]->Fit(pedFits[pixel], "RQ");
			
				pedestal_waveforms[pixel]->Draw("ALP");
			
				pad++;
			}

		} //END conditional to print pedestal plots
		
		///////////////////////////////////////////////////////////////////
		//////// END OF PEDESTAL CODE BLOCK ///////////////////////////////
		///////////////////////////////////////////////////////////////////

		//RESET internal variables

		// numEvents = 0;
		// initialize all pixels as having good events to ensure all data is read
		for(int i = 0; i < 32; i++){
			for(int j = 0; j < 32; j++){
				numGoodEvents[i][j]++;
			}
		}
		
		numHitBitSet = 0;
		numTriggeredEvents = 0;
		
		///////////////////////////////////////////////////////////////////
		//////// BEGINNING OF MAIN DATA CODE BLOCK ///////////////////////////////
		///////////////////////////////////////////////////////////////////
		
		//pixelColumn = 0; pixelRow = 0;

		// counter for each
		int currentEventID = -1;  // Tracks the ID of the event we are currently assembling
		bool firstLineOfFile = true;

		while(!inputCSVFile.eof()){
			// extract ADC values, event number and unpack pixel number from each line of inputcsvfile
			if (getline(inputCSVFile, eventNumber_str, ';') &&
				getline(inputCSVFile, pixelNumber_str, ';') &&
				getline(inputCSVFile, TDC_values_str[7], ';') && getline(inputCSVFile, ADC_values_str[7], ';') && getline(inputCSVFile, hitBits_str[7], ';') &&
				getline(inputCSVFile, TDC_values_str[6], ';') && getline(inputCSVFile, ADC_values_str[6], ';') && getline(inputCSVFile, hitBits_str[6], ';') &&
				getline(inputCSVFile, TDC_values_str[5], ';') && getline(inputCSVFile, ADC_values_str[5], ';') && getline(inputCSVFile, hitBits_str[5], ';') &&
				getline(inputCSVFile, TDC_values_str[4], ';') && getline(inputCSVFile, ADC_values_str[4], ';') && getline(inputCSVFile, hitBits_str[4], ';') &&
				getline(inputCSVFile, TDC_values_str[3], ';') && getline(inputCSVFile, ADC_values_str[3], ';') && getline(inputCSVFile, hitBits_str[3], ';') &&
				getline(inputCSVFile, TDC_values_str[2], ';') && getline(inputCSVFile, ADC_values_str[2], ';') && getline(inputCSVFile, hitBits_str[2], ';') &&
				getline(inputCSVFile, TDC_values_str[1], ';') && getline(inputCSVFile, ADC_values_str[1], ';') && getline(inputCSVFile, hitBits_str[1], ';') &&
				getline(inputCSVFile, TDC_values_str[0], ';') && getline(inputCSVFile, ADC_values_str[0], ';') && getline(inputCSVFile, hitBits_str[0]))
			{
				int lineEventID = atoi(eventNumber_str.c_str());
				int pixel = std::stoi(pixelNumber_str);

				// --- EVENT BOUNDARY TRIGGER ---
				// whenever the event ID changes in the file, mark event as completed
				if (lineEventID != currentEventID && !firstLineOfFile) {
					
					// UNCONDITIONALLY REGISTER EVERY EVENT FRAME CAUGHT
					numTriggeredEvents++;
					numEvents++;
					
					if (numTriggeredEvents >= numOfAnalyzedEvents) {
						cout << "WARNING: Reached target array storage limit (" << numOfAnalyzedEvents << ")." << endl;
						break;
					}
				}

				if (firstLineOfFile) {
					currentEventID = lineEventID;
					firstLineOfFile = false;
				}
				currentEventID = lineEventID;

				// loop through all 8 time bins and fill the graphs
				for (int tBin = 0; tBin < 8; tBin++) {
					int ADC_val = atoi(ADC_values_str[tBin].c_str());

					if (pixel < 1024 && numTriggeredEvents < numOfAnalyzedEvents) {
						adc_event_buffer[numTriggeredEvents][pixel][tBin] = ADC_val;
						h_noise_distributions[pixel]->Fill(ADC_val);
						if (pixel <= 1) {
							// only fill 2 pixels of interest for timebin graphs to avoid wasting memory
							h_noise_per_timebin[pixel][tBin]->Fill(ADC_val);
						}
					}
				}
			}
		} // END of while loop
		
		cout << "Number of events total recorded : " << numEvents << endl;

		inputCSVFile.close();

		if (numTriggeredEvents < numOfAnalyzedEvents) {
			numTriggeredEvents++;
			numEvents++;
		}

		// =========================================================================
		// CALCULATE MULTI-EVENT AVERAGE PROFILE AND TIME-BIN RMS
		// =========================================================================
		cout << "Calculating multi-event average waveforms and RMS noise profiles..." << endl;

		// loop over each of your 16 active plotting channels

		for (int pixel_id = 0; pixel_id < 16; pixel_id++) {
			int pixel = pixels_of_interest[pixel_id];

			double pixel_total_sum = 0.0;
			int pixel_total_count = 0;

			for (int tBin = 0; tBin < 8; tBin++) {
				for (int ev = 0; ev < numTriggeredEvents; ev++) {
					pixel_total_sum += adc_event_buffer[ev][pixel][tBin];
					pixel_total_count++;
				}
			}

			if (pixel_total_count > 0) {
				double pixel_overall_mean = pixel_total_sum / pixel_total_count;
				int col = floor(pixel / 32);
				int row = pixel - (32 * col);
				
				// Fill the 2D map precisely (X axis = Col, Y axis = Row)
				ADC_mean_map->SetBinContent(col + 1, row + 1, pixel_overall_mean);
				
			}

			// loop over each of the 8 sequential time slices
			for (int tBin = 0; tBin < 8; tBin++) {
				
				double sum = 0.0;
				double sum_squares = 0.0;
				int count = 0;
				double global_sum = 0.0;
				int global_count = 0;

				// populate buffer matrix with adc values, per event, pixel and timebin
				for (int ev = 0; ev < numTriggeredEvents; ev++) {
					double adc_val = adc_event_buffer[ev][pixel][tBin];
					sum += adc_val;
					sum_squares += (adc_val * adc_val);
					count++;
					global_sum += adc_event_buffer[ev][pixel][tBin];
					global_count++;
				}

				double mean = 0.0;
				// if statements compute rms
				double rms_uncertainty = 0.0;
				double global_pixel_mean = (global_count > 0) ? (global_sum / global_count) : 0.0;
				
				if (count > 0) {
					mean = sum / count;
					
					if (count > 1) {
						// calculate variance across events
						double variance = (sum_squares / count) - (mean * mean);
						if (variance > 0) {
							rms_uncertainty = TMath::Sqrt(variance);
							pixelColumn = floor(pixel/32);
							pixelRow = pixel - (32*pixelColumn);
							ADC_mean_map->Fill(pixelColumn + 1, pixelRow + 1, mean);
						}
					}
				}

				// calculate mean per timebin for the 16 scatter plots of interest
				for (int tBin = 0; tBin < 8; tBin++) {
					double timebin_sum = 0.0;
					int timebin_count = 0;
					for (int ev = 0; ev < numTriggeredEvents; ev++) {
						timebin_sum += adc_event_buffer[ev][pixel][tBin];
						timebin_count++;
					}
					double timebin_mean = (timebin_count > 0) ? (timebin_sum / timebin_count) : 0.0;

					// calculate the relative delta
					double delta = timebin_mean - global_pixel_mean;

					// add the coordinate point (X = Time Bin, Y = Delta) to the pixel's graph
					g_pixel_drift[pixel_id]->SetPoint(tBin, tBin, delta);
				}
				// set a single average data point per time bin for the pixel graph
				int pointIdx = adc_mean_distributions[pixel]->GetN();
				adc_mean_distributions[pixel]->SetPoint(pointIdx, tBin, mean);
				
				// vertical error bar now represents the standard dev across frames
				adc_mean_distributions[pixel]->SetPointError(pointIdx, 0.0, rms_uncertainty);

				// set individual pixels for the layered scatter plot
				g_pixel_drift[pixel_id]->SetMarkerStyle(20 + (pixel_id % 4)); // Varies marker shape (circle, square, triangle, etc.)
				g_pixel_drift[pixel_id]->SetMarkerSize(1.2);
				g_pixel_drift[pixel_id]->SetMarkerColor(pixel_colors[pixel_id]);
				g_pixel_drift[pixel_id]->SetLineColor(pixel_colors[pixel_id]);

				// add this pixel's graph to the multigraph container
				mg_drift->Add(g_pixel_drift[pixel_id]);
			}
		}
		cout << "Averages and graph error properties successfully calculated!" << endl;

	} //END for loop over file indices
	
	cout << Form("Number of events triggered in pixel %d : ", trigger_pixel) << numTriggeredEvents << endl;

	///////////////////////////////////////////////////////////////////
	//////// END OF MAIN DATA CODE BLOCK ///////////////////////////////
	///////////////////////////////////////////////////////////////////
		/**************
	
	These sections are used for doing s-curve calculations, which I now do in a different macro. However, it can be brought back here, as well.
	
	To do an S-Curve analysis, you need to use the "good" event concept, but on a pixel by pixel basis, and as a function of global threshold.
	The idea is to see how many events in that pixel, for a given threshold, had a hitBit == 1, compared to the number of events you sampled. 
	The ratio of those two numbers, as a function of threshold, is the S-curve.
	
	These can be fit with Logistical functions to extract information about the behavior, and perform calibrations.
		****************/
		// temporarily manually setting number of good events
/*
		cout << "\n --- number of good events in pixels --- " << endl;
		cout << " THRESHOLD = " << 100 << endl;
		for(int i = 0; i < 32; i++){
			for(int j = 0; j < 32; j++){
				numGoodEvents[i][j] = 10;
			
				int pixel = getPixelIndex(i, j);
				
				if(pixel != -1){
					s_curve[pixel]->AddPoint(100, numGoodEvents[i][j]/numEvents);
				}
			
				cout << "pixel( " << i << ", " << j << ") = " << numGoodEvents[i][j]/numEvents << endl;
				
			}
		}

	} //end of file loop
	
	// s curves not working yet
	TGraph s_curve = new TGraph(threshold_values.size(), &threshold_values[0], &efficiency[0]);
	
	TCanvas * sCurveCanvas = new TCanvas("canv1", "canv1", 1600, 1600);
	sCurveCanvas->Divide(4,4);
	
	double thresholdMin = 300.0;
	double thresholdMax = 0.0;
	
	int thresh_50_percent_values[16];
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
	
			int pixel = getPixelIndex(i, j);
			int canvasBin = getPixelIndex(j, i);
	
			if(pixel != -1){
				
				sCurveCanvas->cd(canvasBin+1);
				
				s_curve[pixel]->SetTitle(Form("S-Curve, pixel (%d,%d), Q_{injected} = %d DACu", i, j, charge_DACu));
				s_curve[pixel]->GetXaxis()->SetTitle("Threshold [DACu]");
				s_curve[pixel]->GetYaxis()->SetTitle("Efficiency");
				s_curve[pixel]->SetMarkerStyle(markerStyle[pixel]);
				s_curve[pixel]->SetMarkerColor(markerColor[pixel]);
				s_curve[pixel]->SetLineColor(markerColor[pixel]);
				
				double thresh = 0;
				double eff = 0;
				
				for(int point = 0; point < s_curve[pixel]->GetN(); point++){
					
					thresh = s_curve[pixel]->GetPointX(point);
					eff = s_curve[pixel]->GetPointY(point);
					
					if(eff > 0.3 && eff < 0.7 && thresh > 370){ break; }
				}
					
				
				//TF1 * sCurveFit = new TF1("scurvefit", "-1*[0]*TMath::Erf((x-[1])/[2])", thresh - 100, thresh + 150);
				//TF1 * sCurveFit = new TF1("scurvefit", "[0] / (1 + exp(−1*[1]*(x−[2])))", thresh - 100, thresh + 150);
				//sCurveFit->SetParLimits(0, 0.99, 1.01);
				//sCurveFit->SetParLimits(1, thresh - 25, thresh + 25);
				//sCurveFit->SetParLimits(2, 10, 40);
				
				TF1 *fLogistic = new TF1("fLogistic", "[0] / (1.0 + TMath::Exp(-[1] * (x - [2]))) + [3]", thresh - 70, thresh + 250);

				fLogistic->SetParNames("Amplitude", "Steepness", "Midpoint", "Offset");
				fLogistic->SetParameters(1.0, 0.2, thresh, 0.0);
				
				s_curve[pixel]->Fit(fLogistic, "RQN");
				
				TF1 *fInverse = new TF1("fInverse", "(-1/[0]) * TMath::Log(([1] - x)/x) + [2]", 0.0, 1.0);
				fInverse->SetParameter(0, fLogistic->GetParameter(1));
				fInverse->SetParameter(1, fLogistic->GetParameter(0));
				fInverse->SetParameter(2, fLogistic->GetParameter(2));
				
				cout << "pixel " << pixel << "  50 percent threshold = " << fInverse->Eval(0.5) << endl;
				
				thresh_50_percent_values[pixel] = (int)fInverse->Eval(0.5);
				
				if(fInverse->Eval(0.5) < thresholdMin) { thresholdMin = fInverse->Eval(0.5); }
				if(fInverse->Eval(0.5) > thresholdMax) { thresholdMax = fInverse->Eval(0.5); }
				
				TLine * line = new TLine(fInverse->Eval(0.5), 0, fInverse->Eval(0.5), 1.0);
				line->SetLineWidth(3);
				line->SetLineColor(markerColor[pixel]);
				
				s_curve[pixel]->Draw("ALP");
				line->Draw("SAME");
				
			}
		}
	}

	
	cout << "minimum threshold = " << thresholdMin << endl;
	cout << "maximum threshold = " << thresholdMax << endl;
	cout << "\n";
	
	for(int pixel = 0; pixel < 16; pixel++){
		
		
		//int offset_value =  (int)thresholdMax - thresh_50_percent_values[pixel] ;
		int offset_value =  467 - thresh_50_percent_values[pixel] ;
		
		cout << "Offset for pixel " << pixel << " = " << offset_value;
		cout << " --> " << std::bitset<7>(offset_value) << endl;
		
	}
			
	TCanvas * sCurvesOneCanvas = new TCanvas("canv2", "canv2", 800, 800);
	sCurvesOneCanvas->cd();
	
	TLine * line_50_percent = new TLine(350, 0.5, 650, 0.5);
	line_50_percent->SetLineWidth(3);
	line_50_percent->SetLineColor(kRed);
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
	
			int pixel = getPixelIndex(i, j);
			int canvasBin = getPixelIndex(j, i);
	
			if(pixel != -1){
				
				//if(pixel == 0  || pixel == 1  || pixel == 2  || pixel == 3 ){continue;}
				//if(pixel == 4  || pixel == 5  || pixel == 6  || pixel == 7 ){continue;}
				//if(pixel == 8  || pixel == 9  || pixel == 10 || pixel == 11 ){continue;}
				//if(pixel == 12 || pixel == 13 || pixel == 14 || pixel == 15 ){continue;}
				
				//s_curve[pixel]->GetXaxis()->SetRangeUser(350, 650);
				
				
				
				if(pixel == 0){ 
					s_curve[pixel]->Draw("ALP");
					line_50_percent->Draw("SAME");
				}
				else s_curve[pixel]->Draw("SAME LP");
				
			}
		}
	}
		
		
	int pad = 1;
		
	TCanvas * tdcCan = new TCanvas("canv3", "canv3", 1600, 1600);
	tdcCan->Divide(4,4);

	pad = 1;
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
		
			tdcCan->cd(pad);
			tdc_distributions[i][j]->Draw();
			pad++;
		}
	}	
	
	
	pad = 1;
	
	
	TCanvas * adcCan = new TCanvas("canv4", "canv4", 1600, 1600);
	adcCan->Divide(4,4);

	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
		
			adcCan->cd(pad);
			adc_distributions[i][j]->Draw();
			pad++;
		}
	}

	pad = 1;
	
	TCanvas * adcRAWCan = new TCanvas("canv5", "canv5", 1600, 1600);
	adcRAWCan->Divide(4,4);

	
	for(int tBin = 0; tBin < 8; tBin++){
		
		adcRAWCan->cd(pad);
		adc_RAW_distributions[30][tBin][0]->Draw();
		pad++;
		
	}
	
	//calculate ADC mean in each time bin here
	
	for(int thresh = 0; thresh < numOfThresholdsScanned; thresh++){
		for(int pixel = 0; pixel < 16; pixel++){
			for(int tBin = 0; tBin < 8; tBin++){
	
				double meanADC     = adc_RAW_distributions[thresh][tBin][pixel]->GetMean();
				double meanADC_err = adc_RAW_distributions[thresh][tBin][pixel]->GetMeanError();
				
				adc_mean_distributions[thresh][pixel]->SetBinContent(8 - tBin, meanADC);
				adc_mean_distributions[thresh][pixel]->SetBinError(8 - tBin, meanADC_err);
			}
		}
	}	
	*/

	///////////////////////////////////////
	////// BEGIN PLOT DRAWING BLOCK ///////
	///////////////////////////////////////

	// SCATTER NOISE PLOT FOR ALL 16 PIXELS OF INTEREST
	pad = 1;

	TCanvas *c_scatter = new TCanvas("c_scatter", "16 Scatter Plot", 1200, 800);
	c_scatter->SetGrid(); // Adds a background grid to make centering at 0 easy to verify visually

	// Draw Option "AP": 
	// 'A' draws the overall bounding axes frame
	// 'P' draws the points/markers (creates the scatter plot)
	// If you also want faint trend lines connecting the dots, use "ALP" instead
	mg_drift->Draw("AP");

	// force y axis limits to zoom in around 0
	mg_drift->GetYaxis()->SetRangeUser(-2, 2); 

	// legend
	TLegend *legend = new TLegend(0.85, 0.15, 0.98, 0.85);
	legend->SetHeader("Pixel No.","C");
	for (int pixel_id = 0; pixel_id < 16; pixel_id++) {
		legend->AddEntry(g_pixel_drift[pixel_id], Form("Pixel %d", pixels_of_interest[pixel_id]), "p");
	}
	legend->Draw();

	c_scatter->SaveAs("ADC_Timebin_Drift_Scatter.png");
	
	// 2D MAP OF ADC MEAN CANVAS
	TCanvas * adcMeanCan = new TCanvas("ADC mean distribution per timebin", "ADC mean distribution per timebin", 1600, 1600);
	adcMeanCan->Divide(4,4);	
	
	bool firstHistogram = true;
	
	double peak_ADC[1024];
	int pixel_peak_ADC = -1;
	
	// DRAW SINGLE AVERAGE WAVEFORMS INTO A 4x4 CANVAS
	cout << "\nDrawing 16 averaged pixel tracks onto 4x4 canvas layout..." << endl;
	adcMeanCan->Clear();
	adcMeanCan->Divide(4, 4);

	for (int pixel_id = 0; pixel_id < 16; pixel_id++) {
			int pixel = pixels_of_interest[pixel_id];
			
			// 1. Calculate the average height (baseline offset) of this specific graph
			double total_y = 0;
			int n_points = adc_mean_distributions[pixel]->GetN();
			
			for (int p = 0; p < n_points; p++) {
				total_y += adc_mean_distributions[pixel]->GetPointY(p);
			}
			double graph_baseline = total_y / n_points;

			// 2. Loop through again to subtract that baseline from every single point coordinate
			for (int p = 0; p < n_points; p++) {
				double current_x = adc_mean_distributions[pixel]->GetPointX(p);
				double current_y = adc_mean_distributions[pixel]->GetPointY(p);
				
				// Update the point to sit perfectly relative to 0
				adc_mean_distributions[pixel]->SetPoint(p, current_x, current_y - graph_baseline);
			}

			// Now navigate to the pad and draw the cleanly centered line
			adcMeanCan->cd(pixel_id + 1);
			adc_mean_distributions[pixel]->SetMarkerSize(4); // Sets marker to 1.5 times the default size
			adc_mean_distributions[pixel]->SetMinimum(-0.6);
			adc_mean_distributions[pixel]->SetMaximum(0.6);
			adc_mean_distributions[pixel]->Draw("APX");
	}

	adcMeanCan->SaveAs(Form("ADC_AVERAGE_distributions_trigger_pixel_%d.png", trigger_pixel));
	cout << "SUCCESS: Saved all 16 average pixel paths!" << endl;

	// NOISE HISTOGRAMS: GLOBAL NOISE
	TCanvas *noiseCan = new TCanvas("noiseCan", "Global Noise Spectrum", 1600, 1600);
	noiseCan->Divide(4, 4); // Divide into a 4x4 grid

	for (int pixel_id = 0; pixel_id < 16; pixel_id++) {
		int pixel = pixels_of_interest[pixel_id];
		
		// Navigate to the respective pad slot (1 to 16)
		noiseCan->cd(pixel_id + 1);

		// Style the histogram
		Double_t mean = h_noise_distributions[pixel]->GetMean();

		// Retrieve current axis bounds
		TAxis *xAxis = h_noise_distributions[pixel]->GetXaxis();
		Double_t xMin = xAxis->GetXmin();
		Double_t xMax = xAxis->GetXmax();

		// Shift the axis limits left by the mean value
		xAxis->SetLimits(xMin - mean, xMax - mean);

		h_noise_distributions[pixel]->SetLineColor(kBlue);
		h_noise_distributions[pixel]->SetFillColorAlpha(kBlue, 0.15); 
		h_noise_distributions[pixel]->SetMarkerStyle(20);
		h_noise_distributions[pixel]->SetMarkerSize(0.5);

		// Fit a Gaussian profile to the peak to extract true mean and sigma (RMS noise)
		h_noise_distributions[pixel]->Fit("gaus", "Q");

		// Enable the stats box automatically on the canvas pad
		gStyle->SetOptFit(1111); 
		
		// Draw the structured noise graph
		h_noise_distributions[pixel]->Draw("E HIST"); 
	}
	
	noiseCan->SaveAs("Global_ADC_Noise_Distribution.png");
	cout << "SUCCESS: Saved global noise distributions grid to Global_ADC_Noise_Distribution.png" << endl;


	// DRAW NOISE HISTOGRAMS BY TIMEBIN FOR 2 PIXELS OF INTEREST
	TCanvas *timebinNoiseCan = new TCanvas("timebinNoiseCan", "Timebin Noise Spectrum", 1000, 2400);
	timebinNoiseCan->Divide(2,8);

	for (int pixel = 0; pixel < 2; pixel++) {
		for (int tbin = 0; tbin < 8; tbin++) {
			timebinNoiseCan->cd((pixel * 8) + tbin + 1);
			// Style the histogram
			Double_t mean = h_noise_per_timebin[pixel][tbin]->GetMean();

			// Retrieve current axis bounds
			TAxis *xAxis = h_noise_per_timebin[pixel][tbin]->GetXaxis();
			Double_t xMin = xAxis->GetXmin();
			Double_t xMax = xAxis->GetXmax();

			// Shift the axis limits left by the mean value
			xAxis->SetLimits(xMin - mean, xMax - mean);

			h_noise_per_timebin[pixel][tbin]->SetLineColor(kBlue);
			h_noise_per_timebin[pixel][tbin]->SetFillColorAlpha(kBlue, 0.15); 
			h_noise_per_timebin[pixel][tbin]->SetMarkerStyle(20);
			h_noise_per_timebin[pixel][tbin]->SetMarkerSize(0.5);

			// Fit a Gaussian profile to the peak to extract true mean and sigma (RMS noise)
			h_noise_per_timebin[pixel][tbin]->Fit("gaus", "Q");

			// Enable the stats box automatically on the canvas pad
			gStyle->SetOptFit(1111); 
			
			// Draw the structured noise graph
			h_noise_per_timebin[pixel][tbin]->Draw("E HIST"); 
		}
	}
	
	
	cout << "Extracted ADC sum from 2D histogram = " << trigger_bin_ADC_sum << endl;
	
	if (trigger_bin_ADC_sum != 0.0) {
		double scaleFactor = 1.0 / trigger_bin_ADC_sum;
		ADC_mean_map->Scale(scaleFactor);
	} else {
		cout << "WARNING: trigger_bin_ADC_sum is zero! Skipping scaling operation to prevent Infs." << endl;
	}
	
	// DRAW ADC MAP OF WHOLE CHIP
	TCanvas * ADC_map_whole_chip = new TCanvas("ADC mean map", "ADC mean map", 1600, 800);
	
	double trigger_pixel_total = ADC_mean_map->GetBinContent(2, 3);
	
	cout << "trigger_pixel_total = " << trigger_pixel_total << endl;

	ADC_map_whole_chip->cd(1);
	
	ADC_mean_map->SetStats(0);

	// set text size
	gStyle->SetPaintTextFormat("4.2f"); // Optional: Formats numbers to 2 decimal places so they fit better
	gStyle->SetTextSize(0.025);
	
	ADC_mean_map->Draw("COLZ TEXT");
	
	ADC_map_whole_chip->cd(2);
	
	ADC_map_whole_chip->SaveAs(Form("ADC_MAP_ALL_PIXELS_trigger_pixel_%d.png", trigger_pixel));
	
	
	return;

}

// helper functions to retrieve indices for pixels and canvases
int getPixelIndex(int col, int row){
	return (col*32) + row;
}

int getPixelCanvasIndex(int i){
	return floor(i/4);
	
}
