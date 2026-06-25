//Lightweight EICROC digital output analyzer
//
//
// Author: Alex Jentsch


// Event number (16 of them, pixels (0,0) ... (3,3)
//          
//                                                                                                           TDC ADC  HB
//57876376839;   0;  68;   0;  0;  79;  0;  0;  74;  0;  0;  63;  0;  0;  57;  0;  0;  59;  0;  0;  63;  0;  0;  73;  0

//Data format is read form the end to the beginning (it's stored backwards).
//Format (right to left) is hit bit, ADC, and then TDC.
//each event is 25ns, and there are 8 time bins (the triplets stored back to front).
//There are 16 "copies" of the event number, where each one is the various pixels.

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

After that, the code will begin processing the Sr-90 (or other) data file, byt reading the data lines 
for each pixel, and for each unique event ID into a buffer array, and then storing information in 
ROOT histograms for analysis. 

This part can easily be modified to do things differently.

*********************/


#include <string>

using namespace std;

int getPixelIndex(int i, int j);
int getPixelCanvasIndex(int i);

//Color_t markerColor[16] = {kBlack, kRed+1, kRed, kGreen, kGreen+2, kBlue+1, kBlue+3, kMagenta, kBlue, kBlack, kBlue+1, kBlue+2, kGreen+2, kGreen+3, kMagenta+2, kMagenta+3};
Color_t markerColor[16] = {kGray, kGray+1, kGray+2, kGray+3, kGreen, kGreen+1, kGreen+2, kGreen+3, kBlue, kBlue+1, kBlue+2, kBlue+3, kRed, kRed+1, kRed+2, kRed+3};
//Color_t markerColor[16] = {kBlack, kGreen, kBlue, kRed, kBlack+1, kGreen+1, kBlue+1, kRed+1, kBlack+2, kGreen+2, kBlue+2, kRed+2, kBlack+3, kGreen+3, kBlue+3, kRed+3  };

int markerStyle[16] = {20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 29, 29, 29, 29};

void analyzeCSVData_radSource(TString inputFileName = "", TString inputFileName_pedestal = ""){

	ifstream inputCSVFile;
	ifstream inputCSVFile_pedestal;
	
	bool subtractPedestals = false;
	bool usePedestalFile = false;
	
	bool subtractADC_baseline = false; //this is to set the middle of the ADC distributions to 0 to try and mimic Arzoo's analysis
	
	double thresholdStepSize   = 10;
	const int numOfAnalyzedEvents = 5000; //This sets the size of buffer array for the event by event analysis - meaning, number of event user deems "good"
	
	int trigger_pixel = 5; //This sets the trigger pixel to focus on, but this logic can be disabled in the main code block, if desired
	
	TGraph * s_curve[16]; // storage graphs for s_curves -- not used in this code
	for(int pixel = 0; pixel < 16; pixel++){s_curve[pixel] = new TGraph();}
	
	std::vector<double> threshold_values;
	std::vector<double> efficiency;
	
	int pad = 1; //Sets a pad number that you can use later for drawing histograms on a divide canvas
	
	//various histograms for storing information for a 4x4 pixel array - can be modified to handle 32x32 pixel array
	
	TH1D * adc_distributions[4][4];
	TH1D * tdc_distributions[4][4];
	
	TH1D * adc_distributions_PEDESTAL[4][4];

	// number of thresholds, number of time bins, number of pixels
	TH1D * adc_RAW_distributions[numOfAnalyzedEvents][8][16];
	TH1D * adc_RAW_distributions_pedestal[8][16];
	
	TH1D * hitPixel = new TH1D("hit_pixel", "hit_pixel", 16, 0.0, 16.0);
	TH2D * ADC_mean_map = new TH2D("ADC_mean_map", "ADC_mean_map", 4, 0, 4, 4, 0, 4);
	
	TH2D * hitBitMap_trigger_pixel = new TH2D(Form("hitBitMap_trigger_pixel_%d", trigger_pixel), Form("hitBitMap_trigger_pixel_%d", trigger_pixel), 4, 0, 4, 4, 0, 4);
	
	TH1D * peak_ADC_trigger_pixel = new TH1D("peak_ADC_trigger_pixel", "peak_ADC_trigger_pixel", 8, 0.0, 8.0);
	
	// number of thresholds, AVERAGED OVER TIME BINS, number of pixels
	TGraph * adc_mean_distributions[numOfAnalyzedEvents][16];
	TGraph * adc_mean_distributions_PEDESTAL[numOfAnalyzedEvents][16];

	TGraphErrors * pedestal_waveforms[16];

	for(int tBin = 0; tBin < 8; tBin++){
		for(int pixel = 0; pixel < 16; pixel++){
			
			TString title;
			title.Form("adc_max_distribution_pixel_%d_timeBin_%d_PEDESTAL", pixel, tBin);

			adc_RAW_distributions_pedestal[tBin][pixel] = new TH1D(title, "; ADC value [DACu]; counts", 256, 0, 255);
			adc_RAW_distributions_pedestal[tBin][pixel]->SetTitle(title);
				
		}
	}
	

	for(int ev = 0; ev < numOfAnalyzedEvents; ev++){
		for(int tBin = 0; tBin < 8; tBin++){
			for(int pixel = 0; pixel < 16; pixel++){

				TString title;
				title.Form("adc_max_distribution_thresh_%d_pixel_%d_timeBin_%d", ev, pixel, tBin);

				adc_RAW_distributions[ev][tBin][pixel] = new TH1D(title, "; ADC value [DACu]; counts", 256, 0, 255);
				adc_RAW_distributions[ev][tBin][pixel]->SetTitle(title);
				
	
			}
		}
	}	
	
	/*
	
	//used to get the ADC mean distributions to draw  the ADC waveforms --> good for charge injection, not for Sr-90 source
	
	for(int ev = 0; ev < numOfAnalyzedEvents; ev++){
		for(int pixel = 0; pixel < 16; pixel++){

			TString title;
			title.Form("adc_MEAN_distribution_thresh_%d_pixel_%d", ev, pixel);

			adc_mean_distributions[ev][pixel] = new TGraph();
			
			adc_mean_distributions[ev][pixel]->GetXaxis()->SetTitle("time bin [time bin = 25ns"); 
			adc_mean_distributions[ev][pixel]->GetYaxis()->SetTitle("ADC value [DACu]");
			adc_mean_distributions[ev][pixel]->SetTitle(title);
	
		}
	}	
	*/		//adc_mean_distributions[61][16] = new TH1D(Form("adc_max_distribution_pixel_%d%d", i, j), "counts; ADC value [DACu]", 256, 0, 255);

	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
		
			adc_distributions[i][j] = new TH1D(Form("adc_max_distribution_pixel_%d%d", i, j), "counts; ADC value [DACu]", 256, 0, 255);
			tdc_distributions[i][j] = new TH1D(Form("tdc_distribution_pixel_%d%d", i, j), "counts; TDC value [DACu]", 1024, 0, 1023);
				
			adc_distributions_PEDESTAL[i][j] = new TH1D(Form("adc_max_distribution_pixel_%d%d_PEDESTAL", i, j), "counts; ADC value [DACu]", 256, 0, 255);
		}
	}
	
	double numFiles = 0;
	
	int charge_DACu = 45;
	
	int numTriggeredEvents = 0;
	
	//older data structure to hold pre-calculated pedstals -- not used right now, but left here for comparison testing, if needed.
	
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
	 
	if you do this, you will need to make sure your file name is parsed properly
	
	Like this --> inputFileName.Form("[path_to_data_files]/[my_file_name]_%.0f_1.csv", fileIdx);
	
	then you can loop over files from the same run if you broke them up into smaller pieces, or did a scan, etc.
	
	***********************/
	
	for(int fileIdx = 0; fileIdx < 1; fileIdx++){ //BEGIN for loop over file indices
	
		//These are just older files and left here for quick tests and comparisons
		//
		//
	
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_LPC_recentered_source_5_27_26_0_1.csv");
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_LPC_recentered_source_06_02_26_1_1.csv");
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_LPC_recentered_source_06_02_26_2_1.csv");
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_LPC_recentered_source_longer_delay_06_17_2026_1.csv");
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_LPC_recentered_source__longer_delay_no_CMD_06_17_2026_1.csv");
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_LPC_recentered_source_06_10_26_2_1.csv");
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_LPC_recentered_source_6_16_2026_1.csv");
			//inputFileName_pedestal.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_LPC_recentered_source_ped_06_02_26_0_1.csv");
			//inputFileName_pedestal.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_LPC_recentered_source_PEDESTAL_using_charge_inj_approach_06052026_1.csv");
			//inputFileName_pedestal.Form("radioactive_source_digital_data/sr90_unbiased_330threshold_LPC_recentered_source_PEDESTAL_using_charge_inj_approach_06_10_26_0_1.csv");
		
		
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_TRGOUT_delay_6_timebins_06_17_2026_1.csv");
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_TRGOUT_delay_4_timebins_06_17_2026_1.csv");
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_TRGOUT_delay_3_timebins_06_17_2026_1.csv");
		
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_TRGOUT_delay_2_timebins_CALIBRATIONS_06_23_2026_1.csv");
		
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_TRGOUT_delay_2_timebins_CALIBRATIONS_no_clk_gating_06_23_2026_1.csv");
		
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_TRGOUT_delay_1_timebins_06_18_2026_1.csv");
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_TRGOUT_delay_1_timebins_NO_Vref_CORR_06_18_2026_1.csv");
		
			//inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_TRGOUT_delay_4_timebins_NO_Vref_CORR_06_18_2026_1.csv");
			//inputFileName_pedestal.Form("radioactive_source_digital_data/sr90_170Vbias_330threshold_NO_CMD_PULSE_PEDESTAL_06_17_2026_1.csv");
		
		
		//Input files go here --> one is for a pedestal file, the other is for the main data file
		//inputFileName_pedestal.Form("radioactive_source_digital_data/sr90_unbiased_330threshold_NO_CMD_PULSE_PEDESTAL_06_17_2026_1.csv");
		inputFileName.Form("radioactive_source_digital_data/sr90_170Vbias_400threshold_TRGOUT_delay_2_timebins_CALIBRATIONS_no_clk_gating_06_24_2026_1.csv");
		
		cout << inputFileName << endl;
		cout << inputFileName_pedestal << endl;
	
		inputCSVFile.open(inputFileName.Data());
		inputCSVFile_pedestal.open(inputFileName_pedestal.Data());

		if(!inputCSVFile){ 
			cout << "File does not exist, skipping...." << endl; 
			return;
		}
		
		if(!inputCSVFile_pedestal){ 
			cout << "PEDESTAL File does not exist, using internal pedestal map..." << endl; 
			usePedestalFile = false;
		}

		numFiles++;

		string eventNumber[4][4];
		double numGoodEvents[4][4];
		double ADC_values[8];
		int TDC_values[8];
		int hitBits[8];
		
		double ADC_values_event[16][8];
		int TDC_values_event[16][8];
		int hitBits_event[16][8];

		double ADC_sum_all_timebins = 0.0;

		double ADC_max_pixel[4][4];
		double ADC_max_whole_ASIC = 0;
		int pixel_with_max_ADC_ASIC = 0;

		string eventNumber_str;
		string ADC_values_str[8];
		string TDC_values_str[8];
		string hitBits_str[8];

		double tmp_ADC_MEAN_pedestal[8];
		double tmp_ADC_MEAN_err_pedestal[8];
		
		double tmp_ADC_MEAN_pedestal_BINS[8]      = {0, 1, 2, 3, 4, 5, 6, 7};
		double tmp_ADC_MEAN_err_pedestal_x_err[8] = {0, 0, 0, 0, 0, 0, 0, 0};
		
		double avg_ADC_offset_pedestal_pixel[16];
		double max_ADC_pedestal_pixel[16];
		
		int pixel_has_hit[16];

		int timeBin     = 7;
		int pixelRow    = 0;
		int pixelColumn = 0;
		
		int ADC_max_time_bin = 0;

		string tmp;

		int numEvents = 0;
		
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				numGoodEvents[i][j] = 0;
			}
		}
		
		int numHitBitSet = 0;
		bool hitBitInTriggerPixel = false;
		bool hitBitInCentralTimeBins = false;
		bool hitBitPixelFifteen = false;
		
		bool pedestalNotSynched = true;

		int hitBitPixel[4][4];

		///////////////////////////////////////////////////////////////////
		//////// BEGINNING OF PEDESTAL CODE BLOCK ///////////////////////////////
		///////////////////////////////////////////////////////////////////

		while(!inputCSVFile_pedestal.eof() ){ //BEGIN while loop over PEDESTAL file
	
			if (getline(inputCSVFile_pedestal, eventNumber_str, ';') &&
				getline(inputCSVFile_pedestal, TDC_values_str[7], ';') && getline(inputCSVFile_pedestal, ADC_values_str[7], ';') && getline(inputCSVFile_pedestal, hitBits_str[7], ';') &&
				getline(inputCSVFile_pedestal, TDC_values_str[6], ';') && getline(inputCSVFile_pedestal, ADC_values_str[6], ';') && getline(inputCSVFile_pedestal, hitBits_str[6], ';') &&
				getline(inputCSVFile_pedestal, TDC_values_str[5], ';') && getline(inputCSVFile_pedestal, ADC_values_str[5], ';') && getline(inputCSVFile_pedestal, hitBits_str[5], ';') &&
				getline(inputCSVFile_pedestal, TDC_values_str[4], ';') && getline(inputCSVFile_pedestal, ADC_values_str[4], ';') && getline(inputCSVFile_pedestal, hitBits_str[4], ';') &&
				getline(inputCSVFile_pedestal, TDC_values_str[3], ';') && getline(inputCSVFile_pedestal, ADC_values_str[3], ';') && getline(inputCSVFile_pedestal, hitBits_str[3], ';') &&
				getline(inputCSVFile_pedestal, TDC_values_str[2], ';') && getline(inputCSVFile_pedestal, ADC_values_str[2], ';') && getline(inputCSVFile_pedestal, hitBits_str[2], ';') &&
				getline(inputCSVFile_pedestal, TDC_values_str[1], ';') && getline(inputCSVFile_pedestal, ADC_values_str[1], ';') && getline(inputCSVFile_pedestal, hitBits_str[1], ';') &&
				getline(inputCSVFile_pedestal, TDC_values_str[0], ';') && getline(inputCSVFile_pedestal, ADC_values_str[0], ';') && getline(inputCSVFile_pedestal, hitBits_str[0]))
	
			{

				eventNumber[pixelColumn][pixelRow] = eventNumber_str; //stoi(eventNumber_str);

				//////////////////////////////////////////////////////////////////////////////////////////
				//loop through the events and buffer ADC and TDC values, and count up the hitBits == 1
				//////////////////////////////////////////////////////////////////////////////////////////

				for(int tBin = 0; tBin < 8; tBin++){ //loop over 8 slices in one event's time bin

					ADC_values[tBin] = stoi(ADC_values_str[tBin]);
	    			TDC_values[tBin] = stoi(TDC_values_str[tBin]);
	    			hitBits[tBin]    = stoi(hitBits_str[tBin]);
				
					int pixel = getPixelIndex(pixelColumn, pixelRow);
				
					ADC_values_event[pixel][tBin] = ADC_values[tBin];
					TDC_values_event[pixel][tBin] = TDC_values[tBin];
					hitBits_event[pixel][tBin]    = hitBits[tBin];
								
					if(hitBits[tBin] == 1){ 
						numHitBitSet++; 
						
					}
				
				}
					
				//////////////////////////////////////////////////////////////////////////////////////////
				//once all 16 pixels looped, check for NO hitbits in event, if true, use for pedestal
				//////////////////////////////////////////////////////////////////////////////////////////
	
				pixelRow++;
				if(pixelRow == 4 && pixelColumn < 3) { pixelColumn++; pixelRow = 0;}
				if(pixelRow == 4 && pixelColumn == 3){ pixelColumn++;}
				if(pixelRow == 4 && pixelColumn == 4){ 
					
					if(numHitBitSet == 0 ){
						
						numTriggeredEvents++;
						
						for(int i = 0; i < 4; i++){
							for(int j = 0; j < 4; j++){
								
								int pixel = getPixelIndex(i, j);
									
								//cout << eventNumber[i][j] << "\t";
																
								for(int tBin = 0; tBin < 8; tBin++){
																
									if(pixel != -1){
										
										//cout << TDC_values_event[pixel][tBin] << ", " << ADC_values_event[pixel][tBin] << ", " << hitBits_event[pixel][tBin] << ", \t";
										
										adc_RAW_distributions_pedestal[tBin][pixel]->Fill(ADC_values_event[pixel][tBin]);
										
										//if(tBin == 7){ cout << "\n"; }
									}
								}
							}
						}
						
						//cout << "\n\n";
					}
					
					pixelRow = 0; pixelColumn = 0; numHitBitSet = 0;
					
				    numEvents++;

				}
			}

		} //END while loop over PEDESTAL file

		inputCSVFile_pedestal.close();

		cout << "number of pedestal events --> " << numTriggeredEvents << endl;
		cout << "\n\n";
		
		if(usePedestalFile){ //BEGIN conditional to print pedestal plots --> only used if you had a pedestal file, and set the flag

			TCanvas * adcRAWCan[8]; 
		
			for(int tBin = 0; tBin < 8; tBin++){
				adcRAWCan[tBin] = new TCanvas(Form("can_tBin_%d", tBin), Form("can_tBin_%d", tBin), 2000, 1600);
				adcRAWCan[tBin]->Divide(4,4);
			
			}
		
			pad = 1;
		
			adcRAWCan[0]->cd(1);
		
			for(int pixel = 0; pixel < 16; pixel++){
			
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
		
		
			TF1 * pedFits[16];
		
			for(int pixel = 0; pixel < 16; pixel++){
			
			
				int idx = getPixelCanvasIndex(pad);
			
				pedestal_fit_canvas->cd(idx);
			
				//pedestal_waveforms[pixel]->GetXaxis()->SetRangeUser(-1, 8);
			
				pedFits[pixel] = new TF1(Form("pixel_%d_ped_fit", pixel), "[0] + [1]*sin([2]*x + [3])", 0, 7);
			
				double amplitude =  max_ADC_pedestal_pixel[pixel] - avg_ADC_offset_pedestal_pixel[pixel];
			
				cout << "pixel " << pixel << " avg ADC = " << avg_ADC_offset_pedestal_pixel[pixel] << " and max ADC = " << max_ADC_pedestal_pixel[pixel] << " possible amp = " << amplitude << endl;
			
			
				pedFits[pixel]->SetParLimits(1, 25, 35);
				pedFits[pixel]->SetParLimits(0, 35, 45);
			
				//pedFits[pixel]->SetParLimits(1, amplitude-5, 2*amplitude);
				//pedFits[pixel]->SetParLimits(0, avg_ADC_offset_pedestal_pixel[pixel]-5, avg_ADC_offset_pedestal_pixel[pixel]+5);
			
				//pedFits[pixel]->SetParameter(0, avg_ADC_offset_pedestal_pixel[pixel]);
				//pedFits[pixel]->SetParameter(1, amplitude);
			
			
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

		numEvents = 0;
		
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				numGoodEvents[i][j] = 0;
			}
		}
		
		
		
		
		numHitBitSet = 0;
		numTriggeredEvents = 0;
			
		
		///////////////////////////////////////////////////////////////////
		//////// BEGINNING OF MAIN DATA CODE BLOCK ///////////////////////////////
		///////////////////////////////////////////////////////////////////
		


		while(!inputCSVFile.eof()){//BEGIN while loop over main data file
	
			if (getline(inputCSVFile, eventNumber_str, ';') &&
				getline(inputCSVFile, TDC_values_str[7], ';') && getline(inputCSVFile, ADC_values_str[7], ';') && getline(inputCSVFile, hitBits_str[7], ';') &&
				getline(inputCSVFile, TDC_values_str[6], ';') && getline(inputCSVFile, ADC_values_str[6], ';') && getline(inputCSVFile, hitBits_str[6], ';') &&
				getline(inputCSVFile, TDC_values_str[5], ';') && getline(inputCSVFile, ADC_values_str[5], ';') && getline(inputCSVFile, hitBits_str[5], ';') &&
				getline(inputCSVFile, TDC_values_str[4], ';') && getline(inputCSVFile, ADC_values_str[4], ';') && getline(inputCSVFile, hitBits_str[4], ';') &&
				getline(inputCSVFile, TDC_values_str[3], ';') && getline(inputCSVFile, ADC_values_str[3], ';') && getline(inputCSVFile, hitBits_str[3], ';') &&
				getline(inputCSVFile, TDC_values_str[2], ';') && getline(inputCSVFile, ADC_values_str[2], ';') && getline(inputCSVFile, hitBits_str[2], ';') &&
				getline(inputCSVFile, TDC_values_str[1], ';') && getline(inputCSVFile, ADC_values_str[1], ';') && getline(inputCSVFile, hitBits_str[1], ';') &&
				getline(inputCSVFile, TDC_values_str[0], ';') && getline(inputCSVFile, ADC_values_str[0], ';') && getline(inputCSVFile, hitBits_str[0]))
	
			{
	
				eventNumber[pixelColumn][pixelRow] = eventNumber_str;

				//int numHitBitSet = 0;
				
				double ADC_value_at_timebin = 0.0;
				
				ADC_sum_all_timebins = 0.0;
				
				ADC_max_pixel[pixelColumn][pixelRow] = 0;
				hitBitPixel[pixelColumn][pixelRow] = 0;

				for(int tBin = 0; tBin < 8; tBin++){ //loop over 8 slices in one event's time bin

					ADC_values[tBin] = stoi(ADC_values_str[tBin]); //converts strings to integer numbers 
	    			TDC_values[tBin] = stoi(TDC_values_str[tBin]);
	    			hitBits[tBin]    = stoi(hitBits_str[tBin]);
				
					adc_distributions[pixelColumn][pixelRow]->Fill(ADC_values[tBin]); //fills the 
				
					int pixel = getPixelIndex(pixelColumn, pixelRow);
				
					ADC_values_event[pixel][tBin] = ADC_values[tBin];
					TDC_values_event[pixel][tBin] = TDC_values[tBin];
					hitBits_event[pixel][tBin]    = hitBits[tBin];
				
					ADC_sum_all_timebins += ADC_values[tBin];
				
					if(tBin == 0) { ADC_value_at_timebin = ADC_values[tBin]; }
				
					if(ADC_values[tBin] > ADC_max_pixel[pixelColumn][pixelRow]) { 
						ADC_max_time_bin = tBin;
						ADC_max_pixel[pixelColumn][pixelRow] = ADC_values[tBin];
					}
				
					if(hitBits[tBin] == 1){ 
						numHitBitSet++; 
						if(pixel == trigger_pixel)  { hitBitInTriggerPixel = true; hitBitInCentralTimeBins = true; }
						if(pixel == 15) { hitBitPixelFifteen = true; }
						
						hitBitPixel[pixelColumn][pixelRow]++;// = 1;
						
						hitPixel->Fill(pixel + 0.5);
						
					}
				
				}
					
				//double ADC_baseline = ADC_sum_all_timebins/8.0;
				double ADC_baseline = ADC_value_at_timebin;
			
				if(hitBitPixel[pixelColumn][pixelRow] > 1){ //This is quickly see if some events/pixels are noisy -- meaning they have #hitBits == 1 > 1
					
					int pixel = getPixelIndex(pixelColumn, pixelRow);
					
					cout << "numHitBits in pixel "<< pixel << " = " << hitBitPixel[pixelColumn][pixelRow] <<  endl;
				}
			
				if(subtractADC_baseline){ //This is to subtract an ADC baseline to make all ADC distributions have the same baseline -- it's not necessary, so there is a flag
					
					int pixel = getPixelIndex(pixelColumn, pixelRow);
				
					for(int tBin = 0; tBin < 8; tBin++){ 
						ADC_values_event[pixel][tBin] = ADC_values_event[pixel][tBin] - ADC_baseline;
					}
				}
			
	
				pixelRow++;
				
				//once you iterate through the last (row, column) the analysis for the event will begin
				
				if(pixelRow == 4 && pixelColumn < 3) { pixelColumn++; pixelRow = 0;}
				if(pixelRow == 4 && pixelColumn == 3){ pixelColumn++;}
				if(pixelRow == 4 && pixelColumn == 4){ 
					
					double ADC_far_pixel_one[8];
					double ADC_far_pixel_two[8];
					double ADC_far_pixel_three[8];
					double ADC_far_pixel_four[8];
					double ADC_far_pixel_five[8];
					double ADC_far_pixel_six[8];
					double ADC_far_pixel_seven[8];
					double avg_ADC_pedestal[8];
					
					/*********************
					
					This block is used to try and perform "far" pixel pedestal substraction
					This is done because the noise is driven by the clock, and you need to ensure the noise
					pedestal has the same phase or the pedestal subtration will behave poorly
					Alternatively, you can perform pedestal subtraction from a file, especially useful
					for charge injections where you will always have the same phase
					
					*******************/
					
					for(int i = 0; i < 4; i++){
						for(int j = 0; j < 4; j++){
							
							int pixel = getPixelIndex(i, j);
							
							for(int tBin = 0; tBin < 8; tBin++){
								
								
								if(pixel == 3){ ADC_far_pixel_one[tBin]   = ADC_values_event[pixel][tBin]; }
								if(pixel == 7){ ADC_far_pixel_two[tBin]   = ADC_values_event[pixel][tBin]; }
								if(pixel == 11){ ADC_far_pixel_three[tBin]   = ADC_values_event[pixel][tBin]; }
								if(pixel == 12){ ADC_far_pixel_four[tBin]   = ADC_values_event[pixel][tBin]; }
								if(pixel == 13){ ADC_far_pixel_five[tBin]   = ADC_values_event[pixel][tBin]; }
								if(pixel == 14){ ADC_far_pixel_six[tBin]   = ADC_values_event[pixel][tBin]; }
								if(pixel == 15){ 
									ADC_far_pixel_seven[tBin] = ADC_values_event[pixel][tBin];
									
									if(trigger_pixel == 5){
										avg_ADC_pedestal[tBin] = (ADC_values_event[3][tBin] + ADC_values_event[7][tBin] + ADC_values_event[11][tBin]) +
										                     (ADC_values_event[12][tBin] + ADC_values_event[13][tBin]);
										avg_ADC_pedestal[tBin] = avg_ADC_pedestal[tBin]/5.0;
									}
									
									/*
									if(trigger_pixel == 5){
										avg_ADC_pedestal[tBin] = ADC_values_event[15][tBin];
									}
									*/
									if(trigger_pixel == 2){
										avg_ADC_pedestal[tBin] = ADC_values_event[12][tBin];
									}
									/*
									if(trigger_pixel == 5){
										avg_ADC_pedestal[tBin] = (ADC_values_event[7][tBin] +ADC_values_event[13][tBin]);
										avg_ADC_pedestal[tBin] = avg_ADC_pedestal[tBin]/2.0;
									}
									*/
									if(trigger_pixel == 6){
										avg_ADC_pedestal[tBin] = ADC_values_event[0][tBin];
										
									}
									if(trigger_pixel == 10){
										avg_ADC_pedestal[tBin] = ADC_values_event[0][tBin];
										
									}
									if(trigger_pixel == 12){
										avg_ADC_pedestal[tBin] = ADC_values_event[3][tBin];
										
									}
									//else if(trigger_pixel != 5 && trigger_pixel != 6 && trigger_pixel != 10){avg_ADC_pedestal[tBin] = 0.0;}
									
								}
							}
						}
					}
					
					
					
					/*********************
					
					This block is where the pedestal subtration is actually performed. 
					Flags are set at the top of the code to pick which method you want to use, 
					or perform no pedestal subtraction.
					
					*******************/
					for(int i = 0; i < 4; i++){
						for(int j = 0; j < 4; j++){
							
							int pixel = getPixelIndex(i, j);
							
							for(int tBin = 0; tBin < 8; tBin++){
					
								
								int pedastalSubtracted = ADC_values_event[pixel][tBin];
								
								//if(usePedestalFile){ pedastalSubtracted = ADC_values_event[pixel][tBin] - adc_RAW_distributions_pedestal[tBin][pixel]->GetMean(); }
					
								if(subtractPedestals && !usePedestalFile){ pedastalSubtracted = ADC_values_event[pixel][tBin] - (int)avg_ADC_pedestal[tBin]; }
								if(subtractPedestals && usePedestalFile){ pedastalSubtracted = ADC_values_event[pixel][tBin] - adc_RAW_distributions_pedestal[tBin][pixel]->GetMean(); }
					
								if(pedastalSubtracted > ADC_max_whole_ASIC){
									ADC_max_whole_ASIC = pedastalSubtracted;
									pixel_with_max_ADC_ASIC = getPixelIndex(i, j);
				
								}
							}
						}
					}
					
					if(numHitBitSet == 1 && hitBitInTriggerPixel ){ peak_ADC_trigger_pixel->Fill(ADC_max_time_bin+0.5); } //This just fills a QA histogram with the peak ADC values
					
					/***************
					
					This code block is where we decide if an event is "good" - which is based on your needs.
					The various comments are there for different options used previously by me.
					
					***************/
					
					//if(numHitBitSet > 0 &&  hitBitInTriggerPixel && !hitBitPixelFifteen && hitBitInCentralTimeBins && numTriggeredEvents < numOfAnalyzedEvents){
					//if(numHitBitSet > 0 &&  hitBitInTriggerPixel && !hitBitPixelFifteen && numTriggeredEvents < numOfAnalyzedEvents){
					//if(numHitBitSet == 1 &&  hitBitInTriggerPixel && numTriggeredEvents < numOfAnalyzedEvents && pixel_with_max_ADC_ASIC == trigger_pixel){
					//if(numHitBitSet == 1 &&  hitBitInTriggerPixel && numTriggeredEvents < numOfAnalyzedEvents && ADC_max_time_bin == 1){
					//if(numHitBitSet == 1 &&  hitBitInTriggerPixel && numTriggeredEvents < numOfAnalyzedEvents){
					
					if(numHitBitSet > 0 &&  numTriggeredEvents < numOfAnalyzedEvents ){// BEGIN good triggered event
						
						//counter for the number of events which eventually passed your cuts -- remember, you set the upper-bound at the top of the code in numOfAnalyzedEvents
						numTriggeredEvents++; 
						
						for(int i = 0; i < 4; i++){
							for(int j = 0; j < 4; j++){
								
								TF1 * dataFit = new TF1("data_fit", "[0] + [1]*sin([2]*x + [3])", 0, 7);
								
								int pixel = getPixelIndex(i, j);
								
								
								
								if(hitBitPixel[i][j] == 1){ hitBitMap_trigger_pixel->Fill(i+0.5, (3-j)+0.5); }
								
								TString title;
								title.Form("adc_MEAN_distribution_event_%d_pixel_%d", numTriggeredEvents-1, pixel);

								adc_mean_distributions[numTriggeredEvents-1][pixel] = new TGraph();
			
								adc_mean_distributions[numTriggeredEvents-1][pixel]->GetXaxis()->SetTitle("time bin [time bin = 25ns"); 
								adc_mean_distributions[numTriggeredEvents-1][pixel]->GetYaxis()->SetTitle("ADC value [DACu]");
								adc_mean_distributions[numTriggeredEvents-1][pixel]->SetTitle(title);
								
								int ADC_minimum = 0;
								
								for(int tBin = 0; tBin < 8; tBin++){
																
									if(pixel != -1){
										
										//if(tBin == 0) { cout << eventNumber[i][j] << "\t"; }
										//cout << TDC_values_event[pixel][tBin] << ", "; 
										
										double ADC_tmp = ADC_values_event[pixel][tBin];
										double pedestal_tmp = ADC_tmp - pedestal_map[pixel][tBin];
										
										int pedastalSubtracted = ADC_values_event[pixel][tBin];
										
										
										if(subtractPedestals && !usePedestalFile){ pedastalSubtracted = ADC_values_event[pixel][tBin] - (int)avg_ADC_pedestal[tBin]; }
										if(subtractPedestals && usePedestalFile){ pedastalSubtracted = ADC_values_event[pixel][tBin] - adc_RAW_distributions_pedestal[tBin][pixel]->GetMean(); }
										
										//cout << pedastalSubtracted << ", " << hitBits_event[pixel][tBin] << ", \t";
										
										if(pedastalSubtracted < ADC_minimum){ ADC_minimum = pedastalSubtracted;}
										
										adc_mean_distributions[numTriggeredEvents-1][pixel]->AddPoint(tBin, pedastalSubtracted);
										
										if(tBin == 7){ 
											adc_mean_distributions[numTriggeredEvents-1][pixel]->SetMinimum(ADC_minimum); 
											//cout << "ADC_minimum = " << ADC_minimum << endl;
											//cout << "\n";
										
										}
										
										/*
										dataFit->SetParLimits(1, 25, 35);
										dataFit->SetParLimits(0, 35, 45);
			
										dataFit->SetParLimits(2, 1.0, 2.0);
										dataFit->SetParLimits(3, -TMath::PiOver2(), TMath::PiOver2());
			
										adc_mean_distributions[numTriggeredEvents-1][pixel]->Fit(dataFit, "RQ");
										adc_mean_distributions[numTriggeredEvents-1][pixel]->Draw("ALP");
										*/
										pad++;
										
									}
								}
							}
						}
						
						//cout << "\n\n";
						
					}// END good triggered event
					
					//reset internal counter and flags
					pixelRow = 0; pixelColumn = 0; numHitBitSet = 0; hitBitInTriggerPixel = false; hitBitPixelFifteen = false;
					
					numEvents++; //global event counter, different from the "good" event counter

				}
			}

		} //END while loop over main data file

		cout << "Number of events total recorded : " << numEvents << endl;

		inputCSVFile.close();

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
		
		
		//cout << "total number of events = " << numEvents << endl;

		cout << "\n --- number of good events in pixels --- " << endl;
		cout << " THRESHOLD = " << 100 + fileIdx*thresholdStepSize << endl;
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
			
			
				int pixel = getPixelIndex(i, j);
				
				if(pixel != -1){	
					s_curve[pixel]->AddPoint(100 + fileIdx*thresholdStepSize, numGoodEvents[i][j]/numEvents);
				}
			
				cout << "pixel( " << i << ", " << j << ") = " << numGoodEvents[i][j]/numEvents << endl;
				
			}
		}

	} //end of file loop
	
	//s_curve = new TGraph(threshold_values.size(), &threshold_values[0], &efficiency[0]);
	
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
	
	
	
	**********/
	
	/***********
	
	histogramming and plotting of data and information from Sr-90 scans
	
	************/
	
	
	
	
	TCanvas * can1 = new TCanvas("canv1", "canv1", 500, 500);
	
	hitPixel->Draw();
	
	pad = 1;
	
	TCanvas * adcMeanCan = new TCanvas("canv6", "canv6", 1600, 1600);
	adcMeanCan->Divide(4,4);			
	
	
	bool firstHistogram = true;
	
	double peak_ADC[16];
	int pixel_peak_ADC = -1;
	
	for(int ev = 0; ev < numTriggeredEvents; ev++){
		
		//cout << "\n\n";
		//cout << "triggered event " << ev << endl;
		
		//cout << "pixel 5 minimum = " << adc_mean_distributions[ev][5]->GetMinimum() << endl;
		
		//if(adc_mean_distributions[ev][5]->GetMinimum() > -30) {continue;}
		
		int peak_ADC_whole_sensor = 0;
		bool peak_in_pixel_five = false;
		
		for(int pixel = 0; pixel < 16; pixel++){
			for(int tBin = 0; tBin < 8; tBin++){
			
				if(adc_mean_distributions[ev][pixel]->GetPointY(tBin) < peak_ADC_whole_sensor){ 
					peak_ADC_whole_sensor = adc_mean_distributions[ev][pixel]->GetPointY(tBin);
					if(pixel == trigger_pixel){peak_in_pixel_five = true;}
					else peak_in_pixel_five = false;
				}
			}
		}
		
		//if(!peak_in_pixel_five ){continue;}
		//if(!peak_in_pixel_five ){continue;}
		
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
		
				int pixel = getPixelIndex(i, j);
				int canvasBin = getPixelIndex(j, i);
				
				peak_ADC[pixel] = 0.0;
		
				for(int tBin = 0; tBin < 8; tBin++){
					
					if( adc_mean_distributions[ev][pixel]->GetPointY(tBin) > peak_ADC[pixel]){peak_ADC[pixel] = adc_mean_distributions[ev][pixel]->GetPointY(tBin);}
					
				}
			
				//cout << "pixel " << pixel << "peak ADC = " << peak_ADC[pixel] << endl;
			}
		}
		
		double tmp_max = peak_ADC[0];
		
		for(int pixel = 0; pixel < 16; pixel++){
			
			if(peak_ADC[pixel] > tmp_max){
				
				tmp_max = peak_ADC[pixel];
				pixel_peak_ADC = pixel;
				
			}
			
		}
		
		//if(subtractPedestals && (pixel_peak_ADC != trigger_pixel) ){continue;}
		
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
		
				int pixel = getPixelIndex(i, j);
				int canvasBin = getPixelIndex(j, i);
			
				
				
				
				ADC_mean_map->Fill(i+0.5, (3-j)+0.5, -1*peak_ADC[pixel]);
				
				if(ev < 10){
					
					cout << "Event " << ev << " --> pixel 0 baseline mean = " << adc_mean_distributions[ev][pixel]->GetMean(2) << endl;
					
					
				}
			
		   	 	adcMeanCan->cd(canvasBin+1);
				if(firstHistogram){
					adc_mean_distributions[ev][pixel]->GetYaxis()->SetRangeUser(-120, 256);
					adc_mean_distributions[ev][pixel]->Draw("ALP");
				}
				else
				{
					adc_mean_distributions[ev][pixel]->GetYaxis()->SetRangeUser(-120, 256);
					adc_mean_distributions[ev][pixel]->SetLineColor(markerColor[pixel]);
					adc_mean_distributions[ev][pixel]->SetMarkerColor(markerColor[pixel]);
					adc_mean_distributions[ev][pixel]->Draw("LP SAME");
				}
			}
	    }
		
		
		
		
		
		firstHistogram = false;
	}
	
	adcMeanCan->SaveAs(Form("ADC_distributions_ALL_PIXELS_trigger_pixel_%d.png", trigger_pixel));
	
	double trigger_bin_ADC_sum = 0.0;
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
	
			int pixel = getPixelIndex(i, j);
			
			if(pixel == trigger_pixel){
				
				trigger_bin_ADC_sum = ADC_mean_map->GetBinContent(i+1, (3 - j) + 1);
				
				cout << "trigger pixel bin (i,j) = " << i << ", " << j << endl;
				
			}
		}
	}
	
	cout << "Extracted ADC sum from 2D histogram = " << trigger_bin_ADC_sum << endl;
	
	double scaleFactor = 1.0/trigger_bin_ADC_sum;
	
	ADC_mean_map->Scale(scaleFactor);
	
	TCanvas * ADC_map_whole_chip = new TCanvas("canv7", "canv7", 1600, 800);
	ADC_map_whole_chip->Divide(2,1);
	
	double trigger_pixel_total = ADC_mean_map->GetBinContent(2, 3);
	
	cout << "trigger_pixel_total = " << trigger_pixel_total << endl;
	
	TPaveText * neighbor_1_text = new TPaveText(0.53, 0.52, 0.6, 0.55, "NB NDC");
	neighbor_1_text->SetFillColor(0);
	neighbor_1_text->SetTextColor(kRed);
	neighbor_1_text->AddText(Form("%.3f", ADC_mean_map->GetBinContent(3,3)/trigger_pixel_total));
	
	TPaveText * neighbor_2_text = new TPaveText(0.19, 0.52, 0.25, 0.55, "NB NDC");
	neighbor_2_text->SetFillColor(0);
	neighbor_2_text->SetTextColor(kRed);
	neighbor_2_text->AddText(Form("%.3f", ADC_mean_map->GetBinContent(1,3)/trigger_pixel_total));
	
	for(int i = 1; i < 5; i++){ ADC_mean_map->GetXaxis()->SetBinLabel(i, Form("%d", i-1)); hitBitMap_trigger_pixel->GetXaxis()->SetBinLabel(i, Form("%d", i-1));}
	for(int j = 1; j < 5; j++){ ADC_mean_map->GetYaxis()->SetBinLabel(j, Form("%d", 4-j)); hitBitMap_trigger_pixel->GetYaxis()->SetBinLabel(j, Form("%d", 4-j));}
	
	ADC_map_whole_chip->cd(1);
	
	ADC_mean_map->SetStats(0);
	
	ADC_mean_map->Draw("COLZ TEXT");
	//neighbor_1_text->Draw("SAME");
	//neighbor_2_text->Draw("SAME");
	
	ADC_map_whole_chip->cd(2);
	
	hitBitMap_trigger_pixel->Draw("COLZ TEXT");
	
	ADC_map_whole_chip->SaveAs(Form("ADC_MAP_ALL_PIXELS_trigger_pixel_%d.png", trigger_pixel));
	
	TCanvas * ADC_trg_pixel = new TCanvas("canv8", "canv8", 800, 800);
	
	peak_ADC_trigger_pixel->SetMinimum(0);
	peak_ADC_trigger_pixel->Draw();
	
	
	return;

}

int getPixelIndex(int i, int j){
	
	if(i == 0 && j == 0){ return 0;}
	if(i == 0 && j == 1){ return 1;}
	if(i == 0 && j == 2){ return 2;}
	if(i == 0 && j == 3){ return 3;}
	if(i == 1 && j == 0){ return 4;}
	if(i == 1 && j == 1){ return 5;}
	if(i == 1 && j == 2){ return 6;}
	if(i == 1 && j == 3){ return 7;}
	if(i == 2 && j == 0){ return 8;}
	if(i == 2 && j == 1){ return 9;}
	if(i == 2 && j == 2){ return 10;}
	if(i == 2 && j == 3){ return 11;}
	if(i == 3 && j == 0){ return 12;}
	if(i == 3 && j == 1){ return 13;}
	if(i == 3 && j == 2){ return 14;}
	if(i == 3 && j == 3){ return 15;}
	
	else return -1;
	
	
}

int getPixelCanvasIndex(int i){
	
	if(i == 1 ){ return 1;}
	if(i == 2 ){ return 5;}
	if(i == 3 ){ return 9;}
	if(i == 4 ){ return 13;}
	if(i == 5 ){ return 2;}
	if(i == 6 ){ return 6;}
	if(i == 7 ){ return 10;}
	if(i == 8 ){ return 14;}
	if(i == 9 ){ return 3;}
	if(i == 10 ){ return 7;}
	if(i == 11 ){ return 11;}
	if(i == 12 ){ return 15;}
	if(i == 13 ){ return 4;}
	if(i == 14 ){ return 8;}
	if(i == 15 ){ return 12;}
	if(i == 16 ){ return 16;}
	
	else return -1;
	
	
}
