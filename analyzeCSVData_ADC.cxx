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

#include <string>

using namespace std;

int getPixelIndex(int i, int j);
int getPixelCanvasIndex(int i);

Color_t markerColor[16] = {kGray, kGray+1, kGray+2, kGray+3, kGreen, kGreen+1, kGreen+2, kGreen+3, kBlue, kBlue+1, kBlue+2, kBlue+3, kRed, kRed+1, kRed+2, kRed+3};
//Color_t markerColor[16] = {kBlack, kGreen, kBlue, kRed, kBlack+1, kGreen+1, kBlue+1, kRed+1, kBlack+2, kGreen+2, kBlue+2, kRed+2, kBlack+3, kGreen+3, kBlue+3, kRed+3  };

int markerStyle[16] = {20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 29, 29, 29, 29};


//double ADC_offsets_coarse[16] = {120, 119, 127, 111, 116, 118, 114, 115, 119, 118, 116, 113, 114, 125, 118, 122};
//double ADC_offsets_coarse[16] = {92, 91, 99, 83, 88, 90, 86, 87, 91, 90, 88, 85, 86, 97, 90, 94};

double ADC_offsets_coarse[16] = {86, 89, 96, 84, 79, 86, 86, 89, 86, 91, 93, 90, 91, 102, 93, 98};


void analyzeCSVData_ADC(TString inputFileName = ""){

	ifstream inputCSVFile;
	
	double chargeStepSize             = 3;
	const int numOfChargesScanned     = 22;
	double startingValue              = 0; 
	
	TGraph * s_curve[16];
	TGraph * adc_MAX_vs_charge[16];
	TGraphErrors * TDC_sigma_vs_charge[16];
	
	for(int pixel = 0; pixel < 16; pixel++){
		
		s_curve[pixel]             = new TGraph();
		adc_MAX_vs_charge[pixel]   = new TGraph();
		TDC_sigma_vs_charge[pixel] = new TGraphErrors();
	
	}
	
	//double chargeBinold_values[startingValue];
	//double efficiency[startingValue];
	
	std::vector<double> charge_values;
	//std::vector<double> efficiency;
	
	TH1D * adc_distributions[4][4];
	TH1D * tdc_distributions[4][4];

	// number of thresholds, number of time bins, number of pixels
	TH1D * adc_RAW_distributions[numOfChargesScanned][8][16];
	
	// number of thresholds, AVERAGED OVER TIME BINS, number of pixels
	TH1D * adc_mean_distributions[numOfChargesScanned][16];
	TH1D * adc_mean_distributions_pedestal[16];
	
	// number of charges,  number of pixels
	TH1D * TDC_distributions[numOfChargesScanned][16];
	
	// number of charges,  number of pixels
	

	for(int chargeBin = 0; chargeBin < numOfChargesScanned; chargeBin++){
		for(int tBin = 0; tBin < 8; tBin++){
			for(int pixel = 0; pixel < 16; pixel++){

				TString title;
				title.Form("adc_RAW_distribution_charge_%.0f_pixel_%d_timeBin_%d", startingValue + chargeBin*chargeStepSize, pixel, tBin);

				adc_RAW_distributions[chargeBin][tBin][pixel] = new TH1D(title, "; ADC value [DACu]; counts", 256, 0, 255);
				adc_RAW_distributions[chargeBin][tBin][pixel]->SetTitle(title);
				
				
	
			}
		}
	}	
	
	for(int chargeBin = 0; chargeBin < numOfChargesScanned; chargeBin++){
		for(int pixel = 0; pixel < 16; pixel++){

			TString title;
			title.Form("adc_MEAN_distribution_charge_%.0f_pixel_%d", startingValue + chargeBin*chargeStepSize, pixel);

			adc_mean_distributions[chargeBin][pixel] = new TH1D(title, "; ADC bin [time bin = 25ns]; counts ", 8, 0, 8);
			adc_mean_distributions[chargeBin][pixel]->SetTitle(title);
			
			title.Form("TDC_distribution_charge_%.0f_pixel_%d", startingValue + chargeBin*chargeStepSize, pixel);

			TDC_distributions[chargeBin][pixel] = new TH1D(title, "; TDC bin [time bin = 25ns/1024]; counts ", 1024, 0, 1023);
			TDC_distributions[chargeBin][pixel]->SetTitle(title);
			
			if(chargeBin == 0){
			
				title.Form("adc_MEAN_distribution_charge_PEDESTAL_%.0f_pixel_%d", startingValue + chargeBin*chargeStepSize, pixel);
				adc_mean_distributions_pedestal[pixel] = new TH1D(title, "; ADC bin [time bin = 25ns]; counts ", 8, 0, 8);
				adc_mean_distributions_pedestal[pixel]->SetTitle(title);
			}
		}
	}	
			//adc_mean_distributions[61][16] = new TH1D(Form("adc_max_distribution_pixel_%d%d", i, j), "counts; ADC value [DACu]", 256, 0, 255);

	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
		
			adc_distributions[i][j] = new TH1D(Form("adc_max_distribution_pixel_%d%d", i, j), "counts; ADC value [DACu]", 256, 0, 255);
			tdc_distributions[i][j] = new TH1D(Form("tdc_distribution_pixel_%d%d", i, j), "counts; TDC value [DACu]", 1024, 0, 1023);
				
		}
	}
	
	double numFiles = 0;
	
	for(int fileIdx = 0; fileIdx < numOfChargesScanned; fileIdx++){ //61
	
		//inputFileName.Form("Maya_Feb_12_2026/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("Alex_Feb_20_2026_ZU706_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("Maya_Mar_4_2026_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("Ashik_Mar5_2026/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("Maya_Mar_9_2026_ZU706_Biased170_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("Alex_scan_newish_firmware_April_2026/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("output_charge_scan_april_22_2026/eic_data_charge_%.0f_dacp%.0f_th300_cDelay0_1.csv", 0 + fileIdx*chargeStepSize, 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("ChargeScan_zc706_4_24_26/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("chargeScan_zc706_4_28_26/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		
		//inputFileName.Form("output_charge_scan_old_firmware_4_29_2026/eic_data_charge_%.0f_dacp%.0f_th350_cDelay0_1.csv", 0 + fileIdx*chargeStepSize, 0 + fileIdx*chargeStepSize);
		
		//new firmware with updated delay parameters
		
		inputFileName.Form("new_firmware_new_tests_good_delays_parameters/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//GOOD GOOD GOOD GOOD GOOD
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/output_charge_scan_biased_shifted_pulse_may1_2026/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/output_charge_scan_biased_shifted_pulse_noscope_may1_2026/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/output_charge_scan_Unbiased_shifted_pulse_noscope_may1_2026/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/output_charge_scan_biased_may1_2026/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/output_charge_scan_pixel_1_1_april30_2026/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//ZCU106
		
		//inputFileName.Form("ZCU106_January_2025_firmware_good_timing_parameters/charge_scan_threshold_300_unbiased_same_timing_ZC706/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_100VbiasedSensor_zc706_5_15_26_0/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_100VbiasedSensor_zc706_5_15_26_1/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);	
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_100VbiasedSensor_zc706_5_18_26_0/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_100VbiasedSensor_zc706_5_18_26_1/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_170VbiasedSensor_zc706_5_29_26_0/output_charge_scan_biased/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T650_170VbiasedSensor_zc706_5_29_26_1/output_charge_scan_biased/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
	    
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_170VbiasedSensor_zc706_5_29_26_1/output_charge_scan_biased/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_170VbiasedSensor_zc706_5_29_26_DigProbeOff/output_charge_scan_biased/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_01_26_0/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_01_26_1/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_01_26_2/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_01_26_3/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_01_26_4/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		
		
	    //inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_02_26_0/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_02_26_1/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_02_26_3/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_02_26_4/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_02_26_5/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_02_26_6/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_zc706_06_02_26_7/output_charge_scan/eic_data_charge_%.0f_1.csv", 0 + fileIdx*chargeStepSize);
	
		cout << inputFileName << endl;
	
		inputCSVFile.open(inputFileName.Data());

		if(!inputCSVFile){ 
			cout << "File does not exist, skipping...." << endl; 
			//threshold_values[fileIdx] = 0;
			//efficiency[fileIdx] = 0;
			continue;
		}

		numFiles++;

		string eventNumber[4][4];
		double numGoodEvents[4][4];
		int ADC_values[8];
		int TDC_values[8];
		int hitBits[8];

		string eventNumber_str;
		string ADC_values_str[8];
		string TDC_values_str[8];
		string hitBits_str[8];

		if(!inputCSVFile){ cout << "ERROR - no input file found!!!! Exiting...." << endl; return; }

		int timeBin     = 7;
		int pixelRow    = 0;
		int pixelColumn = 0;

		string tmp;

		double numEvents = 0;
		
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				numGoodEvents[i][j] = 0;
			}
		}
		

		while(!inputCSVFile.eof()){
	
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
	
				eventNumber[pixelColumn][pixelRow] = eventNumber_str; //stoi(eventNumber_str);

				int numHitBitSet = 0;

				for(int tBin = 0; tBin < 8; tBin++){ //loop over 8 slices in one event's time bin

					ADC_values[tBin] = stoi(ADC_values_str[tBin]);
	    			TDC_values[tBin] = stoi(TDC_values_str[tBin]);
	    			hitBits[tBin]    = stoi(hitBits_str[tBin]);
				
					int pixelNum = getPixelIndex(pixelColumn, pixelRow);
				
					adc_distributions[pixelColumn][pixelRow]->Fill(ADC_values[tBin]);
					adc_RAW_distributions[fileIdx][tBin][pixelNum]->Fill(ADC_values[tBin]);
					/*
					if(startingValue + fileIdx*chargeStepSize == 350){
						if(tBin == 0){
							cout << eventNumber[pixelColumn][pixelRow] << " ; " << TDC_values_str[7-tBin] << ", " << ADC_values_str[7-tBin] << ", " << hitBits_str[7-tBin];
						}
						if(tBin > 0 && tBin != 7){
							cout <<" ; " << TDC_values_str[7-tBin] << ", " << ADC_values_str[7-tBin] << ", " << hitBits_str[7-tBin];
						}
						if(tBin == 7){
							cout <<" ; " << TDC_values_str[7-tBin] << ", " << ADC_values_str[7-tBin] << ", " << hitBits_str[7-tBin] << endl;
						}
					}
					*/
					if(hitBits[tBin] == 1){ numHitBitSet++; }
					
				}

				//cout << "pixel (" << pixelColumn << ", " << pixelRow << ")" << endl;
	
				bool foundTDC = false;
	
				if(numHitBitSet == 1){
				//if(numHitBitSet > 0 && numHitBitSet < 5){
				//if(numHitBitSet > 0 ){
					
					for(int tBin = 0; tBin < 8; tBin++){
						if(TDC_values[tBin] > 1 && TDC_values[tBin] < 1023 && !foundTDC){ 
							
							foundTDC = true;
							
							
							//tdc_distributions[pixelColumn][pixelRow]->Fill(TDC_values[tBin]); 
							
							int pixelNum = getPixelIndex(pixelColumn, pixelRow);
							
							if(fileIdx == 10 && pixelNum == 5){ cout << "Event " << numEvents << " :  Q_inj = 30 DACu --> TDC = " << TDC_values[tBin] << endl; }
							
							TDC_distributions[fileIdx][pixelNum]->Fill(TDC_values[tBin]);
							
						
						} 
					}
					
					//adc_distributions[pixelColumn][pixelRow]->Fill(ADC_values[tBin]);
					//tdc_distributions[pixelColumn][pixelRow]->Fill(TDC_values[tBin]);
					numGoodEvents[pixelColumn][pixelRow]++;
			
				}
			
	
				pixelRow++;
				if(pixelRow == 4 && pixelColumn < 3) { pixelColumn++; pixelRow = 0;}
				if(pixelRow == 4 && pixelColumn == 3){ pixelColumn++;}
				if(pixelRow == 4 && pixelColumn == 4){ 
					pixelRow = 0; pixelColumn = 0; 

					//cout << eventNumber[0][0] << " ";		

					numEvents++;

					//for(int tBin = 0; tBin < 8; tBin++){
			
						//cout << TDC_values[7-tBin] << " " << ADC_values[7-tBin] << " " << hitBits[7-tBin] << " ";
						//cout << TDC_values_str[7-tBin] << " " << ADC_values_str[7-tBin] << " " << hitBits_str[7-tBin] << " ";
		
					//}
					//cout << "\n";		

				}
			}

		} //end while loop over file

		inputCSVFile.close();

		

		/*
		
		*/
		
		
		//cout << "total number of events = " << numEvents << endl;
		/*
		cout << "\n --- number of good events in pixels --- " << endl;
		cout << " THRESHOLD = " << startingValue + fileIdx*chargeStepSize << endl;
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
			
			
				int pixel = getPixelIndex(i, j);
				
				if(pixel != -1){	
					s_curve[pixel]->AddPoint(startingValue + fileIdx*chargeStepSize, numGoodEvents[i][j]/numEvents);
				}
			
				cout << "pixel( " << i << ", " << j << ") = " << numGoodEvents[i][j]/numEvents << endl;
				
			}
		}
		*/

	} //end of file loop
	
	//s_curve = new TGraph(threshold_values.size(), &threshold_values[0], &efficiency[0]);
	/*
	TCanvas * sCurveCanvas = new TCanvas("canv1", "canv1", 1600, 1600);
	sCurveCanvas->Divide(4,4);
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
	
			int pixel = getPixelIndex(i, j);
			int canvasBin = getPixelIndex(j, i);
	
			if(pixel != -1){
				
				sCurveCanvas->cd(canvasBin+1);
				
				s_curve[pixel]->SetTitle(Form("S-Curve, pixel (%d,%d), Q_{injected} = 63 DACu", i, j));
				s_curve[pixel]->GetXaxis()->SetTitle("Threshold [DACu]");
				s_curve[pixel]->GetYaxis()->SetTitle("Efficiency");
				s_curve[pixel]->SetMarkerStyle(20);
				s_curve[pixel]->SetMarkerColor(markerColor[pixel]);
				s_curve[pixel]->SetLineColor(markerColor[pixel]);
				s_curve[pixel]->Draw("ALP");
				
			}
		}
	}
			
	TCanvas * sCurvesOneCanvas = new TCanvas("canv2", "canv2", 800, 800);
	sCurvesOneCanvas->cd();
	*/
	/*
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
	
			int pixel = getPixelIndex(i, j);
			int canvasBin = getPixelIndex(j, i);
	
			if(pixel != -1){
				
				if(pixel == 0){ s_curve[pixel]->Draw("ALP");}
				else s_curve[pixel]->Draw("SAME LP");
				
			}
		}
	}
	*/	
		
	int pad = 1;
	/*	
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
	*/
	pad = 1;
	/*
	TCanvas * adcCan = new TCanvas("canv4", "canv4", 1600, 1600);
	adcCan->Divide(4,4);

	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
		
			adcCan->cd(pad);
			adc_distributions[i][j]->Draw();
			pad++;
		}
	}
	*/
	pad = 1;
	
	TCanvas * adcRAWCan = new TCanvas("canv5", "canv5", 1600, 1600);
	adcRAWCan->Divide(4,4);

	//for(int thresh = 0; thresh < numOfChargesScanned; chargeBin++){
	for(int chargeBin = 0; chargeBin < numOfChargesScanned; chargeBin++){
		//pad = 1;
		
		for(int pixel = 0; pixel < 16; pixel++){
			
			int idx = getPixelCanvasIndex(pad);
			adcRAWCan->cd(idx);
			for(int tBin = 0; tBin < 8; tBin++){
		
				adc_RAW_distributions[chargeBin][tBin][pixel]->SetLineColor(markerColor[tBin]);
				
				if(pixel == 0){ 
					if(tBin == 0){adc_RAW_distributions[chargeBin][tBin][pixel]->Draw();}
					else adc_RAW_distributions[chargeBin][tBin][pixel]->Draw("SAME");
					
					if(tBin == 7){pad++;}
				}
				
		
			}
			
		}
	}
	
	
	
	//calculate ADC mean in each time bin here
	
	for(int chargeBin = 0; chargeBin < numOfChargesScanned; chargeBin++){
		for(int pixel = 0; pixel < 16; pixel++){
			for(int tBin = 0; tBin < 8; tBin++){
	
				double meanADC     = adc_RAW_distributions[chargeBin][tBin][pixel]->GetMean();
				double meanADC_err = adc_RAW_distributions[chargeBin][tBin][pixel]->GetMeanError();
				
				adc_mean_distributions[chargeBin][pixel]->SetBinContent(8 - tBin, meanADC);
				adc_mean_distributions[chargeBin][pixel]->SetBinError(8 - tBin, meanADC_err);
				
				if(chargeBin == 0){
					adc_mean_distributions_pedestal[pixel]->SetBinContent(8 - tBin, meanADC);
					adc_mean_distributions_pedestal[pixel]->SetBinError(8 - tBin, meanADC_err);
				}
				
				//adc_mean_distributions[chargeBin][pixel]->Add(adc_mean_distributions_pedestal[pixel], -1);
			}
		}
	}	
	
	
	
	pad = 1;
	
	TCanvas * adcMeanCan = new TCanvas("canv6", "canv6", 1600, 1600);
	adcMeanCan->Divide(4,4);			
	
	//for(int chargeBin = 0; chargeBin < numOfChargesScanned; chargeBin++){
	for(int chargeBin = 0; chargeBin < numOfChargesScanned; chargeBin++){
		
		//pad = 1;
		
		for(int pixel = 0; pixel < 16; pixel++){
		
			//if(chargeBin != 0){ adc_mean_distributions[chargeBin][pixel]->Add(adc_mean_distributions[0][pixel], -1); }
		
			int idx = getPixelCanvasIndex(pad);
			adcMeanCan->cd(idx);
			if(pixel == 0) {
				
				//adc_mean_distributions[chargeBin][pixel]->Add(adc_mean_distributions[0][pixel], -1);
				adc_mean_distributions[chargeBin][pixel]->SetMaximum(160);
				adc_mean_distributions[chargeBin][pixel]->SetMinimum(40);
				adc_mean_distributions[chargeBin][pixel]->Draw(); 
				pad++;
			}
			//pad++;
		}
	}
	
	pad = 1;
	
	TCanvas * adcMaxCan = new TCanvas("canv7", "canv7", 1600, 1600);
	adcMaxCan->Divide(4,4);
	
	double TDC_pixel_zero[numOfChargesScanned];

	int numPoints[numOfChargesScanned][16];
	
	for(int chargeBin = 0; chargeBin < numOfChargesScanned; chargeBin++){
		for(int pixel = 0; pixel < 16; pixel++){
			numPoints[chargeBin][pixel] = 0;
		}
	}
	
	for(int chargeBin = 0; chargeBin < numOfChargesScanned; chargeBin++){
		
		double charge_value = startingValue + chargeBin*chargeStepSize;
	
		
		for(int pixel = 0; pixel < 16; pixel++){
	
			int adc_max_bin      = adc_mean_distributions[chargeBin][pixel]->GetMaximumBin();
			double adc_max_value = adc_mean_distributions[chargeBin][pixel]->GetBinContent(adc_max_bin);
				
			adc_MAX_vs_charge[pixel]->AddPoint(charge_value, adc_max_value);
			
			//double scale_factor = (25.0/1024.0) * 1000.0;
			
			//TDC_distributions[chargeBin][pixel]->Scale(scale_factor);
			
			double TDC_sigma_raw      = TDC_distributions[chargeBin][pixel]->GetRMS(); //using RMS for now to extract the sigma value
			double TDC_sigma_raw_err  = TDC_distributions[chargeBin][pixel]->GetRMSError(); //using RMS for now to extract the sigma value
			
			if(pixel == 5){
				cout << "chargeBin = " << chargeBin << " TDC_sigma_raw = " << TDC_sigma_raw << endl;
			}
			
			//if(pixel > 0)  {TDC_sigma_raw = TMath::Abs(TDC_sigma_raw - TDC_distributions[chargeBin][pixel-1]->GetRMS() );}
			//if(pixel == 0) {TDC_sigma_raw = TMath::Abs(TDC_sigma_raw - TDC_distributions[chargeBin][pixel+1]->GetRMS() );}
			
			double TDC_sigma_time       = (25.0/1024.0) * TDC_sigma_raw * 1000.0; //1000 is to convert from ns to ps
			double TDC_sigma_time_error = (25.0/1024.0) * TDC_sigma_raw_err * 1000.0; //1000 is to convert from ns to ps
			
			//double TDC_sigma_time       = TDC_sigma_raw; //1000 is to convert from ns to ps
			//double TDC_sigma_time_error = TDC_sigma_raw_err; //1000 is to convert from ns to ps
			
			if(TDC_sigma_time_error/TDC_sigma_time < 0.04 && TDC_sigma_time != 0.0){
			//if(TDC_sigma_time != 0.0){
				TDC_sigma_vs_charge[pixel]->AddPoint(charge_value, TDC_sigma_time); 
				TDC_sigma_vs_charge[pixel]->SetPointError(numPoints[chargeBin][pixel], 0.0, TDC_sigma_time_error);
				numPoints[chargeBin][pixel]++;
			}
			
		}
	}
	
	
	
	///////////////////////////////////////////////////////////////////////////////////////
	//calculate ADC correction offsets here
	////////////
	//////////////////////////////////////////////////////////////////////////
	
	
	pad = 1;
	
	TLine * ADC_45_line[16];
	
	for(int pixel = 0; pixel < 16; pixel++){
	
		adc_MAX_vs_charge[pixel]->SetTitle(Form("ADC_vs_charge_pixel_%d", pixel));
		
		adc_MAX_vs_charge[pixel]->GetYaxis()->SetTitle("ADC_{max} [DACu] (no pedestal subtraction)");
		adc_MAX_vs_charge[pixel]->GetXaxis()->SetTitle("Injected Charge [DACu]");
	
		ADC_45_line[pixel] = new TLine(41, adc_MAX_vs_charge[pixel]->Eval(45), 49, adc_MAX_vs_charge[pixel]->Eval(45) );
		ADC_45_line[pixel]->SetLineWidth(2);
		ADC_45_line[pixel]->SetLineColor(markerColor[pixel]);
	
		int idx = getPixelCanvasIndex(pad);
		adcMaxCan->cd(idx);
		adc_MAX_vs_charge[pixel]->SetMarkerStyle(20);
		adc_MAX_vs_charge[pixel]->SetMarkerSize(0.5);
		adc_MAX_vs_charge[pixel]->Draw("ALP");
		//ADC_45_line[pixel]->Draw("SAME");
		pad++;
		
	}
	
	
	cout << "\n\n";
	cout << "------------------------------ calculate ADC offsets using new approach ---------------------------" << "\n\n";
	
	TCanvas * can8 = new TCanvas("canv8", "canv8", 1600, 600);
	can8->Divide(2,1);
	
	can8->cd(1);
	
	double min_ADC_all_channels = 256;
	double max_ADC_all_channels = 0;
	double ADC_at_q_equal_zero = 256;
	double ADC_at_q_equal_45   = 256;
	
	double average_ADC = 0.0;
	double tmp_average_ADC = 0.0;
	
	TLine * pivot_charge_line = new TLine(45, 0, 45, 256);
	pivot_charge_line->SetLineWidth(2);
	pivot_charge_line->SetLineColor(kRed);
	
	TH1D * adc_spread_at_charge = new TH1D("ADC_spread", "ADC_spread", 256, 0, 256);
	
	for(int pixel = 0; pixel < 16; pixel++){
	
		tmp_average_ADC = 0.0;
		
		ADC_at_q_equal_zero = adc_MAX_vs_charge[pixel]->Eval(45);
		
		adc_spread_at_charge->Fill(adc_MAX_vs_charge[pixel]->Eval(45));
		
		for(int q = 0; q < 5; q++){
			tmp_average_ADC += adc_MAX_vs_charge[pixel]->Eval(q + 2*q);
		}
		
		average_ADC += ADC_at_q_equal_zero; //tmp_average_ADC/5.0;
		
		if(ADC_at_q_equal_zero < min_ADC_all_channels){min_ADC_all_channels = ADC_at_q_equal_zero; }
		if(ADC_at_q_equal_zero > max_ADC_all_channels){max_ADC_all_channels = ADC_at_q_equal_zero; }
		
	}
	
	//find median value
	
	double min_max_average = 0.5*(max_ADC_all_channels - min_ADC_all_channels);
	double delta_min_max_average = min_max_average;
	
	int mid_most_pixel = 0;
	
	for(int pixel = 0; pixel < 16; pixel++){
	
		
		ADC_at_q_equal_zero = adc_MAX_vs_charge[pixel]->Eval(45);
		
		double tmp = TMath::Abs(ADC_at_q_equal_zero - min_max_average);
		
		if(tmp < delta_min_max_average){ delta_min_max_average = tmp; mid_most_pixel = pixel;}
		
	}
	
	for(int pixel = 0; pixel < 16; pixel++){
	
		adc_MAX_vs_charge[pixel]->SetMinimum(0);
		adc_MAX_vs_charge[pixel]->SetMaximum(256);
		//adc_MAX_vs_charge[pixel]->SetMaximum(81);
		
		adc_MAX_vs_charge[pixel]->GetXaxis()->SetRangeUser(0, 63);
	
		adc_MAX_vs_charge[pixel]->SetMarkerColor(markerColor[pixel]);
		adc_MAX_vs_charge[pixel]->SetMarkerStyle(markerStyle[pixel]);
		adc_MAX_vs_charge[pixel]->SetMarkerSize(0.5);
		
		//if(pixel > 3) { continue; }
		if(pixel == 0){adc_MAX_vs_charge[pixel]->Draw("ALP"); pivot_charge_line->Draw("SAME");}// ADC_45_line[pixel]->Draw("SAME");}
		else adc_MAX_vs_charge[pixel]->Draw("SAME LP"); //ADC_45_line[pixel]->Draw("SAME"); 		
		
		
	}
	can8->cd(2);
	adc_spread_at_charge->Draw();
	
	average_ADC = average_ADC/16.0;
	
	cout << "Minimum ADC = " << min_ADC_all_channels << endl;
	cout << "Maximum ADC = " << max_ADC_all_channels << endl;
	cout << "Average ADC = " << average_ADC << endl;
	cout << "Median ADC  = " << adc_MAX_vs_charge[mid_most_pixel]->Eval(45) << endl;
	
	TLine * average_line = new TLine(0, average_ADC, 63, average_ADC);
	average_line->SetLineWidth(2);
	average_line->SetLineColor(kRed);
	
	//average_line->Draw("SAME");
		
	bool usingCoarseCorrectionRun = false;
	if( max_ADC_all_channels > 150 && min_ADC_all_channels > 100){ usingCoarseCorrectionRun = true; }
	
	double correction_ADC = 0.0; 
	
	if(usingCoarseCorrectionRun){ correction_ADC = max_ADC_all_channels - 127;}
	
	if(!usingCoarseCorrectionRun){ //in this case, this is the second pass where we are now doing the find adjustment
		
		//correction_ADC = max_ADC_all_channels;
		correction_ADC = average_ADC;// adc_MAX_vs_charge[mid_most_pixel]->Eval(45);
		
	}
	
	
	
	//First, let's print the array for the coarse offsets
	
	for(int pixel = 0; pixel < 16; pixel++){
	
		double current_ADC = adc_MAX_vs_charge[pixel]->Eval(45);
		double offset = current_ADC - correction_ADC;
		
		int offset_per_channel = static_cast<int>(std::lround(offset));
		
		
		if(pixel == 0 ){ cout << "{" << offset_per_channel << ", "; }
		else if(pixel == 15){ cout <<  offset_per_channel << "};" << endl; }
		else { cout << offset_per_channel << ", "; }
	}
	
	for(int i = 0; i < 16; i++){
	
		int pixel = getPixelCanvasIndex(i+1) - 1;
	
		double current_ADC = adc_MAX_vs_charge[pixel]->Eval(45);
		double offset = current_ADC - correction_ADC;
		
		if(!usingCoarseCorrectionRun){
			
			offset = offset + ADC_offsets_coarse[pixel];
				
		}
		
		int offset_per_channel = static_cast<int>(std::lround(offset));
		
		cout << "ADC V_ref offset, pixel " << pixel << "\t--->\t" << offset_per_channel << "  (real = "<< offset << ")\t--->\t" << std::bitset<7>(offset_per_channel) << endl; 
		
		
	
	}
	
	for(int pixel = 0; pixel < 16; pixel++){
	
		double current_ADC = adc_MAX_vs_charge[pixel]->Eval(45);
		double offset = current_ADC - correction_ADC;
		
		if(!usingCoarseCorrectionRun){
			
			offset = offset + ADC_offsets_coarse[pixel];
				
		}
		
		int offset_per_channel = static_cast<int>(std::lround(offset));
		
		if(pixel == 0 ){ cout << "{" << offset_per_channel << ", "; }
		else if(pixel == 15){ cout <<  offset_per_channel << "};" << endl; }
		else { cout << offset_per_channel << ", "; }
		
	
	}
	
	
	
	/*
	pad = 1;
	
	TCanvas * TDCCan = new TCanvas("canv9", "canv9", 1600, 1600);
	TDCCan->Divide(4,4);
	
	for(int pixel = 0; pixel < 16; pixel++){
	
		int idx = getPixelCanvasIndex(pad);
		TDCCan->cd(idx);
	
		TDC_distributions[4][pixel]->GetXaxis()->SetRange(400, 700); // 4 = 12, 15 = 45
		TDC_distributions[4][pixel]->Draw();
	
		pad++;
	}
	
	pad = 1;
	
	TCanvas * TDC_vs_charge_Can = new TCanvas("canv10", "canv10", 1600, 1600);
	TDC_vs_charge_Can->Divide(4,4);
	
	for(int pixel = 0; pixel < 16; pixel++){
	
		TDC_sigma_vs_charge[pixel]->RemovePoint(0);
	
		TDC_sigma_vs_charge[pixel]->SetTitle(Form("TDC_sigma_vs_charge_pixel_%d", pixel));
		
		TDC_sigma_vs_charge[pixel]->GetYaxis()->SetTitle("#sigma_{TDC} [ps]");
		TDC_sigma_vs_charge[pixel]->GetXaxis()->SetTitle("Injected Charge [DACu]");
	
		int idx = getPixelCanvasIndex(pad);
		TDC_vs_charge_Can->cd(idx);
		TDC_sigma_vs_charge[pixel]->SetMarkerStyle(20);
		TDC_sigma_vs_charge[pixel]->SetMarkerSize(0.5);
		TDC_sigma_vs_charge[pixel]->Draw("ALP");
		pad++;
	}
	*/
	
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
