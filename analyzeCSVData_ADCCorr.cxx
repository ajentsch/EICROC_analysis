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

Color_t markerColor[16] = {kBlack, kRed+1, kRed, kGreen, kGreen+2, kBlue+1, kBlue+3, kMagenta, kBlue, kBlack, kBlue+1, kBlue+2, kGreen+2, kGreen+3, kMagenta+2, kMagenta+3};

void analyzeCSVData_ADCCorr(TString inputFileName = ""){

	ifstream inputCSVFile;
	
	double adcStepSize                = 1;
	const int numOfADCScanned         = 128;
	double startingValue              = 0; 
	
	TGraph * adc_pedestal_vs_adc_vref[16];
	TGraph * adc_pedestal_vs_adc_vref_INVERSE[16];
	TGraphErrors * TDC_sigma_vs_charge[16];
	
	for(int pixel = 0; pixel < 16; pixel++){
		
		adc_pedestal_vs_adc_vref[pixel]   = new TGraph();
		adc_pedestal_vs_adc_vref_INVERSE[pixel]   = new TGraph();
		TDC_sigma_vs_charge[pixel] = new TGraphErrors();
	
	}
	
	//double chargeBinold_values[startingValue];
	//double efficiency[startingValue];
	
	std::vector<double> charge_values;
	//std::vector<double> efficiency;
	
	TH1D * adc_distributions[4][4];
	TH1D * tdc_distributions[4][4];

	// number of thresholds, number of time bins, number of pixels
	TH1D * adc_RAW_distributions[numOfADCScanned][8][16];
	
	// number of thresholds, AVERAGED OVER TIME BINS, number of pixels
	TH1D * adc_mean_distributions[numOfADCScanned][16];
	
	// number of charges,  number of pixels
	TH1D * TDC_distributions[numOfADCScanned][16];
	
	// number of charges,  number of pixels
	

	for(int adcBin = 0; adcBin < numOfADCScanned; adcBin++){
		for(int tBin = 0; tBin < 8; tBin++){
			for(int pixel = 0; pixel < 16; pixel++){

				TString title;
				title.Form("adc_RAW_distribution_charge_%.0f_pixel_%d_timeBin_%d", startingValue + adcBin*adcStepSize, pixel, tBin);

				adc_RAW_distributions[adcBin][tBin][pixel] = new TH1D(title, "; ADC value [DACu]; counts", 256, 0, 255);
				adc_RAW_distributions[adcBin][tBin][pixel]->SetTitle(title);
				
				
	
			}
		}
	}	
	
	for(int adcBin = 0; adcBin < numOfADCScanned; adcBin++){
		for(int pixel = 0; pixel < 16; pixel++){

			TString title;
			title.Form("adc_MEAN_distribution_charge_%.0f_pixel_%d", startingValue + adcBin*adcStepSize, pixel);

			adc_mean_distributions[adcBin][pixel] = new TH1D(title, "; ADC bin [time bin = 25ns]; counts ", 8, 0, 8);
			adc_mean_distributions[adcBin][pixel]->SetTitle(title);
			
			title.Form("TDC_distribution_charge_%.0f_pixel_%d", startingValue + adcBin*adcStepSize, pixel);

			TDC_distributions[adcBin][pixel] = new TH1D(title, "; TDC bin [time bin = 25ns/1024]; counts ", 1024, 0, 1023);
			TDC_distributions[adcBin][pixel]->SetTitle(title);
			
	
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
	
	for(int fileIdx = 0; fileIdx < numOfADCScanned; fileIdx++){ //61
	
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_Vref_63DACuGlobal_100VbiasedSensor_zc706_5_15_26_0/output_charge_scan_Vref_cor/eic_data_charge_%.0f_1.csv", 0 + fileIdx*adcStepSize);
		inputFileName.Form("new_firmware_new_tests_good_delays_parameters/charge_scan_T330_unbiasedSensor_Vref_Cor_zc706_06_01_26_0/output_charge_scan_Vref_cor/eic_data_charge_%.0f_1.csv", 0 + fileIdx*adcStepSize);
		
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

						

				}
			}

		} //end while loop over file

		inputCSVFile.close();

		

		

	} //end of file loop
	
	
	
		
	int pad = 1;
	
	pad = 1;
	
	TCanvas * adcRAWCan = new TCanvas("canv5", "canv5", 1600, 1600);
	adcRAWCan->Divide(4,4);

	//for(int thresh = 0; thresh < numOfChargesScanned; chargeBin++){
	for(int adcBin = 0; adcBin < numOfADCScanned; adcBin++){
		//pad = 1;
		
		for(int pixel = 0; pixel < 16; pixel++){
			
			int idx = getPixelCanvasIndex(pad);
			adcRAWCan->cd(idx);
			for(int tBin = 0; tBin < 8; tBin++){
		
				adc_RAW_distributions[adcBin][tBin][pixel]->SetLineColor(markerColor[tBin]);
				
				if(pixel == 0){ 
					if(tBin == 0){adc_RAW_distributions[adcBin][tBin][pixel]->Draw();}
					else adc_RAW_distributions[adcBin][tBin][pixel]->Draw("SAME");
					
					if(tBin == 7){pad++;}
				}
				
		
			}
			
		}
	}
	
	
	
	//calculate ADC mean in each time bin here
	
	for(int adcBin = 0; adcBin < numOfADCScanned; adcBin++){
		for(int pixel = 0; pixel < 16; pixel++){
			for(int tBin = 0; tBin < 8; tBin++){
	
				double meanADC     = adc_RAW_distributions[adcBin][tBin][pixel]->GetMean();
				double meanADC_err = adc_RAW_distributions[adcBin][tBin][pixel]->GetRMS();
				
				adc_mean_distributions[adcBin][pixel]->SetBinContent(8 - tBin, meanADC);
				adc_mean_distributions[adcBin][pixel]->SetBinError(8 - tBin, meanADC_err);
			}
		}
	}	
	
	pad = 1;
	
	TCanvas * adcMeanCan = new TCanvas("canv6", "canv6", 1600, 1600);
	adcMeanCan->Divide(4,4);			
	
	//for(int chargeBin = 0; chargeBin < numOfChargesScanned; chargeBin++){
	for(int adcBin = 0; adcBin < numOfADCScanned; adcBin++){
		
		//pad = 1;
		
		for(int pixel = 0; pixel < 16; pixel++){
		
		
			int idx = getPixelCanvasIndex(pad);
			adcMeanCan->cd(idx);
			if(pixel == 0) {
				
				//adc_mean_distributions[chargeBin][pixel]->Add(adc_mean_distributions[0][pixel], -1);
				adc_mean_distributions[adcBin][pixel]->Draw(); 
				pad++;
			}
			//pad++;
		}
	}
	
	pad = 1;
	
	TCanvas * adcMaxCan = new TCanvas("canv7", "canv7", 1600, 1600);
	adcMaxCan->Divide(4,4);
	
	double TDC_pixel_zero[numOfADCScanned];

	int numPoints[numOfADCScanned][16];
	
	for(int adcBin = 0; adcBin < numOfADCScanned; adcBin++){
		for(int pixel = 0; pixel < 16; pixel++){
			numPoints[adcBin][pixel] = 0;
		}
	}
	
	int adc_max_val_at_max_vref = 10;
	int worst_pixel = 0;
	
	for(int adcBin = 0; adcBin < numOfADCScanned; adcBin++){
		
		double adc_vref_value = startingValue + adcBin*adcStepSize;
	
		
		for(int pixel = 0; pixel < 16; pixel++){
	
			int adc_max_bin      = adc_mean_distributions[adcBin][pixel]->GetMaximumBin();
			double adc_pedestal_value = adc_mean_distributions[adcBin][pixel]->GetBinContent(adc_max_bin); //get max ADC
				
			adc_pedestal_vs_adc_vref[pixel]->AddPoint(adc_vref_value, adc_pedestal_value);
			adc_pedestal_vs_adc_vref_INVERSE[pixel]->AddPoint(adc_pedestal_value, adc_vref_value);
			
			if(adcBin == (numOfADCScanned - 1) && (int)adc_pedestal_value >  adc_max_val_at_max_vref ){
				
				adc_max_val_at_max_vref = (int)adc_pedestal_value;
				worst_pixel = pixel;
				
				cout << "pixel " << pixel << "   ADC value at maximum offset = " << adc_pedestal_value << endl;
			}
			
		}
	}
	
	cout << "At maximum offset, the worst pixel is " << worst_pixel << endl;
	
	TLine * adc_worst_line = new TLine(0, adc_max_val_at_max_vref, 127, adc_max_val_at_max_vref);
	adc_worst_line->SetLineColor(kRed);
	adc_worst_line->SetLineWidth(2);
	
	for(int pixel = 0; pixel < 16; pixel++){
	
		adc_pedestal_vs_adc_vref[pixel]->SetTitle(Form("ADC_vs_charge_pixel_%d", pixel));
		
		adc_pedestal_vs_adc_vref[pixel]->GetYaxis()->SetTitle("ADC_{pedestal} [DACu]");
		adc_pedestal_vs_adc_vref[pixel]->GetXaxis()->SetTitle("ADC Vref [DACu]");
		
		adc_pedestal_vs_adc_vref[pixel]->GetYaxis()->SetRangeUser(0, 220);
		//adc_pedestal_vs_adc_vref[pixel]->GetXaxis()->SetTitle("ADC Vref [DACu]");
		
		double real_offset = adc_pedestal_vs_adc_vref_INVERSE[pixel]->Eval(adc_max_val_at_max_vref);
		
		int offset_per_channel = static_cast<int>(std::lround(adc_pedestal_vs_adc_vref_INVERSE[pixel]->Eval(adc_max_val_at_max_vref)));
	
		TLine * offset_line_per_channel = new TLine(offset_per_channel, 0, offset_per_channel, 220);
		offset_line_per_channel->SetLineColor(kGreen+2);
		offset_line_per_channel->SetLineWidth(2);
		
		cout << "Vref offset value for pixel " << pixel << " = " << offset_per_channel << "(real = "<< real_offset << ") ---> " << std::bitset<7>(offset_per_channel) << endl; 
		
		int idx = getPixelCanvasIndex(pad);
		adcMaxCan->cd(idx);
		adc_pedestal_vs_adc_vref[pixel]->SetMarkerStyle(20);
		adc_pedestal_vs_adc_vref[pixel]->SetMarkerSize(0.5);
		adc_pedestal_vs_adc_vref[pixel]->Draw("ALP");
		adc_worst_line->Draw("SAME");
		offset_line_per_channel->Draw("SAME");
		pad++;
		
		
		//int adc_pedestal = adc_MAX_vs_charge[pixel]->Eval(0);
		
		//cout << "pixel " << pixel << " ADC at minimum charge = " << adc_pedestal << "  ---> " << std::bitset<7>(adc_pedestal) << endl;
		
	}
	
	pad = 1;
	/*
	TCanvas * TDCCan = new TCanvas("canv8", "canv8", 1600, 1600);
	TDCCan->Divide(4,4);
	
	for(int pixel = 0; pixel < 16; pixel++){
	
		int idx = getPixelCanvasIndex(pad);
		TDCCan->cd(idx);
	
		TDC_distributions[4][pixel]->GetXaxis()->SetRange(400, 700); // 4 = 12, 15 = 45
		TDC_distributions[4][pixel]->Draw();
	
		pad++;
	}
	
	pad = 1;
	
	TCanvas * TDC_vs_charge_Can = new TCanvas("canv9", "canv9", 1600, 1600);
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
