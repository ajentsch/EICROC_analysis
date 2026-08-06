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

Color_t markerColor[16] = {kGray, kGray+1, kGray+2, kGray+3, kGreen, kGreen+1, kGreen+2, kGreen+3, kBlue, kBlue+1, kBlue+2, kBlue+3, kRed, kRed+1, kRed+2, kRed+3};
//Color_t markerColor[16] = {kBlack, kGreen, kBlue, kRed, kBlack+1, kGreen+1, kBlue+1, kRed+1, kBlack+2, kGreen+2, kBlue+2, kRed+2, kBlack+3, kGreen+3, kBlue+3, kRed+3  };

int markerStyle[16] = {20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 22, 22, 29, 29, 29, 29};

void analyzeCSVData_thresholdCorr(TString inputFileName = ""){

	ifstream inputCSVFile;
	
	double thresholdStepSize   = 1;
	const int numOfThresholdsScanned = 127;
	
	int minThreshold = 0;
	
	TGraph * s_curve[16];
	
	for(int pixel = 0; pixel < 16; pixel++){s_curve[pixel] = new TGraph();}
	
	//double threshold_values[100];
	//double efficiency[100];
	
	std::vector<double> threshold_values;
	std::vector<double> efficiency;
	
	TH1D * adc_distributions[4][4];
	TH1D * tdc_distributions[4][4];

	// number of thresholds, number of time bins, number of pixels
	TH1D * adc_RAW_distributions[numOfThresholdsScanned][8][16];
	
	// number of thresholds, AVERAGED OVER TIME BINS, number of pixels
	TH1D * adc_mean_distributions[numOfThresholdsScanned][16];

	for(int thresh = 0; thresh < numOfThresholdsScanned; thresh++){
		for(int tBin = 0; tBin < 8; tBin++){
			for(int pixel = 0; pixel < 16; pixel++){

				TString title;
				title.Form("adc_max_distribution_thresh_%.0f_pixel_%d_timeBin_%d", minThreshold + thresh*thresholdStepSize, pixel, tBin);

				adc_RAW_distributions[thresh][tBin][pixel] = new TH1D(title, "; ADC value [DACu]; counts", 256, 0, 255);
				adc_RAW_distributions[thresh][tBin][pixel]->SetTitle(title);
	
			}
		}
	}	
	
	for(int thresh = 0; thresh < numOfThresholdsScanned; thresh++){
		for(int pixel = 0; pixel < 16; pixel++){

			TString title;
			title.Form("adc_MEAN_distribution_thresh_%.0f_pixel_%d", minThreshold + thresh*thresholdStepSize, pixel);

			adc_mean_distributions[thresh][pixel] = new TH1D(title, "; ADC bin [time bin = 25ns]; counts ", 8, 0, 8);
			adc_mean_distributions[thresh][pixel]->SetTitle(title);
	
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
	
	int charge_DACu = 45;
	
	for(int fileIdx = 0; fileIdx < numOfThresholdsScanned; fileIdx++){ //61
	
		
	
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/correction_thresholdScan_Q45_5_7_26_zc706_unbaisedSensor/output_threshold_scan_Q45/eic_data_ThreshHold_Q45_%.0f_1.csv", minThreshold + fileIdx*thresholdStepSize);
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/correction_threshold_Q45_unbiasedSensor_zc706_5_8_26_0/output_threshold_scan_Q45/eic_data_ThreshHold_Q45_%.0f_1.csv", minThreshold + fileIdx*thresholdStepSize);
		
		//special combination from two tests [250, 417]
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/correction_threshold_Q45_unbiasedSensor_zc706_COMBINED_DATA_v0/output_threshold_scan_Q45/eic_data_ThreshHold_Q45_%.0f_1.csv", minThreshold + fileIdx*thresholdStepSize);
		
		//special combination from two tests [250, 417]
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/correction_threshold_Q45_unbiasedSensor_zc706_COMBINED_DATA_v1/output_threshold_scan_Q45/eic_data_ThreshHold_Q45_%.0f_1.csv", minThreshold + fileIdx*thresholdStepSize);
	
		//single scan broad range [250, 417]
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/correction_threshold_Q45_unbiasedSensor_zc706_5_8_26_1/output_threshold_scan_Q45/eic_data_ThreshHold_Q45_%.0f_1.csv", minThreshold + fileIdx*thresholdStepSize);
		
		//single scan 280 to 477, 197 total
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/correction_threshold_Q45_biasedSensor_zc706_5_11_26_0/output_threshold_scan_Q45/eic_data_ThreshHold_Q45_%.0f_1.csv", minThreshold + fileIdx*thresholdStepSize);
	
		//single scan 316 to 556. with 64 DACu offset
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/correction_threshold_Q45_64DACuOffset_100VbiasedSensor_zc706_5_12_26_0/output_threshold_scan_Q45/eic_data_ThreshHold_Q45_%.0f_1.csv", minThreshold + fileIdx*thresholdStepSize);
	
		//VTh_cor scan using edited main code with direct register editing
		
		//inputFileName.Form("new_firmware_new_tests_good_delays_parameters/Vth_cor_threshold_Q45_100VbiasedSensor_zc706_5_14_26_0/output_threshold_scan_Q45_Vth_cor/eic_data_ThreshHold_Vth_corc_Q45_%.0f_1.csv", minThreshold + fileIdx*thresholdStepSize);
	
		inputFileName.Form("new_firmware_new_tests_good_delays_parameters/Vth_cor_threshold_Q45_100VbiasedSensor_zc706_5_15_26_0/output_threshold_scan_Q45_Vth_cor/eic_data_ThreshHold_Vth_corc_Q45_%.0f_1.csv", minThreshold + fileIdx*thresholdStepSize);
	
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
				
					adc_distributions[pixelColumn][pixelRow]->Fill(ADC_values[tBin]);
				
					if(minThreshold + fileIdx*thresholdStepSize == 350){
						if(tBin == 0){
							cout << eventNumber[pixelColumn][pixelRow] << " ; " << TDC_values_str[7-tBin] << ", " << ADC_values_str[7-tBin] << ", " << hitBits_str[7-tBin];
						}
						if(tBin > 0 && tBin != 7){
							cout <<" ; " << TDC_values_str[7-tBin] << ", " << ADC_values_str[7-tBin] << ", " << hitBits_str[7-tBin];
						}
						if(tBin == 7){
							cout <<" ; " << TDC_values_str[7-tBin] << ", " << ADC_values_str[7-tBin] << ", " << hitBits_str[7-tBin] << endl;
						}
						
						int pixelNum = getPixelIndex(pixelColumn, pixelRow);
						
						adc_RAW_distributions[fileIdx][tBin][pixelNum]->Fill(ADC_values[tBin]);
						
					}
					
					if(hitBits[tBin] == 1){ numHitBitSet++; }
					
				}

				//cout << "pixel (" << pixelColumn << ", " << pixelRow << ")" << endl;
	
				if(numHitBitSet == 1){
				//if(numHitBitSet > 0 && numHitBitSet < 5){
				//if(numHitBitSet > 0 ){
					
					for(int tBin = 0; tBin < 8; tBin++){
						if(TDC_values[tBin] > 0){ tdc_distributions[pixelColumn][pixelRow]->Fill(TDC_values[tBin]); } 
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

		cout << "\n --- number of good events in pixels --- " << endl;
		cout << " THRESHOLD = " << minThreshold + fileIdx*thresholdStepSize << endl;
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
			
			
				int pixel = getPixelIndex(i, j);
				
				if(pixel != -1){	
					s_curve[pixel]->AddPoint(fileIdx*thresholdStepSize, numGoodEvents[i][j]/numEvents);
				}
			
				cout << "pixel( " << i << ", " << j << ") = " << numGoodEvents[i][j]/numEvents << endl;
				
			}
		}

	} //end of file loop
	
	//s_curve = new TGraph(threshold_values.size(), &threshold_values[0], &efficiency[0]);
	
	TCanvas * sCurveCanvas = new TCanvas("canv1", "canv1", 1600, 1600);
	sCurveCanvas->Divide(4,4);
	
	double thresholdMin = 0.0;
	double thresholdMax = 127.0;
	
	int thresh_50_percent_values[16];
	
	TGraph * fitted_50_percent_values = new TGraph();
	
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
					
					if(eff > 0.3 && eff < 0.7 && thresh > 40 ){ break; }
				}
					
				
				//TF1 * sCurveFit = new TF1("scurvefit", "-1*[0]*TMath::Erf((x-[1])/[2])", thresh - 100, thresh + 150);
				//TF1 * sCurveFit = new TF1("scurvefit", "[0] / (1 + exp(−1*[1]*(x−[2])))", thresh - 100, thresh + 150);
				
				//if( pixel == 10){ continue; }
				
				//sCurveFit->SetParLimits(0, 0.99, 1.01);
				//sCurveFit->SetParLimits(1, thresh - 25, thresh + 25);
				//sCurveFit->SetParLimits(2, 10, 40);
				
				TF1 *fLogistic = new TF1("fLogistic", "[0] / (1.0 + TMath::Exp(-[1] * (x - [2]))) + [3]", thresh - 100, thresh + 100);

				fLogistic->SetParNames("Amplitude", "Steepness", "Midpoint", "Offset");
				fLogistic->SetParameters(1.0, 0.2, thresh, 0.0);
				
				s_curve[pixel]->Fit(fLogistic, "RQ");
				
				TF1 *fInverse = new TF1("fInverse", "(-1/[0]) * TMath::Log(([1] - x)/x) + [2]", 0.0, 1.0);
				fInverse->SetParameter(0, fLogistic->GetParameter(1));
				fInverse->SetParameter(1, fLogistic->GetParameter(0));
				fInverse->SetParameter(2, fLogistic->GetParameter(2));
				
				cout << "pixel " << pixel << "  50 percent threshold = " << fInverse->Eval(0.5) << endl;
				
				thresh_50_percent_values[pixel] = (int)fInverse->Eval(0.5);
				
				
				
				if(fInverse->Eval(0.5) < thresholdMin || pixel == 0) { thresholdMin = fInverse->Eval(0.5); }
				if(fInverse->Eval(0.5) > thresholdMax || pixel == 0) { thresholdMax = fInverse->Eval(0.5); }
				
				fitted_50_percent_values->AddPoint(fInverse->Eval(0.5), 0.5);
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
	
	cout << "MEAN threshold = " << fitted_50_percent_values->GetMean() << endl;
	cout << "\n";
	
	//TLine * mean_50_percent = new TLine(fitted_50_percent_values->GetMean(), 0.0, fitted_50_percent_values->GetMean(), 1.0);
	//TLine * mean_50_percent = new TLine(thresholdMax, 0.0, thresholdMax, 1.0);
	TLine * mean_50_percent = new TLine(thresholdMin, 0.0, thresholdMin, 1.0);
	mean_50_percent->SetLineWidth(3);
	mean_50_percent->SetLineColor(kBlack);
	
	TLine * threhold_50_percent_line = new TLine(thresholdMin, 0.5, thresholdMax, 0.5);
	threhold_50_percent_line->SetLineWidth(3);
	threhold_50_percent_line->SetLineColor(kBlack);
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
	
			int pixel = getPixelIndex(i, j);
			int canvasBin = getPixelIndex(j, i);
	
			if(pixel != -1){
				
				sCurveCanvas->cd(canvasBin+1);
				mean_50_percent->Draw("SAME");
				threhold_50_percent_line->Draw("SAME");
			}
		}
	}
	
	for(int pixel = 0; pixel < 16; pixel++){
		
		
		//int offset_value = ( (int)thresholdMax - thresh_50_percent_values[pixel] );
		int offset_value = ( thresh_50_percent_values[pixel] - (int)thresholdMin);
		
		//if offset_value is positive, then need to REDUCE VTH_corr to ADD voltage
		//if(offset_value > 0){ offset_value = 64 - offset_value*0.5; }
		//if(offset_value < 0){ offset_value = 64 + offset_value*0.5; }
		
		cout << "Offset for pixel " << pixel << " = " << offset_value;
		cout << " --> " << std::bitset<7>(offset_value) << endl;
		
	}
			
	TCanvas * sCurvesOneCanvas = new TCanvas("canv2", "canv2", 1600, 600);
	sCurvesOneCanvas->Divide(2,1);
	
	sCurvesOneCanvas->cd(1);
	
	
	
	for(int i = 0; i < 4; i++){
		for(int j = 0; j < 4; j++){
	
			int pixel = getPixelIndex(i, j);
			int canvasBin = getPixelIndex(j, i);
	
			if(pixel != -1){
				
				//if(pixel == 1 || pixel == 5 || pixel == 7 || pixel == 10 || pixel == 12){continue;}
				
				if(pixel == 0){ 
					s_curve[pixel]->GetYaxis()->SetRangeUser(0.0, 1.2);
					//s_curve[pixel]->GetXaxis()->SetRangeUser(380, 500);
					s_curve[pixel]->Draw("ALP");
					mean_50_percent->Draw("SAME");
				}
				else{ 
					s_curve[pixel]->Draw("SAME LP");
					//mean_50_percent->Draw("SAME");
				}
				
			}
		}
	}
		
	sCurvesOneCanvas->cd(2);	
	fitted_50_percent_values->Draw("ALP");
		
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
	
	pad = 1;
	
	TCanvas * adcMeanCan = new TCanvas("canv6", "canv6", 1600, 1600);
	adcMeanCan->Divide(4,4);			
	
	for(int pixel = 0; pixel < 8; pixel++){
		
		adcMeanCan->cd(pad);
		adc_mean_distributions[30][pixel]->Draw();
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
