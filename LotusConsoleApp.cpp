// LotusConsoleApp.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <windows.h>
#include <iostream>
#include <vector>
#include <deque>
#include <string>
#include <bitset>
#include "..\iwxLotusAPI\iwxLotusAPI.h"

//#include <vld.h>   // visual leak detector. 


using namespace std;
const int string_buf_size = 256;
enum {
	SUCCESS,
};

#define CHECK_BIT(var,pos) ((var) & (1<<(pos)))

bool is_video_on = false;
bool manual_stimulate = false;
bool manual_stimulation_stop = false;
bool update_stim_amplitude = false;
bool communication_error = false;
float stim1_amp = 2;
// when there is an error in communications the Error Callback function is called.
void __stdcall ErrorNotification(int p_param)
{
	if(p_param == CONTROL_HARDWARE_NOT_FOUND)	wcout << "Error Communicating with Control device";
	if (p_param == USB_COMMUNICATION_ERROR) { // During acquisition if there is a error in USB communications with the stimulator or acquisition boards
		// or during EMC test, this error function will be called.
		wcout << "Error Communicating with Acquisition or Stimulator";
		// if recording was in progress. Stop the recording and restart the recording. 
		// software should remember what the settings were and resend them down. 
		communication_error = true;		
	}

	if (p_param == SYNC_SIGNAL_ERROR) { // stimulator board is not getting the sync signal from the acquisition board
		wcout << "Error Stimulator not receiving sync pulses";
	}
}

// When the user makes a change to the control module, This callback function is called by the API
void __stdcall ControlParamChangeNotification(uint8_t* p_param)
{
	if(p_param[SW_Stim1_5_RIGHT] )  wcout << "Control Stim 5 Right Switch\n ";
	if (p_param[SW_Stim1_5_LEFT]) {
		wcout << "Control Stim 5 Left Switch\n ";
		// set the stimulator mode to : LED_GRN_Stim1_5_mode to CC50
	}
	if (p_param[SW_Stim1_5_ON_OFF]) {
		wcout << "Control Stim 5 On Off Switch\n ";
		manual_stimulate = true;
	}
	if (p_param[SW_Stim2_6_RIGHT] )  wcout << "Control Stim 6 Right Switch\n ";
	if (p_param[SW_Stim2_6_LEFT] )  wcout << "Control Stim 6 Left Switch\n ";
	if (p_param[SW_Stim2_6_ON_OFF] )  wcout << "Control Stim 6 On Off Switch\n ";

	if (p_param[SW_VIDEO_ON_OFF])  wcout << "Control Video On-Off Switch\n ";
	if (p_param[SW_SCREENSHOT] )  wcout << "Control Screenshot Switch\n ";
	/*
	   The change inthe Volume and the Stim1 and Stim2 amplitude  will be reported.  Say the user turns the volume knob positive by 3 clicks
	   The value reported back will be 103, ie corresponding to a positive change of 3. 
	   The software can use this information to calculate the required change. 
	   This way the software, can keep track of what the amplitude should be,
	   and not have to tell the control board about it. 
	   It give the software the freedom to set the step size of each change, ie the user may want each click to be worth 1mA or 0.1mA or 5mA
	*/
	if (p_param[Volume_AMP] != 100) { 
		wcout << "Control Volume Amp :"; wcout << ((int)p_param[Volume_AMP] - 100); wcout << "\n ";
	}
	if (p_param[STIM1_AMP] != 100) { // use this to change the amplitude on the fly
		int amp_change = (int)p_param[STIM1_AMP] - 100;
		wcout << "Control STIM1 Amp :"; wcout << (amp_change); wcout << "\n ";
		stim1_amp += (float)amp_change / 10.0;
		update_stim_amplitude = true;
	}
	if (p_param[STIM2_AMP] != 100) {
		wcout << "Control STIM2 Amp :"; wcout << ((int)p_param[STIM2_AMP] - 100); wcout << "\n ";
	}
//	if (p_param[STIM1_MODE] )  wcout << "Control Stim 1 Mode\n ";
//	if (p_param[STIM2_MODE] )  wcout << "Control Stim 2 Mode\n ";
}


const int NUM_SAMPLES_PER_LOOP = 5000;
class CData {
public:
	CData() {
		// the data buffer has to be big enough to hold this data, 
// size of data buffer should be greater than num_samples_per_ch*num_channels_recorded
		m_data.resize(16 + 6); // we have upto 16 acq channels + 6 stimulator;
		m_filtered_data.resize(16 + 6); // we have upto 16 acq channels + 6 stimulator;
		number_of_samples_to_read_per_loop = 0; // no memory allocated yet.
		acqdata = stimdata = NULL;
	};
	~CData() {
		FreeMemory();
	}
	bool AllocateMemory(int samples_per_loop) {
		FreeMemory();
		acqdata = (float*) malloc(samples_per_loop*16 * sizeof(float));
		stimdata = (float*)malloc(samples_per_loop * 6 * sizeof(float));
		number_of_samples_to_read_per_loop = samples_per_loop;
		if (acqdata && stimdata) return true;
		number_of_samples_to_read_per_loop = 0;
		FreeMemory();//  error free the memory allocated. 
		return false;
	}
	bool FreeMemory() {
		if (acqdata) {
			free(acqdata);
			acqdata = NULL;
		}
		if (stimdata) {
			free(stimdata);
			stimdata = NULL;
		}
		return true;
	}
	std::vector<std::vector<float>> m_data;
	std::vector<std::vector<float>> m_filtered_data;
	float* acqdata;// [NUM_SAMPLES_PER_LOOP * 16] ;// greater than num_acq_channels_recorded*number_of_samples_to_read_per_loop
	float* stimdata;// [NUM_SAMPLES_PER_LOOP * 6] ;
	int number_of_samples_to_read_per_loop;
	int m_sampling_speed;
};

class ChannelImpedance {
public:
	float m_pos_impedance;
	float m_neg_impedance;
	bool m_isvalid;
};

class ChannelTestValues {

	double ch_mean_val;
	double ch_max_val;
	double ch_min_val;
public:	
	int count;

	void Reset() {
		ch_mean_val = count =  0;
		ch_max_val = -100000;
		ch_min_val = +100000;
	}
	void AddValue(double val) {
		if (ch_max_val < val)ch_max_val = val;
		if (ch_min_val > val)ch_min_val = val;
		ch_mean_val += val;
		++count;
	}

	int GetStat(double& max, double& min, double& mean) {
		max = ch_max_val;
		min = ch_min_val;
		if (count > 0) {
			mean = ch_mean_val / (double)count;
		}
		return count;
	}
	double GetMaxMin() {
		return ch_max_val - ch_min_val;
	}
};




void ShowErrorMesage(int err) {
	const wchar_t* err_msg = GetErrorMessage(err);
	wcout << err_msg;
	wcout << L"\n";
};

int GetStats(std::vector<ChannelTestValues>& ch_stat, std::vector<std::vector<float>>& data, int num_acq_channels_recorded, int total_num_samples_per_ch, int pts_to_skip, bool print_stats, FILE* fout ) {
	int k = 0;
	for (k = 0; k < num_acq_channels_recorded; ++k) {
		ch_stat[k].Reset();
	}
	for (int j = 0; j < total_num_samples_per_ch; ++j) {
		for (k = 0; k < num_acq_channels_recorded; ++k) {
			if (j > pts_to_skip) ch_stat[k].AddValue(data[k].at(j)); 
		}
	}
	if( print_stats){
		if (fout) {
			fprintf(fout, "Printing Stats for the Recording\n\n");
			fprintf(fout, "Channel, Mean, Max, Min, Max-Min, RMS \n");
		}
		double max, min, mean;
		for (k = 0; k < num_acq_channels_recorded; ++k) {
			ch_stat[k].GetStat(max, min, mean);
			// calculate RMS value
			double rms = 0;
			for (int j = pts_to_skip; j < total_num_samples_per_ch; ++j) {
				double diff = (data[k].at(j) - mean);
				rms += diff * diff;
			}
			rms = sqrt(rms / (total_num_samples_per_ch - pts_to_skip));

			wcout << "Channel "; wcout << k + 1;
			wcout << " Mean = "; wcout << mean;
			wcout << " Max = "; wcout << max;
			wcout << " Min = "; wcout << min;
			wcout << " Max-Min = "; wcout << max - min;
			wcout << " RMS = "; wcout << rms;
			wcout << "\n";
			if (fout) {
				fprintf(fout, "%d, %g, %g, %g, %g, %g \n", k+1, mean,max, min, max-min,  rms);
			}
		}
	}
	return SUCCESS;
}

int RecordData(CData& emgdata, int& total_num_samples_per_ch_recorded, int total_number_of_samples_per_ch_to_record, int num_acq_channels_recorded, int num_stim_channels, bool show_leadoff) {
	 // total_number_of_samples_per_ch_to_record  is the number of samples per channel you want to record
	//total_num_samples_per_ch_recorded is the number that is actually recorded. 
	int iRet;
	total_num_samples_per_ch_recorded = 0;
	int total_num_stim_samples_per_ch_recorded = 0;
	int num_acq_samples_per_ch = 0;
	int j, m, k;
	Sleep(10);
	StartAcq();
	for (k = 0; k < num_acq_channels_recorded; ++k) {
		emgdata.m_data[k].clear();
	}

	Sleep(10);
	uint32_t ch_leadoff;
	
	uint16_t  enable_online_notch_filter = 0x0D; // 1101b  we are only applying the notch filter to Channel 1, 3, 4 
	if (enable_online_notch_filter) {
		OnlineNotchFilterSetup(emgdata.m_sampling_speed, 50, false);
	}


	while (total_num_samples_per_ch_recorded < total_number_of_samples_per_ch_to_record) { //continue recording data until stop condition is reached.
		// if you want to user control, run this function in a separate thread and have it look for a varioble such as m_is_recording_data
		// then check the value of m_is_recording_data  in the above while loop.
		num_acq_samples_per_ch = 0; // we need to set this to zero so we read all the data available 
		//depending on the speed of the application there is a time delay between subsequent calls of ReadDataFromAcquisitionDevice in the while loop
		// during this time the API is collecting data and storing it in a buffer.
		// makes sure that emgdata.number_of_samples_to_read_per_loop is always greater than the num_acq_samples_per_ch returned by this function
		iRet = ReadDataFromAcquisitionDevice(num_acq_samples_per_ch, emgdata.acqdata, emgdata.number_of_samples_to_read_per_loop);
		// check ch_leadoff to make sure leads were connected.
		if (iRet == 0) { // acquired data copied to data.  pts_recorded = num_samples_per_ch*num_channels_recorded;
			// check if leadoff was ok
			if (show_leadoff) {
				GetAcquisitionChannelLeadOff(ch_leadoff);  //  0 - 15 ch 
				if (ch_leadoff > 0) { // display leadoff if channel disconnected.
					wcout << "LeadOff Channels : ";
					wcout << ch_leadoff; wcout << "  Num acq Samples: "; wcout << num_acq_samples_per_ch;  wcout << "\n";
				}
			}
			m = 0;
			for (j = 0; j < num_acq_samples_per_ch; ++j) {
				for (k = 0; k < num_acq_channels_recorded; ++k) {
					if (CHECK_BIT(enable_online_notch_filter,k)) { // do an online notch filter on the data.
						float data = emgdata.acqdata[m++];
						emgdata.m_data[k].push_back(OnlineNotchFilterData(k, data));
					}
					else {
						emgdata.m_data[k].push_back(emgdata.acqdata[m++]);
					}
				}
			}
			total_num_samples_per_ch_recorded += num_acq_samples_per_ch;
			wcout << '\r';
			wcout << "Number of samples Recorded : ";  wcout << total_num_samples_per_ch_recorded; // wcout << "\n";
			wcout << std::flush;
		}
		if (iRet > 0) { // Error condition 
			StopAcq();// Stop Acquisition
			// Show the error message if needed
			ShowErrorMesage(iRet);
			
			return iRet;
		}
		if (iRet < 0) { // warning, 
			//typical warning -3 ==> that no new data is available wait longer
		}
		if (num_stim_channels > 0) {
			//num_acq_samples_per_ch will be what is returned from the ReadDataFromAcquisitionDevice, 
			//we want to read the same number of data points from the stimulator device
			iRet = ReadDataFromStimulatorDevice(num_acq_samples_per_ch, emgdata.stimdata, emgdata.number_of_samples_to_read_per_loop);
			if (iRet == 0) { // acquired data copied to data.  pts_recorded = num_samples_per_ch*num_channels_recorded;
				m = 0;
				for (j = 0; j < num_acq_samples_per_ch; ++j) {
					for (k = num_acq_channels_recorded; k < num_acq_channels_recorded+ num_stim_channels; ++k) { // append the stimulator channels after the acquisition channels
						emgdata.m_data[k].push_back(emgdata.stimdata[m++]);
					}
				}
			}
			if (iRet > 0) { // Error condition 
				StopAcq();// Stop Acquisition
				
				return iRet;
			}
			if (iRet < 0) { // warning, 
							//typical warning -3 ==> that no new data is available wait longer
			}
			total_num_stim_samples_per_ch_recorded += num_acq_samples_per_ch;
		}
	}
	if (total_num_samples_per_ch_recorded > total_number_of_samples_per_ch_to_record) total_num_samples_per_ch_recorded = total_number_of_samples_per_ch_to_record; // only interested in total_number_of_samples_per_ch_to_record drop the rest of the data
	wcout << "\n";
	StopAcq();// Stop Acquisition
	
	return iRet;
}

int PrintDatatoFile(FILE* fout, std::vector<std::vector<float>>& data, int total_num_samples_per_ch, int num_acq_channels_recorded) {
	int k = 0;
	for (k = 0; k < num_acq_channels_recorded; ++k) {
		if (total_num_samples_per_ch > data[k].size()) total_num_samples_per_ch = data[k].size();
	}
	if(fout){
		for (int j = 0; j < total_num_samples_per_ch; ++j) {
			for (int k = 0; k < num_acq_channels_recorded; ++k) {
				fprintf(fout, "%g,", data[k].at(j));
			}
			fprintf(fout, "\n");
		}
	}
	return 0;
}

int RecordTEMG(CData& emgdata, FILE* fout)
{
	//  now lets record 100 sweeps, that are 100msec ( 0.1sec)  in duration and then average them to get one average sweep
	// we will be recording all the 100 sweeps continously at 2000 samples per sec
	// if we want to filter the data it would make things easy if we recorded 2 extra sweeps one at the begining and one at the end, since the sweep duration is so small.

	int number_of_sweeps_in_group = 5; // a group of sweeps will have the same stimulation amplitude
	int number_of_groups_of_sweeps = 1;
	stim1_amp = 10; // this is the starting stim amplitude, this can be changed with the control module using function ControlParamChangeNotification

	float m_sweep_length_in_sec = 0.5;

	//	wcout << "Set Sampling speed: \n";
	int speed = 16000;
	//	wcin >> speed;

	emgdata.m_sampling_speed = speed;

	int samples_per_sweep = emgdata.m_sampling_speed * m_sweep_length_in_sec;
	int total_number_of_samples_to_record = number_of_sweeps_in_group * emgdata.m_sampling_speed * m_sweep_length_in_sec;
	// setup the stimulator

	SetStimHVPowerSupply(true); // it takes about 300msec for the HV power supply to settle.

	CStimulationTrain stim_param;
	//	wcout << "Set Delay time in msec: \n";
	int delaytime = 110;
	//	wcin >> delaytime;

	stim_param.m_delay = delaytime;//time in msec wait for 10msec
	stim_param.m_num_pulses = 1; // 1 pulse
	stim_param.m_bipolar = false;
	stim_param.m_pulse_width = 0.25; // 0.25 msec pulse
	stim_param.m_pulse_off_time = 1.75; // 29msec off time . m_pulse_off_time + m_pulse_width = pulse Period 
	stim_param.m_num_trains = number_of_sweeps_in_group;  // set to 0 for continous trains, we have 100 sweeps that we are recording continously
	// since we want one set of pulses for every train,
	// each train duration has to match the sweep length, that way all the trains will align for all the sweeps
	// sweep length - delay - m_num_pulses*(m_pulse_width + m_pulse_off_time)
	stim_param.m_intertrain_duration = m_sweep_length_in_sec * 1000 - stim_param.m_num_pulses * (stim_param.m_pulse_width + stim_param.m_pulse_off_time);
	stim_param.m_pulse_amplitude = stim1_amp; // 10 V, for CV20 ; 10mA if CC20
	stim_param.m_stim_mode = CC20;//  CC20;
	stim_param.m_start_with_recording = true;

//	wcout << "Stimulator Channel: \n";
	int stim_number = 1;
//	wcin >> stim_number;



	int ret = SetStimulationParameters(stim_number - 1, stim_param);
	if (ret > SUCCESS) { // error in stimulation parameters
		ShowErrorMesage(ret);
	}

	// setup the acquisition
	int num_acq_channels_recorded = 1;
	int num_stim_channels_recorded = 1;  // this is the number of stimulator channels to record
	unsigned int channels_to_record = 0x01;// record channel 1 and channel 3 , ch 2 , channel 4-16 are skipped
	unsigned int stim_channels_to_record = 0x01; // because we are recording data from Stim5 only 
	if (stim_number == 2)stim_channels_to_record = 0x02;
	/*  if we want to record only from stim6 then we would set  unsigned int stim_channels_to_record = 0x02;   and num_stim_channels_recorded = 1
	* since 0x02  sets the second bit, which corresponds to the second stimulator channel
	* If you want both stim 5 and stim 6, then set  unsigned int stim_channels_to_record = 0x03; and num_stim_channels_recorded = 2
	*since 0x03  sets the first and second bit, which corresponds to the first and second stimulator channel
	*/
	//	int num_stim_ch_to_record = 1;

	SetAcquisitionParameters(emgdata.m_sampling_speed, LOTUS_RECORD, channels_to_record, stim_channels_to_record);

	wcout << "EMG Sweeps with Stimulation\n";
	if (fout) fprintf(fout, "EMG Sweeps with Stimulation\n");
	// size of data buffer should be greater than num_samples_per_ch*num_channels_recorded 

	int number_of_samples_to_read_per_loop = 1000;
	float lowcutoff = 100;
	float highcutoff = 500;
	bool use_bandpass = true;
	bool use_notchfilter = false;
	int order = 51;
	int i, j, k, m;
	
	for (int s = 0; s < number_of_groups_of_sweeps; ++s) {
		wcout << "Recording Sweep Group: "; wcout << s + 1; wcout << "\n";
		int total_num_samples_per_ch = 0;
		int iRet = RecordData(emgdata, total_num_samples_per_ch, total_number_of_samples_to_record, num_acq_channels_recorded, num_stim_channels_recorded, false);
		if (iRet > 0) { // Error condition 
			ShowErrorMesage(iRet);
		}
		else {
			// make sure that we have enough data
			for (i = 0; i < num_acq_channels_recorded + stim_channels_to_record; ++i) {
				emgdata.m_filtered_data[i].resize(total_num_samples_per_ch);
				emgdata.m_data[i].resize(total_num_samples_per_ch, 0); // incase eonough data was not recorded, fill data with zeros. 
				if (i < num_acq_channels_recorded) { // only filter acquisition channels
					if (use_bandpass) {
						// filter the channel data with a band pass filter 
						FilterData(&(emgdata.m_data[i][0]), &(emgdata.m_filtered_data[i][0]), total_num_samples_per_ch, emgdata.m_sampling_speed, lowcutoff, highcutoff, order);
					}
					if (use_notchfilter) {
						NotchFilterData(&(emgdata.m_data[i][0]), &(emgdata.m_filtered_data[i][0]), total_num_samples_per_ch, emgdata.m_sampling_speed, 50, false, true);
					}
				}
				else {// just copy the stimulation data
					emgdata.m_filtered_data[i] = emgdata.m_data[i];
				}
			}
			if (fout) {
				fprintf(fout, "Raw Data\n");
				PrintDatatoFile(fout, emgdata.m_data, total_num_samples_per_ch, num_acq_channels_recorded + num_stim_channels_recorded);
				fprintf(fout, "Filtered Data\n");
				PrintDatatoFile(fout, emgdata.m_filtered_data, samples_per_sweep, num_acq_channels_recorded + num_stim_channels_recorded);
			}
		}
	}
	SetStimHVPowerSupply(false); // Turn off the HV power supply.

	return SUCCESS;
}

int RecordSweepsWithStimulation(CData& emgdata, FILE* fout)
{
//  now lets record 100 sweeps, that are 100msec ( 0.1sec)  in duration and then average them to get one average sweep
// we will be recording all the 100 sweeps continously at 2000 samples per sec
// if we want to filter the data it would make things easy if we recorded 2 extra sweeps one at the begining and one at the end, since the sweep duration is so small.

	int number_of_sweeps_to_average_in_group = 20; // a group of sweeps will have the same stimulation amplitude
	int number_of_groups_of_sweeps = 100; 
	stim1_amp = 1; // this is the starting stim amplitude, this can be changed with the control module using function ControlParamChangeNotification

	float m_sweep_length_in_sec = 0.5;

//	wcout << "Set Sampling speed: \n";
	int speed = 16000;
//	wcin >> speed;

	emgdata.m_sampling_speed = speed;

	int samples_per_sweep = emgdata.m_sampling_speed * m_sweep_length_in_sec;
	int total_number_of_samples_to_record = number_of_sweeps_to_average_in_group * emgdata.m_sampling_speed * m_sweep_length_in_sec;
	// setup the stimulator

	SetStimHVPowerSupply(true); // it takes about 300msec for the HV power supply to settle.

	CStimulationTrain stim_param;
//	wcout << "Set Delay time in msec: \n";
	int delaytime = 110;
//	wcin >> delaytime;

	stim_param.m_delay = delaytime;//time in msec wait for 10msec
	stim_param.m_num_pulses = 5; // 1 pulse
	stim_param.m_bipolar = false;
	stim_param.m_pulse_width = 0.5; // 0.25 msec pulse
	stim_param.m_pulse_off_time = 4.5; // 29msec off time . m_pulse_off_time + m_pulse_width = pulse Period 
	stim_param.m_num_trains = number_of_sweeps_to_average_in_group;  // set to 0 for continous trains, we have 100 sweeps that we are recording continously
	// since we want one set of pulses for every train,
	// each train duration has to match the sweep length, that way all the trains will align for all the sweeps
	// sweep length - delay - m_num_pulses*(m_pulse_width + m_pulse_off_time)
	stim_param.m_intertrain_duration = m_sweep_length_in_sec * 1000 - stim_param.m_num_pulses * (stim_param.m_pulse_width + stim_param.m_pulse_off_time);
	stim_param.m_pulse_amplitude = stim1_amp; // 10 V, for CV20 ; 10mA if CC20
	stim_param.m_stim_mode = CC20;//  CC20;
	stim_param.m_start_with_recording = true;

//	wcout << "Stimulator Channel: \n";
	int stim_number = 1;
//	wcin >> stim_number;



	int ret = SetStimulationParameters(stim_number-1, stim_param);
	if (ret > SUCCESS) { // error in stimulation parameters
		ShowErrorMesage(ret);
	}

	// setup the acquisition
	int num_acq_channels_recorded = 2;
	int num_stim_channels_recorded = 1;  // this is the number of stimulator channels to record
	unsigned int channels_to_record = 0x05;// record channel 1 and channel 3 , ch 2 , channel 4-16 are skipped
	unsigned int stim_channels_to_record = 0x01; // because we are recording data from Stim5 only 
	if(stim_number == 2 )stim_channels_to_record = 0x02;
	/*  if we want to record only from stim6 then we would set  unsigned int stim_channels_to_record = 0x02;   and num_stim_channels_recorded = 1
	* since 0x02  sets the second bit, which corresponds to the second stimulator channel
	* If you want both stim 5 and stim 6, then set  unsigned int stim_channels_to_record = 0x03; and num_stim_channels_recorded = 2
	*since 0x03  sets the first and second bit, which corresponds to the first and second stimulator channel
	*/
//	int num_stim_ch_to_record = 1;

	SetAcquisitionParameters(emgdata.m_sampling_speed, LOTUS_RECORD, channels_to_record, stim_channels_to_record);
	
	wcout << "EMG Sweeps with Stimulation\n";
	if(fout) fprintf(fout, "EMG Sweeps with Stimulation\n");
	// size of data buffer should be greater than num_samples_per_ch*num_channels_recorded 
//	int number_of_samples_to_read_per_loop = 1000;

	int total_num_samples_per_ch = 0;
	int iRet = RecordData(emgdata, total_num_samples_per_ch, total_number_of_samples_to_record, num_acq_channels_recorded, num_stim_channels_recorded, false);
	SetStimHVPowerSupply(false); // Turn off the HV power supply.
	if (iRet > 0) { // Error condition 
		ShowErrorMesage(iRet);
	}
	else {
		float lowcutoff = 100;
		float highcutoff = 500;
		bool use_bandpass = true;
		bool use_notchfilter = false;
		int order = 51;
		int i, j, k, m;
		// make sure that we have enough data

		for (i = 0; i < num_acq_channels_recorded + stim_channels_to_record; ++i) {
			emgdata.m_filtered_data[i].resize(total_num_samples_per_ch);
			int size_of_array = emgdata.m_data[i].size();
			if(size_of_array < total_num_samples_per_ch)   emgdata.m_data[i].resize(total_num_samples_per_ch, 0); // incase eonough data was not recorded, fill data with zeros. 
			if (i < num_acq_channels_recorded) { // only filter acquisition channels
				if (use_bandpass) {
					// filter the channel data with a band pass filter 
					FilterData(&(emgdata.m_data[i][0]), &(emgdata.m_filtered_data[i][0]), total_num_samples_per_ch, emgdata.m_sampling_speed, lowcutoff, highcutoff, order);
				}
				if (use_notchfilter) {
					NotchFilterData(&(emgdata.m_data[i][0]), &(emgdata.m_filtered_data[i][0]), total_num_samples_per_ch, emgdata.m_sampling_speed, 50, false, true);
				}
			}
			else {// just copy the stimulation data
				emgdata.m_filtered_data[i] = emgdata.m_data[i];
			}

			// now average the sweeps skipping the first sweep
			//for (j = 0; j < samples_per_sweep; ++j) {
			//	(emgdata.m_filtered_data[i])[j] = 0;// the first sweep was extra 
			//}
			//m = samples_per_sweep; // lets start at the begining of the second sweep
			//for (k = 0; k < num_sweep; ++k) {
			//	for (j = 0; j < samples_per_sweep; ++j) {
			//		(emgdata.m_filtered_data[i])[j] += (emgdata.m_filtered_data[i])[m++];
			//	}
			//}
			// now find the average
			//for (j = 0; j < samples_per_sweep; ++j) {
			//	(emgdata.m_filtered_data[i])[j] = (emgdata.m_filtered_data[i])[j] / num_sweep;
			//}
		}
		if (fout) {
			fprintf(fout, "Raw Data\n");
			PrintDatatoFile(fout, emgdata.m_data, total_num_samples_per_ch, num_acq_channels_recorded + num_stim_channels_recorded);
			fprintf(fout, "Filtered Data\n");
			PrintDatatoFile(fout, emgdata.m_filtered_data, samples_per_sweep, num_acq_channels_recorded + num_stim_channels_recorded);
		}
	}

	return SUCCESS;
}
int TestDataLoss(CData& emgdata, FILE* fout)
{

	return SUCCESS;
}

int RecordSweepsWithStimulationAmplitudeChange(CData& emgdata, FILE* fout)
{
	// if you want to average 5 sweeps, they will have the same amplitude, these sweeps will be in a group.
	// and there can be a a number of these groups,  You can change the amplitude of the stimulation between groups

	int number_of_sweeps_to_average_in_group = 5; // a group of sweeps will have the same stimulation amplitude  this can be 1
	int number_of_groups_of_sweeps = 100;

	stim1_amp = 1; // this is the starting stim amplitude, this can be changed with the control module using function ControlParamChangeNotification

	float m_sweep_length_in_sec = 1;

	emgdata.m_sampling_speed = 16000;

	int samples_per_sweep = emgdata.m_sampling_speed * m_sweep_length_in_sec;
	int total_number_of_samples_to_record = number_of_sweeps_to_average_in_group * emgdata.m_sampling_speed * m_sweep_length_in_sec;
	// setup the stimulator

	SetStimHVPowerSupply(true); // it takes about 300msec for the HV power supply to settle.

	CStimulationTrain stim_param;
	wcout << "Set Delay time in msec: \n";
	int delaytime = 10;
	wcin >> delaytime;

	stim_param.m_delay = delaytime;//time in msec wait for 10msec
	stim_param.m_num_pulses = 2; // 1 pulse
	stim_param.m_bipolar = false;
	stim_param.m_pulse_width = 0.25; // 0.25 msec pulse
	stim_param.m_pulse_off_time = 9.75; // 29msec off time . m_pulse_off_time + m_pulse_width = pulse Period 
	stim_param.m_num_trains = number_of_sweeps_to_average_in_group;  // set to 0 for continous trains, we have 100 sweeps that we are recording continously
	// since we want one set of pulses for every train,
	// each train duration has to match the sweep length, that way all the trains will align for all the sweeps
	// sweep length - delay - m_num_pulses*(m_pulse_width + m_pulse_off_time)
	stim_param.m_intertrain_duration = m_sweep_length_in_sec * 1000 - stim_param.m_num_pulses * (stim_param.m_pulse_width + stim_param.m_pulse_off_time);
	stim_param.m_pulse_amplitude = stim1_amp; // 10 V, for CV20 ; 10mA if CC20
	stim_param.m_stim_mode = CC20;//  CC20;
	stim_param.m_start_with_recording = true;

	int ret = SetStimulationParameters(0, stim_param);
	if (ret > SUCCESS) { // error in stimulation parameters
		ShowErrorMesage(ret);
	}

	// setup the acquisition
	int num_acq_channels_recorded = 2;
	int num_stim_channels_recorded = 1;  // this is the number of stimulator channels to record
	unsigned int channels_to_record = 0x05;// record channel 1 and channel 3 , ch 2 , channel 4-16 are skipped
	unsigned int stim_channels_to_record = 0x01; // because we are recording data from Stim5 only 
	/*  if we want to record only from stim6 then we would set  unsigned int stim_channels_to_record = 0x02;   and num_stim_channels_recorded = 1
	* since 0x02  sets the second bit, which corresponds to the second stimulator channel
	* If you want both stim 5 and stim 6, then set  unsigned int stim_channels_to_record = 0x03; and num_stim_channels_recorded = 2
	*since 0x03  sets the first and second bit, which corresponds to the first and second stimulator channel
	*/
	int num_stim_ch_to_record = 1;

	SetAcquisitionParameters(emgdata.m_sampling_speed, LOTUS_RECORD, channels_to_record, stim_channels_to_record);

	wcout << "EMG Sweeps with Stimulation\n";
	if (fout) fprintf(fout, "EMG Sweeps with Stimulation\n");
	// size of data buffer should be greater than num_samples_per_ch*num_channels_recorded 

	int iRet = SUCCESS;
	for (int g = 0; g < number_of_groups_of_sweeps && iRet == SUCCESS; ++g) {
		stim_param.m_pulse_amplitude = stim1_amp; // 10 V, for CV20 ; 10mA if CC20
		iRet = SetStimulationParameters(0, stim_param);
		int total_num_samples_per_ch = 0;
		iRet = RecordData(emgdata, total_num_samples_per_ch, total_number_of_samples_to_record, num_acq_channels_recorded, num_stim_ch_to_record, false);
		if (fout) {
			fprintf(fout, "Raw Data\n");
			PrintDatatoFile(fout, emgdata.m_data, total_num_samples_per_ch, num_acq_channels_recorded + num_stim_channels_recorded);
		}
	}
	SetStimHVPowerSupply(false); // Turn off the HV power supply.
	if (iRet > 0) { // Error condition 
		ShowErrorMesage(iRet);
	}
	return SUCCESS;
}

void TestControlBoard() {
	
	uint8_t p_control_parameters[NUM_CTRL_LED_PARAM];
	wstring p_control_name[NUM_CTRL_LED_PARAM];

	p_control_name[LED_VIDEO] = L"Video: ";
	p_control_name[LED_RED_Stim2_6_Alarm] = L"Stim2 Alarm: ";
	p_control_name[LED_GRN_Stim2_6_Selected] = L"Stim2 Selected: ";
	p_control_name[LED_RED_Stim1_5_Alarm] = L"Stim1 Alarm: ";
	p_control_name[LED_GRN_Stim1_5_Selected] = L"Stim1 Selected: ";
	p_control_name[LED_GRN_Stim1_5_On] = L"Stim1 On: ";
	p_control_name[LED_GRN_Stim2_6_On] = L"Stim2 On: ";
	p_control_name[LED_SCREENSHOT] = L"Screenshot: ";
	p_control_name[LED_GRN_Stim2_6_mode] = L"Stim2 mode: ";
	p_control_name[LED_GRN_Stim2_6_Dynamic1] = L"Stim2 Dynamic1: ";
	p_control_name[LED_GRN_Stim2_6_Dynamic2] = L"Stim2 Dynamic2: ";
	p_control_name[LED_GRN_Stim2_6_Dynamic3] = L"Stim2 Dynamic3: ";
	p_control_name[LED_GRN_Stim1_5_mode] = L"Stim1 mode: ";
	p_control_name[LED_GRN_Stim1_5_Dynamic1] = L"Stim1 Dynamic1: ";
	p_control_name[LED_GRN_Stim1_5_Dynamic2] = L"Stim1 Dynamic2: ";
	p_control_name[LED_GRN_Stim1_5_Dynamic3] = L"Stim1 Dynamic3: ";

	wcout << "Testing Control Board\n ";
	wcout << "Setting Control Parameters\n Press Q to Quit\n";
	wchar_t ch;
	wcin >> ch;
	while (!(ch == 'Q' || ch == 'q')) {
		const uint8_t* p_ctrlparam_read = GetControlStimulatorParameters();
		// save the parameters 
		
		for (int i = 0; i < NUM_CTRL_LED_PARAM; ++i) {
			for (int j = 0; j < NUM_CTRL_LED_PARAM; ++j) {
				p_control_parameters[j] = 0;
			}
			p_control_parameters[i] = 1;
			SetControlStimulatorParameters(p_control_parameters);
			wcout << "Setting Control Parameters ";
			wcout << p_control_name[i];
			wcout << "\n";
			Sleep(2000); // this is so that we can see the change on the control board visually to check that the correct LED is turnign on
		}
		// test stim mode leds
		for (int j = 0; j < NUM_CTRL_LED_PARAM; ++j) {
			p_control_parameters[j] = 0;
		}
		for (int j = 0; j < NUM_STIM_MODES; ++j) {
			p_control_parameters[LED_GRN_Stim1_5_mode] = j;
			SetControlStimulatorParameters(p_control_parameters);
			wcout << "Setting Stim1 mode ";
			wcout << j;
			wcout << "\n";
			Sleep(2000);// this is so that we can see the change on the control board visually to check that the correct LED is turnign on
		}

		// turn the Stim Mode to Off if the stimulator is not being used:	
		p_control_parameters[LED_GRN_Stim1_5_mode] = S_OFF;
		SetControlStimulatorParameters(p_control_parameters);

		for (int j = 0; j < NUM_CTRL_LED_PARAM; ++j) {
			p_control_parameters[j] = 0;
		}
		for (int j = 0; j < NUM_STIM_MODES; ++j) {
			p_control_parameters[LED_GRN_Stim2_6_mode] = j;
			SetControlStimulatorParameters(p_control_parameters);
			wcout << "Setting Stim2 mode ";
			wcout << j; 
			wcout << "\n";
			Sleep(2000);// this is so that we can see the change on the control board visually to check that the correct LED is turnign on
		}
		// turn the Stim Mode to Off if the stimulator is not being used:	
		p_control_parameters[LED_GRN_Stim2_6_mode] = S_OFF;
		SetControlStimulatorParameters(p_control_parameters);

		wcin >> ch;
	}		// change the parameters 	
}

int TestFilterTimeOffset(CData& emgdata) {
	emgdata.m_sampling_speed = 16000;
	int num_points = 16000;
	int k = 0;
	for ( k = 0; k < 1; ++k) {
		emgdata.m_data[k].clear();
		emgdata.m_data[k].resize(num_points, 0);
		emgdata.m_filtered_data[k].clear();
		emgdata.m_filtered_data[k].resize(num_points, 0);
	}
	int i = 0;
	for (i = 4000; i < 4010; ++i) {
		emgdata.m_data[0][i] = 1;
	}
	int filter_order = 1501;
	k = 0;
	for (int hp = 10; hp < 100; hp = hp + 50) {
		filter_order = 1501;
		int reported_offset = OnlinebandPassFilterSetup(-1, emgdata.m_sampling_speed, hp, 2000, filter_order); // setting all filters the same.
		int pulse_start_index = -1;
		for (i = 0; i < num_points; ++i) {
			emgdata.m_filtered_data[0][i] = OnlineBandpassFilterData(0, emgdata.m_data[0][i]);
			if (pulse_start_index < 0 && emgdata.m_filtered_data[0][i] > 0.5)pulse_start_index = i;
		}
		wcout << "Filter Order : "; wcout << filter_order;
		wcout << "\tHigh Pass: "; wcout << hp;
		wcout << "\tTime Offset: "; wcout << (pulse_start_index - 4000);
		wcout << "\tReported Time Offset: "; wcout << reported_offset;
		
		wcout << "\n";
	}


	return SUCCESS;
}

int Record_1ch_SEMG(CData& emgdata, FILE* acqfout, int recording_time) {
	//S - EMG：16 channel, 16Ksps, 
	//Notch filter is enabled, disable second harmonic，disable third harmonic；
	//bandpass filter is enabled, 16Ksps, order is odd(11~1001) adjustable; high - cut and low - cut can be adjustable；
//Display : dynamtically display the waves;
		// turn off all stimulators 
	SetStimulatorOff(0);
	SetStimulatorOff(1);
	SetStimulatorOff(2);
	SetStimulatorOff(3);
	SetStimulatorOff(4);
	SetStimulatorOff(5);

	// the controller board should also be updated to show that the stimulators are off.
	uint8_t p_control_parameters[NUM_CTRL_LED_PARAM];
	for (int j = 0; j < NUM_CTRL_LED_PARAM; ++j) {
		p_control_parameters[j] = 0;
	}
	SetControlStimulatorParameters(p_control_parameters);

	// Record 16 channels 
	int num_acq_channels_recorded = 1;
	int num_stim_channels_recorded = 0;
	unsigned int channels_to_record = 0x01;// we are recording 1 channel
	unsigned int stim_channels_to_record = 0x00; // we are not recording any stimualtion channels 
	int j, m, k;
	emgdata.m_sampling_speed = 16000; // 1k sampling

	// NOTE : if we are filtering the data at 1khz there is not much point in sampling the data at 16Khz 

	OnlineNotchFilterSetup(emgdata.m_sampling_speed, 50, true); // if using the LOTUS_NOTCH_FILTER_50 do not use the OnlineNotchFilterSetup

	// to bandpass filter the data between 10hz and 1kHz
	int filter_order = 1601;
	//OnlinebandPassFilterSetup will change the filter order to match the High and Low Pass freq requested. 
	// for eg. if you set the filter order to 1001, and ask for a 10Hz HP, the filter order will be changed to 1601 since that is the min order required
	// if you set it to 2001 then the order will not be changed since it is greater than 1601
	for (k = 0; k < 16; ++k) {
		OnlinebandPassFilterSetup(k, emgdata.m_sampling_speed, 10, 1000, filter_order); // setting all filters the same.
	}


	SetAcquisitionParameters(emgdata.m_sampling_speed, LOTUS_RECORD_WITH_LEADOFF, channels_to_record, stim_channels_to_record);
	
	int iRet;
	int total_num_samples_per_ch_recorded = 0;
	int total_number_of_samples_per_ch_to_record = recording_time* emgdata.m_sampling_speed;
	int num_acq_samples_per_ch = 0;
	//clear the data buffers
	for (k = 0; k < num_acq_channels_recorded; ++k) {
		emgdata.m_data[k].clear();
		emgdata.m_filtered_data[k].clear();
		emgdata.m_data[k].reserve(total_number_of_samples_per_ch_to_record + 100);
		emgdata.m_filtered_data[k].reserve(total_number_of_samples_per_ch_to_record + 100);
	}
	Sleep(10);
	StartAcq();
	Sleep(10);
	uint32_t ch_leadoff;
	
	while (total_num_samples_per_ch_recorded < total_number_of_samples_per_ch_to_record) { //continue recording data until stop condition is reached.
		// if you want to user control, run this function in a separate thread and have it look for a varioble such as m_is_recording_data
		// then check the value of m_is_recording_data  in the above while loop.
		iRet = ReadDataFromAcquisitionDevice(num_acq_samples_per_ch, emgdata.acqdata, emgdata.number_of_samples_to_read_per_loop); //25600  num_acq_samples_per_ch=160
		// check if leadoff was ok
		GetAcquisitionChannelLeadOff(ch_leadoff);  //  0 - 15 ch   if ch_leadoff == 0 then all recorded channels are connected
		if (ch_leadoff > 0) { // display leadoff if channel disconnected. ch_leadoff has information about the disconnected channel
			wcout << "LeadOff Channels : ";
			wcout << ch_leadoff; wcout << "\n";
		}
		// check ch_leadoff to make sure leads were connected.
		if (iRet == 0) { // acquired data copied to data.  pts_recorded = num_samples_per_ch*num_channels_recorded;
			m = 0;
			for (j = 0; j < num_acq_samples_per_ch; ++j) {
				for (k = 0; k < num_acq_channels_recorded; ++k) {
					float data = emgdata.acqdata[m++];
					emgdata.m_data[k].push_back(data);
					float filtered_data = OnlineBandpassFilterData(k, data);
					filtered_data = OnlineNotchFilterData(k, filtered_data);
					emgdata.m_filtered_data[k].push_back(filtered_data);
				}
			}
			total_num_samples_per_ch_recorded += num_acq_samples_per_ch;
			wcout << "Number of samples Recorded : ";  wcout << total_num_samples_per_ch_recorded;  wcout << "\n";
		}
		if (iRet > 0) { // Error condition 
			StopAcq();// Stop Acquisition
			// Show the error message if needed
			ShowErrorMesage(iRet);

			return iRet;
		}
	}
	StopAcq();// Stop Acquisition
//	PrintDatatoFile(acqfout, emgdata.m_data, total_num_samples_per_ch_recorded, num_acq_channels_recorded + stim_channels_to_record);

	return iRet;
}

int Record_SEMG(CData& emgdata, FILE* acqfout, bool use_online_filters, bool cal_ch1 = false) {
	//S - EMG：16 channel, 16Ksps, 
	//Notch filter is enabled, disable second harmonic，disable third harmonic；
	//bandpass filter is enabled, 16Ksps, order is odd(11~1001) adjustable; high - cut and low - cut can be adjustable；
//Display : dynamtically display the waves;
		// turn off all stimulators 
	SetStimulatorOff(0);
	SetStimulatorOff(1);
	SetStimulatorOff(2);
	SetStimulatorOff(3);
	SetStimulatorOff(4);
	SetStimulatorOff(5);
	int i, j, m, k;
	// the controller board should also be updated to show that the stimulators are off.
	uint8_t p_control_parameters[NUM_CTRL_LED_PARAM];
	for (j = 0; j < NUM_CTRL_LED_PARAM; ++j) {
		p_control_parameters[j] = 0;
	}
	SetControlStimulatorParameters(p_control_parameters);

	// Record 16 channels 
	int num_acq_channels_recorded = 15;
	int num_stim_channels_recorded = 0;
	unsigned int channels_to_record = 0x0FFFE;// we are recording 16 channels 
	unsigned int stim_channels_to_record = 0x00; // we are not recording any stimualtion channels 
	
	emgdata.m_sampling_speed = 16000; // 16k sampling

	// setup the notch filter
	bool use_usb_notch_filter = false; // we cannot use this with the calibration signal on 2 channels
	int time_offset[16];

	float lowcutoff = 10;
	float highcutoff = 1000;
	int filter_order = 1601;

	for (k = 0; k < 16; ++k) time_offset[k] = 0;
	if (use_online_filters) {
		if (!use_usb_notch_filter) {
			OnlineNotchFilterSetup(emgdata.m_sampling_speed, 50, true); // if using the LOTUS_NOTCH_FILTER_50 do not use the OnlineNotchFilterSetup
			// add the time offset for the notch filter 
		}

		// to bandpass filter the data between 10hz and 1kHz
		//OnlinebandPassFilterSetup will change the filter order to match the High and Low Pass freq requested. 
		// for eg. if you set the filter order to 1001, and ask for a 10Hz HP, the filter order will be changed to 1601 since that is the min order required
		// if you set it to 2001 then the order will not be changed since it is greater than 1601

		for (k = 0; k < 16; ++k) {
			time_offset[k] += OnlinebandPassFilterSetup(k, emgdata.m_sampling_speed, lowcutoff, highcutoff, filter_order); // setting all filters the same.
		}
		if (filter_order > 2001) {
			// the required filter order was greater than 2001 and was modified.
		//	If the filter order is used in any calculations, update those calculations as needed
		}
	}

	// NOTE : if we are filtering the data at 1khz there is not much point in sampling the data at 16Khz 

	if (use_usb_notch_filter) {
		SetAcquisitionParameters(emgdata.m_sampling_speed, LOTUS_NOTCH_FILTER_50, channels_to_record, stim_channels_to_record);
	}
	else {
		FlipLoff(); // toggle the LOFF current direction to discharge the capacitance. 
		if (cal_ch1) {
			SetAcquisitionParameters(emgdata.m_sampling_speed, LOTUS_CAL_SIGNAL_CH1, channels_to_record, stim_channels_to_record);
		}
		else {
			SetAcquisitionParameters(emgdata.m_sampling_speed, LOTUS_RECORD_WITH_LEADOFF, channels_to_record, stim_channels_to_record);

		}
	}
	int iRet;
	int total_num_samples_per_ch_recorded = 0;
	int total_number_of_samples_per_ch_to_record = 960000;
	int num_acq_samples_per_ch = 0;
	//clear the data buffers
	for (k = 0; k < num_acq_channels_recorded; ++k) {
		emgdata.m_data[k].clear();
		emgdata.m_filtered_data[k].clear();
		emgdata.m_data[k].reserve(total_number_of_samples_per_ch_to_record + 100);
		emgdata.m_filtered_data[k].reserve(total_number_of_samples_per_ch_to_record + 100);
	}
	
	Sleep(10);
	StartAcq();
	Sleep(10);
	uint32_t ch_leadoff;
	while (total_num_samples_per_ch_recorded < total_number_of_samples_per_ch_to_record) { //continue recording data until stop condition is reached.
		// if you want to user control, run this function in a separate thread and have it look for a varioble such as m_is_recording_data
		// then check the value of m_is_recording_data  in the above while loop.
		iRet = ReadDataFromAcquisitionDevice(num_acq_samples_per_ch, emgdata.acqdata, emgdata.number_of_samples_to_read_per_loop); //25600  num_acq_samples_per_ch=160
		// check if leadoff was ok
		GetAcquisitionChannelLeadOff(ch_leadoff);  //  0 - 15 ch 
		if (ch_leadoff > 0) { // display leadoff if channel disconnected.
			wcout << "LeadOff Channels : ";
			wcout << ch_leadoff; wcout << "\n";
		}
		// check ch_leadoff to make sure leads were connected.
		if (iRet == 0) { // acquired data copied to data.  pts_recorded = num_samples_per_ch*num_channels_recorded;
			m = 0;
			for (j = 0; j < num_acq_samples_per_ch; ++j) {
				for (k = 0; k < num_acq_channels_recorded; ++k) {
					float data = emgdata.acqdata[m++];
					emgdata.m_data[k].push_back(data);
					if (use_online_filters) {
						float filtered_data = OnlineBandpassFilterData(k, data);
						if (!use_usb_notch_filter) {
							filtered_data = OnlineNotchFilterData(k, filtered_data);
						}
						//						if (total_num_samples_per_ch_recorded > time_offset[k]) { // offset the filtered data by the bandpass filter time offset. 
						emgdata.m_filtered_data[k].push_back(filtered_data);
						//						}
					}
				}
				++total_num_samples_per_ch_recorded;
			}
			wcout << "Total of samples: ";  wcout << total_num_samples_per_ch_recorded;  
			wcout << ", Sample per Loop:  ";  wcout << num_acq_samples_per_ch;
			wcout << "\n";
			//			wcout << std::flush;
		}
		if (iRet > 0) { // Error condition 
			StopAcq();// Stop Acquisition
			// Show the error message if needed
			return iRet;
		}
		if (communication_error) {
			communication_error = false;
			Sleep(100);
			// restart the acquisition 
			int conn_status = IsAcqAndStimConnected();
			conn_status = conn_status & 0x01;
			while (conn_status) {
				wcout << "Waiting for Acquisition board\n ";
				Sleep(100);
				conn_status = IsAcqAndStimConnected() & 0x01;
			}
			StartAcq();
		}
	}
	StopAcq();// Stop Acquisition
	if (!use_online_filters) {// filter the data after acquisition.
		for (i = 0; i < num_acq_channels_recorded; ++i) {
			FilterData(&(emgdata.m_data[i][0]), &(emgdata.m_filtered_data[i][0]), total_num_samples_per_ch_recorded, emgdata.m_sampling_speed, lowcutoff, highcutoff, filter_order);
		}
	}
	PrintDatatoFile(acqfout, emgdata.m_filtered_data, total_num_samples_per_ch_recorded, num_acq_channels_recorded + stim_channels_to_record);
//	PrintDatatoFile(acqfout, emgdata.m_data, total_num_samples_per_ch_recorded, num_acq_channels_recorded + stim_channels_to_record);

	wchar_t ch;
	wcout << "Continue?\n ";
	wcin >> ch;
	return iRet;
}

int FindAcquisitionChannelImpedance(FILE* acqfout, CData& emgdata) {

	// read the channel impedance offset values
	float ch_impedance_offset[16];
	GetAcquisitionImpedanceOffset(ch_impedance_offset);
	wcout << "Channel Impedance Offsets \n";
	int num_acq_channels_recorded = 16;
	int num_stim_channels_recorded = 0;
	unsigned int channels_to_record = 0x0FFFF;// we are recording 16 channels
	unsigned int stim_channels_to_record = 0x00; // we are not recording any stimualtion channels 
	std::vector<ChannelTestValues> m_ch_stat;
	m_ch_stat.resize(16);
	std::vector<ChannelImpedance> m_acquisition_impedance;
	m_acquisition_impedance.resize(16);
	int n = 0;
	bool show_leadoff = false;
	double max;
	double min;
	double mean;
	int total_number_of_samples_to_record = 2000; // record 2 seconds of data.
	int num_samples_to_report = total_number_of_samples_to_record - 100; // drop the last 100 samples
	int samples_start_index = 500;
	int num_samples_to_analyze = 936;
	int samples_stop_index = 500 + num_samples_to_analyze;
	bool reference_connected = true;
	int total_num_samples_per_ch_recorded = 0;
	int iRet;
	double ch_baseline_val[16];
	double ch_abs_mean_val[16];
	double ch_abs_max_val[16];

	LOFF_CURRENT_Values current = LOTUS_LOFF_CURRENT_6nA;
	AcqSetLoffCurrentFreq(LOTUS_LOFF_DC_FREQ, current);
	float Impedance_gain = 1;
	switch (current) {
	case LOTUS_LOFF_CURRENT_6nA://!< Lead Off 6nA excitation current
		Impedance_gain = 4.0 ;
		break;
	case LOTUS_LOFF_CURRENT_24nA://!<Lead Off 24nA excitation current
		Impedance_gain = 1.0;
		break;
	case LOTUS_LOFF_CURRENT_6uA://!<Lead Off 6uA excitation current
		Impedance_gain = 1.0/250;
		break;
	case LOTUS_LOFF_CURRENT_24uA://!<Lead Off 24uA excitation current
		Impedance_gain = 1.0/1000;
		break;
	}
	int i = 0;
	int j = 0;
	SetAcquisitionParameters(1000, LOTUS_LEADOFFPOS, channels_to_record, stim_channels_to_record);
	wcout << "Channel Test with LeadOffPos\n";
	if (acqfout)fprintf(acqfout, "Channel Test with LeadOffPos\n");
	total_num_samples_per_ch_recorded = 0;
	iRet = RecordData(emgdata, total_num_samples_per_ch_recorded, total_number_of_samples_to_record, num_acq_channels_recorded, 0, show_leadoff);
	if (iRet > 0) { // Error condition 
		ShowErrorMesage(iRet);
	}
	else {
		for ( i = 0; i < num_acq_channels_recorded; ++i) {
			ch_baseline_val[i] = 0;
			ch_abs_mean_val[i] = 0;
			ch_abs_max_val[i] = 0;
			emgdata.m_filtered_data[i].clear();
			for ( j = samples_start_index; j < samples_stop_index; ++j) {
				// find mean for eah channel;
				ch_baseline_val[i] + emgdata.m_data[i].at(j);
			}
			ch_baseline_val[i] = ch_baseline_val[i] / (samples_stop_index - samples_start_index);
			// subtract mean and find abs value;
			for ( j = samples_start_index; j < samples_stop_index; ++j) {
				emgdata.m_filtered_data[i].push_back( emgdata.m_data[i].at(j) - ch_baseline_val[i]);
			}
			for ( j = samples_start_index; j < samples_stop_index; ++j) {
				emgdata.m_filtered_data[i].push_back(emgdata.m_data[i].at(j) - ch_baseline_val[i]);
			}
			for (j = 0; j < num_samples_to_analyze; ++j) {
				double abs_val = emgdata.m_filtered_data[i].at(j);
				ch_abs_mean_val[i] += abs_val;
				if (ch_abs_max_val[i] < abs_val)ch_abs_max_val[i] = abs_val;
			}
		}
		GetStats(m_ch_stat, emgdata.m_filtered_data, num_acq_channels_recorded, num_samples_to_analyze, 0, true, acqfout);
		PrintDatatoFile(acqfout, emgdata.m_filtered_data, num_samples_to_analyze, num_acq_channels_recorded);
	}


	SetAcquisitionParameters(1000, LOTUS_LEADOFFNEG, channels_to_record, stim_channels_to_record);
	wcout << "Channel Test with LeadOffNeg\n";
	if (acqfout)fprintf(acqfout, "Channel Test with LeadOffNeg\n");
	total_num_samples_per_ch_recorded = 0;
	 iRet = RecordData(emgdata, total_num_samples_per_ch_recorded, total_number_of_samples_to_record, num_acq_channels_recorded, 0, show_leadoff);
	if (iRet > 0) { // Error condition 
		ShowErrorMesage(iRet);
	}
	else {
//		GetStats(m_ch_stat, emgdata.m_data, num_acq_channels_recorded, num_samples_to_report, 1000, true, acqfout);
		PrintDatatoFile(acqfout, emgdata.m_data, num_samples_to_report, num_acq_channels_recorded);

	}

	return SUCCESS;
}

int main()
{
	bool control_device_ok = false;
	bool acq_device_ok = false;
	bool stim_device_ok = false;
	int speed = 1000;


	// Set Demo mode if required
	wchar_t ch = 'n';
//	wcout << "Start Demo Mode?: \n";
//	wcin >> ch;
	bool demo_mode = false;
	bool test_notch_filter = false;
	if (ch == 'Y' || ch == 'y') {
		demo_mode = true;	
		wcout << "Test Notch Filter?: \n";
		wcin >> ch;
		if (ch == 'Y' || ch == 'y') {
			wcout << "Wait: Demo mode , Loading Data in notch_test.txt\n";
			SetDemoMode(L"playbackfiles/notch_test.txt", NULL);
			test_notch_filter = true;
			speed = 10000;
		}
		else {
			wcout << "Wait:  Demo Mode Loading Data in 2ch_emg_1k.csv\n";
			SetDemoMode(L"playbackfiles/2ch_emg_1k.csv", NULL);
		}
		
		// if you have the data as an array already loaded you can use the following function
		//SetDemoModeWithArray(const float* acq_data, int num_acq_point, const float* stim_data, int num_stim_points);
	}

	wcout << "Opening Devices\n";
	OpenDevices(L"LotusLog.txt");//	Open Devices calls FindHardware() to find all the devices.
	// Get info about the acquisition device
	wchar_t acq_manufacturer_id_buffer[100] = L"None";
	wchar_t acq_sn_buffer[100] = L"None";
	

	int iRet = FindAcquisitionHardware(acq_manufacturer_id_buffer, acq_sn_buffer, 100);

	uint8_t acq_device_info[8];
	if (iRet != SUCCESS)ShowErrorMesage(iRet);
	else {
		int iRet = ReportAcquisitionSelfCheck(acq_device_info);//self check see LotusTechnicalRequirementsDocument.ods for info
		wcout << " Acquisition Self Check";
		for (int i = 0; i < 8; ++i) {
			wcout << std::hex << acq_device_info[i]; wcout << ", ";
		}
		wcout << std::dec;
		wcout << "\n";
		// description of the self check feedback 

		if( acq_device_info[0] & 0x01 )  wcout << "ACQ USB bod_stat Set\n";
		if (acq_device_info[0] & 0x02)  wcout << "ACQ USB por_reset Set\n";
		if (acq_device_info[0] & 0x04)  wcout << "ACQ USB pin_reset Set\n";
		if (acq_device_info[0] & 0x08)  wcout << "ACQ USB bod_reset Set\n";
		if (acq_device_info[0] & 0x10)  wcout << "ACQ USB wdt_reset Set\n";
		if (acq_device_info[0] & 0x20)  wcout << "ACQ USB USB_SettledAt480Mb   Set\n";
		if (acq_device_info[0] & 0x40)  wcout << "ACQ USB UART_Active  Set\n";
		if (acq_device_info[0] & 0x80)  wcout << "ACQ USB UART_DetectedOnPowerUp \n";

		if (acq_device_info[2] & 0x02)  wcout << "ACQ USB Led Power Set\n";
		if (acq_device_info[2] & 0x04)  wcout << "ACQ USB Sync Clock Set\n";
		if (acq_device_info[2] & 0x08)  wcout << "ACQ USB Sync Mode Set\n";

		if (acq_device_info[4] & 0x01)  wcout << "ACQ ISO bod_stat Set\n";
		if (acq_device_info[4] & 0x02)  wcout << "ACQ ISO por_reset Set\n";
		if (acq_device_info[4] & 0x04)  wcout << "ACQ ISO pin_reset Set\n";
		if (acq_device_info[4] & 0x08)  wcout << "ACQ ISO wdt_reset Set\n";
		if (acq_device_info[4] & 0x20)  wcout << "ACQ ISO Vcc Fault  Set\n";
		if (acq_device_info[4] & 0x40)  wcout << "ACQ ISO Vdd Fault  Set\n";

		if (acq_device_info[6] & 0x01)  wcout << "ACQ ISO Led Power Set\n";
		if (acq_device_info[6] & 0x02)  wcout << "ACQ ISO AFE 1 powerdown Set\n";
		if (acq_device_info[6] & 0x04)  wcout << "ACQ ISO AFE 2 powerdown Set\n";
		if (acq_device_info[6] & 0x20)  wcout << "ACQ ISO AFE nReset State\n";
		if (acq_device_info[6] & 0x40)  wcout << "ACQ ISO +5.3 Vcc nShutdown Set\n";

		wcout << "Acquisition Hardware: ";
		wcout << acq_manufacturer_id_buffer;
		wcout << "   Serial Number:";
		wcout << acq_sn_buffer;
		wcout << "\n";
		acq_device_ok = true;

		// set acquisition impedance offset 
		ch = 'y';
		wcout << "Set the Acquisition Impedance Offset? ";
		wcin >> ch;
		if (ch == 'y') {
			ch = 'n';
			wcout << "Short all the inputs to the Acquisition board and press Y ";
			wcin >> ch;
			if (ch == 'y') {
				wcout << "Performing Impedance offset measurements... ";
				SetAcquisitionImpedanceOffset();
			}
		}
	}

	// find the stimulation device
	wchar_t stim_manufacturer_id_buffer[100] = L"None";
	wchar_t stim_sn_buffer[100] = L"None";

	uint8_t stim_device_info[8];
	uint8_t ch_id_info[512];
	int i = 0;
	iRet = FindStimulationHardware(stim_manufacturer_id_buffer, stim_sn_buffer, 100);
	if (iRet != SUCCESS)ShowErrorMesage(iRet);
	else {
		iRet = ReportStimulatorSelfCheck(stim_device_info);
		bool isolated_section_is_powered = false;
		for ( i = 4; i < 8; ++i) {
			if (stim_device_info[i] > 0) isolated_section_is_powered = true;
		}
		if (!isolated_section_is_powered) {
			Sleep(200); // wait 200msec and try again;
			iRet = ReportStimulatorSelfCheck(stim_device_info);
		}
			//self check failed
		wcout << "Stimulator Self Check: ";
		wcout << std::hex;
		for ( i = 0; i < 8; ++i) {
			wcout <<  stim_device_info[i]; wcout << ", ";
		}
		wcout << std::dec;
		wcout << "\n";

		if (stim_device_info[0] & 0x01)  wcout << "STIM USB bod_stat Set\n";
		if (stim_device_info[0] & 0x02)  wcout << "STIM USB por_reset Set\n";
		if (stim_device_info[0] & 0x04)  wcout << "STIM USB pin_reset Set\n";
		if (stim_device_info[0] & 0x08)  wcout << "STIM USB bod_reset Set\n";
		if (stim_device_info[0] & 0x10)  wcout << "STIM USB wdt_reset Set\n";
		if (stim_device_info[0] & 0x20)  wcout << "STIM USB USB_SettledAt480Mb   Set\n";
		if (stim_device_info[0] & 0x40)  wcout << "STIM USB UART_Active  Set\n";
		if (stim_device_info[0] & 0x80)  wcout << "STIM USB UART_DetectedOnPowerUp \n";

		if (stim_device_info[2] & 0x01)  wcout << "STIM USB DCDC_Enable Set\n";
		if (stim_device_info[2] & 0x02)  wcout << "STIM USB Led Power Set\n";
		if (stim_device_info[2] & 0x04)  wcout << "STIM USB Sync Clock Set\n";
		if (stim_device_info[2] & 0x08)  wcout << "STIM USB Sync Mode Set\n";

		if (stim_device_info[4] & 0x01)  wcout << "STIM ISO bod_stat Set\n";
		if (stim_device_info[4] & 0x02)  wcout << "STIM ISO por_reset Set\n";
		if (stim_device_info[4] & 0x04)  wcout << "STIM ISO pin_reset Set\n";
		if (stim_device_info[4] & 0x08)  wcout << "STIM ISO bod_reset Set\n";
		if (stim_device_info[4] & 0x10)  wcout << "STIM ISO wdt_reset  Set\n";
		if (stim_device_info[4] & 0x20)  wcout << "STIM ISO Vcc Fault  Set\n";
		if (stim_device_info[4] & 0x40)  wcout << "STIM ISO Vdd Fault  Set\n";
		if (stim_device_info[4] & 0x80)  wcout << "STIM ISO Avdd Fault  Set\n";


		if (stim_device_info[5] & 0x01)  wcout << "STIM ISO HVS fault Set\n";
		if (stim_device_info[5] & 0x02)  wcout << "STIM ISO 1 wire fault Set\n";
		if (stim_device_info[5] & 0x04)  wcout << "STIM ISO I2C init failed Set\n";
		if (stim_device_info[5] & 0x08)  wcout << "STIM ISO I2C ExcessNack  Faled\n";
		if (stim_device_info[5] & 0x10)  wcout << "STIM ISO I2C Timeout Failed \n";
		if (stim_device_info[5] & 0x20)  wcout << "STIM ISO I2C Undefined  Failed\n";
		if (stim_device_info[5] & 0x40)  wcout << "STIM ISO DS2484  Failed\n";


		if (stim_device_info[6] & 0x01)  wcout << "STIM ISO Led Power\n";
		if (stim_device_info[6] & 0x02)  wcout << "STIM ISO Sync Clock Set\n";
		if (stim_device_info[6] & 0x04)  wcout << "STIM ISO HVS En\n";
		if (stim_device_info[6] & 0x20)  wcout << "STIM ISO Unused\n";
		if (stim_device_info[6] & 0x40)  wcout << "STIM ISO I2C Address Ack Set\n";

		wcout << "Stimulation Hardware: ";
		wcout << stim_manufacturer_id_buffer;
		wcout << "   Serial Number:";
		wcout << stim_sn_buffer;
		wcout << "\n";
		stim_device_ok = true;
		int n = 0;
		int waittime = 1;
		wcout << "Time "; wcout << (waittime); wcout << "\n";
		GetStimChannelInfo(ch_id_info, waittime); // this waittime is for debugging only, and will be removed in final API. 

		// Display the Channel data
		for (n = 0; n < 10; ++n) {
			wcout << "Channel "; wcout << (n + 1); wcout << "Lotus Serial ID : ";
			for (int m = 0; m < 17; ++m) {
				ch = ch_id_info[n * 25 + m];
				wcout << ch;
			}
			wcout << "\n";
		}

		// End display channel data
		wcout << "Calibrate Stim Impedance? ";
		wcin >> ch;
		if (ch == 'y') {
			wcout << "Unplug any connectors from the stimulation channels and press y to set the stimulator offset";
			wcin >> ch;
			if (ch == 'y') {
				SetStimulatorImpedanceOffset();
			}
		}
		wcout << "Calibrate Stim Pulse Amplitude Offset? ";
		wcin >> ch;
		if (ch == 'y') {
			int ch1_offset, ch2_offset;
			wcout << "Connect a 10 kOhm resistor across Channel 1. Connect a Voltmeter across the Load  and press y to calculate the offset";
			wcin >> ch;
			if (ch == 'y') {
				StartStimIdleAmpOffsetMeasurement(0);
				wcout << "Enter the voltmeter value in mV";
				wcin >> ch1_offset;

				wcout << "Connect a 10 kOhm resistor across Channel 2.  Connect a Voltmeter across the Load and press y to calculate the offset";
				wcin >> ch;
				if (ch == 'y') {
					StartStimIdleAmpOffsetMeasurement(1);
					wcout << "Enter the voltmeter value in mV";
					wcin >> ch2_offset;
					SaveStimIdleAmpOffsetToFlash(ch1_offset, ch2_offset);
				}
			}
		}
		
		// The GetStimImpedanceValues is used by the API 
		const int16_t* stim_offset_values = GetStimImpedanceOffsetValues();
		wcout << "Stim Impedance Offset: ";
		// the API  uses the Stim Impedance offset values in the calculation of the impedance. 
		// The application does not need to do anything with this. 
		for (n = 0; n < 4; ++n) {
			wcout << stim_offset_values[n]; wcout << ", ";
		}
		wcout << "\n";
		wcout << "Stim Channel Offset: ";
		// the channel offset is used by the API to adjust the pulse amplitude.
		for (n = 0; n < 2; ++n) {
			wcout << stim_offset_values[n + 12]; wcout << ", ";
		}
		wcout << "\n";



		ch = 'y';
		wcout << "Test Stim Impedance? ";
		wcin >> ch;
		while(ch == 'y'){
			for (n = 0; n < 2; ++n) {
				int stim_ch_impedance = 0;
				GetStimChannelImpedance(n, stim_ch_impedance); // get the stimulator channel impedance for channel 1.
				wcout << "Channel "; wcout << (n + 1); wcout << "Impedance : ";
				wcout << stim_ch_impedance;
				wcout << "\n";
			}
			wcout << "ReTest Stim Impedance? ";
			wcin >> ch;
		}

	
	}
	// find the control device
	wchar_t ctrl_manufacturer_id_buffer[100] = L"None";
	wchar_t ctrl_sn_buffer[100] = L"None";

	uint8_t ctrl_device_info[8];
	iRet = FindControlHardware(ctrl_manufacturer_id_buffer, ctrl_sn_buffer, 100);
	if (iRet != SUCCESS)ShowErrorMesage(iRet);
	else {
		iRet = ReportControlSelfCheck(ctrl_device_info);
			//self check failed
			wcout << " Control Self Check";
			for (int i = 0; i < 8; ++i) {
				wcout << std::hex << ctrl_device_info[i]; wcout << ", ";
			}
			wcout << std::dec;
			wcout << "\n";

			if (ctrl_device_info[0] & 0x01)  wcout << "CTRL USB bod_stat Set\n";
			if (ctrl_device_info[0] & 0x02)  wcout << "CTRL USB por_reset Set\n";
			if (ctrl_device_info[0] & 0x04)  wcout << "CTRL USB pin_reset Set\n";
			if (ctrl_device_info[0] & 0x08)  wcout << "CTRL USB bod_reset Set\n";
			if (ctrl_device_info[0] & 0x10)  wcout << "CTRL USB wdt_reset Set\n";
			if (ctrl_device_info[0] & 0x20)  wcout << "CTRL USB USB_SettledAt480Mb   Set\n";

			wcout << "Control Hardware: ";
			wcout << ctrl_manufacturer_id_buffer;
			wcout << "   Serial Number:";
			wcout << ctrl_sn_buffer;
			wcout << "\n";
			control_device_ok = true;
			SetControlParamChangeCallback(&ControlParamChangeNotification);  // This only needs to be called once after the control device has been opened
			// The device class will use the call back function to notify the main program of any changes in the controller board.
			// if the control device is closed then the call back function is set to NULL.  
			//The Application will then need to reopen the control device and call this function again
		
	}
	// After all the hardware has been detected
	SetLotusErrorCallback(ErrorNotification);
	if(control_device_ok) TestControlBoard();

	std::vector<ChannelTestValues> m_ch_stat;
	m_ch_stat.resize(16);
	

//	if (control_device_ok && acq_device_ok && stim_device_ok) {
	if (acq_device_ok ) {
		// setup output files to save the data.
		FILE* acqfout = NULL;
		int err = _tfopen_s(&acqfout, L"output_acq_file.csv", L"w");
		if (err) {
			printf("Error Opening output_acq_file.csv");
		}

		CData  emgdata;
		emgdata.AllocateMemory(10000);


		TestFilterTimeOffset(emgdata);
		// Unlimited data can be recorded from thehardware, 
		//each time the hardware is asked for data it will send the acquired data
		int total_num_samples_per_ch_recorded = 0; // variable to hold the total number of samples per channel that have been recorded.
		int total_number_of_samples_to_record = 10000; // How much data we want to record.
		int number_of_samples_to_read_per_loop = 10000; // this is the max number of samples that we can read per each call to ReadDataFromAcquisitionDevice
		// we want to read 1000 samples for the 16 channels;
		
		int num_acq_channels_recorded = 16;
		int num_stim_channels_recorded = 0;
		unsigned int channels_to_record = 0x0FFFF;// we are recording 16 channels
		unsigned int stim_channels_to_record = 0x00; // we are not recording any stimualtion channels 



		int num_acq_samples_per_ch = 0; // this will hold the number of acquisition samples per channel that have been recorded since the last time ReadDataFromAcquisitionDevice is called
		int num_stim_samples_per_ch = 0;// this will hold the number of stimulation samples per channel that have been recorded since the last time ReadDataFromStimulationDevice is called

		// turn off all stimulators 
		SetStimulatorOff(0);
		SetStimulatorOff(1);
		SetStimulatorOff(2);
		SetStimulatorOff(3);
		SetStimulatorOff(4);
		SetStimulatorOff(5);
		bool noerror = true;
		int i, k, j, m;
		double max, min, mean;
		
		emgdata.m_sampling_speed = speed;

		if (demo_mode) {
			wcout << "Record Demo EMG: Press 1 \n 50Hz Notch Filter: Press 2 \n Record Filtered EMG 100Hz to 10KHz: Press 3:\n Record Filtered EMG 10Hz to 1KHz: Press 4:\n  ";
			wcin >> ch;
			
			num_acq_channels_recorded = 2;// 2ch_emg_1k  data file only has 2 channels of data
			unsigned int channels_to_record = 0x00003;// we are recording 2 channels
			int buffer_time = 5; // in sec , time of the sweep, for scope mode or the time that the software may be away for chartmode 
			int numberofbytes_in_buffer = num_acq_channels_recorded * 4 * emgdata.m_sampling_speed * buffer_time;
			SetAcquisitionParameters(1000, LOTUS_RECORD, channels_to_record, stim_channels_to_record, numberofbytes_in_buffer);
			int filter_order = 1001;
			bool notchfilter = false;
			bool bandpass = false;
			if (ch == '2') {
				// setup the notch filter
				OnlineNotchFilterSetup(emgdata.m_sampling_speed, 50, true);
				notchfilter = true;
			}
			else if (ch == '3') {
				bandpass = true;
				OnlinebandPassFilterSetup(-1, emgdata.m_sampling_speed, 100, 10000, filter_order); // setting all filters the same.
			}
			else if (ch == '4') {
				bandpass = true;
				OnlinebandPassFilterSetup(-1, emgdata.m_sampling_speed, 10, 1000, filter_order); // setting all filters the same.
			}
			else {
			}
			wcout << "Recording Demo Data \n";
			if (acqfout)fprintf(acqfout, "Recording Demo Data\n");
			StartAcq();
			for (k = 0; k < num_acq_channels_recorded; ++k) {
				emgdata.m_data[k].clear();
			}
			int total_number_of_samples_per_ch_to_record = 400000;

			while (total_num_samples_per_ch_recorded < total_number_of_samples_per_ch_to_record) { //continue recording data until stop condition is reached.
	// if you want to user control, run this function in a separate thread and have it look for a variable such as m_is_recording_data
	// then check the value of m_is_recording_data  in the above while loop.
				iRet = ReadDataFromAcquisitionDevice(num_acq_samples_per_ch, emgdata.acqdata, emgdata.number_of_samples_to_read_per_loop);
				// check if leadoff was ok
//				GetAcquisitionChannelLeadOff(ch_leadoff);  //  0 - 15 ch positive and 16 - 31 negative 
				// check ch_leadoff to make sure leads were connected.
				if (iRet == 0) { // acquired data copied to data.  pts_recorded = num_samples_per_ch*num_channels_recorded;
					m = 0;
					for (j = 0; j < num_acq_samples_per_ch; ++j) {
						for (k = 0; k < num_acq_channels_recorded; ++k) {
							float data = emgdata.acqdata[m++];
							if (bandpass) {
								data = OnlineBandpassFilterData(k, data);
							}
							if (notchfilter) { // do an online notch filter on the data.
								data = OnlineNotchFilterData(k, data);
							}
							emgdata.m_data[k].push_back(data);
						}
					}
					total_num_samples_per_ch_recorded += num_acq_samples_per_ch;
				}
				if (iRet > 0) { // Error condition 
					StopAcq();// Stop Acquisition
					// Show the error message if needed
					ShowErrorMesage(iRet);
					return iRet;
				}
			}
			StopAcq();// Stop Acquisition			if (iRet > 0) { // Error condition 
			PrintDatatoFile(acqfout, emgdata.m_data, total_num_samples_per_ch_recorded, num_acq_channels_recorded);
		}
		else {
			ch = 'n';
			wcout << "Find Acquisition Channel Impedance?";
			wcin >> ch;
			while (ch == 'y') {
				FindAcquisitionChannelImpedance(acqfout, emgdata);
				wcout << "ReTest Acquisition Impedance? ";
				wcin >> ch;
			}

			// Check input shorted 
			ch = '1';
			wcout << "Data in micro Volts\n";
			if (acqfout)fprintf(acqfout, "Channel Test with Inputs Shorted\n");
			while ((ch >= '1' && ch <= '9') || (ch >= 'a' && ch <= 'd')) {
				wcout << "Perform Input Shorted Test: Press 1 \n Channel Test with Cal Signal: Press 2:\n Channel Test with LeadOff:Press 3:\nChannel Test with LOTUS_LEADOFFPOS:Press 4\n Channel Test with LOTUS_LEADOFFNEG: Press 5\n Record SEMG: Press 6\n Record SEMG with LeadOff: Press 7\n Record LEADOFF with Gain x4: Press 8\nRecord LEADOFFPOS with Gain x4: Press 9\nRecord LEADOFFNEG with Gain x4: Press a\nTo Quit: Press any key\n";
				wcin >> ch;
				bool show_leadoff = false;
				if (ch == '1') {
					SetAcquisitionParameters(1000, LOTUS_INPUT_SHORTED, channels_to_record, stim_channels_to_record);
					wcout << "Channel Test with Inputs Shorted: Data in micro Volts\n\n";
					if (acqfout)fprintf(acqfout, "Channel Test with Inputs Shorted: Data in micro Volts\n");
				}
				if (ch == '2') {
					SetAcquisitionParameters(1000, LOTUS_CAL_SIGNAL, channels_to_record, stim_channels_to_record);
					wcout << "Channel Test with Cal Signal\n The Max-Min should be about 3750: Data in micro Volts\n";
					if (acqfout)fprintf(acqfout, "Channel Test with Cal Signal\n The Max-Min should be about 3750: Data in micro Volts\n");
				}
				if (ch == '3') {
					SetAcquisitionParameters(1000, LOTUS_LEADOFF, channels_to_record, stim_channels_to_record);
					wcout << "Channel Test with LeadOff\n";
					if (acqfout)fprintf(acqfout, "Channel Test with LeadOff\n");
				}
				if (ch == '4') {
					AcqSetLoffCurrentFreq(LOTUS_LOFF_312_FREQ, LOTUS_LOFF_CURRENT_24nA);
					SetAcquisitionParameters(1000, LOTUS_GAIN4_LEADOFFPOS, channels_to_record, stim_channels_to_record);
					wcout << "Channel Test with LOTUS_GAIN4_LEADOFFPOS\n";
					if (acqfout)fprintf(acqfout, "Channel Test with LOTUS_GAIN4_LEADOFFPOS\n");
				}
				if (ch == '5') {
					AcqSetLoffCurrentFreq(LOTUS_LOFF_312_FREQ, LOTUS_LOFF_CURRENT_24nA);
					SetAcquisitionParameters(1000, LOTUS_GAIN4_LEADOFFNEG, channels_to_record, stim_channels_to_record);
					wcout << "Channel Test with LOTUS_GAIN4_LEADOFFNEG\n";
					if (acqfout)fprintf(acqfout, "Channel Test with LOTUS_GAIN4_LEADOFFNEG\n");
				}
				if (ch == '6') {
					SetAcquisitionParameters(1000, LOTUS_RECORD, channels_to_record, stim_channels_to_record);
					wcout << "Record EMG\n";
					if (acqfout)fprintf(acqfout, "Record EMG\n");
				}
				if (ch == '7') {
					SetAcquisitionParameters(1000, LOTUS_RECORD_WITH_LEADOFF, channels_to_record, stim_channels_to_record);
					wcout << "Record  EMG with LeadOff\n";
					show_leadoff = true;
					if (acqfout)fprintf(acqfout, "Record EMG with LeadOff\n");
				}
				if (ch == '8') {
					SetAcquisitionParameters(1000, LOTUS_GAIN4_WITH_LEADOFF, channels_to_record, stim_channels_to_record);
					wcout << "Record EMG with Gain x4  \n";
					if (acqfout)fprintf(acqfout, "Record EMG with Gain x4 \n");
				}

				if (ch == '9') {
					SetAcquisitionParameters(1000, LOTUS_GAIN4_LEADOFFPOS, channels_to_record, stim_channels_to_record);
					wcout << "Record EMG with Gain x4  \n";
					if (acqfout)fprintf(acqfout, "Record EMG with Gain x4 \n");
				}

				if (ch == 'a') {
					SetAcquisitionParameters(1000, LOTUS_GAIN4_LEADOFFNEG, channels_to_record, stim_channels_to_record);
					wcout << "Record EMG with Gain x4  \n";
					if (acqfout)fprintf(acqfout, "Record EMG with Gain x4 \n");
				}

				total_number_of_samples_to_record = 2000; // record 2 seconds of data.
				if ((ch >= '1' && ch <= '9') || (ch >= 'a' && ch <= 'd')) {
					iRet = RecordData(emgdata, total_num_samples_per_ch_recorded, total_number_of_samples_to_record, num_acq_channels_recorded, 0, show_leadoff);
					if (iRet > 0) { // Error condition 
						ShowErrorMesage(iRet);
					}
					else {
						GetStats(m_ch_stat, emgdata.m_data, num_acq_channels_recorded, total_num_samples_per_ch_recorded, 500, true, acqfout);
						PrintDatatoFile(acqfout, emgdata.m_data, total_num_samples_per_ch_recorded, num_acq_channels_recorded);
					}
				}
			}
		}// if(demo mode)
	/*
		*/
		// Record 16 channels of SEMG at 16k
		wcout << "Record 15 channels of SEMG at 16k with 2 calibration channels without filters: ";
		wcin >> ch;
		if (ch == 'y' || ch == 'Y') {
//			Record_SEMG(emgdata, acqfout, true, false); // use usb notch filter 
			Record_SEMG(emgdata, acqfout, false, true); // use 2 channel cal signal with no online filter
		}
		// Record 15 channels of SEMG at 16k
		wcout << "Record 15 channels of SEMG at 16k with 2 calibration channels with filters: ";
		wcin >> ch;
		if (ch == 'y' || ch == 'Y') {
			//	Record_SEMG(emgdata, acqfout, true, false); // use usb notch filter 
			Record_SEMG(emgdata, acqfout, true, true); // use 2 channel cal signal with no online filter
		}
		// Record 15 channels of SEMG at 16k
		wcout << "Record 15 channels of SEMG at 16k without calibration channels with filters: ";
		wcin >> ch;
		if (ch == 'y' || ch == 'Y') {
			Record_SEMG(emgdata, acqfout, true, false); // use 2 channel cal signal with no online filter
		}
		// Record 15 channels of SEMG at 16k
		wcout << "Record 15 channels of SEMG at 16k without calibration channels without filters: ";
		wcin >> ch;
		if (ch == 'y' || ch == 'Y') {
			Record_SEMG(emgdata, acqfout, false, false); // use 2 channel cal signal with no online filter
		}

		if (!stim_device_ok) { // if no stimulation device stop here
			if (acqfout)fclose(acqfout);
			CloseDevices(); // to delete the devices created and clear memory.
			emgdata.FreeMemory();
			return 0;		// there is no stimulator board so we are stopping here.
		}
		// Record TEMG
		wcout << "Perform Leadoff Timing Test: ";
		wcin >> ch;
		if (ch == 'y' || ch == 'Y') {
			Record_1ch_SEMG(emgdata, acqfout, 10); // change to Record_SEMG or change the number of channels recorded in Record_1ch_SEMG
			RecordTEMG(emgdata, acqfout);
			Record_1ch_SEMG(emgdata, acqfout, 100);
		}
		
		// Record Sweep with Stimulation
		wcout << "Perform RecordSweepsWithStimulation Test: ";
		wcin >> ch;
		if (ch == 'y' || ch == 'Y') {
			RecordSweepsWithStimulation(emgdata, acqfout);
		}	

		wcout << "Perform RecordSweepsWithStimulationAmplitudeChange Test: ";
		wcin >> ch;
		if (ch == 'y' || ch == 'Y') {
			RecordSweepsWithStimulationAmplitudeChange(emgdata, acqfout);
		}

		wcout << "Test Data loss : ";
		wcin >> ch;
		if (ch == 'y' || ch == 'Y') {
			TestDataLoss(emgdata, acqfout);
		}

		//Record continous Data with with manual stimulation
		wcout << "Perform Record continous Data with with manual stimulation Test: ";

		wcin >> ch;
		while (ch == 'y' || ch == 'Y') {
			int ret = IsAcqAndStimConnected();
			if (ret > 0) wcout << "Error communicating with device:"; wcout << ret; wcout << "\n";
			wcout << "Set Amplitude: \n";
			float amplitude = 10;
			wcin >> amplitude;
			emgdata.m_sampling_speed = 16000;
			total_number_of_samples_to_record = 10 * emgdata.m_sampling_speed; // record 10 seconds of data
			// setup the stimulator
			// we are using the HV stimulator so we need to turn on the HV supply
			SetStimHVPowerSupply(true); // it takes about 100msec for the HV power supply to settle.

			CStimulationTrain stim_param;
			stim_param.m_delay = 10;//time in msec wait for 10msec
			stim_param.m_num_pulses = 1; // 4 pulses
			stim_param.m_bipolar = true;
			stim_param.m_pulse_width = 0.25; // 1 msec pulse
			stim_param.m_pulse_off_time = 19; // 29msec off time . m_pulse_off_time + m_pulse_width = pulse Period 
			stim_param.m_num_trains = 0;  // set to 0 for continous trains, 
			// since we want one set of pulses for every train,
			// each train duration has to match the sweep length, that way all the trains will align for all the sweeps
			// sweep length - delay - m_num_pulses*(m_pulse_width + m_pulse_off_time)
			float train_duration = 30; // in msec;
			stim_param.m_intertrain_duration = train_duration - stim_param.m_delay - stim_param.m_num_pulses * (stim_param.m_pulse_width + stim_param.m_pulse_off_time);
			stim_param.m_pulse_amplitude = stim1_amp = amplitude; // 5mA  // we are using the stim1_amp later to change the stimulation amplitude. 
			stim_param.m_stim_mode = CC50;
			stim_param.m_start_with_recording = true;
			SetStimulationParameters(0, stim_param);

			// setup the acquisition
			num_acq_channels_recorded = 2;
			num_stim_channels_recorded = 1;
			channels_to_record = 0x0A;// record channel 2 and channel 4 , ch 1 , Ch 3, channel 5-16 are skipped
			stim_channels_to_record = 0x01;

			SetAcquisitionParameters(emgdata.m_sampling_speed, LOTUS_RECORD, channels_to_record, stim_channels_to_record);
			wcout << "EMG with Manual Stimulation\n";
			fprintf(acqfout, "EMG with Manual Stimulation\n");

			// Record Data
			total_num_samples_per_ch_recorded = 0;
			num_acq_samples_per_ch = 0;
			//		int j, m, k;
			Sleep(10);
			StartAcq();
			for (k = 0; k < num_acq_channels_recorded; ++k) {
				emgdata.m_data[k].clear();
			}
			Sleep(10);
			uint32_t ch_leadoff;
			// To enable manual stimulation, you need to run the acquisition in a separate thread so that the user can press the stimulate button 
			// to emulate this we are going to assume that after 2000 samples are recorded the user presses the stimulate button. 

			while (total_num_samples_per_ch_recorded < total_number_of_samples_to_record) {
				wcout << '\r';
				wcout << "Number of samples Recorded : ";  wcout << total_num_samples_per_ch_recorded; 
				wcout << std::flush;
				num_acq_samples_per_ch = 0;
				iRet = ReadDataFromAcquisitionDevice(num_acq_samples_per_ch, emgdata.acqdata, emgdata.number_of_samples_to_read_per_loop);
				// check if leadoff was ok
				GetAcquisitionChannelLeadOff(ch_leadoff);  
				// check ch_leadoff to make sure leads were connected.
				if (iRet == 0) { // acquired data copied to data.  pts_recorded = num_samples_per_ch*num_channels_recorded;
					m = 0;
					for (j = 0; j < num_acq_samples_per_ch; ++j) {
						for (k = 0; k < num_acq_channels_recorded; ++k) {
							emgdata.m_data[k].push_back(emgdata.acqdata[m++]);
						}
					}
					total_num_samples_per_ch_recorded += num_acq_samples_per_ch;
				}
				if (iRet > 0) { // Error condition 
					StopAcq();// Stop Acquisition
					wcout << "Acquisition Board Error :  Recording Stopped\n";
					return iRet;
				}
				if (iRet < 0) { // warning, 
								//typical warning -3 ==> that no new data is available wait longer
				}
				if (num_stim_channels_recorded > 0) {
					num_acq_samples_per_ch = 0;
					iRet = ReadDataFromStimulatorDevice(num_acq_samples_per_ch, emgdata.stimdata, emgdata.number_of_samples_to_read_per_loop);
					if (iRet == 0) { // acquired data copied to data.  pts_recorded = num_samples_per_ch*num_channels_recorded;
						m = 0;
						for (j = 0; j < num_acq_samples_per_ch; ++j) {
							for (k = num_acq_channels_recorded; k < num_acq_channels_recorded + num_stim_channels_recorded; ++k) { // append the stimulator channels after the acquisition channels
								emgdata.m_data[k].push_back(emgdata.stimdata[m++]);
							}
						}
//						total_num_samples_per_ch_recorded += num_acq_samples_per_ch; // we took care of this in the acquisition section.
					}
					if (iRet > 0) { // Error condition 
						StopAcq();// Stop Acquisition
						wcout << "Stimulation Board Error :  Recording Stopped\n";
						return iRet;
					}
					if (iRet < 0) { // warning, 
									//typical warning -3 ==> that no new data is available wait longer
					}
				}
				if (manual_stimulate) { // use if you want to restart the stimulation protocol 
					manual_stimulate = false;
					FireStimulator(0x01); // we are firing the first stimulator.
					wcout << "Stimulation Started\n";
				}
				if (manual_stimulation_stop) {
					StopStimulator(0x01); // stop the stimulator 
					wcout << "Stimulation Stopped\n";
				}
				if (update_stim_amplitude) { // During continous stimulation use this to change the amplitude.  
					update_stim_amplitude = false;
					SetStimulatorAmplitude(0, stim1_amp);
				}
				// if the control parameters have changed handle the change here.
			}
			SetStimHVPowerSupply(false); // Turn off the HV power supply.
			StopAcq();// Stop Acquisition
			// End Record Data

			if (iRet > 0) { // Error condition 
				ShowErrorMesage(iRet);
			}
			else {
				float lowcutoff = 100;
				float highcutoff = 500;
				int order = 51;
				for (i = 0; i < num_acq_channels_recorded; ++i) {
					emgdata.m_filtered_data[i].resize(total_num_samples_per_ch_recorded);
					// filter the channel data
					FilterData(&(emgdata.m_data[i][0]), &(emgdata.m_filtered_data[i][0]), total_num_samples_per_ch_recorded, emgdata.m_sampling_speed, lowcutoff, highcutoff, order);
				}
				for (i = num_acq_channels_recorded; i < num_acq_channels_recorded+ stim_channels_to_record; ++i) {
					emgdata.m_filtered_data[i] = emgdata.m_data[i];
					// do not filter the stimulation data. since we want to preserve the location of the pulses
				}
				//emgdata.m_filtered_data has the filtered data but the first "order" number of data points and the last "order" number of data points are not valid.
				PrintDatatoFile(acqfout, emgdata.m_filtered_data, total_num_samples_per_ch_recorded, num_acq_channels_recorded + stim_channels_to_record);
			}
			wcout << " Repeat Record continous Data with with manual stimulation Test: ";

			wcin >> ch;
		}
		// close the files
		if(acqfout)fclose(acqfout);
	}

	/*
	wcout << "Testing Acq and Stim Disconnect: for 100 seconds\n";
	int count = 1000;
	// use a Timer to call the IsAcqAndStimConnected periodically
	while (--count > 0 && (acq_device_ok || stim_device_ok)) {
		int ret = IsAcqAndStimConnected();
		if (acq_device_ok && ret & 0x01) {
			wcout << "Acquisition Hardware disconnected\n";
			acq_device_ok = false;
		}
		if (stim_device_ok && ret & 0x02) {
			wcout << "Stimulator Hardware disconnected\n";
			stim_device_ok = false;
		}
		Sleep(100);
	}


	if (control_device_ok) {
		wcout << "Testing Control device: Press q to quit\n";
		wcin >> ch;
	}
	*/
	
	CloseDevices(); // to delete the devices created and clear memory.
    return 0;
}



