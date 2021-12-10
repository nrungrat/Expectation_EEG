function eegtemp = rbdf_stroopica(filename)
%This function reads in .bdf files and organizes the file into something sean-style.
%Uses the left and right mastoids as the reference electrodes
%
%inputs:
%	filename 		- raw file name of bdf file
%outputs:
%	eeg.filenames 		- raw file name of bdf file
%	eeg.samplingrate 	- recording sampling rate
%	eeg.data		- 64 channel data, average referenced
%	eeg.xdata		- 8 external channels (e.g., mastoid, EOG, etc.)
%	eeg.chanlabels		- channel labels of the 64 channels
%	eeg.xchanlabels		- channel labels of the 7 external channels
%	eeg.triggers.type	- event type
%	eeg.triggers.time	- sample index corresponding to event in 'type'
%
%


eegtemp = pop_biosig(filename); % load bdf 

%eeg.filename = filename;
%eeg.samplingrate = eegtemp.srate; 

%bandpass filter 0.25 hz - 58 hz 
[b,a] = butter(3, .25/(eegtemp.srate/2),'high'); % high pass
[b1,a1] = butter(3, 58/(eegtemp.srate/2),'low'); % low pass

%reference to mastoids
refdata = squeeze(mean(eegtemp.data([65 66],:))); % take mean pf left and right mastriods 
eegtemp.data = eegtemp.data-repmat(refdata,[eegtemp.nbchan 1]);

eegtemp.data = filtfilt(b,a,double(eegtemp.data(:,:))')'; % 1:64 chan filter data here
eegtemp.data = filtfilt(b1,a1,double(eegtemp.data(:,:))')';


for k = 1:length(eegtemp.event) % 
eegtemp.event(k).type = abs(2^16-eegtemp.event(k).type-256); % convert bdf code to binary code 

end;

