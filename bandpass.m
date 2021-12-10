eegtemp.srate = 512;

[b,a] = butter(3,4/(eegtemp.srate/2),'high'); % high pass
[b1,a1] = butter(3,8/(eegtemp.srate/2),'low'); % low pass

for i = 1:1625
    i
filt.data(i, :, :) = filtfilt(b,a,double(squeeze(epoch.data(:,:, i)))')'; % 1:64 chan filter data here
filt.data(i, :, :) = filtfilt(b1,a1,double(squeeze(filt.data(i,:, :)))')';
end