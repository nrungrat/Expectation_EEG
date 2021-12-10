%updated
%031318: for ssvep; we will be using the erp locked to ntg so
%epochrmntg.data not the .datatg one
%091517: run ssvep 30 and ssvep 50
clear all; close all; 

for sub = [1 2 4:9 11:14 15 16 18 19 20 21];
cd (['data/sbj' num2str(sub)]);

% %the data here were epoched [-1.5 1.5] tg-locked
%load (['eegEpochArtFree_rmbase_sbj' num2str(sub) '.mat']);

%08/09/17
%the data here were epoched [-1.5 2.5] tg-locked
load(['eegEpochArtFree_rmbasentg_bigepoch_sbj' num2str(sub) '.mat']);

numtrial= size(epochrmntg.data, 3);
numtrialtg= size(epochrmntg.datatg, 3);
timeLength = size(epochrmntg.data,2);
timeLengthtg = size(epochrmntg.datatg,2);
chanNo = 72; 
freq = 1:58;
%fractional_bandwidth = [.2 .2 .2 .05 .05]; default
fractional_bandwidth = [.5 .5]; %09/15/17

sampRate = 512;

%for big epoch
sname = {['sbj_' num2str(sub) 'ssvep33_big'] ; ... 
  ['sbj_' num2str(sub) 'ssvep50_big'] };

snametg = {['sbj_' num2str(sub) 'ssvep33tg_big'] ; ... 
  ['sbj_' num2str(sub) 'ssvep50tg_big'] };

time = epochrmntg.time;
timetg = epochrmntg.timetg;

wband{1} = 33.33333333333; %ssvep30
wband{2} = 50; %ssvep50


for floop = 1:2 % loop through each freq
   clear amp
   amp = nan(numtrial, chanNo, timeLength);   
    for trial = 1:numtrial;
        fprintf ([sname{floop} ' stim onset: freq/ trial ' num2str(trial) '\n']);
        
        for chan = 1:72% loop through each electrode
            clear signal;
            signal = squeeze(epochrmntg.data(chan, :, trial));
            
           fcnt = 0;
           clear coEff  
           for fq = wband{floop} 
               fcnt = fcnt+1;

           coEff(fcnt, :) = gaussian_filter_signal_pcl ...
                ('raw_signal', signal, 'sampling_rate', sampRate,...
                'center_frequency',fq,'fractional_bandwidth', fractional_bandwidth(floop));
           end
            
           %if floop <=2
            %   amp(trial, chan, :) = mean(abs(coEff));
           %else %ssvep
               amp(trial, chan, :) = abs(coEff);
           %end
        end
        
    end
    
    save(sname{floop}, 'amp');
end

for floop = 1:2 % loop through each freq
    clear amp
   amp = nan(numtrialtg, chanNo, timeLengthtg);   
    for trial = 1:numtrialtg;
        
        fprintf ([sname{floop} ' tg onset: freq/ trial ' num2str(trial) '\n']);
        
        for chan = 1:72% loop through each electrode

            clear signal;
            signal = squeeze(epochrmntg.datatg(chan, :, trial));
            
             fcnt = 0;
           clear coEff  
           for fq = wband{floop} 
               fcnt = fcnt+1;

           coEff(fcnt, :) = gaussian_filter_signal_pcl ...
                ('raw_signal', signal, 'sampling_rate', sampRate,...
                'center_frequency',fq,'fractional_bandwidth', fractional_bandwidth(floop));
           end
            
           if floop <=2
               amp(trial, chan, :) = mean(abs(coEff));
           else
               amp(trial, chan, :) = (coEff);
           end
        end
        
    end
    
    save(snametg{floop}, 'amp');
end


cd ../..
end
    
%compute time domain SD of gaussian wavelet
%frequency_domain_standard_deviation=...
%                varargin{n+1}*(center_frequency/(2*(2*log(2))^.5));
% time_domain_standard_deviation=...
%     (2*pi*frequency_domain_standard_deviation).^-1;
%           when tdsd is in seconds, fdsd is in Hz
% bandwidth = 0.2;
% cnt = 0;
% for iii = 5:12 %freq
%     cnt = cnt+1;
%     fdsd(cnt) = bandwidth*(iii/(2*(2*log(2))^.5));
%     tdsd(cnt) = (2*pi*fdsd(cnt)).^-1;
% end

