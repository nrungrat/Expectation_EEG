%updated: 
%03/20/18: cutting out the first 20 trials
%%
% p.freq = 1 --> red fast 50hz/ blue slow 33.33 Hz
clear all; close all;

%% for the purpose of baseline removal
bandwidth = 0.2;
cnt = 0;
for iii = 4:8 % theta freq
    cnt = cnt+1;
    fdsd(cnt) = bandwidth*(iii/(2*(2*log(2))^.5));
    tdsd(cnt) = ((2*pi*fdsd(cnt)).^-1)*1000; %in ms
end
max_tdsd = round(max(tdsd), 0);

tlength = 200; %200 ms; this is the same length as what is removed from ERPs (-200 to 0 ms)
td_pretgdur = -1*(tlength+max_tdsd): -1*max_tdsd; %time-domain (ms)
fq = 1000/512;

for sub = [1 2 4:9 11:14 16 18 19 20 21]; %17 subjects
    sub
    cd (['data/sbj' num2str(sub)]);
    
    %load(['sbj_' num2str(sub) 'alphatg.mat']);
    %load(['eegEpochArtFree_rmbase_sbj' num2str(sub)]);
    load(['sbj_' num2str(sub) 'thetatg_big.mat']);
    load(['eegEpochArtFree_rmbasentg_bigepoch_sbj' num2str(sub) '.mat']);
    
    tp_pretgdur = nearest(epochrmntg.timetg, td_pretgdur(1)): nearest(epochrmntg.timetg, td_pretgdur(end));
    tp_prerespdur = tp_pretgdur;
    
    close all;
    clear pindex
    [b,a] = butter(3, 10/(512/2), 'low'); % low pass filter
    
    [reorder chanNum chanLabel chanNumFlip goodChan chanCon chanIps poz_mid] = chanMontage;
    
    fstarget = [ 1 2; 2 1]; % set up the converntion of the target color and speed
    
    %% this loop sorts out trial types
    for bias = 1:3 % color bias/ orientation bias/ object-response bias
        for prior = 1:3 % 70 50 30 --- note 50 is from different blocks (1, 2,3 = expected, neutral, unexpected)
            for fastslow = 1:2 % fast and slow flicker rates of the target (50/33.33 Hz)
                for hit = 1:3 % 1 = incorrect+miss / 2 = correct/ 3 all
                    
                    if hit <3 % if incorrect or correct
                        if prior ==1 || prior ==3 % trials that are not neutral
                            pindex = p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == bias & p.hit ==hit-1 ...
                                & p.prior_705030 == prior  ...
                                & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2)));
                            
                        elseif prior ==2
                            pindex = ...
                                p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == 0 & p.hit ==hit-1 ...
                                & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2)));
                        end
                        
                    elseif hit ==3 % all trials correct+ incorrect
                        if prior ==1 || prior ==3
                            pindex =  ...
                                p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == bias  ...
                                & p.prior_705030 == prior  & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2)));
                            
                        elseif prior ==2
                            pindex =  ...
                                p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == 0  ...
                                & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2)));
                        end
                    end
                    
                    %%
                    RT_all = p.RT_fromtgonset(pindex); %p.RT_fromtgonset has reaction times
             
                    RT_indexpool{bias, prior, fastslow, hit} = pindex;
                 
                    eindexpool{bias, prior, fastslow, hit} = find(ismember(epochrmntg.typetg(:, 1), ...
                        RT_indexpool{bias, prior, fastslow, hit})); %eeg trial label
                
                    clear eegtglocked_pool eegtglocked_fast eegtglocked_slow eegresplocked_pool eegresplocked_fast eegresplocked_slow
                    
                    %amp is a matrix of the size [remaining_trial_numx72x1536]
                    eegtglocked_pool = amp(eindexpool{bias, prior, fastslow, hit}, :, :);
                   
                    respdeadline = 1400; %disregard trials where responset >= 1400 ms after tg onset
                    beforeresp = 1500; %ms; how far we want to go before resp onset
                    afterresp = 100; %ms; how far we want to go after resp onset
                    
                    eegresplocked_pool(1, :, :) = nan(1, 72, round((beforeresp + afterresp)./fq) + 1);
                  
                    %% response locked: for all trials (collapsed fast/slow trials)
                    
                    for trial_RT_pool = 1:length(eindexpool{bias, prior, fastslow, hit})
                        clear responset
                        responset = (find(epochrmntg.timetg == 0) -1) + ...
                            round((p.RT_fromtgonset(epochrmntg.typetg(eindexpool{bias, prior, fastslow, hit}...
                            (trial_RT_pool)))*1000)./fq); %find the timepoint for resp onset
                        
                        %1) everything after plus sign is just to convert RT t
                        %to timepoint
                        %2) then the term before the plus sign is just the
                        %timepoint of 0 ms
                        %3) epochrmntgrm.timetg(responset) is RT on that trial
                        

                        if isnan(responset) || responset > nearest(epochrmntg.timetg, respdeadline)
                            %in the case of nan-onset or onset that's too
                            %late (i.e., occured >= 1500-respdeadline ms)
                            eegresplocked_pool(trial_RT_pool, :, :) = nan(1, 72, round((beforeresp + afterresp)./fq) + 1);
                            
                        else
                            %remove baseline: based on tp_prerespdur
                            %defined up top
                            eegresplocked_pool(trial_RT_pool, :, :) =   ... 
                                eegtglocked_pool(trial_RT_pool, :, ...
                                responset - round(beforeresp./fq): responset + round(afterresp./fq)) ...
                                - repmat(mean(eegtglocked_pool(trial_RT_pool, :, tp_prerespdur), 3), ...
                                [1 1 size(responset - round(beforeresp./fq): responset + round(afterresp./fq), 2)]);
                            
                            %don't remove baseline
                            eegresplocked_pool_nope(trial_RT_pool, :, :) =   ...
                                eegtglocked_pool(trial_RT_pool, :, ...
                                responset - round(beforeresp./fq): responset + round(afterresp./fq));
                        end
                    end
                    
                 
                    %%
                    thetatg_pool(bias, prior, fastslow, hit,  :, :) = ... 
                        mean(eegtglocked_pool, 1)- repmat(mean(mean(eegtglocked_pool(:, :, tp_pretgdur), 3), 1), ...
                        [1 1 size(eegtglocked_pool, 3)]);
                    thetatg_pool_nope(bias, prior, fastslow, hit,  :, :) = ...
                        mean(eegtglocked_pool, 1);
       
                    thetaresp_pool(bias, prior, fastslow, hit,  :, :) = ...
                        nanmean(eegresplocked_pool, 1);
                    thetaresp_pool_nope(bias, prior, fastslow, hit,  :, :) = ...
                        nanmean(eegresplocked_pool_nope, 1);
                    
                   
                end
            end
        end
    end
    
    %baseline not removed
    save (['big_theta_tg_nope_cut20_sbj' num2str(sub)],  'thetatg_pool_nope'); 
    save (['big_theta_resp_nope_cut20_sbj' num2str(sub)],  'thetaresp_pool_nope');
    %with baseline being removed up to -1SD
    save (['big_theta_tg_sd_cut20_sbj' num2str(sub)],  'thetatg_pool');
    save (['big_theta_resp_sd_cut20_sbj' num2str(sub)],  'thetaresp_pool'); 
    
    
    cd ../..
    
end

