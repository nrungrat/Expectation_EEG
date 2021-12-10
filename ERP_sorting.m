% This script sorts the EEG trials based on conditions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% p.freq = 1 --> red fast 50hz/ blue slow 33.33 Hz
clear all; close all;

for sub = [1 2 4:9 11:14 16 18 19 20 21];
    %sub = 21;
    cd (['data/sbj' num2str(sub)]);
    
    load(['eegEpochArtFree_rmbase_sbj' num2str(sub)]);
    
    close all;
    clear eindex
    [b,a] = butter(3, 10/(512/2),'low'); % low pass filter
    
    [reorder chanNum chanLabel chanNumFlip goodChan chanCon chanIps poz_mid] = chanMontage;
    
    fstarget = [ 1 2; 2 1]; % set up the converntion of the target color and speed
    % this loop sorts out trial types
    for bias = 1:3 % color bias/ orientation bias/ object-response bias
        for prior = 1:3 % 70 50 30 --- note 50 is from different blocks (1, 2,3 = expected, neutral, unexpected)
            for fastslow = 1:2 % fast and slow flicker rates of the target (50/33.33 Hz)
                for hit = 1:3 % 1 = incorrect+miss / 2 = correct/ 3 all
                    %                             for tgcolor = 1:2 % 1 = red target // 2 = blue
                    %                             for tgori = 1:2 % 1 = vertical // 2 = horizontal
                    %
                    if hit <3 % if incorrect or correct
                        if prior ==1 || prior ==3 % trials that are not neutral
                            %%% find trial index to sort EEG
                            %%% data
                            eindex{bias, prior, fastslow, hit} = find(ismember(epochrm.typetg(:, 1), ...
                                p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == bias & p.hit ==hit-1 ...
                                & p.prior_705030 == prior  ...
                                & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))))==1);
                            
                            %                             ermindex{bias, prior, fastslow, hit} = find(ismember(epochrm.type(:, 1), ...
                            %                                 intersect(p.stimcodecum(p.bias_no_color_ori_resp == bias & p.hit ==hit-1 ...
                            %                                 & p.prior_705030 == prior  ...
                            %                                 & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                            %                                 | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))), epochrm.typetg(:, 1)))==1);
                            %
                            rtindex{bias, prior, fastslow, hit} = ...
                                intersect(p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == bias & p.hit ==hit-1 ...
                                & p.prior_705030 == prior  ...
                                & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))), epochrm.typetg(:, 1));
                            
                            
                        elseif prior ==2
                            eindex{bias, prior, fastslow, hit} = find(ismember(epochrm.typetg(:, 1), ...
                                p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == 0 & p.hit ==hit-1 ...
                                & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))))==1);
                            
                            %                             ermindex{bias, prior, fastslow, hit} = find(ismember(epochrm.type(:, 1), ...
                            %                                 intersect(p.stimcodecum(p.bias_no_color_ori_resp == 0 & p.hit ==hit-1 ...
                            %                                 & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                            %                                 | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))), epochrm.typetg(:, 1)))==1);
                            %
                            rtindex{bias, prior, fastslow, hit} =  ...
                                intersect(p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == 0 & p.hit ==hit-1 ...
                                & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))), epochrm.typetg(:, 1));
                        end
                    elseif hit ==3 % all trials correct+ inccorect
                        if prior ==1 || prior ==3
                            eindex{bias, prior, fastslow, hit} = find(ismember(epochrm.typetg(:, 1), ...
                                p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == bias  ...
                                & p.prior_705030 == prior  & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))))==1);
                            
                            %                             ermindex{bias, prior, fastslow, hit} = find(ismember(epochrm.type(:, 1), ...
                            %                                 intersect(p.stimcodecum(p.bias_no_color_ori_resp == bias  ...
                            %                                 & p.prior_705030 == prior  & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                            %                                 | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))), epochrm.typetg(:, 1)))==1);
                            %
                            rtindex{bias, prior, fastslow, hit} =  ...
                                intersect(p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == bias  ...
                                & p.prior_705030 == prior  & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))), epochrm.typetg(:, 1));
                            
                        elseif prior ==2
                            eindex{bias, prior, fastslow, hit} = find(ismember(epochrm.typetg(:, 1), ...
                                p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == 0  ...
                                & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))))==1);
                            
                            %                             ermindex{bias, prior, fastslow, hit} = find(ismember(epochrm.type(:, 1), ...
                            %                                 intersect(p.stimcodecum(p.bias_no_color_ori_resp == 0  ...
                            %                                 & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                            %                                 | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))), epochrm.typetg(:,1)))==1);
                            %
                            rtindex{bias, prior, fastslow, hit} = ...
                                intersect(p.stimcodecum(p.stimcode > 20 & p.bias_no_color_ori_resp == 0  ...
                                & (p.freq ==1 & p.color_redblue == fstarget(fastslow,1) ...
                                | p.freq ==2 & p.color_redblue == fstarget(fastslow,2))), epochrm.typetg(:,1));
                        end
                    end
                    
                    clear eegresplocked  eegtglocked
                    
                    eegtglocked = epochrm.datatg(:, :, eindex{bias, prior, fastslow, hit});
                    eegresplocked(:, :, 1) = nan(72,820,1);
                    
                    
                    for trial_RT = 1:length(rtindex{bias, prior, fastslow, hit})
                        clear responset
                        trial_RT
                        if isempty(trial_RT)
                            eegresplocked(:, :, trial_RT) = nan(72,820,1);
                        else
                            
                            responset = 768+round((p.RT_fromtgonset(rtindex{bias, prior, fastslow, hit}(trial_RT))*1000)./(1000/512));
                            
                            if isnan(responset) || responset > 1536-51
                                eegresplocked(:, :, trial_RT) = nan(72,820,1);
                            else
                                
                                eegresplocked(:, :, trial_RT) =   eegtglocked(:, responset-768: responset+51, trial_RT);
                                
                            end
                        end
                    end
                    erptg(bias, prior, fastslow, hit,  :, :) = mean(eegtglocked,3) ;
                    erpresp(bias, prior, fastslow, hit,  :, :) = nanmean(eegresplocked, 3);
                    
                    erptgf(bias, prior, fastslow, hit,  :, :) = filtfilt(b,a,double(squeeze(erptg(bias, prior, fastslow, hit, :, :)))')';
                    erprespf(bias, prior, fastslow, hit,  :, :) = filtfilt(b,a,double(squeeze(erpresp(bias, prior, fastslow, hit, :, :)))')';
                    
                end
            end
        end
    end
    
    save (['erp_tglocked_cutfirst20_sbj' num2str(sub)],  'erptg', 'erptgf');
    save(['erp_resplocked_cutfirst20_subj' num2str(sub)], 'erpresp', 'erprespf');
    
    cd ../..
    
end
