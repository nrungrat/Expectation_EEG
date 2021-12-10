% Nuttida last updated on 3 March 2018
% Note: 
% - This script analyzes ERPs for all minus the first 20 trials of each
% block
%--------------------------------------------------------------------------

addpath(genpath('/local/home/serenceslab/serenceslab_toolboxes/genTools/'))
clear all; close all;
sub = [1 2 4:9 11:14 16 18 19 20 21]; %17 subjects

%% stacking all subjects' data
for s = 1:size(sub,2)
    cd (['data/sbj' num2str(sub(s))]);
    clear erp erpf erptg erptgf erp_resp erp_respf
    
    %load (['erp_sbj' num2str(sub(s))]); %this one has both trial-onset locked and tg locked
    load(['erp_tglocked_cutfirst20_sbj' num2str(sub(s))]); %tg locked
    load(['erp_resplocked_cutfirst20_subj' num2str(sub(s))]); %resp locked
    
    cd ../..
    
    %trial onset locked
    %     gerp(s, :, :, :, :, :, :) = erp;
    %     gerpf(s, :, :, :, :, :, :) = erpf; %filtered
    
    %tg locked
    gerp_tg(s, :, :, :, :, :, :) = erptg;
    gerp_tgf(s, :, :, :, :, :, :) = erptgf; %filtered
    
    %resp locked
    gerp_resp(s, :, :, :, :, :, :) = erpresp;
    gerp_respf(s, :, :, :, :, :, :) = erprespf; %filtered
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% for figures plotting 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig_gerp_tg_flicker = gerp_tgf - repmat(nanmean((gerp_tgf), 4), [1, 1, 1, 2, 1, 1, 1]) + ...
    repmat(nanmean(nanmean((gerp_tgf), 4)), [s, 1, 1, 2, 1, 1, 1]);
fig_gerp_resp_flicker = gerp_respf - repmat(nanmean((gerp_respf), 4), [1, 1, 1, 2, 1, 1, 1]) + ...
    repmat(nanmean(nanmean((gerp_respf), 4)), [s, 1, 1, 2, 1, 1, 1]);

fig_gerp_tg_exp = gerp_tgf - repmat(nanmean((gerp_tgf), 3), [1, 1, 3, 1, 1, 1, 1]) + ...
    repmat(nanmean(nanmean((gerp_tgf), 3)), [s, 1, 3, 1, 1, 1, 1]);
fig_gerp_resp_exp = gerp_respf - repmat(nanmean((gerp_respf), 3), [1, 1, 3, 1, 1, 1, 1]) + ...
    repmat(nanmean(nanmean((gerp_respf), 3)), [s, 1, 3, 1, 1, 1, 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plotting stuff 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b = beh_analysis(sub); %load behavioral data
fq = 1000/512; %512 = EEG sampling rate

clear erp erpf

timex.tg =  -1500:fq:1499; %time vector in ms --> 1536 time points for tg locked
timex.resp =  -1500:fq:100; %time vector in ms --> 820 time points for resp locked

condi = {'Color Expectation', 'Orientation Expectation', 'Motor Expectation'}; %condition labels
%blabel = {'RT_Color', 'rt_ori', 'rt_resp'}; %labels for behavioral

chan_label = {'FPz', 'AFz', 'Fz', 'FCz' 'Cz', 'CPz', 'Pz', 'POz', '0z'}; %labels for all electrodes

close all;
eoi = [33; ...
    37; ...
    38 ;...
    47  ; ...
    48 ; ...
    32 ;...
    31 ; ...
    30  ; ...
    29  ;];

close all
rtline = -4:.01:8;

% electrodes used to extract CPP and VN
ch_cpz = 6;
ch_oz = 9;

%% Repeated-Measures 3-way ANOVA
%subject, cond (color/ori/resp), prior, flickerSpeed, hit(2), chan , time point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% tg-locked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
binsize = 50;
%chosentime.tg = -100:binsize:1500; %time duration of interest (timepoint 718-1536)
chosentime.tg = -250:binsize:1500;

for b = 1:numel(chosentime.tg)
    diff.tg = abs(timex.tg - chosentime.tg(b));
    timebin.tg(b) = find(diff.tg == min(diff.tg)); %time points that correspond to the starting time (in ms) of each bin
end

for bb = 1:numel(timebin.tg)-1
    gerp_tgbinned( :, :, :, :, :, :, bb) = mean(gerp_tg( :, :, :, :, :, :, timebin.tg(bb):timebin.tg(bb+1)-1 ), 7);
end

cntch = 0;
for ch = [ch_cpz ch_oz] %CPz and Oz
    cntch = cntch+1;
    cntt = 0;
    for ttt = 1:numel(timebin.tg)-1 %a total of 32 time bins of interest
        cntt = cntt +1;
        condcolor_fast_exp = squeeze(gerp_tgbinned(:, 1, 1, 1, 2, eoi(ch, :), ttt)); %50Hz
        condcolor_fast_unexp = squeeze(gerp_tgbinned(:, 1, 3, 1, 2, eoi(ch, :), ttt));
        condcolor_fast_neu = squeeze(gerp_tgbinned(:, 1, 2, 1, 2, eoi(ch, :), ttt));
        condcolor_slow_exp = squeeze(gerp_tgbinned(:, 1, 1, 2, 2, eoi(ch, :), ttt)); %33Hz
        condcolor_slow_unexp = squeeze(gerp_tgbinned(:, 1, 3, 2, 2, eoi(ch, :), ttt));
        condcolor_slow_neu = squeeze(gerp_tgbinned(:, 1, 2, 2, 2, eoi(ch, :), ttt));
        
        condori_fast_exp = squeeze(gerp_tgbinned(:, 2, 1, 1, 2, eoi(ch, :), ttt)); %50Hz
        condori_fast_unexp = squeeze(gerp_tgbinned(:, 2, 3, 1, 2, eoi(ch, :), ttt));
        condori_fast_neu = squeeze(gerp_tgbinned(:, 2, 2, 1, 2, eoi(ch, :), ttt));
        condori_slow_exp = squeeze(gerp_tgbinned(:, 2, 1, 2, 2, eoi(ch, :), ttt)); %33Hz
        condori_slow_unexp = squeeze(gerp_tgbinned(:, 2, 3, 2, 2, eoi(ch, :), ttt));
        condori_slow_neu = squeeze(gerp_tgbinned(:, 2, 2, 2, 2, eoi(ch, :), ttt));
        
        condresp_fast_exp = squeeze(gerp_tgbinned(:, 3, 1, 1, 2, eoi(ch, :), ttt)); %50Hz
        condresp_fast_unexp = squeeze(gerp_tgbinned(:, 3, 3, 1, 2, eoi(ch, :), ttt));
        condresp_fast_neu = squeeze(gerp_tgbinned(:, 3, 2, 1, 2, eoi(ch, :), ttt));
        condresp_slow_exp = squeeze(gerp_tgbinned(:, 3, 1, 2, 2, eoi(ch, :), ttt)); %33Hz
        condresp_slow_unexp = squeeze(gerp_tgbinned(:, 3, 3, 2, 2, eoi(ch, :), ttt));
        condresp_slow_neu = squeeze(gerp_tgbinned(:, 3, 2, 2, 2, eoi(ch, :), ttt));
        
        DV = [condcolor_fast_exp; condcolor_fast_unexp; condcolor_fast_neu; ...
            condcolor_slow_exp; condcolor_slow_unexp; condcolor_slow_neu; ...
            condori_fast_exp; condori_fast_unexp; condori_fast_neu; ...
            condori_slow_exp; condori_slow_unexp; condori_slow_neu; ...
            condresp_fast_exp; condresp_fast_unexp; condresp_fast_neu; ...
            condresp_slow_exp; condresp_slow_unexp; condresp_slow_neu];
        
        sublist = repmat(1:numel(sub), [1,18])'; % subject list
        IVexp = repmat([ones(1,numel(sub)) ones(1,numel(sub)).*2 ones(1,numel(sub)).*3], [1, 6])'; % independent variable: exp/neutral/unexp
        IVtask = repmat([ones(1,numel(sub)*6) ones(1,numel(sub)*6).*2 ones(1,numel(sub)*6).*3], [1 1])'; %task condition: color/ori/resp
        IVflick= repmat([ones(1,numel(sub)*3) ones(1,numel(sub)*3).*2], [1, 3])'; % independent variable: 50Hz/33Hz
        
        [p, tb1, stats] = anovan(DV, {IVexp, IVtask, IVflick, sublist}, ...
            'model', 'full', 'varnames', {'exp/neu/unexp', 'color/ori/resp', '50Hz/33Hz', 'subj'}, 'random', [4]);
        %'random', [4] means that we are
        %doing repeated measure ANOVA
        %within subjects
        %and we want the 4th matrix (subj
        %in this case) to be random
        %*if we don't specify this, it wont
        %be repeated measure
        
       
        panova_tg(cntch, cntt, :) = p; %associated p-values for all factor comparisons
        %size 2x32x15 (chansxtime binsxnumber of all comparisons)
        %in this case we have 4 things:
        %subj, exp, flicker rates, bias
        %conds; so it's
        %nchoosek(4,4)+(4,3)+(4,2)+(4,1)=15
        
        fvalexp_tg(cntch, cntt, :) = tb1(2, 6); %f-values of expectaion factor
        %pvalexp_tg(cntch, cntt, :) = tb1(2, 7);
        fvalflick_tg(cntch, cntt, :) = tb1(4, 6); %f-values of flicker rate factor
        %pvalflick_tg(cntch, cntt, :) = tb1(4, 7);
        fvalcond_tg(cntch, cntt, :) = tb1(3, 6); %f-values of bias cond factor
        %pvalcond_tg(cntch, cntt, :) = tb1(3, 7);
                %interaction effect of exp x flicker rate
        fvalexpflick_tg(cntch, cntt, :) = tb1(7, 6);
    end
end

%p-values by effects
%main effects
panova_cpz_exp_tg = panova_tg(1, :, 1); %expectation
panova_oz_exp_tg = panova_tg(2, :, 1);
panova_cpz_cond_tg = panova_tg(1, :, 2); %task conditions
panova_oz_cond_tg = panova_tg(2, :, 2);
panova_cpz_flick_tg = panova_tg(1, :, 3); %flicker rates
panova_oz_flick_tg = panova_tg(2, :, 3);

%interaction effects
panova_cpz_exp_cond_tg = panova_tg(1, :, 5); %exp x cond
panova_oz_exp_cond_tg = panova_tg(2, :, 5);
panova_cpz_exp_flick_tg = panova_tg(1, :, 6); %exp x flicker rates
panova_oz_exp_flick_tg = panova_tg(2, :, 6);
panova_cpz_cond_flick_tg = panova_tg(1, :, 8); %cond x flicker rates
panova_oz_cond_flick_tg = panova_tg(2, :, 8);
panova_cpz_exp_cond_flick_tg = panova_tg(1, :, 11); %exp x cond cond x flicker rates
panova_oz_exp_cond_flick_tg = panova_tg(2, :, 11);

fvalexp_cpz_tg = fvalexp_tg(1, :);
fvalexp_oz_tg = fvalexp_tg(2, :);
fvalflick_cpz_tg = fvalflick_tg(1, :);
fvalflick_oz_tg = fvalflick_tg(2, :);
fvalcond_cpz_tg = fvalcond_tg(1, :);
fvalcond_oz_tg = fvalcond_tg(2, :);
%interaction effect of exp x flicker rate
fvalexpflick_cpz_tg = fvalexpflick_tg(1, :);
fvalexpflick_oz_tg = fvalexpflick_tg(2, :);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% resp-locked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear b bb tt ttt diff
binsize = 50;
chosentime.resp = -300:binsize:100; %time duration of interest

for b = 1:numel(chosentime.resp)
    diff.resp = abs(timex.resp - chosentime.resp(b));
    timebin.resp(b) = find(diff.resp == min(diff.resp));
end

for bb = 1:numel(timebin.resp)-1
    gerp_respbinned( :, :, :, :, :, :, bb) = mean(gerp_resp( :, :, :, :, :, :, timebin.resp(bb):timebin.resp(bb+1)-1 ), 7);
end

cntch = 0;
for ch = [ch_cpz ch_oz]
    cntch = cntch+1;
    cntt = 0;
    for ttt = 1:numel(timebin.resp)-1
        cntt = cntt +1;
        Rcondcolor_fast_exp = squeeze(gerp_respbinned(:, 1, 1, 1, 2, eoi(ch, :), ttt)); %50Hz
        Rcondcolor_fast_unexp = squeeze(gerp_respbinned(:, 1, 3, 1, 2, eoi(ch, :), ttt));
        Rcondcolor_fast_neu = squeeze(gerp_respbinned(:, 1, 2, 1, 2, eoi(ch, :), ttt));
        Rcondcolor_slow_exp = squeeze(gerp_respbinned(:, 1, 1, 2, 2, eoi(ch, :), ttt)); %50Hz
        Rcondcolor_slow_unexp = squeeze(gerp_respbinned(:, 1, 3, 2, 2, eoi(ch, :), ttt));
        Rcondcolor_slow_neu = squeeze(gerp_respbinned(:, 1, 2, 2, 2, eoi(ch, :), ttt));
        
        Rcondori_fast_exp = squeeze(gerp_respbinned(:, 2, 1, 1, 2, eoi(ch, :), ttt)); %50Hz
        Rcondori_fast_unexp = squeeze(gerp_respbinned(:, 2, 3, 1, 2, eoi(ch, :), ttt));
        Rcondori_fast_neu = squeeze(gerp_respbinned(:, 2, 2, 1, 2, eoi(ch, :), ttt));
        Rcondori_slow_exp = squeeze(gerp_respbinned(:, 2, 1, 2, 2, eoi(ch, :), ttt)); %50Hz
        Rcondori_slow_unexp = squeeze(gerp_respbinned(:, 2, 3, 2, 2, eoi(ch, :), ttt));
        Rcondori_slow_neu = squeeze(gerp_respbinned(:, 2, 2, 2, 2, eoi(ch, :), ttt));
        
        Rcondresp_fast_exp = squeeze(gerp_respbinned(:, 3, 1, 1, 2, eoi(ch, :), ttt)); %50Hz
        Rcondresp_fast_unexp = squeeze(gerp_respbinned(:, 3, 3, 1, 2, eoi(ch, :), ttt));
        Rcondresp_fast_neu = squeeze(gerp_respbinned(:, 3, 2, 1, 2, eoi(ch, :), ttt));
        Rcondresp_slow_exp = squeeze(gerp_respbinned(:, 3, 1, 2, 2, eoi(ch, :), ttt)); %50Hz
        Rcondresp_slow_unexp = squeeze(gerp_respbinned(:, 3, 3, 2, 2, eoi(ch, :), ttt));
        Rcondresp_slow_neu = squeeze(gerp_respbinned(:, 3, 2, 2, 2, eoi(ch, :), ttt));
        
        RDV = [Rcondcolor_fast_exp; Rcondcolor_fast_unexp; Rcondcolor_fast_neu; ...
            Rcondcolor_slow_exp; Rcondcolor_slow_unexp; Rcondcolor_slow_neu; ...
            Rcondori_fast_exp; Rcondori_fast_unexp; Rcondori_fast_neu; ...
            Rcondori_slow_exp; Rcondori_slow_unexp; Rcondori_slow_neu; ...
            Rcondresp_fast_exp; Rcondresp_fast_unexp; Rcondresp_fast_neu; ...
            Rcondresp_slow_exp; Rcondresp_slow_unexp; Rcondresp_slow_neu];
        
        Rsublist = repmat(1:numel(sub), [1,18])'; % subject list
        RIVexp = repmat([ones(1,numel(sub)) ones(1,numel(sub)).*2 ones(1,numel(sub)).*3], [1, 6])'; % independent variable: exp/neutral/unexp
        RIVtask = repmat([ones(1,numel(sub)*6) ones(1,numel(sub)*6).*2 ones(1,numel(sub)*6).*3], [1 1])'; %task condition: color/ori/resp
        RIVflick= repmat([ones(1,numel(sub)*3) ones(1,numel(sub)*3).*2], [1, 3])'; % independent variable: 50Hz/33Hz
        
        [Rp, Rtb1, Rstats] = anovan(RDV, {RIVexp, RIVtask, RIVflick, Rsublist}, ...
            'model', 'full', 'varnames', {'exp/neu/unexp', 'color/ori/resp', '50Hz/33Hz', 'subj'}, 'random', [4]);
        panova_resp (cntch, cntt, :) = Rp;
        fvalexp_resp (cntch, cntt, :) = Rtb1(2, 6);
        fvalflick_resp (cntch, cntt, :) = Rtb1(4, 6);
        fvalcond_resp (cntch, cntt, :) = Rtb1(3, 6);
        %interaction effect of exp x flicker rate
        fvalexpflick_resp (cntch, cntt, :) = Rtb1(7, 6);
    end
end

%p-values by effects
%main effects
panova_cpz_exp_resp = panova_resp(1, :, 1); %expectation
panova_oz_exp_resp = panova_resp(2, :, 1);
panova_cpz_cond_resp = panova_resp(1, :, 2); %cond conditions
panova_oz_cond_resp = panova_resp(2, :, 2);
panova_cpz_flick_resp = panova_resp(1, :, 3); %flicker rates
panova_oz_flick_resp = panova_resp(2, :, 3);

%interaction effects
panova_cpz_exp_cond_resp = panova_resp(1, :, 5); %exp x cond
panova_oz_exp_cond_resp = panova_resp(2, :, 5);
panova_cpz_exp_flick_resp = panova_resp(1, :, 6); %exp x flicker rates
panova_oz_exp_flick_resp = panova_resp(2, :, 6);
panova_cpz_cond_flick_resp = panova_resp(1, :, 8); %cond x flicker rates
panova_oz_cond_flick_resp = panova_resp(2, :, 8);
panova_cpz_exp_cond_flick_resp = panova_resp(1, :, 11); %exp x cond x flicker rates
panova_oz_exp_cond_flick_resp = panova_resp(2, :, 11);

%f-values
fvalexp_cpz_resp = fvalexp_resp(1, :);
fvalexp_oz_resp = fvalexp_resp(2, :);
fvalflick_cpz_resp = fvalflick_resp(1, :);
fvalflick_oz_resp = fvalflick_resp(2, :);
fvalcond_cpz_resp = fvalcond_resp(1, :);
fvalcond_oz_resp = fvalcond_resp(2, :);
%interaction effect of exp x flicker rate
fvalexpflick_cpz_resp = fvalexpflick_resp(1, :);
fvalexpflick_oz_resp = fvalexpflick_resp(2, :);

%stacking p and f values for visual aesthetics
%cpz
cpz_flick_tg = [panova_cpz_flick_tg; cell2mat(fvalflick_cpz_tg)];
%size = 2x32; first row is p-values & 2nd row is f-values
cpz_exp_tg = [panova_cpz_exp_tg; cell2mat(fvalexp_cpz_tg)];
cpz_cond_tg = [panova_cpz_cond_tg; cell2mat(fvalcond_cpz_tg)];
cpz_flick_resp = [panova_cpz_flick_resp; cell2mat(fvalflick_cpz_resp)];
%size = 2x32; first row is p-values & 2nd row is f-values
cpz_exp_resp = [panova_cpz_exp_resp; cell2mat(fvalexp_cpz_resp)];
cpz_cond_resp = [panova_cpz_cond_resp; cell2mat(fvalcond_cpz_resp)];
 %interaction effect
 cpz_exp_flick_tg = [panova_cpz_exp_flick_tg; cell2mat(fvalexpflick_cpz_tg)];
 cpz_exp_flick_resp = [panova_cpz_exp_flick_resp; cell2mat(fvalexpflick_cpz_resp)];
%oz
oz_flick_tg = [panova_oz_flick_tg; cell2mat(fvalflick_oz_tg)];
%size = 2x32; first row is p-values & 2nd row is f-values
oz_exp_tg = [panova_oz_exp_tg; cell2mat(fvalexp_oz_tg)];
oz_cond_tg = [panova_oz_cond_tg; cell2mat(fvalcond_oz_tg)];
oz_flick_resp = [panova_oz_flick_resp; cell2mat(fvalflick_oz_resp)];
%size = 2x32; first row is p-values & 2nd row is f-values
oz_exp_resp = [panova_oz_exp_resp; cell2mat(fvalexp_oz_resp)];
oz_cond_resp = [panova_oz_cond_resp; cell2mat(fvalcond_oz_resp)];
 %interaction effect
 oz_exp_flick_tg = [panova_oz_exp_flick_tg; cell2mat(fvalexpflick_oz_tg)];
 oz_exp_flick_resp = [panova_oz_exp_flick_resp; cell2mat(fvalexpflick_oz_resp)];

%% FDR Correction 
%combine CPz, Oz of both tg- and resp-locked
%main effects
chosentime.tgresp = [chosentime.tg(1:end-1), chosentime.resp(1:end-1)];
%made time time vector from combining chosen time in
%bins of tg locked and resp locked data; size = 1x40;
%can be used with p_fdr to quickly search sig time bins

[p_fdr_flick, p_masked_flick] = fdr([panova_cpz_flick_tg panova_oz_flick_tg ...
    panova_cpz_flick_resp panova_oz_flick_resp], 0.05, 'parametric');
[p_fdr_exp, p_masked_exp] = fdr([panova_cpz_exp_tg panova_oz_exp_tg ...
    panova_cpz_exp_resp panova_oz_exp_resp], 0.05, 'parametric'); %
[p_fdr_cond, p_masked_cond] = fdr([panova_cpz_cond_tg panova_oz_cond_tg ...
    panova_cpz_cond_resp panova_oz_cond_resp], 0.05, 'parametric'); 

%071317
[p_fdr_flickexp, p_masked_flickexp] = fdr([panova_cpz_flick_tg panova_cpz_flick_resp ...
     panova_cpz_exp_tg panova_cpz_exp_resp panova_oz_flick_tg panova_oz_flick_resp ...
     panova_oz_exp_tg panova_oz_exp_resp], 0.05, 'parametric'); 
[p_fdr_cond, p_masked_cond] = fdr([panova_cpz_cond_tg panova_cpz_cond_resp ...
     panova_oz_cond_tg panova_oz_cond_resp], 0.05, 'parametric'); 

% main_effect_cpz = [[1:32, 1:8]', chosentime.tgresp', [p_masked_flick(1:32),...
%     p_masked_flick(65:72)]', [p_masked_exp(1:32), p_masked_exp(65:72)]', ...
%     [p_masked_cond(1:32), p_masked_cond(65:72)]']; %a summary of main effects of cpz channel
% 
% main_effect_oz = [[1:32, 1:8]', chosentime.tgresp', [p_masked_flick(33:64),...
%     p_masked_flick(73:80)]', [p_masked_exp(33:64), p_masked_exp(73:80)]', ...
%     [p_masked_cond(33:64), p_masked_cond(73:80)]']; %a summary of main effects of cpz channel
% %size = 40x5 b/c we have 32+8 time points from tg locked and
% %resp locked.
% %1st row = ith time bin
% %2nd row = time in ms of the begining of each bin
% %3th, 4th, 5th row = p_masked_flicker, p_masked_exp, p_masked_cond

%same as below--just hard code to double check
% main_effect_cpz = [[1:numel(timebin.tg)-1, 1:numel(timebin.resp)-1]', ... 
%     chosentime.tgresp', [p_masked_flickexp(1:43)]'...
%     , [p_masked_flickexp(44:86)]'...
%     , [p_masked_cond(1:43)]']; 
% main_effect_oz = [[1:numel(timebin.tg)-1, 1:numel(timebin.resp)-1]', ... 
%     chosentime.tgresp', [p_masked_flickexp(87:129)]'...
%     , [p_masked_flickexp(130:172)]'...
%     , [p_masked_cond(44:86)]'];  
main_effect_cpz = [[1:numel(timebin.tg)-1, 1:numel(timebin.resp)-1]', ... 
    chosentime.tgresp', [p_masked_flickexp(1:size(chosentime.tgresp, 2))]'...
    , [p_masked_flickexp(size(chosentime.tgresp, 2)+1:(size(chosentime.tgresp, 2)*2))]'...
    , [p_masked_cond(1:size(chosentime.tgresp, 2))]'];   
main_effect_oz = [[1:numel(timebin.tg)-1, 1:numel(timebin.resp)-1]', ... 
    chosentime.tgresp', [p_masked_flickexp((size(chosentime.tgresp, 2)*2)+1:...
    size(chosentime.tgresp, 2)*3)]'...
    , [p_masked_flickexp((size(chosentime.tgresp, 2)*3)+1:size(chosentime.tgresp, 2)*4)]'...
    , [p_masked_cond(size(chosentime.tgresp, 2)+1:size(chosentime.tgresp, 2)*2)]'];  

% all_p_fdr = [p_fdr_flick, p_fdr_exp, p_fdr_cond];
all_p_fdr = [p_fdr_flickexp, p_fdr_cond];

%crit p to check interaction effects
[p_fdr_exp_cond, p_masked_exp_cond] = fdr([panova_cpz_exp_cond_tg ...
    panova_oz_exp_cond_tg panova_cpz_exp_cond_resp panova_oz_exp_cond_resp], 0.05, 'parametric');
[p_fdr_exp_flick, p_masked_exp_flick] = fdr([panova_cpz_exp_flick_tg ...
    panova_oz_exp_flick_tg panova_cpz_exp_flick_resp panova_oz_exp_flick_resp], 0.05, 'parametric');
[p_fdr_cond_flick, p_masked_cond_flick] = fdr([panova_cpz_cond_flick_tg ...
    panova_oz_cond_flick_tg panova_cpz_cond_flick_resp panova_oz_cond_flick_resp], 0.05, 'parametric');
[p_fdr_exp_cond_flick, p_masked_exp_cond_flick] = fdr([panova_cpz_exp_cond_flick_tg ...
    panova_oz_exp_cond_flick_tg panova_cpz_exp_cond_flick_resp panova_oz_exp_cond_flick_resp], 0.05, 'parametric');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%sig bins
  %flicker rate effect
  %cpz tg-locked
   sigbin_tg_flick_cpz = find(main_effect_cpz(1:size(chosentime.tg, 2)-1, 3)==1);
     sigbin_tg_flick_cpz1 = [sigbin_tg_flick_cpz(1:end-1)];
     sigbin_tg_flick_cpz2 = sigbin_tg_flick_cpz(end);
     sigtime_tg_flick_cpz1 = [chosentime.tg(sigbin_tg_flick_cpz1(1)) ...
         chosentime.tg(sigbin_tg_flick_cpz1(end)+1)];
     sigtime_tg_flick_cpz2 = [chosentime.tg(sigbin_tg_flick_cpz2(1)) ...
         chosentime.tg(sigbin_tg_flick_cpz2(end)+1)];    
  %cpz resp-locked
   sigbin_resp_flick_cpz =  find(main_effect_cpz((size(chosentime.tg, 2):end), 3)==1);
     sigbin_resp_flick_cpz1 = sigbin_resp_flick_cpz;
     sigtime_resp_flick_cpz1 = [chosentime.resp(sigbin_resp_flick_cpz1(1)) ...
         chosentime.resp(sigbin_resp_flick_cpz1(end)+1)];
    
   %oz tg-locked
   sigbin_tg_flick_oz = find(main_effect_oz(1:size(chosentime.tg, 2)-1, 3)==1);
   sigtime_tg_flick_oz = [chosentime.tg(sigbin_tg_flick_oz(1)) ...
         chosentime.tg(sigbin_tg_flick_oz(end)+1)]; 
  
  %expectation effect
  %cpz tg-locked
  sigbin_tg_exp_cpz =  find(main_effect_cpz(1:size(chosentime.tg, 2)-1, 4)==1); 
   sigbin_tg_exp_cpz1 = [sigbin_tg_exp_cpz(1:end-1)];
   sigbin_tg_exp_cpz2 = [sigbin_tg_exp_cpz(end)];
   sigtime_tg_exp_cpz1 = [chosentime.tg(sigbin_tg_exp_cpz1(1)) ...
          chosentime.tg(sigbin_tg_exp_cpz1(end)+1)];
   sigtime_tg_exp_cpz2 = [chosentime.tg(sigbin_tg_exp_cpz2(1)) ...
          chosentime.tg(sigbin_tg_exp_cpz2(end)+1)];  
      
%% CPP & Expectation: Post-Hoc T-Tests
%do this only for time windows with significant expectation main effect to
 %see if what pair of exp conditions significantly differ 

%a) time chunk 1: 950-1100 ms
cpz_tg_exp1 = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 1, :, 2, eoi(ch_cpz, :),...
    sigbin_tg_exp_cpz1), 7), 4), 2));
cpz_tg_neu1 = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 2, :, 2, eoi(ch_cpz, :),...
    sigbin_tg_exp_cpz1), 7), 4), 2));
cpz_tg_unexp1 = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 3, :, 2, eoi(ch_cpz, :),...
    sigbin_tg_exp_cpz1), 7), 4), 2));

%1) CPP: exp vs. unexp
clear h p ci t
[h p ci t] = ttest(cpz_tg_exp1, cpz_tg_unexp1, 'Tail', 'left'); %we're using one-tailed tests as follow-ups
cpz_tg_expunexp1.p = p;
cpz_tg_expunexp1.h = h;
cpz_tg_expunexp1.t = t.tstat;
cpz_tg_expunexp1.df = t.df;

%2) CPP: exp vs. neu
clear h p ci t
[h p ci t] = ttest(cpz_tg_exp1, cpz_tg_neu1, 'Tail', 'left');
cpz_tg_expneu1.p = p;
cpz_tg_expneu1.h = h;
cpz_tg_expneu1.t = t.tstat;
cpz_tg_expneu1.df = t.df;

%3) CPP: neu vs. unexp
clear h p ci t
[h p ci t] = ttest(cpz_tg_neu1, cpz_tg_unexp1, 'Tail', 'left');
cpz_tg_neuunexp1.p = p;
cpz_tg_neuunexp1.h = h;
cpz_tg_neuunexp1.t = t.tstat;
cpz_tg_neuunexp1.df = t.df;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%b) time chunk 1: 1150-1200 ms
cpz_tg_exp2 = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 1, :, 2, eoi(ch_cpz, :),...
    sigbin_tg_exp_cpz2), 7), 4), 2));
cpz_tg_neu2 = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 2, :, 2, eoi(ch_cpz, :),...
    sigbin_tg_exp_cpz2), 7), 4), 2));
cpz_tg_unexp2 = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 3, :, 2, eoi(ch_cpz, :),...
    sigbin_tg_exp_cpz2), 7), 4), 2));

%1) CPP: exp vs. unexp
clear h p ci t
[h p ci t] = ttest(cpz_tg_exp2, cpz_tg_unexp2, 'Tail', 'left'); %we're using one-tailed tests as follow-ups
cpz_tg_expunexp2.p = p;
cpz_tg_expunexp2.h = h;
cpz_tg_expunexp2.t = t.tstat;
cpz_tg_expunexp2.df = t.df;

%2) CPP: exp vs. neu
clear h p ci t
[h p ci t] = ttest(cpz_tg_exp2, cpz_tg_neu2, 'Tail', 'left');
cpz_tg_expneu2.p = p;
cpz_tg_expneu2.h = h;
cpz_tg_expneu2.t = t.tstat;
cpz_tg_expneu2.df = t.df;

%3) CPP: neu vs. unexp
clear h p ci t
[h p ci t] = ttest(cpz_tg_neu2, cpz_tg_unexp2, 'Tail', 'left');
cpz_tg_neuunexp2.p = p;
cpz_tg_neuunexp2.h = h;
cpz_tg_neuunexp2.t = t.tstat;
cpz_tg_neuunexp2.df = t.df;

%% CPP & Expectation: Follow-up one-way ANOVA on individual task condition (sig chunks only)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%do this on data collapsed across sig time chunk separately for each
%individual task condition to see if the main expectation effect 
%(tested via ANOVA on data collapsed across 3 task conditions) 

%tg-locked
%a) time chunk one: 920-1100 ms
for cond = 1:3
    cond_exp1 = squeeze(mean(mean(gerp_tgbinned(:, cond, 1, :, 2, eoi(ch_cpz, :),...
        sigbin_tg_exp_cpz1), 7), 4));
    cond_neu1 = squeeze(mean(mean(gerp_tgbinned(:, cond, 2, :, 2, eoi(ch_cpz, :),...
        sigbin_tg_exp_cpz1), 7), 4));
    cond_unexp1 = squeeze(mean(mean(gerp_tgbinned(:, cond, 3, :, 2, eoi(ch_cpz, :),...
        sigbin_tg_exp_cpz1), 7), 4));
    
    DV1 = [cond_exp1; cond_neu1; cond_unexp1];
    sublist = repmat(1:numel(sub), [1, 3])'; %subject list
    IVexp = repmat([ones(1, numel(sub)) ones(1, numel(sub)).*2 ...
        ones(1, numel(sub)).*3], [1, 1])';
    
    [p, tbl, stats] = anovan(DV1, {IVexp, sublist}, 'model', 'full', 'varnames',...
        {'exp/neu/unexp', 'subj'}, 'random', [2]);
    p_indivcond_exp_tg1(cond, :) = p;
    fvalue_indivcond_exp_tg1(cond, :) = tbl(2, 6);    
end

%b) time chunk one: 1150-1200 ms
for cond = 1:3
    cond_exp2 = squeeze(mean(mean(gerp_tgbinned(:, cond, 1, :, 2, eoi(ch_cpz, :),...
        sigbin_tg_exp_cpz2), 7), 4));
    cond_neu2 = squeeze(mean(mean(gerp_tgbinned(:, cond, 2, :, 2, eoi(ch_cpz, :),...
        sigbin_tg_exp_cpz2), 7), 4));
    cond_unexp2 = squeeze(mean(mean(gerp_tgbinned(:, cond, 3, :, 2, eoi(ch_cpz, :),...
        sigbin_tg_exp_cpz2), 7), 4));
    
    DV2 = [cond_exp2; cond_neu2; cond_unexp2];
    sublist = repmat(1:numel(sub), [1, 3])'; %subject list
    IVexp = repmat([ones(1, numel(sub)) ones(1, numel(sub)).*2 ...
        ones(1, numel(sub)).*3], [1, 1])';
    
    [p, tbl, stats] = anovan(DV2, {IVexp, sublist}, 'model', 'full', 'varnames',...
        {'exp/neu/unexp', 'subj'}, 'random', [2]);
    p_indivcond_exp_tg2(cond, :) = p;
    fvalue_indivcond_exp_tg2(cond, :) = tbl(2, 6);    
end   

%% CPP & Flicker rate: Follow-up t-test on individual task condition (sig chunks only)
%1) tg-locked
for cond = 1:3
    cpz_tg_fast1 = squeeze(mean(mean(gerp_tgbinned(:, cond, [1 2 3], 1, 2, ...
        eoi(ch_cpz, :), sigbin_tg_flick_cpz1), 7), 3)); %time window 1
    cpz_tg_slow1 = squeeze(mean(mean(gerp_tgbinned(:, cond, [1 2 3], 2, 2, ...
        eoi(ch_cpz, :), sigbin_tg_flick_cpz1), 7), 3));
    
    cpz_tg_fast2 = squeeze(mean(mean(gerp_tgbinned(:, cond, [1 2 3], 1, 2, ...
        eoi(ch_cpz, :), sigbin_tg_flick_cpz2), 7), 3)); %time window 2
    cpz_tg_slow2 = squeeze(mean(mean(gerp_tgbinned(:, cond, [1 2 3], 2, 2, ...
        eoi(ch_cpz, :), sigbin_tg_flick_cpz2), 7), 3));
    
    clear h p ci t
    [h p ci t] = ttest(cpz_tg_fast1, cpz_tg_slow1, 'Tail', 'right');
    cpz_tg_flick1.p(cond) = p;
    cpz_tg_flick1.h(cond) = h;
    cpz_tg_flick1.t(cond) = t.tstat;
    cpz_tg_flick1.df(cond) = t.df;
    
    clear h p ci t
    [h p ci t] = ttest(cpz_tg_fast2, cpz_tg_slow2, 'Tail', 'left');
    cpz_tg_flick2.p(cond) = p;
    cpz_tg_flick2.h(cond) = h;
    cpz_tg_flick2.t(cond) = t.tstat;
    cpz_tg_flick2.df(cond) = t.df;
end

%2) resp-locked
for cond = 1:3
    cpz_resp_fast1 = squeeze(mean(mean(gerp_respbinned(:, cond, [1 2 3], 1, 2, ...
        eoi(ch_cpz, :), sigbin_resp_flick_cpz1), 7), 3));
    cpz_resp_slow1 = squeeze(mean(mean(gerp_respbinned(:, cond, [1 2 3], 2, 2, ...
        eoi(ch_cpz, :), sigbin_resp_flick_cpz1), 7), 3));
    
    
    clear h p ci t
    [h p ci t] = ttest(cpz_resp_fast1, cpz_resp_slow1, 'Tail', 'right');
    cpz_resp_flick1.p(cond) = p;
    cpz_resp_flick1.h(cond) = h;
    cpz_resp_flick1.t(cond) = t.tstat;
    cpz_resp_flick1.df(cond) = t.df;
    
end

%% VN & Flicker rate: Follow-up t-test on individual task condition (sig chunks only)
%1) tg-locked
for cond = 1:3
    oz_tg_fast = squeeze(mean(mean(gerp_tgbinned(:, cond, [1 2 3], 1, 2, ...
        eoi(ch_oz, :), sigbin_tg_flick_oz), 7), 3));
    oz_tg_slow = squeeze(mean(mean(gerp_tgbinned(:, cond, [1 2 3], 2, 2, ...
        eoi(ch_oz, :), sigbin_tg_flick_oz), 7), 3));
    
    clear h p ci t
    [h p ci t] = ttest(oz_tg_fast, oz_tg_slow, 'Tail', 'left');
    oz_tg_flick.p(cond) = p;
    oz_tg_flick.h(cond) = h;
    oz_tg_flick.t(cond) = t.tstat;
    oz_tg_flick.df(cond) = t.df;
end

%% stacking all p- and F-values
 % p-value and F-value ranges for the flicker rate ane expectation chunks
 % on collapsed data
 
 %flicker rate effect
 %cpz tg-locked chunk#1
 sig_tg_flick_cpz1 = cpz_flick_tg(:, sigbin_tg_flick_cpz1);
 minmax_tg_flick_cpz1 = [sig_tg_flick_cpz1(:, find(sig_tg_flick_cpz1(1,:) == ...
     min(sig_tg_flick_cpz1(1, :)))), sig_tg_flick_cpz1(:, find(sig_tg_flick_cpz1(1,:)...
     == max(sig_tg_flick_cpz1(1, :))))];
 %cpz tg-locked chunk#2
 sig_tg_flick_cpz2 = cpz_flick_tg(:, sigbin_tg_flick_cpz2);
 minmax_tg_flick_cpz2 = [sig_tg_flick_cpz2(:, find(sig_tg_flick_cpz2(1,:) == ...
     min(sig_tg_flick_cpz2(1, :)))), sig_tg_flick_cpz2(:, find(sig_tg_flick_cpz2(1,:)...
     == max(sig_tg_flick_cpz2(1, :))))];
  %cpz resp-locked chunk#1
 sig_resp_flick_cpz1 = cpz_flick_resp(:, sigbin_resp_flick_cpz1);
 minmax_resp_flick_cpz1 = [sig_resp_flick_cpz1(:, find(sig_resp_flick_cpz1(1,:) == ...
     min(sig_resp_flick_cpz1(1, :)))), sig_resp_flick_cpz1(:, find(sig_resp_flick_cpz1(1,:)...
     == max(sig_resp_flick_cpz1(1, :))))];

  %oz tg-locked 
 sig_tg_flick_oz = cpz_flick_tg(:, sigbin_tg_flick_oz);
 minmax_tg_flick_oz = [sig_tg_flick_oz(:, find(sig_tg_flick_oz(1,:) == ...
     min(sig_tg_flick_oz(1, :)))), sig_tg_flick_oz(:, find(sig_tg_flick_oz(1,:)...
     == max(sig_tg_flick_oz(1, :))))];
 
 %expectation effect
 %cpz tg-locked 
 sig_tg_exp_cpz1 = cpz_exp_tg(:, sigbin_tg_exp_cpz1);
 minmax_tg_exp_cpz1 = [sig_tg_exp_cpz1(:, find(sig_tg_exp_cpz1(1,:) == ...
     min(sig_tg_exp_cpz1(1, :)))), sig_tg_exp_cpz1(:, find(sig_tg_exp_cpz1(1,:)...
     == max(sig_tg_exp_cpz1(1, :))))];
 
 sig_tg_exp_cpz2 = cpz_exp_tg(:, sigbin_tg_exp_cpz2);
 minmax_tg_exp_cpz2 = [sig_tg_exp_cpz2(:, find(sig_tg_exp_cpz2(1,:) == ...
     min(sig_tg_exp_cpz2(1, :)))), sig_tg_exp_cpz2(:, find(sig_tg_exp_cpz2(1,:)...
     == max(sig_tg_exp_cpz2(1, :))))];
 
 %NULL EFFECTS
 %expectation--oz tg-locked
 minmax_tg_exp_oz = [oz_exp_tg(:, find(oz_exp_tg(1,:) == ...
     min(oz_exp_tg(1, :)))), oz_exp_tg(:, find(oz_exp_tg(1,:)...
     == max(oz_exp_tg(1, :))))];
 %expectation--oz resp-locked
%  minmax_resp_exp_oz = [oz_exp_resp(:, find(oz_exp_resp(1,:) == ...
%      min(oz_exp_resp(1, :)))), oz_exp_resp(:, find(oz_exp_resp(1,:)...
%      == max(oz_exp_resp(1, :))))];
 
 %interaction effect: flicker rate x exp
 %cpz tg-locked
 minmax_tg_expflick_cpz = [cpz_exp_flick_tg(:, find(cpz_exp_flick_tg(1,:) == ...
     min(cpz_exp_flick_tg(1, :)))), cpz_exp_flick_tg(:, find(cpz_exp_flick_tg(1,:)...
     == max(cpz_exp_flick_tg(1, :))))];
 %cpz resp-locked
 minmax_resp_expflick_cpz = [cpz_exp_flick_resp(:, find(cpz_exp_flick_resp(1,:) == ...
     min(cpz_exp_flick_resp(1, :)))), cpz_exp_flick_resp(:, find(cpz_exp_flick_resp(1,:)...
     == max(cpz_exp_flick_resp(1, :))))];
 %oz tg-locked
 minmax_tg_expflick_oz = [oz_exp_flick_tg(:, find(oz_exp_flick_tg(1,:) == ...
     min(oz_exp_flick_tg(1, :)))), oz_exp_flick_tg(:, find(oz_exp_flick_tg(1,:)...
     == max(oz_exp_flick_tg(1, :))))];
 %oz resp-locked
%  minmax_resp_expflick_oz = [oz_exp_flick_resp(:, find(oz_exp_flick_resp(1,:) == ...
%      min(oz_exp_flick_resp(1, :)))), oz_exp_flick_resp(:, find(oz_exp_flick_resp(1,:)...
%      == max(oz_exp_flick_resp(1, :))))];
 
%% Figures stuff
close all force;

%set gca stuff
%xlim_tg = [-75 1475];
%xtick_tg = [0:500:1475];
xlim_tg = [-250 1475];
xtick_tg = [-250 0 500 1000 1475];

xlim_resp = [-275 1275];
xtick_resp = [-200 0 75];

ylim_cpz = [-1 10];
ytick_cpz = [0 4 8];
starloc_cpz = 9;

ylim_oz = [-5 5];
ytick_oz = [-4 0 4];
starloc_oz = 3;

flim = [0 50];
ftick = [0:50:50];

plim = [0 0.05];
ptick = [0:0.05:0.05];

%color values for plotting
c_fast = [0 0.4 0]; %green
c_slow = [0.4 0 0.6]; %purple
c_exp = [0 0.4 0.8];
c_neu = [0.5 0.5 0.5]; %grey
c_unexp = [0.8 0 0]; %orange

trans = 0.8; %transparency value for banded error
lwidth = 2.5;
lwidth_sub =1.5;

%% CPP_Fig A left: Flicker rates (collapsed)
close all;
 for fig = 1
    figure(1); clf;
    suptitle('FigA left: CPP & flicker rate (collapsed)')
    
    subplot(3, 2, 1)
    %fast flicker rate
    a1 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, :, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))...
        , 'Color', c_fast, 'LineWidth', lwidth); hold on;
    h1 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, :, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_flicker...
        (:, :, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        a1, trans);   

    %slow flicker rate (33Hz)
    a2 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
         (:, :, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))...
         , 'Color', c_slow, 'LineWidth', lwidth); hold on;
    h2 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
         (:, :, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))', ...
         squeeze(std(mean(mean(mean(fig_gerp_tg_flicker...
        (:, :, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        a2, trans);
    
    xlim (xlim_tg);
    ylim (ylim_cpz);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ytick_cpz)
    %title(chan_label(ch), 'FontSize', 15);
    %title('Tg-locked', 'FontSize', 10);
    ylabel('Amplitude (uv)', 'FontSize', 14)
    legend([a1 a2],'Fast Flicker Rate', 'Slow Flicker Rate');
    
    %star plot
    for tt =1:numel(timebin.tg)-1
        if panova_tg(1, tt, 3) > p_fdr_flickexp
            colorsign = [1 1 1];
%         elseif panova(cnt, tt, 3) <= 0.05 &&  panova(cnt, tt, 3) > 0.01
%             colorsign = [0.6 0.6 0.6];
%         elseif panova(cnt, tt, 3) <= 0.01 &&  panova(cnt, tt, 3) > 0.001
%             colorsign = [0.3 0.3 0.3];
        else
            colorsign = [0 0 0];
        end
        
        if panova_tg(1, tt, 3) <= p_fdr_flickexp
            plot(chosentime.tg(tt)+binsize/2, starloc_cpz, '.', 'color', colorsign, 'LineWidth', 10);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3, 2, 2) %resp-locked
    %fast flicker rate 
    a3 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
           (:, :, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))...
        , 'Color', c_fast, 'LineWidth', lwidth); hold on;
    h3 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
           (:, :, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_flicker...
        (:, :, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        a3, trans); 
    
    %slow flicker rate (33Hz)
    a4 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, :, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))...
        , 'Color', c_slow, 'LineWidth', lwidth); hold on;
    h4 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, :, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_flicker...
        (:, :, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        a4, trans); 
     
    %title(chan_label(ch), 'FontSize', 15);
    xlim (xlim_resp);
    ylim (ylim_cpz);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ytick_cpz)
%     title('Resp-locked', 'FontSize', 10);
    ylabel('Amplitude (uv)', 'FontSize', 14)
    legend([a3 a4],'Fast Flicker Rate', 'Slow Flicker Rate');
    
    %star plot
    for tt =1:numel(timebin.resp)-1
        if panova_resp(1, tt, 3) > p_fdr_flickexp
            colorsign = [1 1 1];
            %         elseif panova(cnt, tt, 3) <= 0.05 &&  panova(cnt, tt, 3) > 0.01
            %             colorsign = [0.6 0.6 0.6];
            %         elseif panova(cnt, tt, 3) <= 0.01 &&  panova(cnt, tt, 3) > 0.001
            %             colorsign = [0.3 0.3 0.3];
        else
            colorsign = [0 0 0];
        end
        
        if panova_resp(1, tt, 3) <=  p_fdr_flick
            
            plot(chosentime.resp(tt)+binsize/2, starloc_cpz, '.', 'color', colorsign, 'LineWidth', 10);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3, 2, 3)
    plot(chosentime.tg(1:end-1)+binsize/2, cell2mat(fvalflick_cpz_tg),...
        'color', 'k', 'Linewidth', lwidth_sub);
    ylabel('F', 'FontSize', 14)
    xlim (xlim_tg);
    ylim(flim)
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ftick)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3, 2, 4)
    plot(chosentime.resp(1:end-1)+binsize/2, cell2mat(fvalflick_cpz_resp),...
        'color', 'k', 'Linewidth', lwidth_sub);
    ylabel('F', 'FontSize', 14)
    xlim (xlim_resp);
    ylim(flim)
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ftick)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3, 2, 5)
    plot(chosentime.tg(1:end-1)+binsize/2, panova_tg(1, :, 3),...
        'color', 'k', 'Linewidth', lwidth_sub); hold on;
    plot(timex.tg, repmat(p_fdr_flickexp, [1, size(timex.tg, 2)]),...
       'color', 'k', 'LineStyle', '--'); hold on;
    ylabel('p', 'FontSize', 14)
    xlim (xlim_tg);
    ylim(plim);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ptick)
    xlabel('Time (ms)', 'FontSize', 14)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3, 2, 6)
    plot(chosentime.resp(1:end-1)+binsize/2, panova_resp(1, :, 3),...
        'color', 'k', 'Linewidth', lwidth_sub); hold on;
    plot(timex.resp, repmat(p_fdr_flickexp, [1, size(timex.resp, 2)]),...
        'color', 'k', 'LineStyle', '--'); hold on;
    ylabel('p', 'FontSize', 14)
    xlim (xlim_resp);
    ylim(plim);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ptick)
    xlabel('Time (ms)', 'FontSize', 14)
end
 %topo plot inset
 for fig = 2
     figure(2); clf;
     subplot(2, 2, 1) 
     %tg-locked
     %1)
     flick_topolength_tg_cpz1 = [sigbin_tg_flick_cpz1', sigbin_tg_flick_cpz1(end)+1]; %duration of sig flick effect 
                                                      %purposes
     flick_topotime_tg_cpz1 = timebin.tg(flick_topolength_tg_cpz1);                                                     
     topoplot(squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
          (:, :, :, 1, 2, :, flick_topotime_tg_cpz1),7),3),2)))...
          -squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
          (:, :, :, 2, 2, :, flick_topotime_tg_cpz1),7),3),2))), 'mychan.loc'); %fast-slow
     colorbar
     caxis([-1 1])
     title('1) Tg-locked CPP: Fast - Slow Flicker Rate');
     
     subplot(2, 2, 2)
     %2)
     flick_topolength_tg_cpz2 = [sigbin_tg_flick_cpz2', sigbin_tg_flick_cpz2(end)+1]; %duration of sig flick effect 
                                                      %purposes
     flick_topotime_tg_cpz2 = timebin.tg(flick_topolength_tg_cpz2);                                                     
     topoplot(squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
          (:, :, :, 1, 2, :, flick_topotime_tg_cpz2),7),3),2)))...
          -squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
          (:, :, :, 2, 2, :, flick_topotime_tg_cpz2),7),3),2))), 'mychan.loc'); %fast-slow
     colorbar
     caxis([-1 1])
     title('2) Tg-locked CPP: Fast - Slow Flicker Rate');
     
     subplot(2, 2, 3) 
     %resp-locked
     %1)
     flick_topolength_resp_cpz1 = [sigbin_resp_flick_cpz1', sigbin_resp_flick_cpz1(end)+1]; %duration of sig flick effect 
                                                      %purposes
     flick_topotime_resp_cpz1 = timebin.resp(flick_topolength_resp_cpz1);                                                     
     topoplot(squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
          (:, :, :, 1, 2, :, flick_topotime_resp_cpz1),7),3),2)))...
          -squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
          (:, :, :, 2, 2, :, flick_topotime_resp_cpz1),7),3),2))), 'mychan.loc'); %fast-slow
     colorbar
     caxis([-1 1])
     title('1) Resp-locked CPP: Fast - Slow Flicker Rate');
     
%      subplot(2, 2, 4)
%      %2)
%      flick_topolength_resp_cpz2 = [sigbin_resp_flick_cpz2', sigbin_resp_flick_cpz2(end)+1]; %duration of sig flick effect 
%                                                       %purposes
%      flick_topotime_resp_cpz2 = timebin.resp(flick_topolength_resp_cpz2);                                                     
%      topoplot(squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
%           (:, :, :, 1, 2, :, flick_topotime_resp_cpz2),7),3),2)))...
%           -squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
%           (:, :, :, 2, 2, :, flick_topotime_resp_cpz2),7),3),2))), 'mychan.loc'); %fast-slow
%      colorbar
%      caxis([-1 1])
%      title('2) Resp-locked CPP: Fast - Slow Flicker Rate');
     
 end
 
%% CPP_Fig A right: Flicker rates (individual conditions)
 for fig = 3
clear tt
figure(3); clf;
for bc = 1:3
    subplot(3, 2, ((bc-1)*2)+1); 
    %target-locked
    %fast flicker rate
    b1 =  plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))...
        , 'Color', c_fast, 'LineWidth', lwidth); hold on;
    k1 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        b1, trans);
    
    %slow flicker rate
    b2 =  plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))...
        , 'Color', c_slow, 'LineWidth', lwidth); hold on;
    k2 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        b2, trans);

    %title(chan_label(ch), 'FontSize', 15);
    ylabel('Amplitude (uv)', 'FontSize', 14);
    xlim (xlim_tg);
    ylim (ylim_cpz);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ytick_cpz)
    
    subplot(3, 2, 2*bc); %resp-locked
    %fast flicker rate
    b3 =  plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))...
        , 'Color', c_fast, 'LineWidth', lwidth); hold on;
    k3 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        b3, trans);
    
    %slow flicker rate
    b4 =  plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))...
        , 'Color', c_slow, 'LineWidth', lwidth); hold on;
    k4 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        b4, trans);    

    xlim (xlim_resp);
    ylim (ylim_cpz);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ytick_cpz)
    xlabel('Time (ms)', 'FontSize', 14)

end
legend([b3 b4],'Fast Flicker Rate', 'Slow Flicker Rate');
suptitle('Fig A right: CPP & Flicker Rates (non-collapsed)')
 end

%% CPP_Fig B left: Expectation (collapsed)
for fig = 4
    figure(4); clf;
    
    suptitle('Fig B left: CPP & expectation (collapsed)')

    subplot(3, 2, 1) %target-locked
    %expected
    c1 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_exp, 'LineWidth', lwidth); hold on;
    m1 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c1, trans);     
    
    %neutral
    c2 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_neu, 'LineWidth', lwidth); hold on;
    m2 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c2, trans);
    
    %unexpected
    c3 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_unexp, 'LineWidth', lwidth); hold on;
    m3 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c3, trans);
    
    %title(chan_label(ch), 'FontSize', 15);
    ylabel('Amplitude (uv)', 'FontSize', 14)
    xlim (xlim_tg);
    ylim (ylim_cpz);
%     title('Tg-locked', 'FontSize', 10);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ytick_cpz)
    
    %star plot
    for tt =1:numel(timebin.tg)-1
        if panova_tg(1, tt, 1) > p_fdr_flickexp
            colorsign = [1 1 1];
%         elseif panova(cnt, tt, 1) <= 0.05 &&  panova(cnt, tt, 1) > 0.01
%             colorsign = [0.6 0.6 0.6];
%         elseif panova(cnt, tt, 1) <= 0.01 &&  panova(cnt, tt, 1) > 0.001
%             colorsign = [0.3 0.3 0.3];
        else
            colorsign = [0 0 0];
        end
        
        if panova_tg(1, tt, 1) <= p_fdr_flickexp
            plot(chosentime.tg(tt)+binsize/2, starloc_cpz, '.', 'color', colorsign, 'LineWidth', 10); 
        end
    end
    
    subplot(3, 2, 2) %response-locked 
    %expected
    c4 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_exp, 'LineWidth', lwidth); hold on;
    m4 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c4, trans);
    
    %neutral
    c5 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_neu, 'LineWidth', lwidth); hold on;
    m5 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c5, trans);
    
    %unexpected
    c6 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_unexp, 'LineWidth', lwidth); hold on;
    m6 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c6, trans);
    
    %title(chan_label(ch), 'FontSize', 15);
    ylabel('Amplitude (uv)', 'FontSize', 14)
%     title('Resp-locked', 'FontSize', 10);
    xlim (xlim_resp);
    ylim (ylim_cpz);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ytick_cpz)
    legend([c4 c5 c6], 'Expected', 'Neutral', 'Unexpected', 'FontSize', 14)
    
    %star plot
    for tt =1:numel(timebin.resp)-1
        if panova_cpz_exp_resp(1, tt) > p_fdr_flickexp
            colorsign = [1 1 1];
%         elseif panova(cnt, tt, 1) <= 0.05 &&  panova(cnt, tt, 1) > 0.01
%             colorsign = [0.6 0.6 0.6];
%         elseif panova(cnt, tt, 1) <= 0.01 &&  panova(cnt, tt, 1) > 0.001
%             colorsign = [0.3 0.3 0.3];
        else
            colorsign = [0 0 0];
        end
        
        if panova_cpz_exp_resp(1, tt) <= p_fdr_flickexp
            plot(chosentime.resp(tt)+binsize/2, starloc_cpz, '.', 'color', colorsign, 'LineWidth', 10);
        end
    end
     
    %F-value plot
    subplot(3, 2, 3)
    plot(chosentime.tg(1:end-1)+binsize/2, cell2mat(fvalexp_cpz_tg), ...
        'color', 'k', 'Linewidth', lwidth_sub);
    ylabel('F', 'FontSize', 14)
    xlim (xlim_tg);
    ylim(flim);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ftick)
    
    
    subplot(3, 2, 4)
    plot(chosentime.resp(1:end-1)+binsize/2, cell2mat(fvalexp_cpz_resp), ...
        'color', 'k', 'Linewidth', lwidth_sub);
    ylabel('F', 'FontSize', 14)
    xlim (xlim_resp);
    ylim(flim);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ftick)
    
    %p-value plot
    subplot(3, 2, 5)
    plot(chosentime.tg(1:end-1)+binsize/2, panova_tg(1, :, 1), ...
        'color', 'k', 'Linewidth', lwidth_sub); hold on;
    plot(timex.tg, repmat(p_fdr_flickexp, [1, size(timex.tg, 2)]), ...
        'color', 'k', 'LineStyle', '--');
    ylabel('p', 'FontSize', 14)
    xlim (xlim_tg);
    ylim(plim);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ptick)
    xlabel('Time (ms)', 'FontSize', 14)
    
    subplot(3, 2, 6)
    plot(chosentime.resp(1:end-1)+binsize/2, panova_cpz_exp_resp(1, :), ...
        'color', 'k', 'Linewidth', lwidth_sub); hold on;
    plot(timex.resp, repmat(p_fdr_flickexp, [1, size(timex.resp, 2)]), ...
        'color', 'k', 'LineStyle', '--');
    ylabel('p', 'FontSize', 14)
    xlim (xlim_resp);
    ylim(plim);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ptick)
    xlabel('Time (ms)', 'FontSize', 14)

end
 %topo plot inset
for fig = 5
    figure(5); clf;
    subplot(2, 1, 1)%tg-locked                                                   
     exp_topolength_tg_cpz1 = [sigbin_tg_exp_cpz1', sigbin_tg_exp_cpz1(end)+1];
     
     exp_topotime_tg_cpz1 = timebin.tg(exp_topolength_tg_cpz1);
     
     
     topoplot(squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
              (:, :, 3, :, 2, :, exp_topotime_tg_cpz1),7),4),2)))...
             -squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
              (:, :, 1, :, 2, :, exp_topotime_tg_cpz1),7),4),2))), 'mychan.loc'); %exp-unexp
     
          subplot(2, 1, 2) %tg-locked chuck #2
          exp_topotime_tg_cpz2 = timebin.tg(exp_topolength_tg_cpz2);
exp_topolength_tg_cpz2 = [sigbin_tg_exp_cpz2', sigbin_tg_exp_cpz2(end)+1];
     
     topoplot(squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
              (:, :, 3, :, 2, :, exp_topotime_tg_cpz2),7),4),2)))...
             -squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
              (:, :, 1, :, 2, :, exp_topotime_tg_cpz2),7),4),2))), 'mychan.loc'); %exp-unexp
         
     
     colorbar
     caxis([-1 1])
     title('Tg-locked CPP: Exp-Unexp'); 
end

%% CPP_Fig B right: Expectation (individual conditions)
for fig = 6
    figure(6); clf;
    suptitle('Fig B right: CPP & expectation (non-collapsed)')
    
    for bc = 1:3
    subplot(3, 2, ((bc-1)*2)+1); %target-locked
    %expected
    d1 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_exp, 'LineWidth', lwidth); hold on;
    n1 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d1, trans);
    
    %neutral
    d2 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_neu, 'LineWidth', lwidth); hold on;
    n2 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d2, trans);    
    
    %unexpected
    d3 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_unexp, 'LineWidth', lwidth); hold on;
    n3 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d3, trans);
    
    %title(chan_label(ch), 'FontSize', 15);
    ylabel('Amplitude (uv)', 'FontSize', 14)
    xlim (xlim_tg);
    ylim (ylim_cpz);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ytick_cpz)
    
    subplot(3, 2, 2*bc) %resp-locked
    %expected
    d4 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_exp, 'LineWidth', lwidth); hold on;
    n4 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d4, trans);
    
    %neutral
    d5 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_neu, 'LineWidth', lwidth); hold on;
    n5 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d5, trans);
    
    %unexpected
    d6 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))...
        , 'Color', c_unexp, 'LineWidth', lwidth); hold on;
    n6 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d6, trans);
    
    %title(chan_label(ch), 'FontSize', 15);
    ylabel('Amplitude (uv)', 'FontSize', 14)
    xlim (xlim_resp);
    ylim (ylim_cpz);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ytick_cpz)   
    xlabel('Time (ms)', 'FontSize', 14)
    end
    legend([d4 d5 d6], 'Expected', 'Neutral', 'Unexpected', 'FontSize', 14)
end

%% VN_Fig A left: Flicker rates (collapsed)
clear ch

for fig = 7
    figure(7); clf;
    suptitle('Fig A left: VN & flicker rate (collapsed)')
    subplot(3, 2, 1)
    %fast flicker rate
    a1 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, :, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))...
        , 'Color', c_fast, 'LineWidth', lwidth); hold on;
    h1 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, :, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_flicker...
        (:, :, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        a1, trans);   

    %slow flicker rate (33Hz)
    a2 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
         (:, :, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))...
         , 'Color', c_slow, 'LineWidth', lwidth); hold on;
    h2 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
         (:, :, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))', ...
         squeeze(std(mean(mean(mean(fig_gerp_tg_flicker...
        (:, :, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        a2, trans);
    
    xlim (xlim_tg);
    ylim (ylim_oz);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ytick_oz)
    %title(chan_label(ch), 'FontSize', 15);
    %title('Tg-locked', 'FontSize', 10);
    ylabel('Amplitude (uv)', 'FontSize', 14)
    legend([a1 a2],'Fast Flicker Rate', 'Slow Flicker Rate');
    
    %star plot
    for tt =1:numel(timebin.tg)-1
        if panova_tg(2, tt, 3) > p_fdr_flickexp
            colorsign = [1 1 1];
%         elseif panova(cnt, tt, 3) <= 0.05 &&  panova(cnt, tt, 3) > 0.01
%             colorsign = [0.6 0.6 0.6];
%         elseif panova(cnt, tt, 3) <= 0.01 &&  panova(cnt, tt, 3) > 0.001
%             colorsign = [0.3 0.3 0.3];
        else
            colorsign = [0 0 0];
        end
        
        if panova_tg(2, tt, 3) <= p_fdr_flickexp
            plot(chosentime.tg(tt)+binsize/2, starloc_oz, '.', 'color', colorsign, 'LineWidth', 10);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3, 2, 2) %resp-locked
    %fast flicker rate 
    a3 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
           (:, :, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))...
        , 'Color', c_fast, 'LineWidth', lwidth); hold on;
    h3 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
           (:, :, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_flicker...
        (:, :, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        a3, trans); 
    
    %slow flicker rate (33Hz)
    a4 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, :, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))...
        , 'Color', c_slow, 'LineWidth', lwidth); hold on;
    h4 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, :, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_flicker...
        (:, :, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        a4, trans); 
     
    %title(chan_label(ch), 'FontSize', 15);
    xlim (xlim_resp);
    ylim (ylim_oz);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ytick_oz)
%     title('Resp-locked', 'FontSize', 10);
    ylabel('Amplitude (uv)', 'FontSize', 14)
    legend([a3 a4],'Fast Flicker Rate', 'Slow Flicker Rate');
    
    %star plot
    for tt =1:numel(timebin.resp)-1
        if panova_resp(2, tt, 3) > p_fdr_flickexp
            colorsign = [1 1 1];
            %         elseif panova(cnt, tt, 3) <= 0.05 &&  panova(cnt, tt, 3) > 0.01
            %             colorsign = [0.6 0.6 0.6];
            %         elseif panova(cnt, tt, 3) <= 0.01 &&  panova(cnt, tt, 3) > 0.001
            %             colorsign = [0.3 0.3 0.3];
        else
            colorsign = [0 0 0];
        end
        
        if panova_resp(2, tt, 3) <=  p_fdr_flickexp
            
            plot(chosentime.resp(tt)+binsize/2, starloc_oz, '.', 'color', colorsign, 'LineWidth', 10);
            
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3, 2, 3)
    plot(chosentime.tg(1:end-1)+binsize/2, cell2mat(fvalflick_oz_tg),...
        'color', 'k', 'Linewidth', lwidth_sub);
    ylabel('F', 'FontSize', 14)
    xlim (xlim_tg);
    ylim(flim)
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ftick)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3, 2, 4)
    plot(chosentime.resp(1:end-1)+binsize/2, cell2mat(fvalflick_oz_resp),...
        'color', 'k', 'Linewidth', lwidth_sub);
    ylabel('F', 'FontSize', 14)
    xlim (xlim_resp);
    ylim(flim)
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ftick)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3, 2, 5)
    plot(chosentime.tg(1:end-1)+binsize/2, panova_tg(2, :, 3),...
        'color', 'k', 'Linewidth', lwidth_sub); hold on;
    plot(timex.tg, repmat(p_fdr_flickexp, [1, size(timex.tg, 2)]),...
       'color', 'k', 'LineStyle', '--'); hold on;
    ylabel('p', 'FontSize', 14)
    xlim (xlim_tg);
    ylim(plim);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ptick)
    xlabel('Time (ms)', 'FontSize', 14)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subplot(3, 2, 6)
    plot(chosentime.resp(1:end-1)+binsize/2, panova_resp(2, :, 3),...
        'color', 'k', 'Linewidth', lwidth_sub); hold on;
    plot(timex.resp, repmat(p_fdr_flickexp, [1, size(timex.resp, 2)]),...
        'color', 'k', 'LineStyle', '--'); hold on;
    ylabel('p', 'FontSize', 14)
    xlim (xlim_resp);
    ylim(plim);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ptick)
    xlabel('Time (ms)', 'FontSize', 14)
end

 %topo plot inset
 for fig = 8
     figure(8); clf;
     subplot(2, 2, 1) 
     %tg-locked
     flick_topolength_tg_oz = [sigbin_tg_flick_oz', sigbin_tg_flick_oz(end)+1]; %duration of sig flick effect 
                                                      %purposes
     flick_topotime_tg_oz = timebin.tg(flick_topolength_tg_oz);                                                     
     topoplot(squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
          (:, :, :, 1, 2, :, flick_topotime_tg_oz),7),3),2)))...
          -squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
          (:, :, :, 2, 2, :, flick_topotime_tg_oz),7),3),2))), 'mychan.loc'); %fast-slow
%      subplot(2, 2, 2)
%      topoplot(squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
%           (:, :, :, 1, 2, :, flick_topotime_tg_oz(1):flick_topotime_tg_oz(end)),7),3),2)))...
%           -squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
%           (:, :, :, 2, 2, :, flick_topotime_tg_oz(1):flick_topotime_tg_oz(end)),7),3),2))), 'mychan.loc'); %fast-slow
%      
%       
     colorbar
     caxis([-1 1])
     title('Tg-locked VN: Fast - Slow Flicker Rate');
 end

%% VN_Fig A right: Flicker rates (individual conditions)
 for fig = 9
clear tt
figure(9); clf;
for bc = 1:3
    subplot(3, 2, ((bc-1)*2)+1); 
    %target-locked
    %fast flicker rate
    b1 =  plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))...
        , 'Color', c_fast, 'LineWidth', lwidth); hold on;
    k1 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        b1, trans);
    
    %slow flicker rate
    b2 =  plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))...
        , 'Color', c_slow, 'LineWidth', lwidth); hold on;
    k2 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        b2, trans);

    %title(chan_label(ch), 'FontSize', 15);
    ylabel('Amplitude (uv)', 'FontSize', 14);
    xlim (xlim_tg);
    ylim (ylim_oz);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ytick_oz)
    
    subplot(3, 2, 2*bc); %resp-locked
    %fast flicker rate
    b3 =  plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))...
        , 'Color', c_fast, 'LineWidth', lwidth); hold on;
    k3 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 1, 2, eoi(ch_oz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        b3, trans);
    
    %slow flicker rate
    b4 =  plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))...
        , 'Color', c_slow, 'LineWidth', lwidth); hold on;
    k4 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_flicker...
        (:, bc, [1 2 3], 2, 2, eoi(ch_oz, :), :),6),3),2)))'./sqrt(numel(sub)), ...
        b4, trans);    

    xlim (xlim_resp);
    ylim (ylim_oz);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ytick_oz)
    xlabel('Time (ms)', 'FontSize', 14)
end
legend([b3 b4],'Fast Flicker Rate', 'Slow Flicker Rate');
suptitle('Fig A right: VN & Flicker Rates (non-collapsed)')
 end

%% VN_Fig B left: Expectation (collapsed)
for fig = 10
    figure(10); clf;
   
    suptitle('Fig B left: VN & expectation (collapsed)')

    subplot(3, 2, 1) %target-locked
    %expected
    c1 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_exp, 'LineWidth', lwidth); hold on;
    m1 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c1, trans);     
    
    %neutral
    c2 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_neu, 'LineWidth', lwidth); hold on;
    m2 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c2, trans);
    
    %unexpected
    c3 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_unexp, 'LineWidth', lwidth); hold on;
    m3 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, :, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c3, trans);
    
    %title(chan_label(ch), 'FontSize', 15);
    ylabel('Amplitude (uv)', 'FontSize', 14)
    xlim (xlim_tg);
    ylim (ylim_oz);
%     title('Tg-locked', 'FontSize', 10);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ytick_oz)
    
    %star plot
    for tt =1:numel(timebin.tg)-1
        if panova_tg(2, tt, 1) > p_fdr_flickexp
            colorsign = [1 1 1];
%         elseif panova(cnt, tt, 1) <= 0.05 &&  panova(cnt, tt, 1) > 0.01
%             colorsign = [0.6 0.6 0.6];
%         elseif panova(cnt, tt, 1) <= 0.01 &&  panova(cnt, tt, 1) > 0.001
%             colorsign = [0.3 0.3 0.3];
        else
            colorsign = [0 0 0];
        end
        
        if panova_tg(2, tt, 1) <= p_fdr_flickexp
            plot(chosentime.tg(tt)+binsize/2, starloc_oz, '.', 'color', colorsign, 'LineWidth', 10); 
        end
    end
    
    subplot(3, 2, 2) %response-locked 
    %expected
    c4 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_exp, 'LineWidth', lwidth); hold on;
    m4 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c4, trans);
    
    %neutral
    c5 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_neu, 'LineWidth', lwidth); hold on;
    m5 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c5, trans);
    
    %unexpected
    c6 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_unexp, 'LineWidth', lwidth); hold on;
    m6 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, :, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        c6, trans);
    
    %title(chan_label(ch), 'FontSize', 15);
    ylabel('Amplitude (uv)', 'FontSize', 14)
%     title('Resp-locked', 'FontSize', 10);
    xlim (xlim_resp);
    ylim (ylim_oz);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ytick_oz)
    legend([c4 c5 c6], 'Expected', 'Neutral', 'Unexpected', 'FontSize', 14)
    
    %star plot
    for tt =1:numel(timebin.resp)-1
        if panova_oz_exp_resp(1, tt) > p_fdr_flickexp
            colorsign = [1 1 1];
%         elseif panova(cnt, tt, 1) <= 0.05 &&  panova(cnt, tt, 1) > 0.01
%             colorsign = [0.6 0.6 0.6];
%         elseif panova(cnt, tt, 1) <= 0.01 &&  panova(cnt, tt, 1) > 0.001
%             colorsign = [0.3 0.3 0.3];
        else
            colorsign = [0 0 0];
        end
        
        if panova_oz_exp_resp(1, tt) <= p_fdr_flickexp
            plot(chosentime.resp(tt)+binsize/2, starloc_oz, '.', 'color', colorsign, 'LineWidth', 10);
        end
    end
     
    %F-value plot
    subplot(3, 2, 3)
    plot(chosentime.tg(1:end-1)+binsize/2, cell2mat(fvalexp_oz_tg), ...
        'color', 'k', 'Linewidth', lwidth_sub);
    ylabel('F', 'FontSize', 14)
    xlim (xlim_tg);
    ylim(flim);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ftick)
    
    
    subplot(3, 2, 4)
    plot(chosentime.resp(1:end-1)+binsize/2, cell2mat(fvalexp_oz_resp), ...
        'color', 'k', 'Linewidth', lwidth_sub);
    ylabel('F', 'FontSize', 14)
    xlim (xlim_resp);
    ylim(flim);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ftick)
    
    %p-value plot
    subplot(3, 2, 5)
    plot(chosentime.tg(1:end-1)+binsize/2, panova_tg(2, :, 1), ...
        'color', 'k', 'Linewidth', lwidth_sub); hold on;
    plot(timex.tg, repmat(p_fdr_flickexp, [1, size(timex.tg, 2)]), ...
        'color', 'k', 'LineStyle', '--');
    ylabel('p', 'FontSize', 14)
    xlim (xlim_tg);
    ylim(plim);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ptick)
    xlabel('Time (ms)', 'FontSize', 14)
    
    subplot(3, 2, 6)
    plot(chosentime.resp(1:end-1)+binsize/2, panova_oz_exp_resp(1, :), ...
        'color', 'k', 'Linewidth', lwidth_sub); hold on;
    plot(timex.resp, repmat(p_fdr_flickexp, [1, size(timex.resp, 2)]), ...
        'color', 'k', 'LineStyle', '--');
    ylabel('p', 'FontSize', 14)
    xlim (xlim_resp);
    ylim(plim);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ptick)
    xlabel('Time (ms)', 'FontSize', 14)
end

%% VN_Fig B right: Expectation (individual conditions)
for fig = 11
    figure(11); clf;
    suptitle('Fig B right: VN & expectation (non-collapsed)')
    
    for bc = 1:3
    subplot(3, 2, ((bc-1)*2)+1); %target-locked
    %expected
    d1 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_exp, 'LineWidth', lwidth); hold on;
    n1 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d1, trans);
    
    %neutral
    d2 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_neu, 'LineWidth', lwidth); hold on;
    n2 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d2, trans);    
    
    %unexpected
    d3 = plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_unexp, 'LineWidth', lwidth); hold on;
    n3 = bandedError(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp...
        (:, bc, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d3, trans);
    
    %title(chan_label(ch), 'FontSize', 15);
    ylabel('Amplitude (uv)', 'FontSize', 14)
    xlim (xlim_tg);
    ylim (ylim_oz);
    set(gca, 'XTick', xtick_tg)
    set(gca, 'YTick', ytick_oz)
    
    subplot(3, 2, 2*bc) %resp-locked
    %expected
    d4 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_exp, 'LineWidth', lwidth); hold on;
    n4 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 1, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d4, trans);
    
    %neutral
    d5 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_neu, 'LineWidth', lwidth); hold on;
    n5 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 2, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d5, trans);
    
    %unexpected
    d6 = plot(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))...
        , 'Color', c_unexp, 'LineWidth', lwidth); hold on;
    n6 = bandedError(timex.resp, squeeze(mean(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))', ...
        squeeze(std(mean(mean(mean(fig_gerp_resp_exp...
        (:, bc, 3, :, 2, eoi(ch_oz, :), :),6),4),2)))'./sqrt(numel(sub)), ...
        d6, trans);
    
    %title(chan_label(ch), 'FontSize', 15);
    ylabel('Amplitude (uv)', 'FontSize', 14)
    xlim (xlim_resp);
    ylim (ylim_oz);
    set(gca, 'XTick', xtick_resp)
    set(gca, 'YTick', ytick_oz)   
    xlabel('Time (ms)', 'FontSize', 14)
    end
    legend([d4 d5 d6], 'Expected', 'Neutral', 'Unexpected', 'FontSize', 14)
end

%% CPP_Slope Town
%buildup rate_tg-locked: 200 ms to 450 ms
t_tg1 = 200; %200 ms
t_tg2 = 550; %450 ms

t_resp1 = -350; %-350 ms
t_resp2 = 0; %-50 ms

for ss = 1:numel(sub)
    %collapsed across 3 expectation task conditions
    %tg-locked
    fast_tg.data = squeeze(mean(mean(gerp_tg(ss, :, [1 2 3], 1, 2, ...
        eoi(ch_cpz, :), nearest(timex.tg, t_tg1):nearest(timex.tg, t_tg2)), 3), 2));
    slow_tg.data = squeeze(mean(mean(gerp_tg(ss, :, [1 2 3], 2, 2, ...
        eoi(ch_cpz, :), nearest(timex.tg, t_tg1):nearest(timex.tg, t_tg2)), 3), 2));
      
    exp_tg.data = squeeze(mean(mean(gerp_tg(ss, :, 1, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.tg, t_tg1):nearest(timex.tg, t_tg2)), 4), 2));
    neu_tg.data = squeeze(mean(mean(gerp_tg(ss, :, 2, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.tg, t_tg1):nearest(timex.tg, t_tg2)), 4), 2));
    unexp_tg.data = squeeze(mean(mean(gerp_tg(ss, :, 3, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.tg, t_tg1):nearest(timex.tg, t_tg2)), 4), 2));
    
    %fit a line
    fast_tg.coeffs(ss, :) = polyfit(1:size(fast_tg.data, 1), fast_tg.data', 1);
    fast_tg.fittedX(ss, :) = linspace(1, round(((t_tg2-t_tg1)/fq)+1), 200);
    fast_tg.fittedY(ss, :) = polyval(fast_tg.coeffs(ss, :), fast_tg.fittedX(ss, :));
    
    slow_tg.coeffs(ss, :) = polyfit(1:size(slow_tg.data, 1), slow_tg.data', 1);
    slow_tg.fittedX(ss, :) = linspace(1, round(((t_tg2-t_tg1)/fq)+1), 200);
    slow_tg.fittedY(ss, :) = polyval(slow_tg.coeffs(ss, :), slow_tg.fittedX(ss, :));
    
    exp_tg.coeffs(ss, :) = polyfit(1:size(exp_tg.data, 1), exp_tg.data', 1);
    exp_tg.fittedX(ss, :) = linspace(1, round(((t_tg2-t_tg1)/fq)+1), 200);
    exp_tg.fittedY(ss, :) = polyval(exp_tg.coeffs(ss, :), exp_tg.fittedX(ss, :));
    
    neu_tg.coeffs(ss, :) = polyfit(1:size(neu_tg.data, 1), neu_tg.data', 1);
    neu_tg.fittedX(ss, :) = linspace(1, round(((t_tg2-t_tg1)/fq)+1), 200);
    neu_tg.fittedY(ss, :) = polyval(neu_tg.coeffs(ss, :), neu_tg.fittedX(ss, :));
    
    unexp_tg.coeffs(ss, :) = polyfit(1:size(unexp_tg.data, 1), unexp_tg.data', 1);
    unexp_tg.fittedX(ss, :) = linspace(1, round(((t_tg2-t_tg1)/fq)+1), 200);
    unexp_tg.fittedY(ss, :) = polyval(unexp_tg.coeffs(ss, :), unexp_tg.fittedX(ss, :));

    %plotting
%     figure();
%     subplot(2, 1, 1)
%     plot(fast_tg.data); hold on;
%     plot(fast_tg.fittedX(ss, :), fast_tg.fittedY(ss, :)); hold on;
%     plot(slow_tg.data); hold on;
%     plot(slow_tg.fittedX(ss, :), slow_tg.fittedY(ss, :));
%     subplot(2, 1, 2)
%     plot(exp_tg.data); hold on;
%     plot(exp_tg.fittedX(ss, :), exp_tg.fittedY(ss, :)); hold on;
%     plot(neu_tg.data); hold on;
%     plot(neu_tg.fittedX(ss, :), neu_tg.fittedY(ss, :)); hold on;
%     plot(unexp_tg.data); hold on;
%     plot(unexp_tg.fittedX(ss, :), unexp_tg.fittedY(ss, :)); hold on;
%     
    %resp-locked
    fast_resp.data = squeeze(mean(mean(gerp_resp(ss, :, [1 2 3], 1, 2, ...
        eoi(ch_cpz, :), nearest(timex.resp, t_resp1):nearest(timex.resp, t_resp2)), 3), 2));
    slow_resp.data = squeeze(mean(mean(gerp_resp(ss, :, [1 2 3], 2, 2, ...
        eoi(ch_cpz, :), nearest(timex.resp, t_resp1):nearest(timex.resp, t_resp2)), 3), 2));
    
    exp_resp.data = squeeze(mean(mean(gerp_resp(ss, :, 1, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.resp, t_resp1):nearest(timex.resp, t_resp2)), 4), 2));
    neu_resp.data = squeeze(mean(mean(gerp_resp(ss, :, 2, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.resp, t_resp1):nearest(timex.resp, t_resp2)), 4), 2));
    unexp_resp.data = squeeze(mean(mean(gerp_resp(ss, :, 3, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.resp, t_resp1):nearest(timex.resp, t_resp2)), 4), 2));
    
    %fit a line
    fast_resp.coeffs(ss, :) = polyfit(1:size(fast_resp.data, 1), fast_resp.data', 1);
    fast_resp.fittedX(ss, :) = linspace(1, round(((t_resp2-t_resp1)/fq)+1), 200);
    fast_resp.fittedY(ss, :) = polyval(fast_resp.coeffs(ss, :), fast_resp.fittedX(ss, :));
    
    slow_resp.coeffs(ss, :) = polyfit(1:size(slow_resp.data, 1), slow_resp.data', 1);
    slow_resp.fittedX(ss, :) = linspace(1, round(((t_resp2-t_resp1)/fq)+1), 200);
    slow_resp.fittedY(ss, :) = polyval(slow_resp.coeffs(ss, :), slow_resp.fittedX(ss, :));
    
    exp_resp.coeffs(ss, :) = polyfit(1:size(exp_resp.data, 1), exp_resp.data', 1);
    exp_resp.fittedX(ss, :) = linspace(1, round(((t_resp2-t_resp1)/fq)+1), 200);
    exp_resp.fittedY(ss, :) = polyval(exp_resp.coeffs(ss, :), exp_resp.fittedX(ss, :));
    
    neu_resp.coeffs(ss, :) = polyfit(1:size(neu_resp.data, 1), neu_resp.data', 1);
    neu_resp.fittedX(ss, :) = linspace(1, round(((t_resp2-t_resp1)/fq)+1), 200);
    neu_resp.fittedY(ss, :) = polyval(neu_resp.coeffs(ss, :), neu_resp.fittedX(ss, :));
    
    unexp_resp.coeffs(ss, :) = polyfit(1:size(unexp_resp.data, 1), unexp_resp.data', 1);
    unexp_resp.fittedX(ss, :) = linspace(1, round(((t_resp2-t_resp1)/fq)+1), 200);
    unexp_resp.fittedY(ss, :) = polyval(unexp_resp.coeffs(ss, :), unexp_resp.fittedX(ss, :));
    
    %plotting
%     figure();
%     subplot(2, 1, 1)
%     plot(fast_resp.data); hold on;
%     plot(fast_resp.fittedX(ss, :), fast_resp.fittedY(ss, :)); hold on;
%     plot(slow_resp.data); hold on;
%     plot(slow_resp.fittedX(ss, :), slow_resp.fittedY(ss, :));
%     subplot(2, 1, 2)
%     plot(exp_resp.data); hold on;
%     plot(exp_resp.fittedX(ss, :), exp_resp.fittedY(ss, :)); hold on;
%     plot(neu_resp.data); hold on;
%     plot(neu_resp.fittedX(ss, :), neu_resp.fittedY(ss, :)); hold on;
%     plot(unexp_resp.data); hold on;
%     plot(unexp_resp.fittedX(ss, :), unexp_resp.fittedY(ss, :)); hold on;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %by individual task conditions
    condcnt = 0;
    for cond = 1:3
        condcnt = condcnt + 1;
    %tg-locked
    fast_tg_cond.data = squeeze(mean(gerp_tg(ss, cond, :, 1, 2, ...
        eoi(ch_cpz, :), nearest(timex.tg, t_tg1):nearest(timex.tg, t_tg2)), 3));
    slow_tg_cond.data = squeeze(mean(gerp_tg(ss, cond, :, 2, 2, ...
        eoi(ch_cpz, :), nearest(timex.tg, t_tg1):nearest(timex.tg, t_tg2)), 3));
    
    exp_tg_cond.data = squeeze(mean(gerp_tg(ss, cond, 1, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.tg, t_tg1):nearest(timex.tg, t_tg2)), 4));
    neu_tg_cond.data = squeeze(mean(gerp_tg(ss, cond, 2, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.tg, t_tg1):nearest(timex.tg, t_tg2)), 4));
    unexp_tg_cond.data = squeeze(mean(gerp_tg(ss, cond, 3, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.tg, t_tg1):nearest(timex.tg, t_tg2)), 4));
    
    %resp-locked
    fast_resp_cond.data = squeeze(mean(gerp_resp(ss, cond, :, 1, 2, ...
        eoi(ch_cpz, :), nearest(timex.resp, t_resp1):nearest(timex.resp, t_resp2)), 3));
    slow_resp_cond.data = squeeze(mean(gerp_resp(ss, cond, :, 2, 2, ...
        eoi(ch_cpz, :), nearest(timex.resp, t_resp1):nearest(timex.resp, t_resp2)), 3));
    
    exp_resp_cond.data = squeeze(mean(gerp_resp(ss, cond, 1, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.resp, t_resp1):nearest(timex.resp, t_resp2)), 4));
    neu_resp_cond.data = squeeze(mean(gerp_resp(ss, cond, 2, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.resp, t_resp1):nearest(timex.resp, t_resp2)), 4));
    unexp_resp_cond.data = squeeze(mean(gerp_resp(ss, cond, 3, :, 2, ...
        eoi(ch_cpz, :), nearest(timex.resp, t_resp1):nearest(timex.resp, t_resp2)), 4));
    
    %fit a line
    fast_tg_cond.coeffs((condcnt-1)*numel(sub)+ss, :) = polyfit(1:size(fast_tg_cond.data, 1), fast_tg_cond.data', 1);
    slow_tg_cond.coeffs((condcnt-1)*numel(sub)+ss, :) = polyfit(1:size(slow_tg_cond.data, 1), slow_tg_cond.data', 1);
    exp_tg_cond.coeffs((condcnt-1)*numel(sub)+ss, :) = polyfit(1:size(exp_tg_cond.data, 1), exp_tg_cond.data', 1);
    neu_tg_cond.coeffs((condcnt-1)*numel(sub)+ss, :) = polyfit(1:size(neu_tg_cond.data, 1), neu_tg_cond.data', 1);
    unexp_tg_cond.coeffs((condcnt-1)*numel(sub)+ss, :) = polyfit(1:size(unexp_tg_cond.data, 1), unexp_tg_cond.data', 1);
    
    fast_resp_cond.coeffs((condcnt-1)*numel(sub)+ss, :) = polyfit(1:size(fast_resp_cond.data, 1), fast_resp_cond.data', 1);
    slow_resp_cond.coeffs((condcnt-1)*numel(sub)+ss, :) = polyfit(1:size(slow_resp_cond.data, 1), slow_resp_cond.data', 1);
    exp_resp_cond.coeffs((condcnt-1)*numel(sub)+ss, :) = polyfit(1:size(exp_resp_cond.data, 1), exp_resp_cond.data', 1);
    neu_resp_cond.coeffs((condcnt-1)*numel(sub)+ss, :) = polyfit(1:size(neu_resp_cond.data, 1), neu_resp_cond.data', 1);
    unexp_resp_cond.coeffs((condcnt-1)*numel(sub)+ss, :) = polyfit(1:size(unexp_resp_cond.data, 1), unexp_resp_cond.data', 1);
    end
end

%% %%%%%%%%%%%%%%%%%%% slope comparison: tg-locked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fast_tg.slope = fast_tg.coeffs(:, 1);
slow_tg.slope = slow_tg.coeffs(:, 1);
exp_tg.slope = exp_tg.coeffs(:, 1);
neu_tg.slope = neu_tg.coeffs(:, 1);
unexp_tg.slope = unexp_tg.coeffs(:, 1);

%mean and SEM
fast_tg.mean = mean(fast_tg.slope);
fast_tg.sem = std(fast_tg.slope)/sqrt(numel(sub));
slow_tg.mean = mean(slow_tg.slope);
slow_tg.sem = std(slow_tg.slope)/sqrt(numel(sub));
exp_tg.mean = mean(exp_tg.slope);
exp_tg.sem = std(exp_tg.slope)/sqrt(numel(sub));
neu_tg.mean = mean(neu_tg.slope);
neu_tg.sem = std(neu_tg.slope)/sqrt(numel(sub));
unexp_tg.mean = mean(unexp_tg.slope);
unexp_tg.sem = std(unexp_tg.slope)/sqrt(numel(sub));

%fast vs slow
[h p ci t] = ttest(fast_tg.slope, slow_tg.slope);
slope_tg_fastslow.p = p;
slope_tg_fastslow.h = h;
slope_tg_fastslow.t = t.tstat;
slope_tg_fastslow.df = t.df;

%one-way anova for exp effect
    DV = [exp_tg.slope; neu_tg.slope; unexp_tg.slope];
    sublist = repmat(1:numel(sub), [1, 3])'; %subject list
    IVexp = repmat([ones(1, numel(sub)) ones(1, numel(sub)).*2 ...
        ones(1, numel(sub)).*3], [1, 1])';
    
    [p, tbl, stats] = anovan(DV, {IVexp, sublist}, 'model', 'full', 'varnames',...
        {'exp/neu/unexp', 'subj'}, 'random', [2]);
    
%exp vs neu vs unexp
%1) exp vs unexp
clear h p ci t
[h p ci t] = ttest(exp_tg.slope, unexp_tg.slope);
slope_tg_expunexp.p = p;
slope_tg_expunexp.h = h;
slope_tg_expunexp.t = t.tstat;
slope_tg_expunexp.df = t.df;

%2) exp vs neu
clear h p ci t
[h p ci t] = ttest(exp_tg.slope, neu_tg.slope);
slope_tg_expneu.p = p;
slope_tg_expneu.h = h;
slope_tg_expneu.t = t.tstat;
slope_tg_expneu.df = t.df;

%3) neu vs unexp
clear h p ci t
[h p ci t] = ttest(neu_tg.slope, unexp_tg.slope);
slope_tg_neuunexp.p = p;
slope_tg_neuunexp.h = h;
slope_tg_neuunexp.t = t.tstat;
slope_tg_neuunexp.df = t.df;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%by individual conditions
condcnt = 0;
for cond = 1:3
    condcnt = condcnt + 1;
fast_tg_cond.mean(cond) = mean(fast_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1), 1);
fast_tg_cond.sem(cond) = std(fast_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1))/sqrt(numel(sub));
slow_tg_cond.mean(cond) = mean(slow_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1), 1);
slow_tg_cond.sem(cond) = std(slow_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1))/sqrt(numel(sub));

exp_tg_cond.mean(cond) = mean(exp_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1), 1);
exp_tg_cond.sem(cond) = std(exp_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1))/sqrt(numel(sub));
neu_tg_cond.mean(cond) = mean(neu_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1), 1);
neu_tg_cond.sem(cond) = std(neu_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1))/sqrt(numel(sub));
unexp_tg_cond.mean(cond) = mean(unexp_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1), 1);
unexp_tg_cond.sem(cond) = std(unexp_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1))/sqrt(numel(sub));

%fast vs slow
[h p ci t] = ttest(fast_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1), slow_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1));
slope_tg_cond_fastslow.p(cond) = p;
slope_tg_cond_fastslow.h(cond) = h;
slope_tg_cond_fastslow.t(cond) = t.tstat;
slope_tg_cond_fastslow.df(cond) = t.df;

%one-way anova for exp effect
    DV = [exp_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1); neu_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1); unexp_tg_cond.coeffs((condcnt-1)*numel(sub)+1:...
    condcnt*numel(sub), 1)];
    sublist = repmat(1:numel(sub), [1, 3])'; %subject list
    IVexp = repmat([ones(1, numel(sub)) ones(1, numel(sub)).*2 ...
        ones(1, numel(sub)).*3], [1, 1])';
    
    [p, tbl, stats] = anovan(DV, {IVexp, sublist}, 'model', 'full', 'varnames',...
        {'exp/neu/unexp', 'subj'}, 'random', [2]);
end

%% %%%%%%%%%%%%%%%%%%% slope comparison: resp-locked %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fast_resp.slope = fast_resp.coeffs(:, 1);
slow_resp.slope = slow_resp.coeffs(:, 1);
exp_resp.slope = exp_resp.coeffs(:, 1);
neu_resp.slope = neu_resp.coeffs(:, 1);
unexp_resp.slope = unexp_resp.coeffs(:, 1);

%fast vs slow
[h p ci t] = ttest(fast_resp.slope, slow_resp.slope);
slope_resp_fastslow.p = p;
slope_resp_fastslow.h = h;
slope_resp_fastslow.t = t.tstat;
slope_resp_fastslow.df = t.df;

%one-way anova for exp effect
    DV = [exp_resp.slope; neu_resp.slope; unexp_resp.slope];
    sublist = repmat(1:numel(sub), [1, 3])'; %subject list
    IVexp = repmat([ones(1, numel(sub)) ones(1, numel(sub)).*2 ...
        ones(1, numel(sub)).*3], [1, 1])';
    
    [p, tbl, stats] = anovan(DV, {IVexp, sublist}, 'model', 'full', 'varnames',...
        {'exp/neu/unexp', 'subj'}, 'random', [2]);

%exp vs neu vs unexp
%1) exp vs unexp
clear h p ci t
[h p ci t] = ttest(exp_resp.slope, unexp_resp.slope);
slope_resp_expunexp.p = p;
slope_resp_expunexp.h = h;
slope_resp_expunexp.t = t.tstat;
slope_resp_expunexp.df = t.df;

%2) exp vs neu
clear h p ci t
[h p ci t] = ttest(exp_resp.slope, neu_resp.slope);
slope_resp_expneu.p = p;
slope_resp_expneu.h = h;
slope_resp_expneu.t = t.tstat;
slope_resp_expneu.df = t.df;

%3) neu vs unexp
clear h p ci t
[h p ci t] = ttest(neu_resp.slope, unexp_resp.slope);
slope_resp_neuunexp.p = p;
slope_resp_neuunexp.h = h;
slope_resp_neuunexp.t = t.tstat;
slope_resp_neuunexp.df = t.df;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CPP Amplitude: Bayes Factor
%this is to test null effect of expectation on CPP amplitude
 %using the time bins where we saw sig effects of flick on cpp
 
%1) 250 to 750 ms tg-locked
%1.1) flicker rates
cpz_tg_fast_null = squeeze(mean(mean(mean(gerp_tgbinned(:, :, :, 1, 2, eoi(ch_cpz, :),...
    sigbin_tg_flick_cpz1), 7), 3), 2));
cpz_tg_slow_null = squeeze(mean(mean(mean(gerp_tgbinned(:, :, :, 2, 2, eoi(ch_cpz, :),...
    sigbin_tg_flick_cpz1), 7), 3), 2));

%first compute t stats for each pair
clear h p ci t
[h p ci t] = ttest(cpz_tg_fast_null, cpz_tg_slow_null);
cpz_tg_fastslow_null.p = p;
cpz_tg_fastslow_null.h = h;
cpz_tg_fastslow_null.t = t.tstat;
cpz_tg_fastslow_null.df = t.df;

r = 0.707; %default scale factor = 0.707
bf10.cpz_tg_amp_fastslow = t1smpbf(cpz_tg_fastslow_null.t, numel(sub), 0.707);
 
%1.2) expectation
cpz_tg_exp_null = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 1, :, 2, eoi(ch_cpz, :),...
    sigbin_tg_flick_cpz1), 7), 4), 2));
cpz_tg_neu_null = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 2, :, 2, eoi(ch_cpz, :),...
    sigbin_tg_flick_cpz1), 7), 4), 2));
cpz_tg_unexp_null = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 3, :, 2, eoi(ch_cpz, :),...
    sigbin_tg_flick_cpz1), 7), 4), 2)); 
%first compute t stats for each pair
clear h p ci t
[h p ci t] = ttest(cpz_tg_exp_null, cpz_tg_unexp_null);
cpz_tg_expunexp_null.p = p;
cpz_tg_expunexp_null.h = h;
cpz_tg_expunexp_null.t = t.tstat;
cpz_tg_expunexp_null.df = t.df;
clear h p ci t
[h p ci t] = ttest(cpz_tg_exp_null, cpz_tg_neu_null);
cpz_tg_expneu_null.p = p;
cpz_tg_expneu_null.h = h;
cpz_tg_expneu_null.t = t.tstat;
cpz_tg_expneu_null.df = t.df;
clear h p ci t
[h p ci t] = ttest(cpz_tg_neu_null, cpz_tg_unexp_null);
cpz_tg_neuunexp_null.p = p;
cpz_tg_neuunexp_null.h = h;
cpz_tg_neuunexp_null.t = t.tstat;
cpz_tg_neuunexp_null.df = t.df;

r = 0.707; %defaul scale factor = 0.707
bf10.cpz_tg_amp_expunexp = t1smpbf(cpz_tg_expunexp_null.t, numel(sub), 0.707);
bf10.cpz_tg_amp_expneu = t1smpbf(cpz_tg_expneu_null.t, numel(sub), 0.707);
bf10.cpz_tg_amp_neuunexp = t1smpbf(cpz_tg_neuunexp_null.t, numel(sub), 0.707);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2) -300 to -200 ms resp-locked
%2.1) flicker rates
cpz_resp_fast_null = squeeze(mean(mean(mean(gerp_respbinned(:, :, :, 1, 2, eoi(ch_cpz, :),...
    sigbin_resp_flick_cpz1), 7), 3), 2));
cpz_resp_slow_null = squeeze(mean(mean(mean(gerp_respbinned(:, :, :, 2, 2, eoi(ch_cpz, :),...
    sigbin_resp_flick_cpz1), 7), 3), 2));

%first compute t stats for each pair
clear h p ci t
[h p ci t] = ttest(cpz_resp_fast_null, cpz_resp_slow_null);
cpz_resp_fastslow_null.p = p;
cpz_resp_fastslow_null.h = h;
cpz_resp_fastslow_null.t = t.tstat;
cpz_resp_fastslow_null.df = t.df;

r = 0.707; %default scale factor = 0.707
bf10.cpz_resp_amp_fastslow = t1smpbf(cpz_resp_fastslow_null.t, numel(sub), 0.707);

%2.2) expectation
cpz_resp_exp_null = squeeze(mean(mean(mean(gerp_respbinned(:, :, 1, :, 2, eoi(ch_cpz, :),...
    sigbin_resp_flick_cpz1), 7), 4), 2));
cpz_resp_neu_null = squeeze(mean(mean(mean(gerp_respbinned(:, :, 2, :, 2, eoi(ch_cpz, :),...
    sigbin_resp_flick_cpz1), 7), 4), 2));
cpz_resp_unexp_null = squeeze(mean(mean(mean(gerp_respbinned(:, :, 3, :, 2, eoi(ch_cpz, :),...
    sigbin_resp_flick_cpz1), 7), 4), 2)); 
%first compute t stats for each pair
clear h p ci t
[h p ci t] = ttest(cpz_resp_exp_null, cpz_resp_unexp_null);
cpz_resp_expunexp_null.p = p;
cpz_resp_expunexp_null.h = h;
cpz_resp_expunexp_null.t = t.tstat;
cpz_resp_expunexp_null.df = t.df;
clear h p ci t
[h p ci t] = ttest(cpz_resp_exp_null, cpz_resp_neu_null);
cpz_resp_expneu_null.p = p;
cpz_resp_expneu_null.h = h;
cpz_resp_expneu_null.t = t.tstat;
cpz_resp_expneu_null.df = t.df;
clear h p ci t
[h p ci t] = ttest(cpz_resp_neu_null, cpz_resp_unexp_null);
cpz_resp_neuunexp_null.p = p;
cpz_resp_neuunexp_null.h = h;
cpz_resp_neuunexp_null.t = t.tstat;
cpz_resp_neuunexp_null.df = t.df;

r = 0.707; %defaul scale factor = 0.707
bf10.cpz_resp_amp_expunexp = t1smpbf(cpz_resp_expunexp_null.t, numel(sub), 0.707);
bf10.cpz_resp_amp_expneu = t1smpbf(cpz_resp_expneu_null.t, numel(sub), 0.707);
bf10.cpz_resp_amp_neuunexp = t1smpbf(cpz_resp_neuunexp_null.t, numel(sub), 0.707);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% VN Amplitude: Bayes Factor
%1) 200 to 300 ms tg-locked
%1.1) flicker rates
oz_tg_fast_null = squeeze(mean(mean(mean(gerp_tgbinned(:, :, :, 1, 2, eoi(ch_oz, :),...
    sigbin_tg_flick_oz), 7), 3), 2));
oz_tg_slow_null = squeeze(mean(mean(mean(gerp_tgbinned(:, :, :, 2, 2, eoi(ch_oz, :),...
    sigbin_tg_flick_oz), 7), 3), 2));

%first compute t stats for each pair
clear h p ci t
[h p ci t] = ttest(oz_tg_fast_null, oz_tg_slow_null);
oz_tg_fastslow_null.p = p;
oz_tg_fastslow_null.h = h;
oz_tg_fastslow_null.t = t.tstat;
oz_tg_fastslow_null.df = t.df;

r = 0.707; %default scale factor = 0.707
bf10.oz_tg_amp_fastslow = t1smpbf(oz_tg_fastslow_null.t, numel(sub), 0.707);

%1.2) expectation
oz_tg_exp_null = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 1, :, 2, eoi(ch_oz, :),...
    sigbin_tg_flick_oz), 7), 4), 2));
oz_tg_neu_null = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 2, :, 2, eoi(ch_oz, :),...
    sigbin_tg_flick_oz), 7), 4), 2));
oz_tg_unexp_null = squeeze(mean(mean(mean(gerp_tgbinned(:, :, 3, :, 2, eoi(ch_oz, :),...
    sigbin_tg_flick_oz), 7), 4), 2)); 
%first compute t stats for each pair
clear h p ci t
[h p ci t] = ttest(oz_tg_exp_null, oz_tg_unexp_null);
oz_tg_expunexp_null.p = p;
oz_tg_expunexp_null.h = h;
oz_tg_expunexp_null.t = t.tstat;
oz_tg_expunexp_null.df = t.df;
clear h p ci t
[h p ci t] = ttest(oz_tg_exp_null, oz_tg_neu_null);
oz_tg_expneu_null.p = p;
oz_tg_expneu_null.h = h;
oz_tg_expneu_null.t = t.tstat;
oz_tg_expneu_null.df = t.df;
clear h p ci t
[h p ci t] = ttest(oz_tg_neu_null, oz_tg_unexp_null);
oz_tg_neuunexp_null.p = p;
oz_tg_neuunexp_null.h = h;
oz_tg_neuunexp_null.t = t.tstat;
oz_tg_neuunexp_null.df = t.df;

r = 0.707; %defaul scale factor = 0.707
bf10.oz_tg_amp_expunexp = t1smpbf(oz_tg_expunexp_null.t, numel(sub), 0.707);
bf10.oz_tg_amp_expneu = t1smpbf(oz_tg_expneu_null.t, numel(sub), 0.707);
bf10.oz_tg_amp_neuunexp = t1smpbf(oz_tg_neuunexp_null.t, numel(sub), 0.707);

%% CPP Slope: Bayes Factor 
r = 0.707; %defaul scale factor = 0.707
%1) flicker rate
bf10.cpz_tg_slope_fastslow = t1smpbf(slope_tg_fastslow.t, numel(sub), 0.707);
%2) expectation
bf10.cpz_tg_slope_expunexp = t1smpbf(slope_tg_expunexp.t, numel(sub), 0.707);
bf10.cpz_tg_slope_expneu = t1smpbf(slope_tg_expneu.t, numel(sub), 0.707);
bf10.cpz_tg_slope_neuunexp = t1smpbf(slope_tg_neuunexp.t, numel(sub), 0.707);

%%
time_marker =[f_onset_t; f_peak_t; f_offset_t];
time_marker =[f_onset_t; f_peak_t; f_offset_t];

ylabel('Amplitude (uv)', 'FontSize', 14)
xlabel('Time (ms)', 'FontSize', 14)
suptitle ('Single Subject CPP')
legend('Fast flicker', 'Slow flicker')

%% t-test to compare slopes as a function of flicker rates
%1) pre-peaked slope
clear h p ci t
[h p ci t] = ttest(slope.f.pre, slope.s.pre);
flick_pre.p = p;
flick_pre.h = h;
flick_pre.t = t.tstat;
flick_pre.df = t.df;

%2) post-peaked slope
clear h p ci t
[h p ci t] = ttest(slope.f.post, slope.s.post);
flick_post.p = p;
flick_post.h = h;
flick_post.t = t.tstat;
flick_post.df = t.df;
%% one-way anova to compare slopes as a function of expectation
%1) pre-peaked slope
DV_pre = [slope.e.pre; slope.n.pre; slope.u.pre];
sublist = repmat(1:numel(sub), [1, 3])'; %subject list
IVexp = repmat([ones(1, numel(sub)) ones(1, numel(sub)).*2 ...
    ones(1, numel(sub)).*3], [1, 1])';

[p_pre, tbl_pre, stats_pre] = anovan(DV_pre, {IVexp, sublist}, 'model', 'full', 'varnames',...
    {'exp/neu/unexp', 'subj'}, 'random', [2]);
p_exp_pre = p;
fvalue_exp_pre = tbl_pre(2, 6);

%2) post-peaked slope
DV_post = [slope.e.post; slope.n.post; slope.u.post];
%     sublist = repmat(1:numel(sub), [1, 3])'; %subject list
%     IVexp = repmat([ones(1, numel(sub)) ones(1, numel(sub)).*2 ...
%         ones(1, numel(sub)).*3], [1, 1])';

[p_post, tbl_post, stats_post] = anovan(DV_post, {IVexp, sublist}, 'model', 'full', 'varnames',...
    {'exp/neu/unexp', 'subj'}, 'random', [2]);
p_exp_post = p_post;
fvalue_exp_post = tbl_post(2, 6);

for cond = 1:3
    cond_exp = squeeze(mean(mean(gerp_tgbinned(:, cond, 1, :, 2, eoi(ch_cpz, :), 22:26), 7), 4));
    cond_neu = squeeze(mean(mean(gerp_tgbinned(:, cond, 2, :, 2, eoi(ch_cpz, :), 22:26), 7), 4));
    cond_unexp = squeeze(mean(mean(gerp_tgbinned(:, cond, 3, :, 2, eoi(ch_cpz, :), 22:26), 7), 4));
    
    DV = [cond_exp; cond_neu; cond_unexp];
    sublist = repmat(1:numel(sub), [1, 3])'; %subject list
    IVexp = repmat([ones(1, numel(sub)) ones(1, numel(sub)).*2 ...
        ones(1, numel(sub)).*3], [1, 1])';
    
    [p, tbl, stats] = anovan(DV, {IVexp, sublist}, 'model', 'full', 'varnames',...
        {'exp/neu/unexp', 'subj'}, 'random', [2]);
    p_indivcond_tg(cond, :) = p;
    fvalue_indivcond_tg(cond, :) = tbl(2, 6);
    
end


%%
for ss = 1%:numel(sub)
    plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_exp(:, :, 1, :, 2, eoi(ch_cpz, :), :),6),4),2))) - ...
        squeeze(std(mean(mean(mean(fig_gerp_tg_exp(:, :, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)))./sqrt(17) ...
        , 'b', 'LineWidth', 1); hold on;
    
    
    fast = squeeze(mean(mean(gerp_tg(ss, :, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2));
    slow = squeeze(mean(mean(mean(gerp_tg(ss, :, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2));
    exp = squeeze(mean(mean(mean(mean(gerp_tg(:, :, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)));
    neu = squeeze(mean(mean(mean(mean(gerp_tg(:, :, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)));
    unexp = squeeze(mean(mean(mean(mean(gerp_tg(:, :, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)));
    
    f_binned_slope = (fast(bin(2:end))-fast(bin(1:end-1)))./fq;
    s_binned_slope = (slow(bin(2:end))-slow(bin(1:end-1)))./fq;
    e_binned_slope = (exp(bin(2:end))-exp(bin(1:end-1)))./fq;
    n_binned_slope = (neu(bin(2:end))-neu(bin(1:end-1)))./fq;
    u_binned_slope = (unexp(bin(2:end))-unexp(bin(1:end-1)))./fq;
    
    f_onset_inds = find(f_binned_slope > slope_crit); %find indexes for all bins that have higher slope than some
    f_onset_ind = f_onset_inds(find(f_onset_inds > 4, 1));
    f_onset = bin(f_onset_ind-1); %find the onset of ERP to calculate pre-peaked slope
    
    s_onset_inds = find(s_binned_slope > slope_crit);
    s_onset_ind = s_onset_inds(find(s_onset_inds > 4, 1));
    s_onset = bin(s_onset_ind-1);
    
    e_onset_inds = find(e_binned_slope > slope_crit);
    e_onset_ind = e_onset_inds(find(e_onset_inds > 4, 1));
    e_onset = bin(e_onset_ind-1);
    
    n_onset_inds = find(n_binned_slope > slope_crit);
    n_onset_ind = n_onset_inds(find(n_onset_inds > 4, 1));
    n_onset = bin(n_onset_ind-1);
    
    u_onset_inds = find(u_binned_slope > slope_crit);
    u_onset_ind = u_onset_inds(find(u_onset_inds > 4, 1));
    u_onset = bin(u_onset_ind-1);
    
    figure;
    subplot(2,2,1)
    plot(timex.tg, squeeze(mean(mean(mean(gerp_tg(ss, :, [1 2 3], 1, 2, eoi(ch_cpz, :), :),6),3),2)),'m', 'LineWidth', 1); hold on; %fast
    plot(timex.tg, squeeze(mean(mean(mean(gerp_tg(ss, :, [1 2 3], 2, 2, eoi(ch_cpz, :), :),6),3),2)),'g', 'LineWidth', 1); %slow
    ylim ([-3 16]);
    xlim ([-100 1475]);
    xlabel('time(ms)')
    subplot(2,2,2)
    plot(timex.tg, squeeze(mean(mean(mean(gerp_tg(ss, :, 1, :, 2, eoi(ch_cpz, :), :),6),4),2)),'b', 'LineWidth', 1); hold on; %exp
    plot(timex.tg, squeeze(mean(mean(mean(gerp_tg(ss, :, 2, :, 2, eoi(ch_cpz, :), :),6),4),2)),'k', 'LineWidth', 1); %neu
    plot(timex.tg, squeeze(mean(mean(mean(gerp_tg(ss, :, 3, :, 2, eoi(ch_cpz, :), :),6),4),2)),'r', 'LineWidth', 1); %unexp
    ylim ([-3 16]);
    xlim ([-100 1475]);
    xlabel('time(ms)')
    subplot(2,2,3)
    scatter(timex.tg(bin(1)+13:26:bin(end)-13), f_binned_slope, 'm'); hold on;
    scatter(timex.tg(bin(1)+13:26:bin(end)-13), s_binned_slope, 'g');
    xlim ([-100 1475]);
    title('Binned Slope: Fast vs Slow Flicker')
    subplot(2,2,4)
    plot(timex.tg(bin(1)+13:26:bin(end)-13), e_binned_slope, 'b'); hold on;
    plot(timex.tg(bin(1)+13:26:bin(end)-13), n_binned_slope, 'k'); hold on;
    plot(timex.tg(bin(1)+13:26:bin(end)-13), u_binned_slope, 'r');
    xlim ([-100 1475]);
    title('Binned Slope: Exp vs Neu vs Unexp')
    
    f_peak = find(fast==max(fast));
    s_peak = find(slow==max(slow));
    e_peak = find(exp==max(exp));
    n_peak = find(neu==max(neu));
    u_peak = find(unexp==max(unexp));
    
    slope.f.pre(ss) = (fast(f_peak)-fast(f_onset))./fq;
    slope.f.post(ss) = (fast(end)-fast(f_peak))./fq;
    
    slope.s.pre(ss) = (slow(s_peak)-slow(s_onset))./fq;
    slope.s.post(ss) = (slow(end)-slow(s_peak))./fq;
    
    slope.e.pre(ss) = (exp(e_peak)-exp(e_onset))./fq;
    slope.e.post(ss) = (exp(end)-exp(e_peak))./fq;
    
    slope.n.pre(ss) = (neu(n_peak)-neu(n_onset))./fq;
    slope.n.post(ss) = (neu(end)-neu(n_peak))./fq;
    
    slope.u.pre(ss) = (unexp(u_peak)-unexp(u_onset))./fq;
    slope.u.post(ss) = (unexp(end)-unexp(u_peak))./fq;
    
    
    %     subplot(6, 3, ss)
    %     plot(timex.tg, squeeze(mean(mean(mean(gerp_tg(ss, :, [1 2 3], 1, 2, eoi(ch, :), :),6),3),2)),'m', 'LineWidth', 1); hold on; %fast
    %     plot(timex.tg, squeeze(mean(mean(mean(gerp_tg(ss, :, [1 2 3], 2, 2, eoi(ch, :), :),6),3),2)),'g', 'LineWidth', 1); %slow
    %
    
    %     title(['sub' num2str(sub(ss))], 'FontSize', 10);
end

ylabel('Amplitude (uv)', 'FontSize', 14)
xlabel('Time (ms)', 'FontSize', 14)
suptitle ('Single Subject CPP')
legend('Fast flicker', 'Slow flicker')


cnt = 0;
for ch = [6] % 9]
    cnt = cnt+1;
    
    %%binned data
    %fast = squeeze(mean(mean(mean(mean(gerp_tgbinned(:, :, [1 2 3], 1, 2, eoi(ch, :), :),6),3),2)));
    %slow = squeeze(mean(mean(mean(mean(gerp_tgbinned(:, :, [1 2 3], 2, 2, eoi(ch, :), :),6),3),2)));
    %exp = squeeze(mean(mean(mean(mean(gerp_tgbinned(:, :, 1, :, 2, eoi(ch, :), :),6),4),2)));
    %neu = squeeze(mean(mean(mean(mean(gerp_tgbinned(:, :, 2, :, 2, eoi(ch, :), :),6),4),2)));
    %unexp = squeeze(mean(mean(mean(mean(gerp_tgbinned(:, :, 3, :, 2, eoi(ch, :), :),6),4),2)));
    
    %unbinned data
    fast = squeeze(mean(mean(mean(mean(gerp_tg(:, :, [1 2 3], 1, 2, eoi(ch, :), :),6),3),2)));
    slow = squeeze(mean(mean(mean(mean(gerp_tg(:, :, [1 2 3], 2, 2, eoi(ch, :), :),6),3),2)));
    exp = squeeze(mean(mean(mean(mean(gerp_tg(:, :, 1, :, 2, eoi(ch, :), :),6),4),2)));
    neu = squeeze(mean(mean(mean(mean(gerp_tg(:, :, 2, :, 2, eoi(ch, :), :),6),4),2)));
    unexp = squeeze(mean(mean(mean(mean(gerp_tg(:, :, 3, :, 2, eoi(ch, :), :),6),4),2)));
    
    fast_peak = find(fast==max(fast));
    slow_peak = find(slow==max(slow));
    exp_peak = find(exp==max(exp));
    neu_peak = find(neu==max(neu));
    unexp_peak = find(unexp==max(unexp));
    
    %pre-peaked slope
    slope.tg.fastpre = (fast(2:fast_peak)-fast(1:fast_peak-1))./fq);
    slope.tg.fastpost = (fast(fast_peak+2:end)-fast(fast_peak+1:end-1))./fq);
    fast_drift_free = fast(fast_peak)-fast(nearest(timex.tg, 0))
    
    
    slope.tg.fast(:, cnt) = (fast(2:end)-fast(1:end-1))./fq);
    slope.tg.slow(:, cnt) = (slow(2:end)-slow(1:end-1))./fq);
    slope.tg.exp(:, cnt) = (exp(2:end)-exp(1:end-1))./fq);
    slope.tg.neu(:, cnt) = (neu(2:end)-neu(1:end-1))./fq);
    slope.tg.unexp(:, cnt) = (unexp(2:end)-unexp(1:end-1))./fq);
end

load('behav.mat')
%hit(subject, prior, task condition, fast/slow target)
meanRT = mean(corRT);
rt.prior = mean(mean(meanRT,4),3);%exp/neu/unexp
rt.exp = rt.prior(1);
rt.neu = rt.prior(2);
rt.unexp = rt.prior(3);

rt.flick = mean(mean(meanRT,3),2); %fast/slow
rt.fast = rt.flick(1);
rt.slow = rt.flick(2);

nearest(timex.tg, rt.fast)

plot(timex.tg)
plot(timex.tg, squeeze(mean(mean(mean(mean(fig_gerp_tg_flicker(:, :, [1 2 3], 1, 2, eoi(ch, :), :),6),3),2))) - ...
    squeeze(std(mean(mean(mean(fig_gerp_tg_flicker(:, :, [1 2 3], 1, 2, eoi(ch, :), :),6),3),2)))./sqrt(17) ...
    , 'm', 'LineWidth', 1); hold on;


%this is for plotting the binned slopes
close all;
figure;
for chancnt = 1:2
    subplot(2,2,chancnt)
    chosentime.tg_slope = chosentime.tg+25;
    scatter(chosentime.tg_slope(1:31), slope.tg.fast(:, chancnt), 'r'); hold on;
    scatter(chosentime.tg_slope(1:31), slope.tg.slow(:, chancnt), 'k'); hold on;
    plot(chosentime.tg_slope(1:31), slope.tg.fast(:, chancnt), 'r'); hold on;
    plot(chosentime.tg_slope(1:31), slope.tg.slow(:, chancnt), 'k');
    legend('fast', 'slow');xlim ([chosentime.tg_slope(1) chosentime.tg_slope(31)]);
    title('flicker rate')
    
    subplot(2,2,chancnt+2)
    scatter(chosentime.tg_slope(1:31), slope.tg.exp(:, chancnt), 'g'); hold on;
    scatter(chosentime.tg_slope(1:31), slope.tg.neu(:, chancnt), 'b'); hold on;
    scatter(chosentime.tg_slope(1:31), slope.tg.unexp(:, chancnt), 'r'); hold on;
    plot(chosentime.tg_slope(1:31), slope.tg.exp(:, chancnt), 'g'); hold on;
    plot(chosentime.tg_slope(1:31), slope.tg.neu(:, chancnt), 'b'); hold on;
    plot(chosentime.tg_slope(1:31), slope.tg.unexp(:, chancnt), 'r'); hold on;
    legend('exp', 'neu', 'unexp');xlim ([chosentime.tg_slope(1) chosentime.tg_slope(31)]);
    title('expectation')
end

%continuous derivative
close all;
figure;
for chancnt = 1:2
    subplot(2,2,chancnt)
    plot(timex.tg(1:end-1), slope.tg.fast(:, chancnt), 'r'); hold on;
    plot(timex.tg(1:end-1), slope.tg.slow(:, chancnt), 'k'); hold on;
    legend('fast', 'slow');
    xlim ([-100 1500]);
    title('flicker rate')
    
    subplot(2,2,chancnt+2)
    plot(timex.tg(1:end-1), slope.tg.exp(:, chancnt), 'g'); hold on;
    plot(timex.tg(1:end-1), slope.tg.neu(:, chancnt), 'b'); hold on;
    plot(timex.tg(1:end-1), slope.tg.unexp(:, chancnt), 'r'); hold on;
    legend('exp', 'neu', 'unexp');
    xlim ([-100 1500]);
    title('expectation')
end
suptitle('point-by-point derivative: Cpz(left) & Oz(right)')

%
fasttg_cpz = slope.tg.fast(:, 1);
fasttg_oz = slope.tg.fast(:, 2);
slowtg_cpz = slope.tg.slow(:, 1);
slowtg_oz = slope.tg.slow(:, 2);
exptg_cpz = slope.tg.exp(:, 1);
exptg_oz = slope.tg.exp(:, 2);
neutg_cpz = slope.tg.neu(:, 1);
neutg_oz = slope.tg.neu(:, 2);
unexptg_cpz = slope.tg.unexp(:, 1);
unexptg_oz = slope.tg.unexp(:, 2);

ttt = 1:25:1536;
bin = 25; %25 time points/bin ~50ms bin
for ii = 1:size(ttt)-1
    med.fasttg_cpz(ii) = median(fasttg_cpz((1+(bin.*(ii-1))):(ii.*bin)));
    %(((ii-1)*bin)+1:ii*bin));
    
end
figure;
plot(ttt(1:end - 1),med.fasttg_cpz);

% % median of each window
% cnt = 0;
% for ch = [6, 9]
%     cnt = cnt + 1;
%     for bb = 1:numel(timebin.tg)-1
%         %flicker rates
%         med.fast(cnt, bb) = median(mean(mean(mean(mean(gerp_tg(:, :, [1 2 3], 1, 2, eoi(ch, :), timebin.tg(bb):timebin.tg(bb+1)-1),6),3),2)));
%         med.slow(cnt, bb) = median(mean(mean(mean(mean(gerp_tg(:, :, [1 2 3], 2, 2, eoi(ch, :), timebin.tg(bb):timebin.tg(bb+1)-1),6),3),2)));
%         %expectation
%         med.exp(cnt, bb) = median(mean(mean(mean(mean(gerp_tg(:, :, 1, :, 2, eoi(ch, :), timebin.tg(bb):timebin.tg(bb+1)-1),6),4),2)));
%         med.neu(cnt, bb) = median(mean(mean(mean(mean(gerp_tg(:, :, 2, :, 2, eoi(ch, :), timebin.tg(bb):timebin.tg(bb+1)-1),6),4),2)));
%         med.unexp(cnt, bb) = median(mean(mean(mean(mean(gerp_tg(:, :, 3, :, 2, eoi(ch, :), timebin.tg(bb):timebin.tg(bb+1)-1),6),4),2)));
%     end
% end