clear all;
close all;
sub = input ('sbj: '); 
cd data
cd (['sbj' num2str(sub)]);
bname = dir(['*trip_subj*mat']); % find behav file
sname = (['eegEpochArtFree_sbj' num2str(sub)  '.mat']); % set save name
snamerm = (['eegEpochArtFree_rmbase_sbj' num2str(sub)  '.mat']); % set save name


%eyeThresh = [60 60 80 80 80 80 100]; % threshold rejection 
%eyeThresh = [100 100 100 100 100 100 150];
eyeThresh = [80 80 100 100 100 100 150]; %standard
% hor left /right  vertical left right 

eeg = pop_loadset(['sbj' num2str(sub) '_epochICA_clean.set']); % load clean ICA epoched EEG data --> this is a big epoch
eegrm = pop_rmbase(eeg, [-200 0]); % rm baseline bf stim onset in ms


eegtg = pop_epoch(eeg, {'6091'}, [-1.5 1.5]); % lock to target onset
eegtgrm = pop_rmbase(eegtg, [-200 0]); %remove baseline 
%eegtgrm = pop_epoch(eegrm, {'6091'}, [-1.5 1.5]); % rm basline

%%%%%%%%%%% convert eventype in stuct to matrix %%%%%%%%%%%%%%%
fq = 1000/512;
prestim = 1500;
clear epoch

for t= 1:size(eeg.data,3) % go through each
    j = 1;
    while j <= size(eeg.epoch(t).eventtype,2)
    typ = 1;
   if str2num(cell2mat(eeg.epoch(t).eventtype(j))) < 6000
   epoch.eventtype(t, typ) = str2num(cell2mat(eeg.epoch(t).eventtype(j))); % stimonset 
   epoch.eventlatency(t, typ) = cell2mat(eeg.epoch(t).eventlatency(j)); 
   epoch.eventtimepoint(t, typ) = cell2mat(eeg.epoch(t).eventlatency(j))/fq+prestim/fq;
   j = 100;
    end
        j = j +1;
    end
end

for t= 1:size(eegtg.data,3)
    j = 1;
    while j <= size(eegtg.epoch(t).eventtype,2)

    typ = 1;
    if str2num(cell2mat(eegtg.epoch(t).eventtype(j))) < 6000
   epochtg.eventtype(t, typ) = str2num(cell2mat(eegtg.epoch(t).eventtype(j))); % stimonset 
   epochtg.eventlatency(t, typ) = cell2mat(eegtg.epoch(t).eventlatency(j)); 
   epochtg.eventtimepoint(t, typ) = cell2mat(eegtg.epoch(t).eventlatency(j))/fq+prestim/fq;
   j = 100;
    end
        j = j +1;
    end
end


epoch.data = eeg.data;
epochtg.data = eegtg.data;
epochrm.data = eegrm.data;
epochtgrm.data = eegtgrm.data;
%%%%%%%%%%% create pcat for the whole behavioral tags %%%%%%%%%%%%%%
condi = {'color_redblue', 'ori_VH', 'prior_705030', ...
    'bias_no_color_ori_resp', 'freq', 'expResp', 'run', 'tooslow', ...
    'hit', 'resp', 'falsealarm' , 'RT_fromStimOnset', ...
    'RT_fromtgonset', 'f_tgOnset', 'stimonset', 'tgonset', ...
    'StimToTgDur', 'tgoffset', 'tgDur', 'stimdur', 'blankonset', ...
    'blankdur', 'feedbackonset' , 'feedbackdur', 'ITIonset', 'ITIdur', 'ITI'};
clear pcat
cnt = 0; % cnt trials
sess = 0;
trialnum = 60;
for r = 1:size(bname,1)
    load (bname(r).name);

    cumtrial = ((r-1)*trialnum )+1:trialnum *r; % cnt trials across all blocks

    
    for ci = 1:27
    pcat.(condi{ci})(cumtrial, :) = p.(condi{ci});
    end
    

    pcat.eachframe(cumtrial, :) = p.fdur; 
    pcat.stimcode(cumtrial, :) = (1:trialnum )';
    pcat.stimcodecum(cumtrial, :) = pcat.stimcode(cumtrial, :)+trialnum *(r-1);



end

%%%%%%%%%%%%%% check the code in epoch struct %%%%%%%%%%%%%%%%
figure; 
subplot(1, 3, 1);
plot(squeeze(epoch.eventtype)');
title ('check event type locked to stim onset');
subplot(1, 3, 2);
plot(squeeze(epoch.eventlatency))
title ('check event latency');
subplot(1, 3, 3);
plot(squeeze(epoch.eventtimepoint))
title ('check event time point');

figure; 
subplot(1, 3, 1);
plot(squeeze(epochtg.eventtype)');
title ('check event type locked to tg onset');
subplot(1, 3, 2);
plot(squeeze(epochtg.eventlatency))
title ('check event latency');
subplot(1, 3, 3);
plot(squeeze(epochtg.eventtimepoint))
title ('check event time point');


%%%%%%%%%%%%%%%%%%% threshold rejection here %%%%%%%%%%%%%%%%%%%%%%
%left/ ver low left/ ver up right/ ver low right

maxdat = squeeze(max(abs(epoch.data(67:end, :, :)), [], 2)); % eye chan 3:8
maxdat(7, :) = max(max(abs(epoch.data(1:64, :, :)), [], 1), [], 2); % all electrodes 


for ch = 1:7
    AboveThreshIndex{ch} = find(maxdat(ch, :) >= eyeThresh(ch));
end

artindex = unique([AboveThreshIndex{1} AboveThreshIndex{2} AboveThreshIndex{3} ...
    AboveThreshIndex{4} AboveThreshIndex{5} AboveThreshIndex{6} ...
    AboveThreshIndex{7} ]);

alltrial = 1:size(epoch.data,3);
artfreeindex = alltrial(~ismember(alltrial, artindex));
eegEpochArtFree.data = epoch.data(:, :, artfreeindex);
eegEpochArtFree.type = epoch.eventtype( artfreeindex,:);
eegEpochArtFree.latency = epoch.eventlatency(artfreeindex, :);
eegEpochArtFree.timepoint = epoch.eventtimepoint(artfreeindex, :);
eegEpochArtFreebad = epoch.data(67:end, :, :);

artfreetg = intersect(epoch.eventtype(artfreeindex, :), epochtg.eventtype);
for ai = 1:numel(artfreetg)
    artfreeindextg(ai) = find(epochtg.eventtype== artfreetg(ai));
end

eegEpochArtFree.datatg = epochtg.data(:, :, artfreeindextg);
eegEpochArtFree.typetg = epochtg.eventtype( artfreeindextg,:);
eegEpochArtFree.latencytg = epochtg.eventlatency(artfreeindextg, :);
eegEpochArtFree.timepointtg = epochtg.eventtimepoint(artfreeindextg, :);


eegEpochArtFreerm.data = epochrm.data(:, :, artfreeindex);
eegEpochArtFreerm.datatg = epochtgrm.data(:, :, artfreeindextg);
eegEpochArtFreerm.type = epoch.eventtype( artfreeindex,:);
eegEpochArtFreerm.latency = epoch.eventlatency(artfreeindex, :);
eegEpochArtFreerm.timepoint = epoch.eventtimepoint(artfreeindex, :);
eegEpochArtFreerm.typetg = epochtg.eventtype( artfreeindextg,:);
eegEpochArtFreerm.latencytg = epochtg.eventlatency(artfreeindextg, :);
eegEpochArtFreerm.timepointtg = epochtg.eventtimepoint(artfreeindextg, :);
%%%% plot threshold %%% 
eyechan = {'hor left', 'hor right', 'ver up left' 'ver down left', 'ver up right', 'ver down right'};
figure;
for ch = 1:6
    subplot(4, 2, ch);
    plot(squeeze(eegEpochArtFreebad(ch,: , :)), 'r');
    hold on;
    plot(squeeze(eegEpochArtFree.data(66+ch, :, :)), 'b');
    plot(1:3000, eyeThresh(ch), 'k-');
    plot(1:3000, -eyeThresh(ch), 'k-');
    title (eyechan{ch});
end

subplot(4, 2, 7);
plot(maxdat(7, :), 'r.')
hold on;
plot(artfreeindex, maxdat(7, artfreeindex), 'b.')
plot(1:4000, eyeThresh(7), 'k-');
title ('all electrode')
ylim([0 500]);


eegEpochArtFree.RejPercent = 100-(size(eegEpochArtFree.data,3)*100./(2*960))
eegEpochArtFreerm.RejPercent = 100-(size(eegEpochArtFree.data,3)*100./(2*960))
%%%%%%%%%%%%%%%

clear p epoch 
epoch = eegEpochArtFree;
epochrm = eegEpochArtFreerm;
epoch.time = eeg.times;
epoch.timetg = eegtg.times;
epochrm.time = eegrm.times;
epochrm.timetg = eegtgrm.times;

p = pcat;


saveornot = input ('save? yes = 1/ no = 0'); 
if saveornot ==1 
save( sname, 'epoch',  'p', 'eyeThresh');
save( snamerm, 'epochrm',  'p', 'eyeThresh');
save(['sbj' num2str(sub) '_eyethreshold_noICA'], 'eyeThresh');
save(['sbj' num2str(sub) '_p_behav'], 'p');
end
cd ../..
