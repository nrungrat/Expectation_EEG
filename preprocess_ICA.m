clear all; close all;
sub = input ('sbj: ');
cd data/
cd (['sbj' num2str(sub)])

fname = dir('*.bdf');

% % then open EEGlab-GUI to reject data b/w blocks for nicer ICA
% % type **eeglab** to open the GUI
for r = 1:size(fname,1)

        clear eegtemp event
        eegtemp = rbdf_trip(fname(r).name); % load bdf/ ref to left and right mastriods / band pass filter/
        if (sub == 3 && r <= 16) || (sub == 4 && r >16) ...% replace bad electrode
        || (sub == 6 && r >16) || (sub == 8) ||  (sub == 9 && r <=16) 
        eegtemp.data(32, :) = eegtemp.data(64+6, :);
        eegtemp.data(64+6, :) = eegtemp.data(64+8, :);
        elseif (sub == 15 && r <=16)
        eegtemp.data(64, :) = eegtemp.data(64+6, :);
        eegtemp.data(64+6, :) = eegtemp.data(64+8, :);    
        end
        event = squeeze(cell2mat(struct2cell(eegtemp.event)));
        onset = find(ismember(event(1, :), 1:60)==1); % find the onset 1:60
        othercode = find(ismember(event(1, :), 61:1000)==1);
        event(1, othercode) = event(1, othercode)+ 6000; % add 6000 for other code
        event(1, onset) = event(1, onset)+ 60*(r-1); % accumalate the order of the onset label
        eegtemp.event =cell2struct(mat2cell(event, ones(1,size(event,1)), ones(1,size(event,2))), ...
            {'type', 'latency', 'urevent'}, 1); % convert back to structure;
        
        if r == 1
            eegcat = eegtemp;
        else
            eegcat = pop_mergeset(eegcat, eegtemp);
        end
        
end


for t = 1:60*size(fname,1)
  eventtype{t} = num2str(t);
end

EEG = pop_epoch(eegcat, eventtype , [-1.5 4]); %epoch data
EEGica = pop_runica(EEG, 'icatype', 'runica'); % run ica
EEG = pop_saveset(EEGica, ['sbj' num2str(sub) '_epochICA']);


cd ../..