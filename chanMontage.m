function [reorder  chanNum chanLabel chanNumFlip goodChan chanCon chanIps poz_mid] = chanMontage
chanLabel= {'CL', 'VLU', 'VLD', 'fp1','fpz','fp2','VRD','VRU','CR', ...
    '', '', 'af7','af3', 'afz', 'af4', 'af8', '', '', ...
    'f7', 'f5', 'f3', 'f1', 'fz', 'f2', 'f4', 'f6', 'f8', ...
    'ft7', 'fc5', 'fc3', 'fc1', 'fcz', 'fc2', 'fc4', 'fc6', 'ft8', ...
    't7', 'c5', 'c3', 'c1', 'cz', 'c2', 'c4', 'c6', 't8', ...
    'tp7', 'cp5', 'cp3', 'cp1', 'cpz', 'cp2', 'cp4', 'cp6', 'tp8', ...
    'p7', 'p5', 'p3', 'p1', 'pz', 'p2', 'p4', 'p6', 'p8' ,...
    'ML', 'p9', 'po7', 'po3', 'poz', 'po4', 'po8', 'p10', 'MR', ...
    '', '', '', 'o1' 'oz' 'o2', '', '', 'Iz'};

chanNum = [67 69 70 1 33 34 72 71 68 ...
    0 0 2 3 37 36 35 0 0 ...
    7 6 5 4 38 39 40 41 42 ...
    8 9 10 11 47 46 45 44 43 ...
    15 14 13 12 48 49 50 51 52 ...
    16 17 18 19 32 56 55 54 53 ...
    23 22 21 20 31 57 58 59 60 ...
    65 24 25 26 30 62 63 64 66 ...
    0   0  0 27 29 61 0  0  28];

chanNumFlip = [ 68 71 72 34 33 1 70 69 67 ...   
    0 0 35 36  37 3 2 0 0 ...
    42 41 40 39 38 4 5 6 7  ...
    43 44 45 46 47 11 10 9 8   ...
    52 51 50 49 48 12 13 14 15  ...
    53 54 55 56 32 19 18 17 16 ...
    60 59 58 57 31 20 21 22 23 ...
    66 64 63 62 30 26 25 24 65 ...
    0 0 0 61 29 27 0 0 28];

for k = 1:72
    reorder(k) = chanNumFlip(find(chanNum==k));
end

chanCon =  [63 59 60]; % peak p1 %[61 62 63];
chanIps = [25 22 23];
poz_mid = [20 31 57];

global ENVIRONMENT % evironment for topo
ENVIRONMENT = 'Biosemi64';
goodChan = [1:64]; % good channel for topo plot