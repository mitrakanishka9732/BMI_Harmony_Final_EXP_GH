clear all
close all
clc
%% Load data

sub_num = 5;
trial_num = 1; 
addpath(genpath('./functions'));
FILENAME = sprintf('./test_data/S%i/Exp4_Sub_%i_trialRun_%i.gdf', sub_num, sub_num, trial_num);

%samples x 
[data, hdr] = sload(FILENAME); 
labels=hdr.EVENT.POS(:,:);
labels_typ=hdr.EVENT.TYP(:,:);

data_1= data(1:labels(end-1),:);


%% Delete Channels 

clc; 
%first 32 electrodes 
new_data_2 = data_1(:,1:32); 
%new_data_2 = new_data - mean(new_data,2); 
%removing specific channels: T7, T8, M1, M2, FP1, FPZ, FP2, O1, Oz, O2
% Keep [F7, F8, P7, P8], remove when feature selection 
new_data_2(:,[1,2,3,13,14,18,19,30,31,32]) = [];

%updating channel label list 
chLabel = hdr.Label(1:32,:); 
chLabel([1,2,3,13,14,18,19,30,31,32],:) = []; 

%load and update topoplot map for only 22 channels 
load('ch32Locations.mat'); 
ch32Locations(:,[1,2,3,13,14,18,19,30,31,32]) = [];

%% Bandpass/CAR Filter the data

%Orset et al., 2021 -> changed from 100hz to 45hz(Sam)
fc1 = 1; % first cutoff frequency in Hz 
fc2 = 45; % second cutoff frequency in Hz
fs = hdr.SampleRate; 

% normalize the frequencies
Wp = [fc1 fc2]*2/fs;
[b,a]=butter(2,Wp,'bandpass');

%filter data - user filter, instead of filtfit, because real-time
filter_sig_1 = filtfilt(b,a,new_data_2); 
%result = filter_sig_1

%Common Average Referrence(CAR) filter
car_filt = filter_sig_1 - mean(filter_sig_1,2); 
filter_sig_2 = car_filt; 
%result = filter_sig_2

%% Test spectrogram function

%electrode: C3(10), Cp1(14), PZ(19)
elec_ch = 10;

spectrogram(filter_sig_2, elec_ch, hdr, sub_num, trial_num, chLabel)


%% Load ALL data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%GRAND%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%AVERAGE%%%%%%%%%%%%%%%%%%%%%%%
C{6,1} = 0; 
L{6,1} = 0; 
L_t{6,1} = 0; 

%all runs
for i = 1:1:6
    FILENAME = sprintf('./test_data/S%i/Exp4_Sub_%i_trialRun_%i.gdf', sub_num, sub_num, i);
    [data1, hdr1] = sload(FILENAME); 
    labels1=hdr1.EVENT.POS(:,:);
    L{i,1}=hdr1.EVENT.POS(:,:);
    L_t{i,1}=hdr1.EVENT.TYP(:,:);
    C{i,1} = data1(1:labels1(end-1),:);
end 

C_1{6,1} = 0; 
%delete channels
for i = 1:1:6
    new_data_2 = C{i,1}(:,1:32); 
    new_data_2(:,[1,2,3,13,14,18,19,30,31,32]) = [];
    C_1{i,1} = new_data_2; 
end     

C_2{6,1} = 0; 
%bandpass filter 
fc1 = 1; % first cutoff frequency in Hz 
fc2 = 45; % second cutoff frequency in Hz
fs = 512; 
Wp = [fc1 fc2]*2/fs;
[b,a]=butter(2,Wp,'bandpass');

for i = 1:1:6
    C_2{i,1} = filtfilt(b,a,C_1{i,1}); 
end 

C_3{6,1} = 0; 
%CAR
for i = 1:1:6
    C_3{i,1} = C_2{i,1} - mean(C_2{i,1},2);
end 

%updating channel label list 
chLabel = hdr1.Label(1:32,:); 
chLabel([1,2,3,13,14,18,19,30,31,32],:) = []; 

%load and update topoplot map for only 22 channels 
load('ch32Locations.mat'); 
ch32Locations(:,[1,2,3,13,14,18,19,30,31,32]) = [];


%% Grand spectrogram 
clc
%electrode: C3(10), Cp1(14), PZ(19)
elec_ch = 10;

grandspec(C_3, elec_ch, L, L_t, sub_num, chLabel)


%% EPOCH EXTRACTION

%1 sec window
wsize = 1; 
%62.5msec overlap
hop = 0.0625; 

%trial counter: number of epochs per class * 20 
x_trial = 1;
y_trial = 1; 
z_trial = 1; 
w_trial = 1; 

%iterate through the whole session 
for i = 3:7:129

    %get total number of epoch in Begin MI - not using now 
    num_epoch_bm = (floor((((labels(i+1) - labels(i))/512) - wsize)/hop)+1);
    
    %Begin MI pos
    lab_x = labels(i);
    %iterate through 2sec = 16 samples
    for j = 1:1:16
        epochs_bm(:,:,x_trial) = (filter_sig_2(lab_x:lab_x+511,:));
        lab_x = lab_x + 32;    %iterate lab_x pos by hop size*512 = 32sam
        x_trial = x_trial + 1;
    end

    %END MI
    %get total number of epoch in end MI
    num_epoch_em = (floor((((labels(i+2) - labels(i+1))/512) - wsize)/hop)+1);

    %shift window by 900ms, 0.9*512=461 
    lab_y = labels(i+1); 
    %2000ms(total slide)/62.5ms = 32/2=16 epochs per trial  
    for j = 1:1:16      
        epochs_em(:,:,y_trial) = (filter_sig_2(lab_y:lab_y+511,:));
        lab_y = lab_y + 32;    %iterate lab_x pos by hop size*512 = 32sam
        y_trial = y_trial + 1;
    end

    %Rest(move)
    %get total number of epoch in rest(still)
    num_epoch_rm = (floor((((labels(i+3) - labels(i+2))/512) - wsize)/hop)+1);

    %rest label position 
    lab_w = labels(i+2); 
    for j = 1:1:16
        epochs_rst_mov(:,:,w_trial) = (filter_sig_2(lab_w:lab_w+511,:));
        lab_w = lab_w + 32;    %iterate lab pos by hop size*512 = 32sam
        w_trial = w_trial+1; 
    end

    %Rest(still)
    %get total number of epoch in rest(still)
    num_epoch_rs = (floor((((labels(i+4) - labels(i+3))/512) - wsize)/hop)+1);

    %rest label position 
    lab_z = labels(i+3); 
    %same number of epochs as begin mi 
    for j = 1:1:12
        epochs_rst(:,:,z_trial) = (filter_sig_2(lab_z:lab_z+511,:));
        lab_z = lab_z + 32;    %iterate lab pos by hop size*512 = 32sam
        z_trial = z_trial+1; 
    end
    
    
end


%% PSD features 

wsize = 1;  %resolution 
fs = hdr.SampleRate;
tot_sze = size(epochs_bm,3)+size(epochs_em,3)+size(epochs_rst,3); 
bm_sze = size(epochs_bm,3); 
em_sze = size(epochs_em,3);
rst_sze = size(epochs_rst,3);

%22 channles x 23 freq = 506 +1(label). 
PSD_epoch(tot_sze,507) = 0;

%temp PSD features 
PSD_bm_temp(23,22) = 0; 
PSD_em_temp(23,22) = 0; 
PSD_rst_temp(23,22) = 0; 

%total PSD counter 
cnt = 1; 
%begin MI 
for i = 1:1:bm_sze
    signalOfInterest = epochs_bm(:,:,i);
    [SOIf, freq]=pwelch(signalOfInterest,wsize*fs, 0.5*wsize*fs, [], fs); 
    PSD_bm_temp = SOIf(9:31,:);
    PSD_epoch(cnt,1:506) = PSD_bm_temp(:)'; 
    PSD_epoch(cnt,507) = 1; 
    PSD_bm_temp = 0;
    cnt = cnt + 1; 
end


%end MI
for i = 1:1:em_sze
    signalOfInterest = epochs_em(:,:,i);
    [SOIf, freq]=pwelch(signalOfInterest,wsize*fs, 0.5*wsize*fs, [], fs); 
    PSD_em_temp = SOIf(9:31,:);
    PSD_epoch(cnt,1:506) = PSD_em_temp(:)'; 
    PSD_epoch(cnt,507) = 2; 
    PSD_em_temp = 0;
    cnt = cnt + 1; 
end

%rest
for i = 1:1:rst_sze
    signalOfInterest = epochs_rst(:,:,i);
    [SOIf, freq]=pwelch(signalOfInterest,wsize*fs, 0.5*wsize*fs, [], fs); 
    PSD_rst_temp = SOIf(9:31,:);
    PSD_epoch(cnt,1:506) = PSD_rst_temp(:)'; 
    PSD_epoch(cnt,507) = 3; 
    PSD_rst_temp = 0;
    cnt = cnt + 1; 
end

freq = freq(9:31,:); 

%% Fisher Score 

%class 1 vs rest
fish_score_1(1,506) = 0;
%class 2 vs rest
fish_score_2(1,506) = 0; 
%class 1 vs class 2
fish_score_3(1,506) = 0; 

PSD_bm_a=0;
PSD_bm_v=0;

PSD_em_a=0;
PSD_em_v=0;

PSD_rst_a=0;
PSD_rst_v=0; 


%iterate through all freq/channel
for i = 1:1:506
    PSD_bm_a = mean(PSD_epoch(1:bm_sze,i));
    PSD_bm_v = std(PSD_epoch(1:bm_sze,i));

    PSD_em_a = mean(PSD_epoch(bm_sze+1:(bm_sze+em_sze),i));
    PSD_em_v = std(PSD_epoch(bm_sze+1:(bm_sze+em_sze),i));

    PSD_rst_a = mean(PSD_epoch((bm_sze+em_sze)+1:(bm_sze+em_sze+rst_sze),i));
    PSD_rst_v = std(PSD_epoch((bm_sze+em_sze)+1:(bm_sze+em_sze+rst_sze),i));

    %Begin MI vs. Rest 
    abav1 = abs(PSD_bm_a - PSD_rst_a); 
    totvar1 = sqrt(PSD_bm_v^2 + PSD_rst_v^2);
    fish_score_1(1,i) = abav1/totvar1; 
    
    %END MI vs. Rest 
    abav2 = abs(PSD_em_a - PSD_rst_a); 
    totvar2 = sqrt(PSD_em_v^2 + PSD_rst_v^2);
    fish_score_2(1,i) = abav2/totvar2; 

    %Begin MI vs. END MI
    abav3 = abs(PSD_bm_a - PSD_em_a); 
    totvar3 = sqrt(PSD_bm_v^2 + PSD_em_v^2);
    fish_score_3(1,i) = abav3/totvar3; 

end 

fish_score_1_dl(23,22) = 0; 
fish_score_2_dl(23,22) = 0; 
fish_score_3_dl(23,22) = 0; 
a = 1; 

for i =1:23:506
    fish_score_1_dl(:,a) = fish_score_1(1,i:i+22)'; 
    fish_score_2_dl(:,a) = fish_score_2(1,i:i+22)'; 
    fish_score_3_dl(:,a) = fish_score_3(1,i:i+22)'; 
    a=a+1; 
end 

%histogram(fish_score_1); 


%% Top 10 Features - Begin MI/Rest

Nmax = 10; % get Nmax biggest entries
[ Avec, Ind ] = sort(fish_score_1(:),1,'descend');
max_values = Avec(1:Nmax);
[ ind_row, ind_col ] = ind2sub(size(fish_score_1),Ind(1:Nmax)); % fetch indices

elec_idx = fix((ind_col(:)-1)/23); 
frq_idx = ind_col(:) - 23*(elec_idx); 

disp("Begin MI/Rest: Top Ten Features")
disp("Band: " + (freq(frq_idx,1)) +"Hz, Channel: "+ ...
    chLabel(elec_idx+1,1) + '(' +ind_col(:)+')'+" Fisher Value: " + max_values(:))


%% Top 10 Features - END MI/Rest

Nmax = 10; % get Nmax biggest entries
[ Avec, Ind ] = sort(fish_score_2(:),1,'descend');
max_values1 = Avec(1:Nmax);
[ ind_row1, ind_col1 ] = ind2sub(size(fish_score_2),Ind(1:Nmax)); % fetch indices

elec_idx1 = fix((ind_col1(:)-1)/23); 
frq_idx1 = ind_col1(:) - 23*(elec_idx1); 

disp("END MI/Rest: Top Ten Features")
disp("Band: " + (freq(frq_idx1,1)) +"Hz, Channel: "+ ...
    chLabel(elec_idx1+1,1) + '(' +ind_col1(:)+')'+" Fisher Value: " + max_values1(:))


%% Top 10 Features - Begin MI/END MI

Nmax = 10; % get Nmax biggest entries
[ Avec, Ind ] = sort(fish_score_3(:),1,'descend');
max_values2 = Avec(1:Nmax);
[ ind_row2, ind_col2 ] = ind2sub(size(fish_score_3),Ind(1:Nmax)); % fetch indices

elec_idx2 = fix((ind_col2(:)-1)/23); 
frq_idx2 = ind_col2(:) - 23*(elec_idx2); 

disp("Begin MI/End MI: Top Ten Features")
disp("Band: " + (freq(frq_idx2,1)) +"Hz, Channel: "+ ...
    chLabel(elec_idx2+1,1) + '(' +ind_col2(:)+')'+" Fisher Value: " + max_values2(:))


%% Extract top PSD features

n_max = 5; 
S1_4_PSD_bm_fs(bm_sze, n_max) = 0; 
S1_4_PSD_bm_rst_fs(rst_sze, n_max) = 0; 

S1_4_PSD_em_fs(em_sze, n_max) = 0; 
S1_4_PSD_em_rst_fs(rst_sze, n_max) = 0; 

S1_4_PSD_bm_em(bm_sze, n_max) = 0; 
S1_4_PSD_em_bm(em_sze, n_max) = 0; 

%extracting the PSD from the average signals 
for i = 1:1:n_max
    %begin MI/rest
    S1_4_PSD_bm_fs(:,i) = PSD_epoch(1:bm_sze,ind_col(i)); 
    S1_4_PSD_bm_rst_fs(:,i) = PSD_epoch((bm_sze+em_sze)+1:(bm_sze+em_sze+rst_sze),ind_col(i));
    %end MI/rest
    S1_4_PSD_em_fs(:,i) = PSD_epoch(bm_sze+1:(bm_sze+em_sze),ind_col1(i)); 
    S1_4_PSD_em_rst_fs(:,i) = PSD_epoch((bm_sze+em_sze)+1:(bm_sze+em_sze+rst_sze),ind_col1(i));
    %begin MI/end MI
    S1_4_PSD_bm_em(:, i) = PSD_epoch(1:bm_sze,ind_col2(i)); 
    S1_4_PSD_em_bm(:, i) = PSD_epoch(bm_sze+1:(bm_sze+em_sze),ind_col2(i)); 

end 


%% Saving Data

%save('S1_4_PSDfeature_exp2.mat', ...
%    'S1_4_PSD_bm_fs','S1_4_PSD_bm_rst_fs', ...
%    'S1_4_PSD_em_fs','S1_4_PSD_em_rst_fs', ...
%    'S1_4_PSD_bm_em','S1_4_PSD_em_bm')

%% Plot topo
figure('units','normalized','Position',[0.1,0.1,0.5,0.5])

A = mean(fish_score_1_dl(1:5,:),1);
A_2 = mean(fish_score_1_dl(6:23,:),1);
B = mean(fish_score_2_dl(1:5,:),1);
B_2 = mean(fish_score_2_dl(6:23,:),1);
C = mean(fish_score_3_dl(1:5,:),1);
C_2 = mean(fish_score_3_dl(6:23,:),1);

tiledlayout(3,3)
title('')
nexttile
x1 = topoplot(A,ch32Locations);
title('Mu(8-12Hz): Onset: Begin MI vs. Rest')
colorbar 
caxis([min(A(:)) max(A(:))])
nexttile
x3 = topoplot(C,ch32Locations);
title('Mu(8-12Hz): Termination: Begin MI vs. End MI')
colorbar
caxis([min(C(:)) max(C(:))])
nexttile
x2 = topoplot(B,ch32Locations);
title('Mu(8-12Hz): End MI vs. Rest')
colorbar
caxis([min(B(:)) max(B(:))])

nexttile
x4 = topoplot(A_2,ch32Locations);
title('Beta(13-30Hz): Onset: Begin MI vs. Rest')
colorbar
caxis([min(A_2(:)) max(A_2(:))])
nexttile
x5 = topoplot(C_2,ch32Locations);
title('Beta(13-30Hz): Termination: Begin MI vs. End MI')
colorbar
caxis([min(C_2(:)) max(C_2(:))])
nexttile
x6 = topoplot(B_2,ch32Locations);
title('Beta(13-30Hz): End MI vs. Rest')
colorbar
caxis([min(B_2(:)) max(B_2(:))])

sgtitle("Final Experiment: Subject "+num2str(sub_num)+": Trial "+num2str(trial_num)+": Topoplots");

