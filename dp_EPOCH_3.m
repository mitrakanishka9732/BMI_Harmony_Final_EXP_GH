clear all
close all
clc
acc = zeros(13,6);
tot_accgraph{13,2} = 0;

%% Load data

inc_sub = [4 5 6 7 8 9 10 11 12 13];
%inc_sub = [4]; 

feature_heat_on = zeros(22,23); 
feature_heat_of = zeros(22,23); 
addpath(genpath('./functions'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% subjects
for subs = 1:1:length(inc_sub)
sub_num = inc_sub(subs);

%%%%%%%%%%%%%%%%%%%%%%import calibration data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FILENAME1 = sprintf('./calib_data/Exp4_Sub_%i_Asynchronous_calibration.gdf', sub_num);
[c_data, c_hdr] = sload(FILENAME1); 
fs = c_hdr.SampleRate; 
c_data1 = c_data(5*fs:95*fs,:); 
c_data1_eeg = c_data1(:, 1:32); 
c_data1_eog = c_data1(:, 33:35); 
%c_data1_eeg(:,[1,2,3,13,14,18,19,30,31,32]) = [];
fc1 = 0.1; % first cutoff frequency in Hz 
fc2 = 45; % second cutoff frequency in Hz
% normalize the frequencies
Wp = [fc1 fc2]*2/fs;
[b,a]=butter(2,Wp,'bandpass');
%filter data - user filter, instead of filtfit, because real-time
c_data1_eeg_filt = filtfilt(b,a,c_data1_eeg); 
c_data1_eog_filt = filtfilt(b,a,c_data1_eog); 
%EOG covariance matrix 
EOGfilter = filterEOG(c_data1_eeg_filt, c_data1_eog_filt);

%%%%%%%%%%%%%%%%%%%%%% end calibration data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('./functions'));
PSD_big{6,1}=0; 

tot_acc1 = zeros(6,2);
tot_acc2 = zeros(6,2); 
tot_accgraph1 = zeros(6,9);
tot_accgraph3 = zeros(6,9); 
tot_chacc = zeros(6,2); 
tot_chacc_on = zeros(10,6); 
tot_chacc_of = zeros(10,6); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%all runs(x)
for x = 1:1:6

    FILENAME = sprintf('./test_data/S%i/Exp4_Sub_%i_trialRun_%i.gdf', sub_num, sub_num, x);
    
    [data, hdr] = sload(FILENAME); 
    labels=hdr.EVENT.POS(:,:);
    labels_typ=hdr.EVENT.TYP(:,:);

    data_1 = data(1:labels(end),:);

    %%%%%%%%%%%%%%%%%%%%% Delete Channels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %get EEG data: first 32 electrodes 
    new_data_2 = data_1(:,1:32); 
    eog_data = data_1(:,33:35); 

    %updating channel label list 
    chLabel = hdr.Label(1:32,:); 
    chLabel([1,2,3,13,14,18,19,30,31,32],:) = []; 

    %load and update topoplot map for only 22 channels 
    load('ch32Locations.mat'); 
    ch32Locations(:,[1,2,3,13,14,18,19,30,31,32]) = [];

    %%%%%%%%%%%%%%%%%%%%% Bandpass/CAR Filter the data %%%%%%%%%%%%%%%%%%%%
    %Orset et al., 2021 -> changed from 100hz to 45hz(Sam)
    fc1 = 0.1; % first cutoff frequency in Hz 
    fc2 = 45; % second cutoff frequency in Hz
    fs = hdr.SampleRate; 

    % normalize the frequencies
    Wp = [fc1 fc2]*2/fs;
    [b,a]=butter(2,Wp,'bandpass');

    %filter data - user filter, instead of filtfit, because real-time
    filter_sig_1 = filtfilt(b,a,new_data_2); 
    eog_filt = filtfilt(b,a,eog_data); 
    %result -> filter_sig_1

    %remove EOG 
    filter_sig_2 = filter_sig_1 - eog_filt*EOGfilter;
    
    %removing specific channels: T7, T8, M1, M2, FP1, FPZ, FP2, O1, Oz, O2
    %removing specific channels: T7, T8, M1, M2, FP1, FPZ, FP2, O1, Oz, O2
    % Keep [F7, F8, P7, P8], remove when feature selection
    filter_sig_2(:,[1,2,3,13,14,18,19,30,31,32]) = [];

    %Common Average Referrence(CAR) filter
    filter_sig_3 = filter_sig_2 - mean(filter_sig_2,2); 
    %result -> filter_sig_3

    %%%%%%%%%%%%%%%%%%%%% Epoch Extraction %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %0.5 sec window
wsize = 0.5; 
%62.5msec overlap
hop = 0.0625; 
%only 9 epochs per trial
nr_w_c1 = 1/hop-7;
nr_w_c2 = 1/hop-7;

%trial counter: number of epochs per class * 20 
x_trial = 1;
y_trial = 1; 
z_trial = 1; 
w_trial = 1; 

ind_bm = labels(labels_typ==100 | labels_typ==101);
ind_MI = labels(labels_typ==150 | labels_typ==151); 
ind_em = labels(labels_typ==500 | labels_typ==501);
ind_rst = labels(labels_typ==950 | labels_typ==951);
ind_cnd = labels(labels_typ==300 | labels_typ==301);

epochs_bm = zeros(round(fs*wsize),length(chLabel),nr_w_c1*length(ind_bm)-nr_w_c1); 
epochs_mi = zeros(round(fs*wsize),length(chLabel),nr_w_c2*length(ind_bm)-nr_w_c2); 
epochs_em = zeros(round(fs*wsize),length(chLabel),nr_w_c2*length(ind_bm)-nr_w_c2); 
epochs_rst = zeros(round(fs*wsize),length(chLabel),nr_w_c1*length(ind_bm)-nr_w_c1); 

%iterate through the whole session 
for i = 1:1:length(ind_bm)-1

    %get total number of epoch in Begin MI - not using now 
    %num_epoch_bm = (floor((((labels(i+1) - labels(i))/512) - wsize)/hop)+1);
    
    %Begin MI label position: T = 0
    lab_x = ind_bm(i);
    for j = 1:nr_w_c1
        epochs_bm(:,:,x_trial) = (filter_sig_3(lab_x:lab_x+(wsize*fs)-1,:));
        lab_x = lab_x + 32;    %iterate lab_x pos by hop size*512 = 32sam
        x_trial = x_trial + 1;
    end

    %MI (mantaining MI)
    %MI label position = T+1 after robot moves 
    lab_w = ind_MI(i)+(2*(wsize*fs)); 
    for j = 1:nr_w_c2
        epochs_mi(:,:,w_trial) = (filter_sig_3(lab_w:lab_w+(wsize*fs)-1,:));
        lab_w = lab_w + 32;    %iterate lab pos by hop size*512 = 32sam
        w_trial = w_trial+1; 
    end

    %END MI
    %end MI label position: T = 0
    lab_y = ind_em(i); 
    for j = 1:nr_w_c2      
        epochs_em(:,:,y_trial) = (filter_sig_3(lab_y:lab_y+(wsize*fs)-1,:));
        lab_y = lab_y + 32;    %iterate lab_x pos by hop size*512 = 32sam
        y_trial = y_trial + 1;
    end

    %Rest(still)
    %rest label position:  
    lab_z = ind_cnd(i)-(2*(wsize*fs)); 
    for j = 1:nr_w_c1
        epochs_rst(:,:,z_trial) = (filter_sig_3(lab_z:lab_z+(wsize*fs)-1,:));
        lab_z = lab_z + 32;    %iterate lab pos by hop size*512 = 32sam
        z_trial = z_trial+1; 
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%% PSD features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = hdr.SampleRate;
tot_sze = size(epochs_bm,3)+size(epochs_mi,3)+size(epochs_em,3)+size(epochs_rst,3); 

bm_sze = size(epochs_bm,3); 
mi_sze = size(epochs_mi,3); 
em_sze = size(epochs_em,3);
rst_sze = size(epochs_rst,3);

%22 channles x 23 freq = 506 +1(label). 
PSD_epoch = zeros(tot_sze,507); 

%temp PSD features 
PSD_bm_temp(23,22) = 0; 
PSD_mi_temp(23,22) = 0; 
PSD_em_temp(23,22) = 0; 
PSD_rst_temp(23,22) = 0; 

%total PSD counter 
cnt = 1; 
%begin MI 
for i = 1:1:bm_sze
    signalOfInterest = epochs_bm(:,:,i);
    [SOIf, freq]=pwelch(signalOfInterest,wsize*fs, round(0.4*wsize*fs), 4:1:40, fs); 
    PSD_bm_temp = SOIf(freq >= 8 & freq <= 30,:);
    PSD_epoch(cnt,1:506) = PSD_bm_temp(:)'; 
    PSD_epoch(cnt,507) = 1; %label (1) for easy of use 
    PSD_bm_temp = 0;
    cnt = cnt + 1; 
end

%maintain MI
for i = 1:1:mi_sze
    signalOfInterest = epochs_mi(:,:,i);
    [SOIf, freq]=pwelch(signalOfInterest,wsize*fs, round(0.4*wsize*fs), 4:1:40, fs);  
    PSD_mi_temp = SOIf(freq >= 8 & freq <= 30,:);
    PSD_epoch(cnt,1:506) = PSD_mi_temp(:)'; 
    PSD_epoch(cnt,507) = 2; %label (2) for easy of use 
    PSD_mi_temp = 0;
    cnt = cnt + 1; 
end

%end MI
for i = 1:1:em_sze
    signalOfInterest = epochs_em(:,:,i);
    [SOIf, freq]=pwelch(signalOfInterest,wsize*fs, round(0.4*wsize*fs), 4:1:40, fs);  
    PSD_em_temp = SOIf(freq >= 8 & freq <= 30,:);
    PSD_epoch(cnt,1:506) = PSD_em_temp(:)'; 
    PSD_epoch(cnt,507) = 3; %label (3) for easy of use 
    PSD_em_temp = 0;
    cnt = cnt + 1; 
end

%rest
for i = 1:1:rst_sze
    signalOfInterest = epochs_rst(:,:,i);
    [SOIf, freq]=pwelch(signalOfInterest,wsize*fs, round(0.4*wsize*fs), 4:1:40, fs);  
    PSD_rst_temp = SOIf(freq >= 8 & freq <= 30,:);
    PSD_epoch(cnt,1:506) = PSD_rst_temp(:)'; 
    PSD_epoch(cnt,507) = 4; %label (4) for easy of use 
    PSD_rst_temp = 0;
    cnt = cnt + 1; 
end

        PSD_big{x,1}=PSD_epoch; 

    
end 

%everything else
nr_tr = [1,2,3,4,5,6]; 

for x = 1:1:6

    %remove testing run
    nr_tr(:,[x]) = [];
    PSD_train = [PSD_big{nr_tr(1),1}; PSD_big{nr_tr(2),1}; PSD_big{nr_tr(3),1}; PSD_big{nr_tr(4),1}; PSD_big{nr_tr(5),1}];
    PSD_test = PSD_big{x,1};
    freq = [8:30]'; 

    % Fisher Score 

    %class 1(begin MI) vs rest
fish_score_1(1,506) = 0;
%class 2(End MI) vs Maintain MI
fish_score_2(1,506) = 0; 

PSD_bm_a=0;
PSD_bm_v=0;

PSD_mi_a=0;
PSD_mi_v=0;

PSD_em_a=0;
PSD_em_v=0;

PSD_rst_a=0;
PSD_rst_v=0; 


%iterate through all freq/channel
for i = 1:1:506
    PSD_bm_a = mean(PSD_train(PSD_train(:,507) == 1,i));
    PSD_bm_v = std(PSD_train(PSD_train(:,507) == 1,i));

    PSD_mi_a = mean(PSD_train(PSD_train(:,507) == 2,i));
    PSD_mi_v = std(PSD_train(PSD_train(:,507) == 2,i));

    PSD_em_a = mean(PSD_train(PSD_train(:,507) == 3,i));
    PSD_em_v = std(PSD_train(PSD_train(:,507) == 3,i));

    PSD_rst_a = mean(PSD_train(PSD_train(:,507) == 4,i));
    PSD_rst_v = std(PSD_train(PSD_train(:,507) == 4,i));

    %Begin MI vs. Rest 
    abav1 = abs(PSD_bm_a - PSD_rst_a); 
    totvar1 = sqrt(PSD_bm_v^2 + PSD_rst_v^2);
    fish_score_1(1,i) = abav1/totvar1; 
    
    %END MI vs. Maintain MI 
    abav2 = abs(PSD_em_a - PSD_mi_a); 
    totvar2 = sqrt(PSD_em_v^2 + PSD_mi_v^2);
    fish_score_2(1,i) = abav2/totvar2; 
end

% Top Features - Begin MI/Rest
Nmax = 10; % get Nmax biggest entries
[ Avec, Ind ] = sort(fish_score_1(:),1,'descend');
max_values = Avec(1:Nmax);
[ ind_row, ind_col ] = ind2sub(size(fish_score_1),Ind(1:Nmax)); % fetch indices

elec_idx = fix((ind_col(:)-1)/23)+1;
frq_idx = 23 - (23*(elec_idx) - ind_col(:));

% Top Features - Maintain MI/End MI
[ Avec, Ind ] = sort(fish_score_2(:),1,'descend');
max_values2 = Avec(1:Nmax);
[ ind_row2, ind_col2 ] = ind2sub(size(fish_score_2),Ind(1:Nmax)); % fetch indices

elec_idx2 = fix((ind_col2(:)-1)/23)+1;
frq_idx2 = 23 - (23*(elec_idx2) - ind_col2(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Feature Heatmap! %%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:1:10
    feature_heat_on(elec_idx(i),frq_idx(i)) = feature_heat_on(elec_idx(i),frq_idx(i)) +1;
    feature_heat_of(elec_idx2(i),frq_idx2(i)) = feature_heat_of(elec_idx2(i),frq_idx2(i)) +1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Top Features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n_max = 10; 
Train_PSD_bm_rst = zeros(length(PSD_train(PSD_train(:,507) == 1)), n_max);
Train_PSD_rst_bm = zeros(length(PSD_train(PSD_train(:,507) == 4)), n_max);
Test_PSD_bm_rst = zeros(length(PSD_test(PSD_test(:,507) == 1)), n_max);
Test_PSD_rst_bm = zeros(length(PSD_test(PSD_test(:,507) == 4)), n_max);

Train_PSD_mi_em = zeros(length(PSD_train(PSD_train(:,507) == 2)), n_max); 
Train_PSD_em_mi = zeros(length(PSD_train(PSD_train(:,507) == 3)), n_max);
Test_PSD_mi_em = zeros(length(PSD_test(PSD_test(:,507) == 2)), n_max); 
Test_PSD_em_mi = zeros(length(PSD_test(PSD_test(:,507) == 3)), n_max); 

%extracting the PSD from the average signals 
for i = 1:1:n_max
    %begin MI/rest
    Train_PSD_bm_rst(:,i) = PSD_train(PSD_train(:,507) == 1,ind_col(i)); 
    Train_PSD_rst_bm(:,i) = PSD_train(PSD_train(:,507) == 4,ind_col(i));
    Test_PSD_bm_rst(:,i) = PSD_test(PSD_test(:,507) == 1, ind_col(i));
    Test_PSD_rst_bm(:,i) = PSD_test(PSD_test(:,507) == 4, ind_col(i));

    %maintain MI/end MI
    Train_PSD_mi_em(:, i) = PSD_train(PSD_train(:,507) == 2,ind_col2(i)); 
    Train_PSD_em_mi(:, i) = PSD_train(PSD_train(:,507) == 3,ind_col2(i));
    Test_PSD_mi_em(:, i) = PSD_test(PSD_test(:,507) == 2, ind_col2(i)); 
    Test_PSD_em_mi(:, i) = PSD_test(PSD_test(:,507) == 3, ind_col2(i));

end 

% Class1(Begin MI) vs Rest dataset
    Tr_Data_Class1vsRest(:,1:10)= [Train_PSD_bm_rst(:,:); Train_PSD_rst_bm(:,:)];
    Tr_Data_Class1vsRest(:,11)= [ones(1,length(Train_PSD_bm_rst(:,:))), 2*ones(1,length(Train_PSD_rst_bm(:,:)))];

    Te_Data_Class1vsRest(:,1:10)= [Test_PSD_bm_rst(:,:); Test_PSD_rst_bm(:,:)];
    Te_Data_Class1vsRest(:,11)= [ones(1,length(Test_PSD_bm_rst(:,:))), 2*ones(1,length(Test_PSD_rst_bm(:,:)))];

% %Class1(Maintain MI) vs Class2 (End MI)
    Tr_Data_Class1vsClass2(:,1:10)= [Train_PSD_mi_em(:,:); Train_PSD_em_mi(:,:)];
    Tr_Data_Class1vsClass2(:,11)= [ones(1,length(Train_PSD_mi_em(:,:))), 2*ones(1,length(Train_PSD_em_mi(:,:)))];

    Te_Data_Class1vsClass2(:,1:10)= [Test_PSD_mi_em(:,:); Test_PSD_em_mi(:,:)];
    Te_Data_Class1vsClass2(:,11)= [ones(1,length(Test_PSD_mi_em(:,:))), 2*ones(1,length(Test_PSD_em_mi(:,:)))];

    %{
    %chance level
    for i = 1:1:100
        [validationAccuracy4] = dLDA_CHANCE(Tr_Data_Class1vsRest, Te_Data_Class1vsRest); 
        tot_chacc_on(i,x) = validationAccuracy4;
        [validationAccuracy5] = dLDA_CHANCE(Tr_Data_Class1vsClass2, Te_Data_Class1vsClass2); 
        tot_chacc_of(i,x) = validationAccuracy5;
    end
    %}

    [trainedClassifier1,validationAccuracy1,trainingAccuracy1, trueLabels1, predLabels1, accgraph1] = dLDA_FINAL(Tr_Data_Class1vsRest, Te_Data_Class1vsRest); 
    [trainedClassifier3,validationAccuracy3,trainingAccuracy3, trueLabels3, predLabels3, accgraph3] = dLDA_FINAL(Tr_Data_Class1vsClass2, Te_Data_Class1vsClass2);

    tot_acc1(x,1) = validationAccuracy1;
    tot_acc2(x,1) = validationAccuracy3; 
    tot_acc1(x,2) = trainingAccuracy1;
    tot_acc2(x,2) = trainingAccuracy3; 

    tot_accgraph1(x,:) = accgraph1; 
    tot_accgraph3(x,:) = accgraph3; 
    
    %chance level- by run 
    %tot_chacc(x,1) = (C1(1,1)+C1(2,2))/(sum(C1,'all')); 
    %tot_chacc(x,2) = (C2(1,1)+C2(2,2))/(sum(C2,'all')); 

    nr_tr = [1,2,3,4,5,6];
    
end 

%Print accuracies 
avg_acc1 = mean(tot_acc1(:,1)); 
acc1 = ['Sub:' num2str(sub_num) ' Accuracy: Onset: Begin MI vs. Rest: ', num2str(avg_acc1*100)];
disp(acc1)
avg_acc2 = mean(tot_acc2(:,1)); 
acc2 = ['Sub:' num2str(sub_num) ' Accuracy: Offset: Maintain MI vs. END MI: ', num2str(avg_acc2*100)];
disp(acc2)

%testing accuracy by subject 
acc(sub_num,1) = avg_acc1*100; 
acc(sub_num,2) = avg_acc2*100; 

%training level by subject 
acc(sub_num,3) = mean(tot_acc1(:,2))*100; 
acc(sub_num,4) = mean(tot_acc2(:,2))*100; 

%chance level
m_CA_on = mean(tot_chacc_on); 
m_CA_of = mean(tot_chacc_of); 
acc(sub_num,5) = (mean(m_CA_on) + (2*std(m_CA_on)))*100; 
acc(sub_num,6) = (mean(m_CA_of) + (2*std(m_CA_of)))*100; 

%acc(sub_num,5) = (mean(tot_chacc(:,1)) + (2*std(tot_chacc(:,1)))) *100; 
%acc(sub_num,6) = (mean(tot_chacc(:,1)) + (2*std(tot_chacc(:,1)))) *100; 

%%%%% for the accuracy graphs 
tot_accgraph{sub_num,1} = tot_accgraph1;
tot_accgraph{sub_num,2} = tot_accgraph3;

end 


%% Creat final Table 

OnsetTraining = acc(4:13,3); 
OnsetTesting = acc(4:13,1); 
OffsetTraining = acc(4:13,4); 
OffsetTesting = acc(4:13,2); 
OnsetChance = acc(4:13,5);
OffsetChance = acc(4:13,6);
results = table(OnsetTraining,OnsetTesting,OffsetTraining, OffsetTesting, OnsetChance, OffsetChance); 

%% grand average acc graph 
gradAVG_accgraph_on = zeros(10,9); 
gradAVG_accgraph_of = zeros(10,9); 

l = 1; 
for i = 4:1:13
    gradAVG_accgraph_on(l,:) = mean(tot_accgraph{i,1},1); 
    gradAVG_accgraph_of(l,:) = mean(tot_accgraph{i,2},1); 
    l = l+1; 
end 

m_gradAVG_accgraph_on = mean(gradAVG_accgraph_on,1); 
v_gradAVG_accgraph_on = std(gradAVG_accgraph_on,1);

m_gradAVG_accgraph_of = mean(gradAVG_accgraph_of,1); 
v_gradAVG_accgraph_of = std(gradAVG_accgraph_of,1);

save acctime_fig_3.mat m_gradAVG_accgraph_on v_gradAVG_accgraph_on m_gradAVG_accgraph_of v_gradAVG_accgraph_of

%% plot acc by epoch 

%x-axis 
x1 = 0.5:0.0625:1;

subplot(2,1,1);
shadedErrorBar(x1,m_gradAVG_accgraph_on,v_gradAVG_accgraph_on,'-b')
xlabel('Time(s)'); ylabel('Accuracy (%)');
title("Grand Average: Onset Accuracy across Time") 

subplot(2,1,2); 
shadedErrorBar(x1,m_gradAVG_accgraph_of,v_gradAVG_accgraph_of, '-r')
xlabel('Time(s)'); ylabel('Accuracy (%)');
title("Grand Average: Offset Accuracy across Time") 

%% Top features heatmap

save heatmap_fig_4.mat feature_heat_on feature_heat_of

%plot control 
subplot(1,2,1); 
imagesc(feature_heat_on);
colorbar
caxis([min(feature_heat_on(:)) max(feature_heat_on(:))])
colormap('jet')
xticks(1:1:23)
xticklabels({'8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'})
xtickangle(0)
yticks(1:1:22)
yticklabels({'F7','F3','FZ','F4','F8','FC5','FC1','FC2','FC6','C3','CZ','C4','CP5','CP1','CP2','CP6','P7','P3','PZ','P4','P8','POZ'})
xlabel('Frequency (Hz)'); ylabel('Channels');
title("Top Features for Onset Classification")

%plot control 
subplot(1,2,2); 
imagesc(feature_heat_of);
colorbar
caxis([min(feature_heat_of(:)) max(feature_heat_of(:))])
colormap('jet')
xticks(1:1:23)
xticklabels({'8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','23','24','25','26','27','28','29','30'})
xtickangle(0)
yticks(1:1:22)
yticklabels({'F7','F3','FZ','F4','F8','FC5','FC1','FC2','FC6','C3','CZ','C4','CP5','CP1','CP2','CP6','P7','P3','PZ','P4','P8','POZ'})
xlabel('Frequency (Hz)'); ylabel('Channels');
title("Top Features for Offset Classification")




