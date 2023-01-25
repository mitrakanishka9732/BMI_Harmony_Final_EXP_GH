function spectrogram(signal,chan, hdr, sub_num, trial_num, chLabel)
%Calculating and plotting the spectrogram 

%1 sec window
wsize = 1; 
%62.5msec overlap
hop = 0.0625; 
%extract labels
labels=hdr.EVENT.POS(:,:);
types=hdr.EVENT.TYP(:,:);
fs = hdr.SampleRate;

%trial counter: number of epochs per class * 20 
x_trial = 1;

%%go through 19 trials   %%%%%%%% need to change !!!!!!!
for i = 2:7:128
    

    % MI TRIALS(300) %%%%%%%%%%%%%%%%%%%%%%%
    if (types(i,1) == 300) 

        %COUNT DOWN 
        lab_x = labels(i);  %label for countdown 
        
        % 8 epochs = 1 sec, 144 epochs = 18 sec 
        for a = 1:1:144
            [SOIF0, freq0] = pwelch(signal(lab_x:lab_x+511,chan),wsize*fs, 0.5*wsize*fs, [], fs);
            epochs_freq(:,a,x_trial) = SOIF0(9:31,:);
            lab_x = lab_x + 32;    %iterate lab_x pos by hop size*512 = 32sam
        end
    

        %BASELINE: Rest(still) - 2sec after Rest(still), 12 samples, 1.5
        %epoch
        lab_w = labels(i+6);
        for l = 1:1:12
            [SOIF1, freq1] = pwelch(signal(lab_w:lab_w+511,chan),wsize*fs, 0.5*wsize*fs, [], fs);
            epochs_freq_b(:,l,x_trial) = SOIF1(9:31,:);
            lab_w = lab_w + 32;    %iterate lab_x pos by hop size*512 = 32sam
        end

        x_trial = x_trial + 1;
    end

end 

%MI TRIALS
mean_epochs_freq = mean(epochs_freq,3); 
PSD_bm_freq_2 = flip(mean_epochs_freq);
%baseline
mean_epochs_freq_b = mean(epochs_freq_b,3); 
PSD_baseline = flip(mean_epochs_freq_b);
PSD_baseline_2 = mean(PSD_baseline,2); 



%LOG!!
%MI: after baseline normalization 
PSD_M(23,144) = 0; 
for i = 1:1:23
    for j = 1:1:144
        PSD_M(i,j) = 10*log(PSD_bm_freq_2(i,j)/PSD_baseline_2(i,1));
    end 
end 


%PLOT!! 
%MI 
imagesc(PSD_M);
colorbar
caxis([min(PSD_M(:)) max(PSD_M(:))])
colormap('jet')
xticks([])
yticks(1:1:23)
xline(0,'-black','Countdown')
xline(24,'-black','Begin MI')
xline(32,'-black','Robot Moves')
xline(48,'-black','End MI')
xline(56,'-black','Robot Stops')
xline(104,'-black','Rest(move)')
xline(128,'-black','Rest(still)')
yticklabels({'30','29','28','27','26','25','24','23','22','21','20','19','18','17','16','15','14','13','12','11','10','9','8','7','6','5','4','3','2','1'})
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title("MI: Spectrogram: Subject "+num2str(sub_num)+": Run "+num2str(trial_num)+": "+chLabel(chan,1))


chan_name = string(chLabel(chan,1)); 
FILENAME = sprintf('./EXP_Figs/Spec_sub_%i_run_%i_ch%s.jpg', sub_num, trial_num, chan_name);
saveas(gcf,FILENAME)

end