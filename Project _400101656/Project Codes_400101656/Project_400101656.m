%% Preprocessing

load('Subject1.mat');
load('Subject2.mat');
subject_1=table2array(subject1);
subject_2=table2array(subject2);
subject_1(:,20)=[];
subject_2(:,20)=[];
save("subject_1" , "subject_1");
save("subject_2" , "subject_2");

%% frequency spectrum of Fz

% bandfirSub1=ALLEEG(1).data;
% channelCzSub1=bandfirSub1(5 , :);
% [y1 , x1]=fftCal(channelCzSub1 , 200);
% plot(x1 , y1);
% title("Frequency Spectrum of Sub1 for channel Fz")

bandfirSub2=ALLEEG(1).data;
channelCzSub2=bandfirSub2(5 , :);
[y2 , x2]=fftCal(channelCzSub2 , 200);
plot(x2 , y2);
title("Frequency Spectrum of Sub2 for channel Fz")

%% step 3

% sub1=ALLEEG(4).data;
% % epoch is a 19*600*120 matrix
% epochedSub1=zeros(19 , 600 , 120);
% % looping on the channels
% for i=1:1:19
%     % looping on the trials
%     for j=1:1:120
%         % skip the first 14 secs , epoch starts from 7th sec and goes on
%         % for 3 secs
%         startIndex=14*200+j*7*200;
%         % 3secs * 200 samples per sec = 600 samples
%         endIndex=startIndex+600-1;
%         % extract epoched data (from the starting index to the ending
%         % index of original data)
%         epochedSub1(i,:,j)=sub1(i,startIndex:endIndex);
%     end
% end
% save("epochedSub1","epochedSub1");

sub2=ALLEEG(4).data;
% epoch is a 19*600*120 matrix
epochedSub2=zeros(19 , 600 , 120);
% looping on the channels
for i=1:1:19
    % looping on the trials
    for j=1:1:120
        % skip the first 16.4 secs , epoch starts from 7th sec and goes on
        % for 3 secs
        startIndex=16.4*200+j*7*200;
        % 3secs * 200 samples per sec = 600 samples
        endIndex=startIndex+600-1;
        % extract epoched data (from the starting index to the ending
        % index of original data)
        epochedSub2(i,:,j)=sub2(i,startIndex:endIndex);
    end
end
save("epochedSub2","epochedSub2");

%% step 4

% noisyTrials = [];
% % looping on the channels
% for i = 1:1:19
%     % p = a matrix of frequency spectrums
%     p = zeros(120 , 600);
%     % looping on the trials
%     for j = 1:1:120
%         % data of the jth trial
%         trialData = squeeze(epochedSub1(:,:,j)); 
%         % frequency spectrum of the trial
%         p(j,:) = abs(fft(trialData(i,:))).^2;  
%     end
%     % finding noisy trials (due to the given condition)
%     vr = sum(nanstd(p,[],2).^2,2);
%     noisy_trials = find(abs(zscore(vr))>3.5);
%     % accumulate over all channels
%     accum_noisy_trials1 = union(noisyTrials, noisy_trials);
% end
% 
% % remove noisy trials
% cleanEpochedSub1 = epochedSub1(:, :, setdiff(1:120,accum_noisy_trials1));
% save("cleanEpochedSub1","cleanEpochedSub1");

noisyTrials = [];
% looping on the channels
for i = 1:1:19
    % p = a matrix of frequency spectrums
    p = zeros(120 , 600);
    % looping on the trials
    for j = 1:1:120
        % data of the jth trial
        trialData = squeeze(epochedSub2(:,:,j)); 
        % frequency spectrum of the trial
        p(j,:) = abs(fft(trialData(i,:))).^2;  
    end
    % finding noisy trials (due to the given condition)
    vr = sum(nanstd(p,[],2).^2,2);
    noisy_trials = find(abs(zscore(vr))>3.5);
    % accumulate over all channels
    accum_noisy_trials2 = union(noisyTrials, noisy_trials);
end

% remove noisy trials
cleanEpochedSub2 = epochedSub2(:, :, setdiff(1:120,accum_noisy_trials2));
save("cleanEpochedSub2","cleanEpochedSub2");

%% step 5

load("Normal.mat")

% keep the data corresponding to the channels: Fp1,Fz,Cz,Pz
% channelSelectedEpochSub1 = cleanEpochedSub1([1 , 5 , 10 , 15] , : , :);
% save("channelSelectedEpochSub1","channelSelectedEpochSub1");
% 
% % Save as a struct
% % create an empty struct
% sub1Struct = struct();
% % add clean epoched data
% sub1Struct.selectedCleanEpoched = channelSelectedEpochSub1;
% % add Normal data's odor
% odor = normal(2).odor;
% sub1Struct.odor = odor; 
% % add the noisy trials
% sub1Struct.noisy = accum_noisy_trials1.'; 
% % save the data
% save("sub1Struct","sub1Struct");

% keep the data corresponding to the channels: Fp1,Fz,Cz,Pz
channelSelectedEpochSub2 = cleanEpochedSub2([1 , 5 , 10 , 15] , : , :);
save("channelSelectedEpochSub2","channelSelectedEpochSub2");
% Save as a struct
% create an empty struct
sub2Struct = struct();
% add clean epoched data
sub2Struct.selectedCleanEpoched = channelSelectedEpochSub2;
% add Normal data's odor
odor = normal(2).odor;
sub2Struct.odor = odor; 
% add the noisy trials
sub2Struct.noisy = accum_noisy_trials2.'; 
% save the data
save("sub2Struct","sub2Struct");

%% Processing

load("AD.mat")
load("Normal.mat")
load("MCI.mat")

%% PLV

% for AD patients
plvFreqMeanAD=zeros(13 , 1);
plvRareMeanAD=zeros(13 , 1);
% looping on the patients
for i=1:13
    % get the epoch field for each patient 
    ADepoch=getfield(AD(1 , i) , "epoch");
    % get the odor field for each patient 
    ADodor=getfield(AD(1 , i) , "odor");
    % find the trial indices at which odor is frequent
    freqIdxAD=find(ADodor(: , 1)==0);
    % find the trial indices at which odor is rare
    rareIdxAD=find(ADodor(: , 1)==1);
    % create an array for PLV of each frequent odor trial 
    plvFreqAD=zeros(size(freqIdxAD , 1),1);
    % create an array for PLV of each rare odor trial 
    plvRareAD=zeros(size(rareIdxAD , 1),1);
    % find the frequent odor trials of channels Fz & Cz
    for j=1:size(freqIdxAD , 1)
        freqChannelFzAD = ADepoch(2 , : , freqIdxAD(j , 1)); % Fz channel
        freqChannelCzAD = ADepoch(3 , : , freqIdxAD(j , 1)); % Cz channel
        % calculate PLV on frequents odors between channels Fz and Cz for each 
        % patient
        plvFreqAD(j) = PLVCal(double(freqChannelFzAD) , double(freqChannelCzAD) , 200 , [35 , 40]);
    end
    % get average on the trials
    plvFreqMeanAD(i) = mean(plvFreqAD);   

    % find the rare odor trials of channels Fz & Cz
    for k=1:size(rareIdxAD , 1)
        rareChannelFzAD = ADepoch(2 , : , rareIdxAD(k , 1)); % Fz channel
        rareChannelCzAD = ADepoch(3 , : , rareIdxAD(k , 1)); % Cz channel
        % calculate PLV on rare odors between channels Fz and Cz for each 
        % patient
        plvRareAD(k) = PLVCal(double(rareChannelFzAD) , double(rareChannelCzAD) , 200 , [35 , 40]);
    end
    % get average on the trials
    plvRareMeanAD(i) = mean(plvRareAD);
end

% for Normal patients
plvFreqMeanNormal=zeros(15 , 1);
plvRareMeanNormal=zeros(15 , 1);
% looping on the patients
for i=1:15
    % get the epoch field for each patient 
    Normalepoch=getfield(normal(1 , i) , "epoch");
    % get the odor field for each patient 
    Normalodor=getfield(normal(1 , i) , "odor");
    % find the trial indices at which odor is frequent
    freqIdxNormal=find(Normalodor(: , 1)==0);
    % find the trial indices at which odor is rare
    rareIdxNormal=find(Normalodor(: , 1)==1);
    % create an array for PLV of each frequent odor trial 
    plvFreqNormal=zeros(size(freqIdxNormal , 1),1);
    % create an array for PLV of each rare odor trial 
    plvRareNormal=zeros(size(rareIdxNormal , 1),1);
    % find the frequent odor trials of channels Fz & Cz
    for j=1:size(freqIdxNormal , 1)
        freqChannelFzNormal = Normalepoch(2 , : , freqIdxNormal(j , 1)); % Fz channel
        freqChannelCzNormal = Normalepoch(3 , : , freqIdxNormal(j , 1)); % Cz channel
        % calculate PLV on frequents odors between channels Fz and Cz for each 
        % patient
        plvFreqNormal(j) = PLVCal(double(freqChannelFzNormal) , double(freqChannelCzNormal) , 200 , [35 , 40]);
    end
    % get average on the trials
    plvFreqMeanNormal(i) = mean(plvFreqNormal);   

    % find the rare odor trials of channels Fz & Cz
    for k=1:size(rareIdxNormal , 1)
        rareChannelFzNormal = Normalepoch(2 , : , rareIdxNormal(k , 1)); % Fz channel
        rareChannelCzNormal = Normalepoch(3 , : , rareIdxNormal(k , 1)); % Cz channel
        % calculate PLV on rare odors between channels Fz and Cz for each 
        % patient
        plvRareNormal(k) = PLVCal(double(rareChannelFzNormal) , double(rareChannelCzNormal) , 200 , [35 , 40]);
    end
    % get average on the trials
    plvRareMeanNormal(i) = mean(plvRareNormal);
end

% draw the box plots
figure
subplot(2 , 2 , 1)
boxplot(plvFreqMeanNormal);
title("Frequent Normal");
subplot(2 , 2 , 2)
boxplot(plvRareMeanNormal);
title("Rare Normal");
subplot(2 , 2 , 3)
boxplot(plvFreqMeanAD);
title("Frequent AD");
subplot(2 , 2 , 4)
boxplot(plvRareMeanAD);
title("Rare AD");

% fit a Gaussian distribution on the PLVs
NormalFreqDist = fitdist(plvFreqMeanNormal ,'Normal');
NormalRareDist = fitdist(plvRareMeanNormal ,'Normal');
ADFreqDist = fitdist(plvFreqMeanAD ,'Normal');
ADRareDist = fitdist(plvRareMeanAD ,'Normal');

% plot the Gaussian PDFs
figure
NormalFreqx = 0.8504-sqrt(0.191814):0.001:0.8504+sqrt(0.191814);
y = pdf(NormalFreqDist,NormalFreqx);
subplot(2 , 2 , 1)
plot(NormalFreqx , y)
title("Frequent Normal");

NormalRarex = 0.85151-sqrt(0.195164):0.001:0.85151+sqrt(0.195164);
y = pdf(NormalRareDist,NormalRarex);
subplot(2 , 2 , 2)
plot(NormalRarex , y)
title("Rare Normal");

ADFreqx = 0.728897-sqrt(0.177983):0.001:0.728897+sqrt(0.177983);
y = pdf(ADFreqDist,ADFreqx);
subplot(2 , 2 , 3)
plot(ADFreqx , y)
title("Frequent AD");

ADRarex = 0.729726-sqrt(0.179864):0.001:0.729726+sqrt(0.179864);
y = pdf(ADRareDist,ADRarex);
subplot(2 , 2 , 4)
plot(ADRarex , y)
title("Rare AD");

% P_values for frequent odor
[freqH,freqP] = ttest2(plvFreqMeanNormal , plvFreqMeanAD);
% P_values for rare odor
[rareH,rareP] = ttest2(plvRareMeanNormal , plvRareMeanAD);

% "P_values for frequent odor = "+freqP
% "Hypothesis test result for frequent odor =  "+freqH
% "P_values for rare odor = "+rareP
% "Hypothesis test result for rare odor =  "+rareH

%% Phase diffrence of the i-th subject

% Phase diffrence of the i-th AD subject(i-th patient)
i=3;
% get the epoch field for the i-th subject
ADepoch1=getfield(AD(1 , i) , "epoch");
% get the odor field for the i-th subject
ADodor1=getfield(AD(1 , i) , "odor");
% find the trial indices at which odor is frequent
freqIdxAD1=find(ADodor1(: , 1)==0);
% create an array for Phase diffrence of each frequent odor trial 
plvFreqAD1=zeros(size(freqIdxAD1 , 1),1);
% create an array for phase difference of each frequent odor trial 
phaseDiffAD1=zeros(size(freqIdxAD1 , 1) , 600);
% find the frequent odor trials of channels Fz & Cz
for j=1:size(freqIdxAD1 , 1)
    freqChannelFzAD1 = ADepoch1(2 , : , freqIdxAD1(j , 1)); % Fz channel
    freqChannelCzAD1 = ADepoch1(3 , : , freqIdxAD1(j , 1)); % Cz channel
    phaseDiffAD1(j , :)=phaseDiffCal(freqChannelFzAD1 , freqChannelCzAD1 , 200 , [35 , 40]);
end
% get average on the trials
phaseDiffAD1Mean = mean(phaseDiffAD1 , 1);

figure
subplot(2 , 1 , 1)
polarhistogram(phaseDiffAD1Mean , 10);
title("AD for subject"+i)

% Phase diffrence of the i-th Normal subject(i-th patient)
% get the epoch field for the i-th subject
Normalepoch1=getfield(normal(1 , i) , "epoch");
% get the odor field for the i-th subject
Normalodor1=getfield(normal(1 , i) , "odor");
% find the trial indices at which odor is frequent
freqIdxNormal1=find(Normalodor1(: , 1)==0);
% create an array for Phase diffrence of each frequent odor trial 
plvFreqNormal1=zeros(size(freqIdxNormal1 , 1),1);
% create an array for phase difference of each frequent odor trial 
phaseDiffNormal1=zeros(size(freqIdxNormal1 , 1) , 600);
% find the frequent odor trials of channel Fz & Cz
for j=1:size(freqIdxNormal1 , 1)
    freqChannelFzNormal1 = Normalepoch1(2 , : , freqIdxNormal1(j , 1)); % Fz channel
    freqChannelCzNormal1 = Normalepoch1(3 , : , freqIdxNormal1(j , 1)); % Cz channel
    phaseDiffNormal1(j , :)=phaseDiffCal(freqChannelFzNormal1 , freqChannelCzNormal1 , 200 , [35 , 40]);
end
% get average on the trials
phaseDiffNormal1Mean = mean(phaseDiffNormal1 , 1);

subplot(2 , 1 , 2)
polarhistogram(phaseDiffNormal1Mean , 10);
title("Normal for subject"+i)

%% Mean value among all the subjects

% Mean value among all the subjects of AD group
phaseDiffFreqMeanAD=zeros(13 , 600);
% looping on the patients
for i=1:13
    % get the epoch field for each patient 
    ADepoch=getfield(AD(1 , i) , "epoch");
    % get the odor field for each patient 
    ADodor=getfield(AD(1 , i) , "odor");
    % find the trial indices at which odor is frequent
    freqIdxAD=find(ADodor(: , 1)==0);
    % create an array for phase difference of each frequent odor trial 
    phaseDiffFreqAD=zeros(size(freqIdxAD , 1) , 600);
    % find the frequent odor trials of channel Fz & Cz
    for j=1:size(freqIdxAD , 1)
        freqChannelFzAD = ADepoch(2 , : , freqIdxAD(j , 1)); % Fz channel
        freqChannelCzAD = ADepoch(3 , : , freqIdxAD(j , 1)); % Cz channel
        % calculate PD on frequents odors between channels Fz and Cz for each 
        % patient
        phaseDiffFreqAD(j , :) = phaseDiffCal(double(freqChannelFzAD) , double(freqChannelCzAD) , 200 , [35 , 40]);
    end
    % get average on the trials
    phaseDiffFreqMeanAD(i , :) = mean(phaseDiffFreqAD , 1);   
end
% get average on the subjects
phaseDiffFreqMeanADRes=mean(phaseDiffFreqMeanAD , 1);

% mean value among all the subjects of Normal group
phaseDiffFreqMeanNormal=zeros(15 , 600);
% looping on the patients
for i=1:15
    % get the epoch field for each patient 
    Normalepoch=getfield(normal(1 , i) , "epoch");
    % get the odor field for each patient 
    Normalodor=getfield(normal(1 , i) , "odor");
    % find the trial indices at which odor is frequent
    freqIdxNormal=find(Normalodor(: , 1)==0);
    % create an array for phase difference of each frequent odor trial 
    phaseDiffFreqNormal=zeros(size(freqIdxNormal , 1) , 600);
    % find the frequent odor trials of channels Fz & Cz
    for j=1:size(freqIdxNormal , 1)
        freqChannelFzNormal = Normalepoch(2 , : , freqIdxNormal(j , 1)); % Fz channel
        freqChannelCzNormal = Normalepoch(3 , : , freqIdxNormal(j , 1)); % Cz channel
        % calculate PD on frequents odors between channels Fz and Cz for each 
        % patient
        phaseDiffFreqNormal(j , :)=phaseDiffCal(double(freqChannelFzNormal) , double(freqChannelCzNormal) , 200 , [35 , 40]);
    end
    % get average on the trials
    phaseDiffFreqMeanNormal(i , :) = mean(phaseDiffFreqNormal , 1); 
end
% get average on the subjects
phaseDiffFreqMeanNormalRes=mean(phaseDiffFreqMeanNormal , 1);

figure
subplot(2 , 1 , 1)
polarhistogram(phaseDiffFreqMeanADRes,10);
title("AD Mean value")
subplot(2 , 1 , 2)
polarhistogram(phaseDiffFreqMeanNormalRes,10);
title("Normal mean value")

%% Heatmap PLVs
% AD
plv_AD_Mean_Freq_Fp1_Fz=zeros(13 , 1);
plv_AD_Mean_Rare_Fp1_Fz=zeros(13 , 1);

plv_AD_Mean_Freq_Fp1_Cz=zeros(13 , 1);
plv_AD_Mean_Rare_Fp1_Cz=zeros(13 , 1);

plv_AD_Mean_Freq_Fp1_Pz=zeros(13 , 1);
plv_AD_Mean_Rare_Fp1_Pz=zeros(13 , 1);

plv_AD_Mean_Freq_Fz_Cz=zeros(13 , 1);
plv_AD_Mean_Rare_Fz_Cz=zeros(13 , 1);

plv_AD_Mean_Freq_Fz_Pz=zeros(13 , 1);
plv_AD_Mean_Rare_Fz_Pz=zeros(13 , 1);

plv_AD_Mean_Freq_Cz_Pz=zeros(13 , 1);
plv_AD_Mean_Rare_Cz_Pz=zeros(13 , 1);

% looping on the patients
for i=1:13
    % get the epoch field for each patient 
    ADepoch=getfield(AD(1 , i) , "epoch");
    % get the odor field for each patient 
    ADodor=getfield(AD(1 , i) , "odor");
    % find the trial indices at which odor is frequent
    freqIdxAD=find(ADodor(: , 1)==0);
    % find the trial indices at which odor is rare
    rareIdxAD=find(ADodor(: , 1)==1);
    % create an array for PLV of each frequent odor trial 
    plv_AD_Freq_Fp1_Fz=zeros(size(freqIdxAD , 1),1);
    plv_AD_Freq_Fp1_Cz=zeros(size(freqIdxAD , 1),1);
    plv_AD_Freq_Fp1_Pz=zeros(size(freqIdxAD , 1),1);
    plv_AD_Freq_Fz_Cz=zeros(size(freqIdxAD , 1),1);
    plv_AD_Freq_Fz_Pz=zeros(size(freqIdxAD , 1),1);
    plv_AD_Freq_Cz_Pz=zeros(size(freqIdxAD , 1),1);
    
    % create an array for PLV of each rare odor trial 
    plv_AD_Rare_Fp1_Fz=zeros(size(rareIdxAD , 1),1);
    plv_AD_Rare_Fp1_Cz=zeros(size(rareIdxAD , 1),1);
    plv_AD_Rare_Fp1_Pz=zeros(size(rareIdxAD , 1),1);
    plv_AD_Rare_Fz_Cz=zeros(size(rareIdxAD , 1),1);
    plv_AD_Rare_Fz_Pz=zeros(size(rareIdxAD , 1),1);
    plv_AD_Rare_Cz_Pz=zeros(size(rareIdxAD , 1),1);

    % find the frequent odor trials of channels Fz & Cz
    for j=1:size(freqIdxAD , 1)
        freqChannelAD_Fp1 = ADepoch(1 , : , freqIdxAD(j , 1)); 
        freqChannelAD_Fz = ADepoch(2 , : , freqIdxAD(j , 1)); 
        freqChannelAD_Cz = ADepoch(3 , : , freqIdxAD(j , 1));
        freqChannelAD_Pz = ADepoch(4 , : , freqIdxAD(j , 1));
        % calculate PLV 
        plv_AD_Freq_Fp1_Fz(j)=PLVCal(double(freqChannelAD_Fp1) , double(freqChannelAD_Fz) , 200 , [35 , 40]);
        plv_AD_Freq_Fp1_Cz(j)=PLVCal(double(freqChannelAD_Fp1) , double(freqChannelAD_Cz) , 200 , [35 , 40]);
        plv_AD_Freq_Fp1_Pz(j)=PLVCal(double(freqChannelAD_Fp1) , double(freqChannelAD_Pz) , 200 , [35 , 40]);
        plv_AD_Freq_Fz_Cz(j)=PLVCal(double(freqChannelAD_Fz) , double(freqChannelAD_Cz) , 200 , [35 , 40]);
        plv_AD_Freq_Fz_Pz(j)=PLVCal(double(freqChannelAD_Fz) , double(freqChannelAD_Pz) , 200 , [35 , 40]);
        plv_AD_Freq_Cz_Pz(j)=PLVCal(double(freqChannelAD_Cz) , double(freqChannelAD_Pz) , 200 , [35 , 40]);
    end
    % get average on the trials
    plv_AD_Mean_Freq_Fp1_Fz(i)=mean(plv_AD_Freq_Fp1_Fz);    
    plv_AD_Mean_Freq_Fp1_Cz(i)=mean(plv_AD_Freq_Fp1_Cz);
    plv_AD_Mean_Freq_Fp1_Pz(i)=mean(plv_AD_Freq_Fp1_Pz);
    plv_AD_Mean_Freq_Fz_Cz(i)=mean(plv_AD_Freq_Fz_Cz);
    plv_AD_Mean_Freq_Fz_Pz(i)=mean(plv_AD_Freq_Fz_Pz);
    plv_AD_Mean_Freq_Cz_Pz(i)=mean(plv_AD_Freq_Cz_Pz);
    
   for j=1:size(rareIdxAD , 1)
        rareChannelAD_Fp1 = ADepoch(1 , : , rareIdxAD(j , 1)); 
        rareChannelAD_Fz = ADepoch(2 , : , rareIdxAD(j , 1)); 
        rareChannelAD_Cz = ADepoch(3 , : , rareIdxAD(j , 1));
        rareChannelAD_Pz = ADepoch(4 , : , rareIdxAD(j , 1));
        % calculate PLV 
        plv_AD_Rare_Fp1_Fz(j)=PLVCal(double(rareChannelAD_Fp1) , double(rareChannelAD_Fz) , 200 , [35 , 40]);
        plv_AD_Rare_Fp1_Cz(j)=PLVCal(double(rareChannelAD_Fp1) , double(rareChannelAD_Cz) , 200 , [35 , 40]);
        plv_AD_Rare_Fp1_Pz(j)=PLVCal(double(rareChannelAD_Fp1) , double(rareChannelAD_Pz) , 200 , [35 , 40]);
        plv_AD_Rare_Fz_Cz(j)=PLVCal(double(rareChannelAD_Fz) , double(rareChannelAD_Cz) , 200 , [35 , 40]);
        plv_AD_Rare_Fz_Pz(j)=PLVCal(double(rareChannelAD_Fz) , double(rareChannelAD_Pz) , 200 , [35 , 40]);
        plv_AD_Rare_Cz_Pz(j)=PLVCal(double(rareChannelAD_Cz) , double(rareChannelAD_Pz) , 200 , [35 , 40]);
    end
    % get average on the trials
    plv_AD_Mean_Rare_Fp1_Fz(i)=mean(plv_AD_Rare_Fp1_Fz);    
    plv_AD_Mean_Rare_Fp1_Cz(i)=mean(plv_AD_Rare_Fp1_Cz);
    plv_AD_Mean_Rare_Fp1_Pz(i)=mean(plv_AD_Rare_Fp1_Pz);
    plv_AD_Mean_Rare_Fz_Cz(i)=mean(plv_AD_Rare_Fz_Cz);
    plv_AD_Mean_Rare_Fz_Pz(i)=mean(plv_AD_Rare_Fz_Pz);
    plv_AD_Mean_Rare_Cz_Pz(i)=mean(plv_AD_Rare_Cz_Pz);
    
end


% Normal
plv_Normal_Mean_Freq_Fp1_Fz=zeros(15 , 1);
plv_Normal_Mean_Rare_Fp1_Fz=zeros(15 , 1);

plv_Normal_Mean_Freq_Fp1_Cz=zeros(15 , 1);
plv_Normal_Mean_Rare_Fp1_Cz=zeros(15 , 1);

plv_Normal_Mean_Freq_Fp1_Pz=zeros(15 , 1);
plv_Normal_Mean_Rare_Fp1_Pz=zeros(15 , 1);

plv_Normal_Mean_Freq_Fz_Cz=zeros(15 , 1);
plv_Normal_Mean_Rare_Fz_Cz=zeros(15 , 1);

plv_Normal_Mean_Freq_Fz_Pz=zeros(15 , 1);
plv_Normal_Mean_Rare_Fz_Pz=zeros(15 , 1);

plv_Normal_Mean_Freq_Cz_Pz=zeros(15 , 1);
plv_Normal_Mean_Rare_Cz_Pz=zeros(15 , 1);

% looping on the patients
for i=1:15
    % get the epoch field for each patient 
    Normalepoch=getfield(normal(1 , i) , "epoch");
    % get the odor field for each patient 
    Normalodor=getfield(normal(1 , i) , "odor");
    % find the trial indices at which odor is frequent
    freqIdxNormal=find(Normalodor(: , 1)==0);
    % find the trial indices at which odor is rare
    rareIdxNormal=find(Normalodor(: , 1)==1);
    % create an array for PLV of each frequent odor trial 
    plv_Normal_Freq_Fp1_Fz=zeros(size(freqIdxNormal , 1),1);
    plv_Normal_Freq_Fp1_Cz=zeros(size(freqIdxNormal , 1),1);
    plv_Normal_Freq_Fp1_Pz=zeros(size(freqIdxNormal , 1),1);
    plv_Normal_Freq_Fz_Cz=zeros(size(freqIdxNormal , 1),1);
    plv_Normal_Freq_Fz_Pz=zeros(size(freqIdxNormal , 1),1);
    plv_Normal_Freq_Cz_Pz=zeros(size(freqIdxNormal , 1),1);

    % create an array for PLV of each rare odor trial 
    plv_Normal_Rare_Fp1_Fz=zeros(size(rareIdxNormal , 1),1);
    plv_Normal_Rare_Fp1_Cz=zeros(size(rareIdxNormal , 1),1);
    plv_Normal_Rare_Fp1_Pz=zeros(size(rareIdxNormal , 1),1);
    plv_Normal_Rare_Fz_Cz=zeros(size(rareIdxNormal , 1),1);
    plv_Normal_Rare_Fz_Pz=zeros(size(rareIdxNormal , 1),1);
    plv_Normal_Rare_Cz_Pz=zeros(size(rareIdxNormal , 1),1);

    % find the frequent odor trials of channels Fz & Cz
    for j=1:size(freqIdxNormal , 1)
        freqChannelNormal_Fp1 = Normalepoch(1 , : , freqIdxNormal(j , 1)); 
        freqChannelNormal_Fz = Normalepoch(2 , : , freqIdxNormal(j , 1)); 
        freqChannelNormal_Cz = Normalepoch(3 , : , freqIdxNormal(j , 1));
        freqChannelNormal_Pz = Normalepoch(4 , : , freqIdxNormal(j , 1));
        % calculate PLV 
        plv_Normal_Freq_Fp1_Fz(j)=PLVCal(double(freqChannelNormal_Fp1) , double(freqChannelNormal_Fz) , 200 , [35 , 40]);
        plv_Normal_Freq_Fp1_Cz(j)=PLVCal(double(freqChannelNormal_Fp1) , double(freqChannelNormal_Cz) , 200 , [35 , 40]);
        plv_Normal_Freq_Fp1_Pz(j)=PLVCal(double(freqChannelNormal_Fp1) , double(freqChannelNormal_Pz) , 200 , [35 , 40]);
        plv_Normal_Freq_Fz_Cz(j)=PLVCal(double(freqChannelNormal_Fz) , double(freqChannelNormal_Cz) , 200 , [35 , 40]);
        plv_Normal_Freq_Fz_Pz(j)=PLVCal(double(freqChannelNormal_Fz) , double(freqChannelNormal_Pz) , 200 , [35 , 40]);
        plv_Normal_Freq_Cz_Pz(j)=PLVCal(double(freqChannelNormal_Cz) , double(freqChannelNormal_Pz) , 200 , [35 , 40]);
    end
    % get average on the trials
    plv_Normal_Mean_Freq_Fp1_Fz(i)=mean(plv_Normal_Freq_Fp1_Fz);    
    plv_Normal_Mean_Freq_Fp1_Cz(i)=mean(plv_Normal_Freq_Fp1_Cz);
    plv_Normal_Mean_Freq_Fp1_Pz(i)=mean(plv_Normal_Freq_Fp1_Pz);
    plv_Normal_Mean_Freq_Fz_Cz(i)=mean(plv_Normal_Freq_Fz_Cz);
    plv_Normal_Mean_Freq_Fz_Pz(i)=mean(plv_Normal_Freq_Fz_Pz);
    plv_Normal_Mean_Freq_Cz_Pz(i)=mean(plv_Normal_Freq_Cz_Pz);
    
   for j=1:size(rareIdxNormal , 1)
        rareChannelNormal_Fp1 = Normalepoch(1 , : , rareIdxNormal(j , 1)); 
        rareChannelNormal_Fz = Normalepoch(2 , : , rareIdxNormal(j , 1)); 
        rareChannelNormal_Cz = Normalepoch(3 , : , rareIdxNormal(j , 1));
        rareChannelNormal_Pz = Normalepoch(4 , : , rareIdxNormal(j , 1));
        % calculate PLV 
        plv_Normal_Rare_Fp1_Fz(j)=PLVCal(double(rareChannelNormal_Fp1) , double(rareChannelNormal_Fz) , 200 , [35 , 40]);
        plv_Normal_Rare_Fp1_Cz(j)=PLVCal(double(rareChannelNormal_Fp1) , double(rareChannelNormal_Cz) , 200 , [35 , 40]);
        plv_Normal_Rare_Fp1_Pz(j)=PLVCal(double(rareChannelNormal_Fp1) , double(rareChannelNormal_Pz) , 200 , [35 , 40]);
        plv_Normal_Rare_Fz_Cz(j)=PLVCal(double(rareChannelNormal_Fz) , double(rareChannelNormal_Cz) , 200 , [35 , 40]);
        plv_Normal_Rare_Fz_Pz(j)=PLVCal(double(rareChannelNormal_Fz) , double(rareChannelNormal_Pz) , 200 , [35 , 40]);
        plv_Normal_Rare_Cz_Pz(j)=PLVCal(double(rareChannelNormal_Cz) , double(rareChannelNormal_Pz) , 200 , [35 , 40]);
    end
    % get average on the trials
    plv_Normal_Mean_Rare_Fp1_Fz(i)=mean(plv_Normal_Rare_Fp1_Fz);    
    plv_Normal_Mean_Rare_Fp1_Cz(i)=mean(plv_Normal_Rare_Fp1_Cz);
    plv_Normal_Mean_Rare_Fp1_Pz(i)=mean(plv_Normal_Rare_Fp1_Pz);
    plv_Normal_Mean_Rare_Fz_Cz(i)=mean(plv_Normal_Rare_Fz_Cz);
    plv_Normal_Mean_Rare_Fz_Pz(i)=mean(plv_Normal_Rare_Fz_Pz);
    plv_Normal_Mean_Rare_Cz_Pz(i)=mean(plv_Normal_Rare_Cz_Pz);
    
end

%% Heatmap plot

figure
subplot(2 , 2 , 1);
data_ad_freq=[mean(plv_AD_Mean_Freq_Fp1_Pz,1) , mean(plv_AD_Mean_Freq_Fz_Pz,1) , mean(plv_AD_Mean_Freq_Cz_Pz,1) ,...
    1 ; mean(plv_AD_Mean_Freq_Fp1_Cz,1) , mean(plv_AD_Mean_Freq_Fz_Cz,1) , 1 , mean(plv_AD_Mean_Freq_Cz_Pz,1) ;...
    mean(plv_AD_Mean_Freq_Fp1_Fz,1) , 1 , mean(plv_AD_Mean_Freq_Fz_Cz,1) , mean(plv_AD_Mean_Freq_Fz_Pz,1) ; ...
    1 , mean(plv_AD_Mean_Freq_Fp1_Fz,1) , mean(plv_AD_Mean_Freq_Fp1_Cz,1) , mean(plv_AD_Mean_Freq_Fp1_Pz,1)];
heatmapPlot(data_ad_freq);
title("AD Freq")

subplot(2 , 2 , 2);
data_ad_rare=[mean(plv_AD_Mean_Rare_Fp1_Pz,1) , mean(plv_AD_Mean_Rare_Fz_Pz,1) , mean(plv_AD_Mean_Rare_Cz_Pz,1) ,...
    1 ; mean(plv_AD_Mean_Rare_Fp1_Cz,1) , mean(plv_AD_Mean_Rare_Fz_Cz,1) , 1 , mean(plv_AD_Mean_Rare_Cz_Pz,1) ;...
    mean(plv_AD_Mean_Rare_Fp1_Fz,1) , 1 , mean(plv_AD_Mean_Rare_Fz_Cz,1) , mean(plv_AD_Mean_Rare_Fz_Pz,1) ; ...
    1 , mean(plv_AD_Mean_Rare_Fp1_Fz,1) , mean(plv_AD_Mean_Rare_Fp1_Cz,1) , mean(plv_AD_Mean_Rare_Fp1_Pz,1)];
heatmapPlot(data_ad_rare);
title("AD Rare")

subplot(2 , 2 , 3);
data_norm_freq=[mean(plv_Normal_Mean_Freq_Fp1_Pz,1) , mean(plv_Normal_Mean_Freq_Fz_Pz,1) , mean(plv_Normal_Mean_Freq_Cz_Pz,1) ,...
    1 ; mean(plv_Normal_Mean_Freq_Fp1_Cz,1) , mean(plv_Normal_Mean_Freq_Fz_Cz,1) , 1 , mean(plv_Normal_Mean_Freq_Cz_Pz,1) ;...
    mean(plv_Normal_Mean_Freq_Fp1_Fz,1) , 1 , mean(plv_Normal_Mean_Freq_Fz_Cz,1) , mean(plv_Normal_Mean_Freq_Fz_Pz,1) ; ...
    1 , mean(plv_Normal_Mean_Freq_Fp1_Fz,1) , mean(plv_Normal_Mean_Freq_Fp1_Cz,1) , mean(plv_Normal_Mean_Freq_Fp1_Pz,1)];
heatmapPlot(data_norm_freq);
title("Normal Freq")

subplot(2 , 2 , 4);
data_norm_rare=[mean(plv_Normal_Mean_Rare_Fp1_Pz,1) , mean(plv_Normal_Mean_Rare_Fz_Pz,1) , mean(plv_Normal_Mean_Rare_Cz_Pz,1) ,...
    1 ; mean(plv_Normal_Mean_Rare_Fp1_Cz,1) , mean(plv_Normal_Mean_Rare_Fz_Cz,1) , 1 , mean(plv_Normal_Mean_Rare_Cz_Pz,1) ;...
    mean(plv_Normal_Mean_Rare_Fp1_Fz,1) , 1 , mean(plv_Normal_Mean_Rare_Fz_Cz,1) , mean(plv_Normal_Mean_Rare_Fz_Pz,1) ; ...
    1 , mean(plv_Normal_Mean_Rare_Fp1_Fz,1) , mean(plv_Normal_Mean_Rare_Fp1_Cz,1) , mean(plv_Normal_Mean_Rare_Fp1_Pz,1)];
heatmapPlot(data_norm_rare);
title("Normal Rare")

%% P_vals

% P_values for frequent odors
[freqH_Fp1_Fz,freqP_Fp1_Fz] = ttest2(plv_Normal_Mean_Freq_Fp1_Fz , plv_AD_Mean_Freq_Fp1_Fz);
[freqH_Fp1_Cz,freqP_Fp1_Cz] = ttest2(plv_Normal_Mean_Freq_Fp1_Cz , plv_AD_Mean_Freq_Fp1_Cz);
[freqH_Fp1_Pz,freqP_Fp1_Pz] = ttest2(plv_Normal_Mean_Freq_Fp1_Pz , plv_AD_Mean_Freq_Fp1_Pz); 
[freqH_Fz_Cz,freqP_Fz_Cz] = ttest2(plv_Normal_Mean_Freq_Fz_Cz , plv_AD_Mean_Freq_Fz_Cz); 
[freqH_Fz_Pz,freqP_Fz_Pz] = ttest2(plv_Normal_Mean_Freq_Fz_Pz , plv_AD_Mean_Freq_Fz_Pz);
[freqH_Cz_Pz,freqP_Cz_Pz] = ttest2(plv_Normal_Mean_Freq_Cz_Pz , plv_AD_Mean_Freq_Cz_Pz); 

% P_values for rare odors
[rareH_Fp1_Fz,rareP_Fp1_Fz] = ttest2(plv_Normal_Mean_Rare_Fp1_Fz , plv_AD_Mean_Rare_Fp1_Fz);
[rareH_Fp1_Cz,rareP_Fp1_Cz] = ttest2(plv_Normal_Mean_Rare_Fp1_Cz , plv_AD_Mean_Rare_Fp1_Cz);
[rareH_Fp1_Pz,rareP_Fp1_Pz] = ttest2(plv_Normal_Mean_Rare_Fp1_Pz , plv_AD_Mean_Rare_Fp1_Pz); 
[rareH_Fz_Cz,rareP_Fz_Cz] = ttest2(plv_Normal_Mean_Rare_Fz_Cz , plv_AD_Mean_Rare_Fz_Cz); 
[rareH_Fz_Pz,rareP_Fz_Pz] = ttest2(plv_Normal_Mean_Rare_Fz_Pz , plv_AD_Mean_Rare_Fz_Pz);
[rareH_Cz_Pz,rareP_Cz_Pz] = ttest2(plv_Normal_Mean_Rare_Cz_Pz , plv_AD_Mean_Rare_Cz_Pz); 

%% Bonus_PLVs

% PLV for MCI patients
plvFreqMeanMCI=zeros(7 , 1);
plvRareMeanMCI=zeros(7 , 1);
% looping on the patients
for i=1:7
    % get the epoch field for each patient 
    MCIepoch=getfield(MCI(1 , i) , "epoch");
    % get the odor field for each patient 
    MCIodor=getfield(MCI(1 , i) , "odor");
    % find the trial indices at which odor is frequent
    freqIdxMCI=find(MCIodor(: , 1)==0);
    % find the trial indices at which odor is rare
    rareIdxMCI=find(MCIodor(: , 1)==1);
    % create an array for PLV of each frequent odor trial 
    plvFreqMCI=zeros(size(freqIdxMCI , 1),1);
    % create an array for PLV of each rare odor trial 
    plvRareMCI=zeros(size(rareIdxMCI , 1),1);
    % find the frequent odor trials of channels Fz & Cz
    for j=1:size(freqIdxMCI , 1)
        freqChannelFzMCI = MCIepoch(2 , : , freqIdxMCI(j , 1)); % Fz channel
        freqChannelCzMCI = MCIepoch(3 , : , freqIdxMCI(j , 1)); % Cz channel
        % calculate PLV on frequents odors between channels Fz and Cz for each 
        % patient
        plvFreqMCI(j) = PLVCal(double(freqChannelFzMCI) , double(freqChannelCzMCI) , 200 , [35 , 40]);
    end
    % get average on the trials
    plvFreqMeanMCI(i) = mean(plvFreqMCI);   

    % find the rare odor trials of channels Fz & Cz
    for k=1:size(rareIdxMCI , 1)
        rareChannelFzMCI = ADepoch(2 , : , rareIdxMCI(k , 1)); % Fz channel
        rareChannelCzMCI = ADepoch(3 , : , rareIdxMCI(k , 1)); % Cz channel
        % calculate PLV on rare odors between channels Fz and Cz for each 
        % patient
        plvRareMCI(k) = PLVCal(double(rareChannelFzMCI) , double(rareChannelCzMCI) , 200 , [35 , 40]);
    end
    % get average on the trials
    plvRareMeanMCI(i) = mean(plvRareMCI);
end

figure
subplot(1 , 2 , 1)
boxplot(plvFreqMeanMCI);
title("Frequent MCI");
subplot(1 , 2 , 2)
boxplot(plvRareMeanMCI);
title("Rare MCI");

%% Bonus_P vals
% between normal & AD
% P_values for frequent odor
[freqH1,freqP1] = ttest2(plvFreqMeanNormal , plvFreqMeanAD);
% P_values for rare odor
[rareH1,rareP1] = ttest2(plvRareMeanNormal , plvRareMeanAD);

% between normal & MCI
% P_values for frequent odor
[freqH2,freqP2] = ttest2(plvFreqMeanNormal , plvFreqMeanMCI);
% P_values for rare odor
[rareH2,rareP2] = ttest2(plvRareMeanNormal , plvRareMeanMCI);

% between AD & MCI
% P_values for frequent odor
[freqH3,freqP3] = ttest2(plvFreqMeanMCI , plvFreqMeanAD);
% P_values for rare odor
[rareH3,rareP3] = ttest2(plvRareMeanMCI , plvRareMeanAD);

%% Functions
% frequency spectrum
function [fftRes,y1] = fftCal(data , fs)
    L = length(data);
    y = fft(data); 
    y2 = abs(y/L);
    y1 = y2(1:L/2+1);
    y1(2:end-1) = 2*y1(2:end-1);
    fftRes = (0:(L/2))*fs/L;
end
% PLV calculator
function plv = PLVCal(data1 , data2 , fs , frequencyRange)
        % calculate plv for each frequent odor trial(looping on the trials)
        % apply bandpass filter to extract the desired frequency range
        lowFreq = frequencyRange(1)/(fs/2);
        highFreq = frequencyRange(2)/(fs/2);
        [b , a] = butter(4 , [lowFreq , highFreq] , 'bandpass');
        filteredData1 = filtfilt(b , a , data1);
        filteredData2 = filtfilt(b , a , data2); 
        % apply Hilbert transform to extract the phases
        analyticSignal1 = hilbert(filteredData1);
        analyticSignal2 = hilbert(filteredData2);
        % calculate the phases
        phase1 = angle(analyticSignal1);
        phase2 = angle(analyticSignal2);
        % the phase differences
        phaseDiff = abs(phase1-phase2);
        % calculate PLV due to the formula
        plv = abs(mean(exp(1i*phaseDiff)));
end
% Phase Diff calculator
function phaseDiff = phaseDiffCal(data1 , data2 , fs , frequencyRange)
        % calculate plv for each frequent odor trial(looping on the trials)
        % apply bandpass filter to extract the desired frequency range
        lowFreq = frequencyRange(1)/(fs/2);
        highFreq = frequencyRange(2)/(fs/2);
        [b , a] = butter(1 , [lowFreq , highFreq] , 'bandpass');
        filteredData1 = filtfilt(b , a , data1);
        filteredData2 = filtfilt(b , a , data2); 
        % apply Hilbert transform to extract the phases
        analyticSignal1 = hilbert(filteredData1);
        analyticSignal2 = hilbert(filteredData2);
        % calculate the phases in deg
        phase1 = rad2deg(angle(analyticSignal1));
        phase2 = rad2deg(angle(analyticSignal2));
        % the phase differences
        phaseDiff = abs(phase1-phase2);
end
% Heatmap plotter
function heatmapPlot(data)
    imagesc(data , [0, 1]);
    xticklabels({'Fp1', 'Fz', 'Cz', 'Pz'});
    yticklabels({'Pz', 'Cz', 'Fz', 'Fp1'});
    [X , Y] = meshgrid(1:size(data, 2), 1:size(data, 1));
    text(X(:), Y(:), num2str(data(:) , '%.2f') , 'HorizontalAlignment', 'center');
end



















