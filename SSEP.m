
%% Calculating the normalized FFT

clear, clc
sub = [82:92,98:124,127:134];
FFT = cell(length(sub),1);

for s = 1:length(sub)

address = ['\N',num2str(sub(s)),'\PreFFT18SecNewClean.mat']; % Load preprocessed data of each subject
load(address);
sig = EEG.data;

    tmp = cell(length(EEG.chanlocs),1);
    for i=1:length(EEG.chanlocs)
        tmp{i} = EEG.chanlocs(i).labels;
    end
    
    if length(tmp) == 98
        load('Chan128.mat') % channel locations
    else
        
        load('Chan64.mat')
    end
    
    chanlocs = cell(1,size(chan,2));
    for i = 1:length(chanlocs)
        chanlocs{i} = chan(i).labels;
    end 

    temp = ismember(chanlocs,tmp);
    tmp1 = find(temp==0);
    tmp2 = find(temp==1);
    NewSig(tmp2,:,:) = sig;
    NewSig(tmp1,:,:) = NaN;
    sig = NewSig;    
    clear NewSig   
    
    sig             = nanmean(sig,3);
    sig             = sig - mean(sig,2);
    freq            = abs(fft(sig')/size(sig,2))';
    freq            = freq(:,1:size(sig,2)/2+1);
    freq(:,2:end-1) = 2*freq(:,2:end-1);
    f               = EEG.srate*(0:(size(sig,2)/2))/size(sig,2);       
    
    temp = freq;
    tmp  = find(f>=.5 & f<=15);
    tmp  = [tmp(1),tmp(end)];
    for i=tmp(1):tmp(2)
            tmp1 = [temp(:,i-6:i-3),temp(:,i+3:i+6)];
            ind = [i-6:i-3,i+3:i+6];
            ind = find(ind==dsearchn(f',1.11) | ind==dsearchn(f',1.667) | ind==dsearchn(f',3.33));
            tmp1(:,ind) = [];
            freq(:,i) = temp(:,i)./mean(tmp1,2); 
    end
        
    FFT{s,1} = freq;
end

%% down-sample 128-electrodes to a 64-electrodes

channelson128 = [1 3 4 6 9 11 13 16 19 22 23 24 27 28 29 30 32 33 34 36 37 41 44 45 46 47 51 52 57 58 60 ...
64 67 62 70 72 75 77 83 85 87 92 95 96 97 98 100 102 103 104 105 108 111 112 114 116 117 122 123 124 125 126 127 128];
for i = 1:length(FFT)
if size(FFT{i,1},1) == 128
FFT{i,1} = FFT{i,1}(channelson128,:);
end
end

%% fieldtrip format

f = linspace(0,256,size(FFT{1,1},2)); 
for k = 1:length(FFT)
    load('Chan64.mat')
    EEG1 = pop_importdata('data',FFT{k,1},'srate',512);
    EEG1.chanlocs = chan;
    fielddata1 = eeglab2fieldtrip(EEG1,'timelockanalysis','coord_transform');
    fielddata1.time = f;                      
    temp = fielddata1.elec.pnt(:,1);
    fielddata1.elec.pnt(:,1) = -fielddata1.elec.pnt(:,2);
    fielddata1.elec.pnt(:,2) = temp;
    all_subs{k,1} = fielddata1;
end

%% To plot the spectrum

WGA = [32 31 33 36 36 35 32 35 32 32 33 31 31 31 32 33 34 34 33 33 33 29 30 30 30 33 30 30 29 30 35 35 35 35 34 33 33 30 33 28 28 30 30 31 29 29]';
age = 7*WGA+[6 4 3 4 1 4 3 2 4 6 2 5 6 5 3 0 1 2 5 5 0 5 0 0 2 0 4 6 2 6 4 4 0 1 2 5 6 1 5 5 7 3 3 3 5 4]';
[sortAge,sortAgeIdx] = sort(age);

FFT = FFT(sortAgeIdx,1);
all_subs = all_subs(sortAgeIdx,1);

% Pick the disierd subjects
FFT = FFT(24:46); % 1:23 for the Younger subjects; 24:46 for the Older subjects
all_subs = all_subs(24:46);

signal = zeros(size(FFT{1,1},1),size(FFT{1,1},2),size(FFT,1));
for s = 1:size(FFT,1)
    signal(:,:,s) = FFT{s,1};
end
temp1 = signal;
signal = nanmean(signal,3);
tmp = temp1;
num = zeros(size(tmp,1),1);
for q = 1:size(tmp,1)
    num(q) = sum(~isnan(sum(tmp(q,:,:),2)),3); 
end
sigstd = nanstd(temp1(:,:,:),[],3)./repmat(sqrt(num),1,size(signal,2));

f1 = dsearchn(f',.5);
f2 = dsearchn(f',15);
figure; plot(f(:,f1:f2), nanmean(signal(:,f1:f2)),'color',[0 0 0],'linewidth',1)



