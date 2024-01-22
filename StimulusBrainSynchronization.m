
clear
close all

fs        = 512;
sub = [82:92,98:124,127:134];

%%
signal = cell(length(sub),1);

for s = 1:length(sub)
    address = ['\N',num2str(sub(s)),'\PreFFT18Sec_.5-15Hz.mat']; % Load preprocessed data of each subject
    load(address);
    signal = EEG.data;

    
    tmp = cell(length(EEG.chanlocs),1);
    for i=1:length(EEG.chanlocs)
        tmp{i} = EEG.chanlocs(i).labels;
    end
    
    if EEG.nbchan ==98
    load('Chan128.mat') % channel locations
    chanlocs = cell(1,128);
    else
    load('Chan64.mat')
    chanlocs = cell(1,64);
    end
    for i=1:length(chanlocs)
        chanlocs{i} = chan(i).labels;
    end
     
    temp = ismember(chanlocs,tmp);
    tmp1 = find(temp==0);
    tmp2 = find(temp==1);
    NewSig(tmp2,:,:) = signal;
    NewSig(tmp1,:,:) = NaN;
    signal = NewSig;    
    clear NewSig
    
    % down-sample 128-electrodes to a 64-electrodes
    if size(signal,1) == 128
        channelson128 = [1 3 4 6 9 11 13 16 19 22 23 24 27 28 29 30 32 33 34 36 37 41 44 45 46 47 51 52 57 58 60 ...
    64 67 62 70 72 75 77 83 85 87 92 95 96 97 98 100 102 103 104 105 108 111 112 114 116 117 122 123 124 125 126 127 128];
    signal = signal(channelson128,:,:);
    end
    
    signal{s} = signal;    
end

for s = 1:length(sub)    
    temp = sum(squeeze(nansum(signal{s},2)));
    signal{s}(:,:,temp==0) = [];        
end

%% narrow-band filter EEG trials around beat and meter frequencies

freq = cell(length(sub),1);

[b, a, fCenterList, nTap, fEdges] = filterBank_prop((1/.3), (1/.3), .1504, 2, 512, 'fir1', .25);  % 1/.3 for Beat, 1/.6 for Duple and 1/.9 for triple
b = b{1};
a = a{1};

for k = 1:length(sub)     
    for j = 1:size(signal{k},3)
        sig = signal{k}(:,:,j);
        pow = sig;
        ind = isnan(sum(pow,2));
        pow(ind,:) = [];
        pow = pow';
        pow = filtfilt(b,a,pow);
        pow = pow';
        Pow(find(~ind),:) = pow;
        Pow(find(ind),:) = NaN;
        pow = Pow;
        clear Pow
        pow = pow(:,round(1.8*fs)+1:end-round(1.8*fs));
        freq{k,1}(:,:,j) = pow;
    end
end

for k=1:length(freq)
    kk = 0;
    index = [];
    for i=1:size(freq{k},3)
        if sum(sum(freq{k}(:,:,i))) == 0 
            kk = kk + 1;
            index(kk) = i;
        end
    end
    freq{k}(:,:,index) = [];
end

%% phase calculation

t = linspace(0,36/2,18432/2);
stimulus = .6*sin(2*pi*(1/.3)*t+pi/2); % 1/.3 for Beat, 1/.6 for Duple and 1/.9 for triple
stimulus = angle(hilbert(stimulus));
stimulus = stimulus(:,round(1.8*fs)+1:end-round(1.8*fs));

ang = cell(length(sub),1);
amp = cell(length(sub),1);
for k=1:length(sub)
    for j=1:size(freq{k,1},3)
        temp = hilbert(freq{k}(:,:,j)');
        temp = angle(temp)';
        
        for jj=1:size(temp,1)
            temp(jj,:) = circ_dist(temp(jj,:),stimulus); 
        end
        temp = mean(exp(1i*temp),2);
        
        ang{k}(:,j) = angle(temp);
        amp{k}(:,j) = abs(temp);
    end
end

% avg over blocks
clear tmp
for i=1:size(ang,1)
    tmp(:,i) = mean(exp(1i*ang{i}),2); 
end
