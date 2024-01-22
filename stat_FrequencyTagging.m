
%% Cluster-based permutation analysis 

clear
load('power.mat')
load('chan.mat')
load('neighbours.mat')

%%% The second dimention of 'power' matrix includes 0.56,1.11,1.67,2.22,2.78,3.33 and averaged of un-related frequencies

sub = 24:46; % Younger: 1:23; Older: 24:46
ff = 6;      % 2: Triple; 3: Duple; 6: Beat
signal1 = cell(length(sub),1);
signal2 = cell(length(sub),1);

for i = 1:length(sub)
    s = sub(i);
    signal1{i,1} = power(:,ff,s); % desired friquency 
    signal2{i,1} = power(:,7,s);  % average of power of unrelated friquencies
end

for i=1:length(sub)
    EEG = pop_importdata('data',signal1{i},'srate',1);
    EEG.chanlocs = chan;
    sig1{i} = eeglab2fieldtrip(EEG,'timelockanalysis','coord_transform');
    sig1{i}.dimord = 'chan_time';
    sig1{i}.time = 1;

    EEG = pop_importdata('data',signal2{i},'srate',1);
    EEG.chanlocs = chan;
    sig2{i} = eeglab2fieldtrip(EEG,'timelockanalysis','coord_transform');
    sig2{i}.dimord = 'chan_time';
    sig2{i}.time = 1;
end

cfg = [];
cfg.channel     = 'all';
cfg.neighbours  = neighbours; 
cfg.latency     = 1;
cfg.avgovertime = 'no';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.clusteralpha= 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 5000;
cfg.minnbchan        = 3; % minimal neighbouring channels 
Nsub = length(sub); 
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
 
stat = ft_timelockstatistics(cfg,sig1{:},sig2{:});
g = find(logical(stat.mask) == 1);

%% t-test with Bonferroni correction

cd 'F:\entrainment\entrainment_statistical'
temp = nanmean(power,1);
group = squeeze(temp(:,2:7,sub));

[h1,p1,sigPairs1] = ttest_bonf(group',[1 6;2 6;3 6;4 6;5 6]);

%% correlation of power and age 

load('power_cluster.mat')
WGA = [32 31 33 36 36 35 32 35 32 32 33 31 31 31 32 33 34 34 33 33 33 29 30 30 30 33 30 30 29 30 35 35 35 35 34 33 33 30 33 28 28 30 30 31 29 29]';
age = 7*WGA+[6 4 3 4 1 4 3 2 4 6 2 5 6 5 3 0 1 2 5 5 0 5 0 0 2 0 4 6 2 6 4 4 0 1 2 5 6 1 5 5 7 3 3 3 5 4]';
[sortAge] = sort(age);

ff = 3; % f = 1 for Beat, 2 for Duple, 3 for Triple
param = power_cluster(:,ff); 
[c, pval] = corr(param, sortAge,'type','Spearman');

