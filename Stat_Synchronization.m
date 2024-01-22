
%% To compare strength values and surrogate data

clear
f = 3;       % 1 for Beat; 2 for duple, 3 for Triple
sub = 24:46; % Younger: 1:23; Older: 24:46

load('strength.mat')
load(['strength_surr',num2str(f),'.mat'])


sig{1} = strength(:,sub,f);
sig{2} = strength_surr{1,1}(:,sub,:);

p = nan(64,1);
s = nan(64,1);
for i=1:64
    a = nanmean(sig{1}(i,:));
    b = nanmean(sig{2}(i,:,:),2);
    M = nanmean(b,3);
    Sigma = squeeze(nanstd(b,[],3));
    s(i) = (abs(a-M))./Sigma;
    p(i) = erfc(s(i)./sqrt(2));
end

load('connmat.mat')
cluster = ft_findcluster(p<.05,connmat,3);

%% Correlation of strength with age 

load('strength_cluster.mat')
WGA = [32 31 33 36 36 35 32 35 32 32 33 31 31 31 32 33 34 34 33 33 33 29 30 30 30 33 30 30 29 30 35 35 35 35 34 33 33 30 33 28 28 30 30 31 29 29]';
age = 7*WGA+[6 4 3 4 1 4 3 2 4 6 2 5 6 5 3 0 1 2 5 5 0 5 0 0 2 0 4 6 2 6 4 4 0 1 2 5 6 1 5 5 7 3 3 3 5 4]';
[sortAge] = sort(age);

param = strength_cluster(:,f); 
[c, pval] = corr(param, sortAge,'type','Spearman');
 
%% Rayleigh test and std calculation 

load('phase_beat.mat')

for ch = 1:64 
    g(ch) = circ_rtest(phase(ch,sub));
end

stdtopo = circ_std(phase(:,sub)'); 

cluster = ft_findcluster(g'<.05,connmat,3);

%% circular linear correlation

for ch = 1:64
    [rcl(ch,1) pcl(ch,1)] = circ_corrcl(phase(ch,:),sortAge);
end
sigcorr = find(pcl < 0.05);
