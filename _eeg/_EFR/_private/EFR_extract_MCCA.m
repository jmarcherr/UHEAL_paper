%% plot EFR results
clear all;close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))

run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
subs = dir('/work3/jonmarc/UHEAL_paper/_eeg/_EFR/_outputs/_derivatives/*.mat')
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat');
%load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data_current_EEG/uheal_data.mat');

%% get data
for s=1:length(subs)
    load(subs(s).name)
    clc

    sub_num(s) = str2num(subs(s).name(3:5));
    extract_efr_data_mcca;

    disp(['sub ' subs(s).name(1:5) ' loaded...'])
end

%% rejections
mean(nr_reject)
std(nr_reject)
 
 %% load clinical measures
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat')

% get age groups
CP = ~uheal_data.CP_new
y = find(age<=25 & CP' & sig_idx_avg);
m = find(age>25 & age<50 & CP' & sig_idx_avg);
o = find(age>=50 & CP' & sig_idx_avg);
nh = find(CP' & sig_idx_avg);
nh_all = find(CP' & ~isnan(TS(:,1,1))')

% get colormap
uheal_colormap;

%% mcca
nchans = 16
x = permute(squeeze(TS(nh_all,1:16,tidx)),[3,2,1]);
xx=x(:,:); % concatenate channelwise
%xx = zscore(xx);
C=xx'*xx;
[A,score,AA]=nt_mcca(C,nchans);
z=xx*A; % common space

%% plot first 8 components
close all
t=0:1/fs(1):length(z)/fs(1)-1/fs(1);
figure('Renderer','painters')
for ii=1:8
    subplot(8,1,ii)
plot(t-0.1,z(:,ii))
subtitle(['SC ' num2str(ii)])
end
set(gcf,'position',[440 65 560 637])
fig = gcf;
saveas(fig,'figs/mcca_comp','epsc')

%%
data_z = nan(size(TS));
FFR_mcca = nan(size(TS,1),1);
SNR_mcca = nan(size(TS,1),1);
sig_mcca = zeros(size(TS,1),1);
SNR_chan_mcca = nan(size(SNR_avg));
t_an = find(t>=0 & t<=0.5);
for ss=1:size(x,3)
    z_sub = squeeze(TS(nh_all(ss),:,:))'*AA{(ss)}(:,:);
    z_sub(:,[1 2 6:end]) = 0;
    data_z(nh_all(ss),:,:) = permute(z_sub*pinv(AA{ss}(:,:)),[2,1]);

        % ------------------------ Analysis ------------------------
    foi = [120]; % AM frequency

    %get fft EFR
    [f_tmp,fft_sub_tmp,f_fft_noise_tmp,FFR_tmp,F_tmp,SNR_tmp,F_crit_tmp]=get_fft_efr(squeeze(data_z(nh_all(ss),:,t_an)),foi,fs(1));

    % selected channels
    chaoi_avg = [find(strcmp(data.chan_labels,'Cz')),...
        find(strcmp(data.chan_labels,'FCz')),...
        find(strcmp(data.chan_labels,'Fz'))];

    % get channel average SNR/EFR over Cz, Fz and FCz
    [f_fft_avg,FFR_avg_tmp,F_avg_tmp,SNR_avg_tmp,F_crit_avg_tmp,sig_idx_avg_tmp,noise_avg_tmp]=get_fft_efr_chaoi(f_tmp,fft_sub_tmp,chaoi_avg,foi);
    

    FFR_mcca(nh_all(ss)) = FFR_avg_tmp;
    SNR_mcca(nh_all(ss)) = SNR_avg_tmp;
    sig_mcca(nh_all(ss)) = sig_idx_avg_tmp;
end
%% compare traces
%nh_all = find(CP);
close all
subplot 121
plot(squeeze(nanmean(data_z(nh_all,chaoi_avg,:)))','k')
subplot 122
plot(squeeze(nanmean(TS(nh_all,chaoi_avg,:)))','k')

%% plot SNR of DSS vs. raw
close all
figure('Renderer','painters')
 subplot 121
 scatter(SNR_avg(CP),SNR_avg(CP))
hold on
plot([-40 60],[-40 60])
ylim([-40 60])
xlim([-40 60])
xlabel('raw');ylabel('raw')
title(['raw: sig = ' num2str(length(find(sig_idx_avg(CP)))) '/' num2str(length(find(CP)))] )
subplot 122
scatter(SNR_mcca(CP),SNR_avg(CP))
xlabel('mcca');ylabel('raw')
hold on
plot([-40 60],[-40 60])
ylim([-40 60])
xlim([-40 60])
title(['mcca: sig = ' num2str(length(find(sig_mcca(CP)))) '/' num2str(length(find(CP)))])

fig = gcf;

saveas(fig,'figs/dssvsmcca_scatter','epsc')

%% time series
%%
close all
figure('renderer','painter');

uheal_colormap;
chanoi = [4,10,12]; % Cz, FCz, Fz


% filter for visualization
y_ts = squeeze(mean(nanmean(data_z(y,[10,4,12],:),1),2));
m_ts = squeeze(mean(nanmean(data_z(m,[10,4,12],:),1),2));
o_ts = squeeze(mean(nanmean(data_z(o,[10,4,12],:),1),2));
ts_all = squeeze(mean(nanmean(data_z(find(sig_mcca),[10,4,12],:),1),2));

filt_coef = [110 130]; % filter around f_stim
fs = fs(1);
filt_def = designfilt('bandpassfir','FilterOrder',40, ...
    'CutoffFrequency1',filt_coef(1),'CutoffFrequency2',filt_coef(2), ...
    'SampleRate',fs);
y_ts    = filtfilt(filt_def,y_ts);
m_ts    = filtfilt(filt_def,m_ts);
o_ts    = filtfilt(filt_def,o_ts);
ts_all  = filtfilt(filt_def,ts_all);

%t = time{1}(find(time{1}>=-0.1));
py=plot(t-0.1,y_ts','color',y_col);
hold on
pm=plot(t-0.1,m_ts','color',m_col);
po=plot(t-0.1,o_ts','color',o_col);

plot([0 0.245],[-3 -3]*1e-3,'-','color',[0.5 0.5 0.5],'linewidth',2)
box off
hleg = legend([py,pm,po],'Young','Middle-aged','Older')
hleg.Box = 'off';
hleg.Position = [0.5518 0.6662 0.4607 0.3260];
xlabel('Time [s]');
ylabel('mV');set(gca,'fontsize',14)
set(gcf,'position',[441 475 369 227])
ylim([-3.5 3]*1e-3)
xlim([-.01 0.5])
title('')
fig = gcf;
saveas(fig,'figs/ts_ymo_mcca','epsc')

%% SNR corr
figure('renderer','painters')
%sig_mcca(isnan(sig_mcca)) = 0;
this_idx = find(sig_mcca' & CP');
non_idx = find(~sig_mcca' & CP');
scatter(age(this_idx),SNR_mcca(this_idx),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
set(gca,'fontsize',12)
ll=lsline
set(ll,'linewidth',2,'color','k')
hold on
scatter(age(non_idx),SNR_mcca(non_idx),'r+')
ylabel('SNR (dB)')

xlabel('Age')
set(gcf,'Position',[228 420 280 209]);
plot([10 80], db(F_crit_avg(1:2)),'k--')
[rho,pval]=corr(age(this_idx)',SNR_mcca(this_idx),'type','spearman')
fig = gcf;
saveas(fig,'figs/age_corr_snr','epsc')
%% time series
close all
uheal_colormap
figure('renderer','painters')
% filter for visualization
chaoi_avg = [10 12 4]; % Cz FCz Fz
y_ts = squeeze(mean(nanmean(TS(YNH_idx,chaoi_avg,:))));
m_ts = squeeze(mean(nanmean(TS(MANH_idx,chaoi_avg,:))));
o_ts = squeeze(mean(nanmean(TS(ONH_idx,chaoi_avg,:))));
ts_all = squeeze(mean(nanmean(TS(find(sig_idx_avg),chaoi_avg,:))));

filt_coef = [110 130];
fs = 4096;
filt_def = designfilt('bandpassfir','FilterOrder',40, ...
    'CutoffFrequency1',filt_coef(1),'CutoffFrequency2',filt_coef(2), ...
    'SampleRate',fs);
y_ts =filtfilt(filt_def,y_ts);
m_ts = filtfilt(filt_def,m_ts);
o_ts = filtfilt(filt_def,o_ts);
ts_all = filtfilt(filt_def,ts_all);
t = time(find(time>=-0.1 & time<=0.5));
py=plot(t,y_ts','color',y_col);
hold on
pm=plot(t,m_ts','color',m_col);
po=plot(t,o_ts','color',o_col);
box off
plot([0 0.3],[-0.08 -0.08],'-','color',[0.5 0.5 0.5],'linewidth',2)
hleg = legend([py,pm,po],'Young','Middle-aged','Older')
hleg.Box = 'off';
hleg.Position = [0.5518 0.6662 0.4607 0.3260];
xlabel('Time [s]');
ylabel('\muV');set(gca,'fontsize',12)
set(gcf,'position',[441 475 369 227])
ylim([-0.15 0.15])
xlim([-0.01 0.5])
title('')
fig = gcf;

%saveas(fig,'figs_paper/ts_ymo','epsc')

% all

figure('renderer','painter')
plot(t,ts_all,'k')
hold on
plot([0 0.3],[-0.08 -0.08],'-','color',[0.5 0.5 0.5],'linewidth',2)
xlabel('Time [s]');
ylabel('\muV');set(gca,'fontsize',12)
set(gcf,'position',[441 475 369 227])
xlim([-0.01 0.5])
ylim([-0.1 0.1])
box off
%title(['all sig. n= ' num2str(length(find(sig_idx(:,10))))])
fig = gcf;
%saveas(fig,'figs_paper/ts_all','epsc')

%% Age EFR
close  all
scatter(age(find(sig_idx_avg==1)),SNR_avg(find(sig_idx_avg==1)))
hold  on
scatter(age(find(sig_idx_cz ==1)),FFR_SNR(find(sig_idx_cz==1),10))
[rho,pval]=corr(age(NH_idx)',FFR_SNR(NH_idx,10))
[rho,pval]=corr(age(NH_idx)',SNR_avg(NH_idx)')


%% noise floor
close all
scatter(age(NH_idx),db(noise_avg(NH_idx)),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
xlabel('Age')
ylabel('Noise floor (dB mV)')
ll=lsline
set(gca,'fontsize',12)
set(ll,'linewidth',2,'color','k')
set(gcf,'Position',[228 420 280 209]);
[rho,pval]=corr(age(NH_idx)',db(noise_avg(NH_idx))')
fig = gcf;
%saveas(fig,'figs_paper/noise_floor','epsc')



%% save to UHEAL_data
load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data_current_EEG/uheal_data.mat');
uheal_data.EFR_SNR = nan(size(uheal_data.subid));
uheal_data.EFR_sig = nan(size(uheal_data.subid));

%%
for s=1:length(SNR_sub)
    % get this subid
    thisID = sub_num(s);
    this_idx = find(uheal_data.subid==thisID);
    
    uheal_data.EFR_SNR(this_idx) = db(F_sub(s,10));
    uheal_data.EFR_sig(this_idx) = sig_idx(s,10)';
end

thisdir = cd;
cd('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data_current_EEG/')
%uheal_table = struct2table(uheal_data)

%writetable(uheal_table,'uheal_data.csv')  
save('uheal_data.mat','uheal_data')

cd(thisdir)








function c=jm_topoplot(var1,zlim,tit_string,coff)
load('/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_func/topo_default.mat');
freq.powspctrm = var1;%nanmean(F_sub(YNH_idx,:))';
cfg = [];
cfg.comment = 'no';
cfg.marker = 'on';
cfg.maarkersymbol = '.';
cfg.layout = 'biosemi64.lay';
cfg.channel = freq.cfg.channel;
cfg.parameter = 'powspctrm';
cfg.style = 'straight';
cfg.zlim = zlim;
ft_topoplotER(cfg,freq);
title(tit_string)
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flip(brewermap(100,'RdBu'))) % change the colormap
%colormap(brewermap(64,'YlOrRd')) % change the colormap
if coff
    c=colorbar;
else
    c=nan;
end
end