
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/zhome/7e/f/64621/Desktop/UHEAL_paper/UHEAL_startup.m')
d = dir('/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR/_outputs/_derivatives/*.mat')
clc
for s=1:length(d)
    load([d(s).folder filesep d(s).name]);
    extract_ffr_data_dss;

    disp([subid{s} ' done...'])
end


 %% load clinical measures
load('/work3/jonmarc/UHEAL_master/UHEAL/uheal_data.mat')

% get age groups
CP = ~uheal_data.CP_new
y = find(age<=25 & CP' & sig_avg);
m = find(age>25 & age<50 & CP' & sig_avg);
o = find(age>=50 & CP' & sig_avg);
nh = find(CP' & sig_avg);
nh_all = find(CP' & ~isnan(TS_sub(:,1,1))')

% get colormap
uheal_colormap;

%% mcca
nchans = 16
x = permute(squeeze(TS_sub(nh_all,1:16,tidx{1})),[3,2,1]);
xx=x(:,:); % concatenate channelwise
%xx = zscore(xx);
C=xx'*xx;
[A,score,AA]=nt_mcca(C,nchans);
z=xx*A; % common space

%% plot first 8 components
close all
t=0:1/128:length(z)/128-1/128;
figure('Renderer','painters')
for ii=1:8
    subplot(8,1,ii)
plot(t-1,z(:,ii))
subtitle(['SC ' num2str(ii)])
end
set(gcf,'position',[440 65 560 637])
fig = gcf;
saveas(fig,'figs/mcca_comp','epsc')

%%
data_z = nan(size(TS_sub));
FFR_mcca = nan(size(TS_sub,1),1);
SNR_mcca = nan(size(TS_sub,1),1);
sig_mcca = nan(size(TS_sub,1),1);
SNR_chan_mcca = nan(size(SNR_sub));
for ss=1:size(x,3)
    z_sub = squeeze(TS_sub(nh_all(ss),:,:))'*AA{(ss)}(:,:);
    z_sub(:,3:end) = 0;
    data_z(nh_all(ss),:,:) = permute(z_sub*pinv(AA{ss}(:,:)),[2,1]);

    % get fft
    foi = [326]; % pure tone frequency

    %get fft FFR
    [f_tmp,fft_sub_tmp,f_fft_noise_tmp,FFR_tmp,F_tmp,SNR_tmp,F_crit_tmp]=get_fft(squeeze(data_z(nh_all(ss),:,tidx{1})),foi,fs(s));
    % selected channels
    chaoi_avg = [find(strcmp(data.chan_labels,'Cz')),...
        find(strcmp(data.chan_labels,'FCz')),...
        find(strcmp(data.chan_labels,'Fz'))];
    SNR_chan_mcca(nh_all(ss),:) = SNR_tmp;

    % get channel average SNR/FFR over Cz, Fz and FCz
    [f_fft_avg,FFR_avg_tmp,F_avg_tmp,SNR_avg_tmp,F_crit_avg_tmp,sig_idx_avg_tmp,noise_avg_tmp]=get_fft_chaoi(f_tmp,fft_sub_tmp,chaoi_avg,foi);
    FFR_mcca(nh_all(ss)) = FFR_avg_tmp;
    SNR_mcca(nh_all(ss)) = SNR_avg_tmp;
    sig_mcca(nh_all(ss)) = sig_idx_avg_tmp;
end
%% compare traces
%nh_all = find(CP);
close all
subplot 121
plot(squeeze(nanmean(data_z(nh_all,1:16,:)))')
subplot 122
plot(squeeze(nanmean(TS_sub(nh_all,1:16,:)))')

%% plot SNR of DSS vs. raw
close all
figure('Renderer','painters')
subplot 121
scatter(SNR_dss(CP),SNR_sub(CP))
hold on
plot([-40 60],[-40 60])
ylim([-40 60])
xlim([-40 60])
xlabel('dss');ylabel('raw')
title(['dss: sig = ' num2str(length(find(sig_dss(CP)))) '/' num2str(length(find(CP)))] )
subplot 122
scatter(SNR_mcca(CP),SNR_sub(CP))
xlabel('mcca');ylabel('raw')
hold on
plot([-40 60],[-40 60])
ylim([-40 60])
xlim([-40 60])
title(['mcca: sig = ' num2str(length(find(sig_mcca(CP)))) '/' num2str(length(find(CP)))])

fig = gcf;

saveas(fig,'figs/dssvsmcca_scatter','epsc')
%% Topoplots of SNR over channels (vertical)
close  all

figure('renderer','painter')
subplot(3,2,1)
zlim = [30 50];
c=jm_topoplot(nanmean(SNR_chan_mcca(y,:))',zlim,'Young',1);
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
c.Ticks = [30 50];

subplot(3,2,3)
c=jm_topoplot(nanmean(SNR_chan_mcca(m,:))',zlim,'Middle-aged',1);
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
c.Ticks = [30 50];

subplot(3,2,5)
c=jm_topoplot(nanmean(SNR_chan_mcca(o,:))',zlim,'Older',1);
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
c.Ticks = [30 50];

set(gcf,'position',[441 318 560 420])

%% Topoplots horizontal

figure('renderer','painter')
subplot(2,3,1)
zlim = [25 40];
c=jm_topoplot(nanmean(SNR_sub(y,:))',zlim,'Young',1);
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
c.Ticks = [30 40];

%set(gca,'Colorbar','off')
subplot(2,3,2)
c=jm_topoplot(nanmean(SNR_sub(m,:))',zlim,'Middle-aged',1);
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
c.Ticks = [30 40];

subplot(2,3,3)
c=jm_topoplot(nanmean(SNR_sub(o,:))',zlim,'Older',1);
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
c.Ticks = [30 40];

set(gcf,'position',[441 457 578 245])
fig = gcf;
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
ts_all = squeeze(mean(nanmean(data_z(find(sig_avg),[10,4,12],:),1),2));

filt_coef = [320 330]; % filter around f_stim
fs = fs(1);
filt_def = designfilt('bandpassfir','FilterOrder',40, ...
    'CutoffFrequency1',filt_coef(1),'CutoffFrequency2',filt_coef(2), ...
    'SampleRate',fs);
y_ts    = filtfilt(filt_def,y_ts);
m_ts    = filtfilt(filt_def,m_ts);
o_ts    = filtfilt(filt_def,o_ts);
ts_all  = filtfilt(filt_def,ts_all);

t = time{1}(find(time{1}>=-0.1));
py=plot(t,y_ts','color',y_col);
hold on
pm=plot(t,m_ts','color',m_col);
po=plot(t,o_ts','color',o_col);

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
sig_mcca(isnan(sig_mcca)) = 0;
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


%% noise floor
figure('renderer','painters')
chanoi = [4,10,12];
scatter(age(nh),db(nanmean(noise_sub(nh,chanoi),2)),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
xlabel('Age')
ylabel('Noise floor (dB mV)')
ll=lsline
set(gca,'fontsize',12)
set(ll,'linewidth',2,'color','k')
set(gcf,'Position',[228 420 280 209]);
[rho,pval]=corr(age(nh)',db(nanmean(noise_sub(nh,chanoi),2)),'type','spearman')
fig = gcf;
%saveas(fig,'figs2/noise_floor','epsc')


%% save to UHEAL_data
%load('/work3/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data.mat')
%uheal_data.FFR_SNR = nan(size(uheal_data.subid));
%uheal_data.FFR_sig = nan(size(uheal_data.subid));

%%
% for s=1:length(SNR_sub)
%     % get this subid
%     thisID = str2double(sub_id{s}(3:5))
%     this_idx = find(uheal_data.subid==thisID);
%     
%     uheal_data.FFR_SNR(this_idx) = SNR_chan(s)';
%     uheal_data.FFR_sig(this_idx) = sig_idx_chan(s)';
% end
% thisdir = cd;
% uheal_data.FFR_SNR = uheal_data.FFR_SNR';
% uheal_data.FFR_sig = uheal_data.FFR_sig';
%cd('/work3/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/')
%uheal_table = struct2table(uheal_data)

%writetable(uheal_table,'uheal_data.csv')  
%save('uheal_data.mat','uheal_data')



%%
function c=jm_topoplot(var1,zlim,tit_string,coff)
load('/work3/jonmarc/UHEAL_master/UHEAL/_EEG/_func/topo_default.mat');
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