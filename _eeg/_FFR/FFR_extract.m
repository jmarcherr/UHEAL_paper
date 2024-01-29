
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work1/jonmarc/UHEAL_master/UHEAL_paper/UHEAL_startup.m')
d = dir('_outputs/_derivatives/*.mat')
clc
for s=1:length(d)
    load([d(s).folder filesep d(s).name]);
    extract_ffr_data;

    disp([subid{s} ' done...'])
end

 %% load clinical measures
load('/work1/jonmarc/UHEAL_master/UHEAL/uheal_data.mat')

% get age groups
CP = ~uheal_data.CP_new
y = find(age<=25 & CP' & sig_avg);
m = find(age>25 & age<50 & CP' & sig_avg);
o = find(age>=50 & CP' & sig_avg);
nh = find(CP' & sig_avg);

% get colormap
uheal_colormap;

%% Topoplots of SNR over channels (vertical)

figure('renderer','painter')
subplot(3,2,1)
zlim = [25 40];
c=jm_topoplot(nanmean(SNR_sub(y,:))',zlim,'Young',1);
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
c.Ticks = [30 40];

subplot(3,2,3)
c=jm_topoplot(nanmean(SNR_sub(m,:))',zlim,'Middle-aged',1);
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
c.Ticks = [30 40];

subplot(3,2,5)
c=jm_topoplot(nanmean(SNR_sub(o,:))',zlim,'Older',1);
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
c.Ticks = [30 40];

set(gcf,'position',[441 318 560 420])
fig = gcf;

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

figure('renderer','painter');

uheal_colormap;
chanoi = [4,10,12]; % Cz, FCz, Fz


% filter for visualization
y_ts = squeeze(mean(nanmean(TS_sub(y,[10,4,12],:),1),2));
m_ts = squeeze(mean(nanmean(TS_sub(m,[10,4,12],:),1),2));
o_ts = squeeze(mean(nanmean(TS_sub(o,[10,4,12],:),1),2));
ts_all = squeeze(mean(nanmean(TS_sub(find(sig_avg),[10,4,12],:),1),2));

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

plot([0 0.245],[-0.128 -0.128],'-','color',[0.5 0.5 0.5],'linewidth',2)
box off
hleg = legend([py,pm,po],'Young','Middle-aged','Older')
hleg.Box = 'off';
hleg.Position = [0.5518 0.6662 0.4607 0.3260];
xlabel('Time [s]');
ylabel('mV');set(gca,'fontsize',14)
set(gcf,'position',[441 475 369 227])
ylim([-0.15 0.15])
xlim([-.01 0.5])
title('')
fig = gcf;
%saveas(fig,'figs_paper/ts_ymo_new','epsc')

%% SNR corr
figure('renderer','painters')
this_idx = find(sig_avg & CP');
non_idx = find(~sig_avg & CP');
scatter(age(this_idx),SNR_avg(this_idx),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
set(gca,'fontsize',12)
ll=lsline
set(ll,'linewidth',2,'color','k')
hold on
scatter(age(non_idx),SNR_avg(non_idx),'r+')
ylabel('SNR (dB)')

xlabel('Age')
set(gcf,'Position',[228 420 280 209]);
plot([10 80], db(F_crit_avg(1:2)),'k--')
[rho,pval]=corr(age(this_idx)',SNR_avg(this_idx)','type','spearman')
fig = gcf;
%saveas(fig,'figs2/age_corr_snr','epsc')


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
%load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data.mat')
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
%cd('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/')
%uheal_table = struct2table(uheal_data)

%writetable(uheal_table,'uheal_data.csv')  
%save('uheal_data.mat','uheal_data')



%%
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