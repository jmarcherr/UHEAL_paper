
clear all
close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
d = dir([rootdir '/_eeg/_FFR/_outputs/_derivatives/*.mat'])
savepath = [rootdir '/_paper_figs/S1/figs/'];
clc
for s=1:length(d)
    load([d(s).folder filesep d(s).name]);
    extract_ffr_data;

    disp([subid{s} ' done...'])
end

 %% load clinical measures
load('/work3/jonmarc//UHEAL_paper/_stats/uheal_data.mat')

% get age groups
CP = uheal_data.CP_SG
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
saveas(fig,[savepath 'FFR_top_vertical'],'svg')

%% alt size
figure('renderer','painter')
subplot(3,2,1)
zlim = [25 40];
c=jm_topoplot(nanmean(SNR_sub(y,:))',zlim,[],0);
%c.Label.String = 'SNR (dB)';
%c.FontSize = 10;
%c.Ticks = [30 40];

subplot(3,2,3)
c=jm_topoplot(nanmean(SNR_sub(m,:))',zlim,[],0);
%c.Label.String = 'SNR (dB)';
%c.FontSize = 10;
%c.Ticks = [30 40];

subplot(3,2,5)
c=jm_topoplot(nanmean(SNR_sub(o,:))',zlim,[],0);
%c.Label.String = 'SNR (dB)';
%c.FontSize = 10;
%c.Ticks = [30 40];

set(gcf,'position',[834 282 396 242])

fig = gcf;
saveas(fig,[savepath 'FFR_top_vert_small'],'svg')

% alt size for even smaller

set(gcf,'position',[680 283 520 420])

fig = gcf;
saveas(fig,[savepath 'FFR_top_vert_smaller'],'svg')

% color bar
figure('Renderer','painters')
subplot(3,2,5)
c=jm_topoplot(nanmean(SNR_sub(o,:))',zlim,[],1);
c.Label.String = 'SNR (dB)';
c.FontSize = 10;
c.Ticks = [25  40];
c.Location = 'southoutside';
set(gcf,'position',[440 187 160 420])
fig = gcf;
saveas(fig,[savepath 'FFR_top_colorbar'],'svg')
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

saveas(fig,[savepath 'FFR_top_horizontal'],'svg')
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
set(gcf,'position',[141 251 369 290])
hleg = legend([py,pm,po],'Young','Middle-aged','Older')
hleg.Box = 'off';
hleg.Position = [0.5819 0.7257 0.3794 0.2034];
xlabel('Time (s)');
ylabel('\muV');set(gca,'fontsize',12,'ytick',[-0.1,0,.1])

ylim([-0.15 0.15])
xlim([-.01 0.5])
title('')
fig = gcf;
saveas(fig,[savepath 'FFR_TS'],'svg')

% small size for combined ABR FFR figure
set(gcf,'position',[293 480 448 215])
hleg = legend([py,pm,po],'Young','Middle-aged','Older')
hleg.Box = 'off';
hleg.Position = [0.6344 0.6763 0.3125 0.2744];
fig = gcf;
saveas(fig,[savepath 'FFR_TS_small'],'svg')

%% FFR age

data_ffr = uheal_data;
data_ffr.CP_new = ~uheal_data.CP_SG;
data_ffr.FFR_SNR(find(data_ffr.FFR_sig==0))=nan;
var = 'FFR_SNR'
labels = 'FFR_{SNR} dB'
lims = [-10 75]
figure

[h,ax1,ax2]=stat_plots_uh_age_pt(data_ffr,var,labels,lims)

scatter(uheal_data.Age(find(uheal_data.FFR_sig==0 & uheal_data.gender==1)),uheal_data.FFR_SNR(find(uheal_data.FFR_sig==0 & uheal_data.gender==1)),'ko','markeredgecolor',[0.3 0.3 0.3],'SizeData',10,'linewidth',1)
scatter(uheal_data.Age(find(uheal_data.FFR_sig==0 & uheal_data.gender==2)),uheal_data.FFR_SNR(find(uheal_data.FFR_sig==0 & uheal_data.gender==2)),'ko','markeredgecolor',[0.3 0.3 0.3],'SizeData',10,'linewidth',1)
%set(gcf,'position',[100 100 297 199]*(96/72))
set(gcf,'position',[483 219 396 320])
%set(hleg,'Visible','off')
fig = gcf;
saveas(fig,[savepath 'FFR_corr_age'],'svg')

% small size for ABR FFR combined figure
set(gcf,'position',[100 100 297 169]*(96/72))
fig = gcf;
saveas(fig,[savepath 'FFR_corr_age_small'],'svg')

%%
function c=jm_topoplot(var1,zlim,tit_string,coff)
load('/work3/jonmarc/UHEAL_master/UHEAL/_EEG/_func/topo_default.mat');
freq.powspctrm = var1;%nanmean(F_sub(YNH_idx,:))';
cfg = [];
cfg.comment = 'no';
cfg.marker = 'on';
cfg.maarkersymbol = '.';
cfg.maskparameter = 0.5;
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