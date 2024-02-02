%% plot EFR results
clear all;close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))

run('/work1/jonmarc/UHEAL_master/UHEAL_paper/UHEAL_startup.m')
subs = dir('_outputs/_derivatives/*.mat')
%load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data.mat');
load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data_current_EEG/uheal_data.mat');

%% get data
for s=1:length(subs)
    
    load(subs(s).name)
    clc
    disp(['sub ' subs(s).name(1:5) ' loaded...'])
    sub_num(s) = str2num(subs(s).name(3:5));
    % get FFR
    if isfield(data,'FFR')
        FFR(s,:) = data.FFR;
        FFR_SNR(s,:) = data.FFR_SNR;
        F(s,:) = data.F;
        F_crit(s,:) = data.F_crit;
        sig_idx(s,:) = F(s,:)>F_crit(s,:);
        sig_idx_cz(s) = sig_idx(s,10);
        subid{s} = data.subid;
        f_fft(s,:,:) = data.f_fft;
        fft_freq(s,:) = data.fft_freq;
        TS(s,:,:) = data.FFR_TS;
        time = data.time;
        tidx = data.tidx;
        fft_noise(s,:) = data.fft_noise;
        noisef(s,:,:) = data.noise.fft_sub_noise;
        subinfo{s} = data.subinfo;

        % 3 channel average
        FFR_avg(s) = data.FFR_avg;
        SNR_avg(s) = data.FFR_SNR_avg;
        F_avg(s) = data.F_avg;
        F_crit_avg(s) = data.F_crit_avg;
        sig_idx_avg(s) = data.sig_idx_avg;
        noise_avg(s) = data.noise_avg;

        if isempty(data.subinfo.age)|isempty(data.subinfo.gender)
            age(s) = nan;
            gender(s) = nan;
        else
        age(s) =data.subinfo.age;
        gender(s) = data.subinfo.gender;
        end

        nr_reject(s) =data.nr_reject;
        chan_labels{s} = data.chan_labels;
        chans{s} = data.channels;
        
    elseif ~isfield(data.subinfo,'age')
        subinfo{s} = data.subinfo;
        age(s) = nan;
    else
        FFR(s,:) = nan(1,16);
        FFR_SNR(s,:) = nan(1,16);
        F(s,:)  = nan(1,16);
        F_crit(s,:) = nan(1,16);
        sig_idx(s,:) = zeros(1,16);
        sig_idx_cz(s) = 0;
        f_fft(s,:,:) = nan(16,1025);
        fft_freq(s,:) = nan(1,1025);
        TS(s,:,:) = nan(16,2458);

        fft_noise(s,:) = nan(1,16);
        noisef(s,:,:) = nan(16,1025);
        subinfo{s} = data.subinfo;

        % 3 channel average
        FFR_avg(s) = nan;
        SNR_avg(s) = nan;
        F_avg(s) = nan;
        F_crit_avg(s) = nan;
        sig_idx_avg(s) = 0;
        noise_avg(s) = nan;
        subinfo{s} = data.subinfo;
        age(s) = nan;
    end
    
        CP(s) = uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));

    
end

%% rejections
mean(nr_reject)
std(nr_reject)
 
%% NH idx & significant
%   get significant responses per channel
%for ii

YNH_idx = find(age<=25 & ~CP & sig_idx_avg); 
MANH_idx = find(age>28 & age<50 & ~CP & sig_idx_avg);
ONH_idx = find(age>=50 & ~CP & sig_idx_avg);
NH_idx = find(~CP & sig_idx_avg);

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