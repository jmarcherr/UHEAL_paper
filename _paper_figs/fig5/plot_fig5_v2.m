%% plot FFR_4Hz results
% plot FFR_4Hz and extract peaks
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
subs = dir([rootdir '/_eeg/_FFR_4Hz/_outputs/_derivatives/*.mat'])
load([rootdir '/_stats/uheal_data.mat']);
savepath = [rootdir '/_paper_figs/fig5/figs/v2/'];
%% get data
for s=1:length(subs)
    
    load([subs(s).folder filesep subs(s).name])
    clc
    disp(['sub ' subs(s).name(1:5) ' loaded...'])
    sub_num(s) = str2num(subs(s).name(3:5));
    chansoi = setdiff(1:16,[5 11]);
    % get FFR
    if isfield(data,'itpc')
        itpc(s,:,:) = data.itpc;
        f = data.f;
        TS_sub(s,:) = nanmean(data.TS(chansoi,:));
        TS_sub_chan(s,:,:)=data.TS;

        time = data.time;
        tidx = data.tidx;
        tidx_TS = data.tidx_TS;
      
        subinfo{s} = data.subinfo;
        if isempty(data.subinfo.age)|isempty(data.subinfo.gender)
            age(s) = nan;
            gender(s) = nan;
        else
        age(s) =data.subinfo.age;
        gender(s) = data.subinfo.gender;
        end

        CP(s) =  uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));%data.subinfo.CP;
        nr_reject(s) =data.nr_reject;
        chan_labels{s} = data.chan_labels;
        chans{s} = data.channels;
        
    else

        itpc(s,:,:) =nan(16,1537);
        TS_sub(s,:,:) = nan(1,5632);
        TS_sub_chan(s,:,:) = nan(16,5632);
        subinfo{s} = data.subinfo;
        age(s) = data.subinfo.age;
        gender(s) = data.subinfo.gender;
    end

    
end
 
%% 
mean(nr_reject)
std(nr_reject)

%% get age groups
% mean spectrum
itpc_spec_mean=squeeze((nanmean(itpc(:,chansoi,:),2)));
% groups
YNH_idx = find(age<=25 & ~CP & ~isnan(itpc_spec_mean(:,1))');
MNH_idx = find(age>25 & age<50 & ~CP & ~isnan(itpc(:,1))')
ONH_idx = find(age>=50 & ~CP & ~isnan(itpc(:,1))');
ages = [17 77];
% colormap
uheal_colormap;
% channels
chansoi  = setdiff(1:16,[5 11]); % all channels but T7 and T8

%% age vs. ITPC frequencies
fid = [2:2:20];
close all

% plot idx
idx = find(~isnan(itpc_spec_mean(:,1)) & ~CP');
% find harmonic frequencies
for ii=1:length(fid)
    fid_idx(ii) = find(f==fid(ii));
    itpc_freq(ii,:) = itpc_spec_mean(:,find(f==fid(ii)))';
end

% Get ratio
itpc_F0 = itpc_freq(1,:);
itpc_F19 = mean(itpc_freq(2:10,:));
itpc_ratio = log10(itpc_freq(1,:)./mean(itpc_freq(2:10,:)));
[rho_ratio,pval_ratio]=corr(age(idx)',itpc_ratio(idx)','type','spearman')


%% ITPC harmonics plot
%close all
figure('renderer','painters')
subplot(1,3,1)
harm = [0:9];
noise_freqs = find(f<=20);
noise_freqs = setdiff(noise_freqs,fid_idx);
itpc_noise= nanmean(itpc_spec_mean(:,noise_freqs));
% plotting init
marks_ag = {'^','sq','o'}
col_ag = {'k','b','r'}
col_ag = {y_col,m_col,o_col}
ag_idx = {YNH_idx,MNH_idx,ONH_idx};
        for ag=1:3
            p_y = plot(harm+1,nanmean(itpc_freq(1:10,ag_idx{ag}),2),'-','color',[col_ag{ag} 0.8])%,%'marker',marks_ag{ag},'markeredgecolor',col_ag{ag},'markerfacecolor','w','linewidth',0.5);
            hold on
            for ll=1:length(harm)

                errorbar(harm(ll)+1,nanmean(itpc_freq(ll,ag_idx{ag}),2),nanstd(itpc_freq(ll,ag_idx{ag}),1)/sqrt(length(ag_idx{ag})),'color',col_ag{ag})
                eby(ag) = plot(harm(ll)+1,nanmean(itpc_freq(ll,ag_idx{ag}),2),'marker',marks_ag{ag},'color',[col_ag{ag} 0.8],'MarkerFaceColor',[col_ag{ag}]);
            end

        end

plot([1 10],[mean(itpc_noise) mean(itpc_noise)],'--','color',[0.75 0.75 0.75],'linewidth',1)
xlabel('Harmonic nr.')
ylabel('ITPC')
set(gca,'xtick',[1:2:10],'xticklabels',[0:2:9],'fontsize',12,'FontName','Arial')

xlim([0 11])
ylim([0.05 0.45])
box off
set(gcf,'position',[316 430 937 261])
[229 226 396 318]
hleg=legend([eby(1),eby(2),eby(3)],'Young','Middle-aged','Older','fontsize',10);hleg.Box = 'off';
hleg.Position = [[0.1726 0.6538 0.1722 0.2605]];



fig = gcf;
saveas(fig,[savepath 'itpc_harmonics_all'],'svg')
%% mean spectrum
figure('renderer','painters')
subplot(1,3,1)
po=plot(f,mean(itpc_spec_mean(idx,:)),'color',y_col,'linewidth',1)
xlim([1 21])
box off
xlabel('Frequency (Hz)')
ylabel('ITPC')
set(gcf,'position',[441 582 443 120])
fig = gcf;

%% Time series alone
% LP for plotting
fs = 1024;
filt_coef = [.5 30]; %4 Hz
filt_def = designfilt('bandpassfir', 'FilterOrder', 40, ...
    'CutoffFrequency1', filt_coef(1), 'CutoffFrequency2', filt_coef(2),...
    'SampleRate', fs);
% bypassing filtfilt fieldtrip
addpath('/appl/matlab/990/toolbox/signal/signal/')

%close all
figure('renderer','painters')
%%%%%%%%%%%%%%%%%%%%%%%
% baseline norm
time_TS = time(tidx_TS);
baseline(:) = mean(TS_sub(:,find(time_TS>-0.1 & time_TS<=0)),2);
baseline_chan=mean(TS_sub_chan(:,1:16,find(time_TS>-0.1 & time_TS<=0)),3);
TS_base = TS_sub-baseline';
%TS_base = filt_def,TS_base;
TS_base_chan = TS_sub_chan-baseline_chan;

ynh_trace = squeeze(nanmean(TS_base(YNH_idx,:)));
ynh_var = squeeze(nanstd(TS_base(YNH_idx,:)))/sqrt(length(YNH_idx));
mnh_trace = squeeze(nanmean(TS_base(MNH_idx,:)));
mnh_var = squeeze(nanmean(TS_base(MNH_idx,:)))/sqrt(length(MNH_idx));
onh_trace = squeeze(nanmean(TS_base(ONH_idx,:)));
onh_var = squeeze(nanstd(TS_base(ONH_idx,:)))/sqrt(length(ONH_idx));

% baseline
ynh_base = ynh_trace-mean(ynh_trace(1,find(time_TS>-0.1 & time_TS<=0)));
mnh_base = mnh_trace-mean(mnh_trace(1,find(time_TS>-0.1 & time_TS<=0)));
onh_base = onh_trace-mean(onh_trace(1,find(time_TS>-0.1 & time_TS<=0)));


eb1=shadedErrorBar(time(tidx_TS),filtfilt(filt_def,ynh_base),filtfilt(filt_def,ynh_var),'lineprops',{'Color',y_col})
hold on
eb2=shadedErrorBar(time(tidx_TS),filtfilt(filt_def,mnh_base),filtfilt(filt_def,mnh_var),'lineprops',{'Color',m_col})
eb3=shadedErrorBar(time(tidx_TS),filtfilt(filt_def,onh_base),filtfilt(filt_def,onh_var),'lineprops',{'Color',o_col})
p1=plot(time(tidx_TS),filtfilt(filt_def,ynh_base),'color',y_col)
p2=plot(time(tidx_TS),filtfilt(filt_def,mnh_base),'color',m_col)
p3=plot(time(tidx_TS),filtfilt(filt_def,onh_base),'color',o_col)

xlim([-.1 3.5])
ylim([-6 3])

box off
ylabel('\muV')
xlabel('Time (s)')
set(gca,'fontsize',12)

set(gcf,'position',[101 464 668 231])
hold on
fs = 1024;
tstim  = 0:1/fs:0.25
ramp = hann(100);win=[ramp(1:end/2-1)' ones(length(tstim)-100,1)' ramp(end/2:end)'];
stim = [sin(2*pi*256*tstim).*win zeros(size(tstim))];
stim_all = [];
for ii=1:6
    stim_all = [stim_all stim]
end
stim_all = [stim_all zeros(size([tstim tstim]))];
stime = 0:1/fs:length(stim_all)/fs-1/fs;
plot(stime,(stim_all*0.5)-5,'color',[0.5 0.5 0.5 0.5])
hleg = legend([p1 p2 p3],'Young','Middle-aged','Older','fontsize',10);
hleg.Box = 'off'
hleg.Position = [0.2297 0.8552 0.6683 0.1037];
hleg.Orientation = 'horizontal'
fig = gcf;
saveas(fig,[savepath 'TS_all'],'svg')

%% scatters 
% 4Hz neg
var = 'Neg_4Hz';
labels = '\muV'
lims = [-3 2];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
%set(gcf,'position',[462 556 297 197])
%set(ax1,'FontSize',11);set(ax2,'FontSize',11)
set(gcf,'position',[229 226 396 318])
fig = gcf;
saveas(fig,[savepath 'Neg_age'],'svg')

%4Hz itpc_ratio
var = 'ITPC_ratio';
labels = 'ITPC ratio'
lims = [-1 1];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
%set(gcf,'position',[462 556 297 197])
%set(ax1,'FontSize',11);set(ax2,'FontSize',11)
set(gcf,'position',[229 226 396 318])
fig = gcf;
saveas(fig,[savepath 'itpc_age'],'svg')
%% Topographies mean neg + ITPC ratio
close all
%figure('renderer','painter','position',[680 381 396 322])

top_idx ={YNH_idx,MNH_idx,ONH_idx};
spidx = [1 4 7 10;2 5 8 11];
%chanoi = setdiff(1:16,[5,6,11,13]);
chanoi = 1:16
itpc_sub = itpc;
%subplot(1,4,1)
itpc_ratio = log10(squeeze(itpc_sub(:,:,fid_idx(1)))./squeeze(nanmean(itpc_sub(:,:,fid_idx(2:end)),3)))
for ii=1:3


    zlim1 = [0.1 0.4]; % fundamental    
    zlim2 = [0.05 0.2]; % harmonics
    zlim3 = [-.8 .8]; % negativity
    zlim4 = [-0.1 0.5]; % itpc ratio
    con=0
    figure(1);set(gcf,'renderer','painters')
    subplot(3,3,ii+(ii-1)*2)
    jm_topoplot(squeeze(nanmean(itpc_ratio(top_idx{ii},chanoi),1))',[zlim4],[''],con);
    figure(2);set(gcf,'renderer','painters')
    subplot(3,3,ii+(ii-1)*2)
    jm_topoplot(nanmean(squeeze(nanmean(TS_base_chan(top_idx{ii},chanoi,find(time<=0 & time<=3)),3)))',zlim3,[''],con);
    

end
    %scale
    figure(2)
    subplot(3,3,3)
    c=jm_topoplot(zeros(16,1),[zlim3],'scale',1)
    c.FontSize =10;
    c.Label.String = 'Ampl. mV'
    c.Location = 'southoutside'
    figure(1)
    subplot(3,3,3)
    c=jm_topoplot(zeros(16,1),[zlim4],'scale',1)
    c.FontSize =10;
    c.Label.String = 'ITPC ratio'
    c.Location = 'southoutside'

figure(2)
set(gcf,'position',[680 283 520 420])
fig = gcf;
saveas(fig,[savepath 'mean_neg_top'],'svg')
figure(1)
set(gcf,'position',[680 283 520 420])
fig = gcf;
saveas(fig,[savepath 'itpc_ratio_top'],'svg')

%% alternative peak find
% P1 = 0.045 -  0.065 s
% N1 = 0.065 -  0.15 s
% P2 = 0.15  -  0.25 s
% N2 = 0.2   -  0.4 s
fs = 1032;
data_all = TS_base;% {TS_base;TS_dss};
% peak latencies

P1_idx =[0.045 0.065];
N1_idx =[0.065  0.15];
P2_idx =[0.15 0.25];
N2_idx = [0.2 0.4];
p_idx = {P1_idx,N1_idx,P2_idx,N2_idx};

% loop over subjects
for ss=1:size(data_all,1)
    % find average peak latency for each subject
    % find mean waveform over 6 tones
    for tt=1%:6
        this_idx = [0+0.5*(tt-1) 0.5+0.5*(tt-1)];
        this_mean(tt,:) = data_all(ss,find(time_TS>=this_idx(1) & time_TS<=this_idx(2)));
    end
    % mean over tones
    this_mean = nanmean(this_mean,1);

    mean_time = time_TS(find(time_TS>=0 & time_TS<=0.5));
    % find mean peaks
    for pp=1:4
        this_p =[p_idx{pp}(1) p_idx{pp}(2)];
        this_t = find(mean_time>=this_p(1) & mean_time<=this_p(2));
        if mod(pp,2)
            [p_fix_s(pp),p_l] = max(this_mean(this_t));
        else
            [p_fix_s(pp),p_l] = min(this_mean(this_t));
        end
        % convert to time
        p_fix_s_l(pp) = mean_time(this_t(p_l));
    end

    % loop over tones
    for ii=1:6
        for pp=1:4
            % find latency for this peak
            this_p =[p_idx{pp}(1)+0.5*(ii-1) p_idx{pp}(2)+0.5*(ii-1)];
            this_t = find(time_TS>=this_p(1) & time_TS<=this_p(2));

            % mean over latencies
            p_mean(ss,ii,pp) = mean(data_all(ss,this_t));

            % max/min
            if mod(pp,2) %negative of positive peaks
                [p_max(ss,ii,pp),p_l] = max(data_all(ss,this_t));
            else
                [p_max(ss,ii,pp),p_l] = min(data_all(ss,this_t));
            end
            % convert index to time
            p_max_l(ss,ii,pp) = time_TS(this_t(p_l));

            % fixed latency based on mean over all tones
            p_fix(ss,ii,pp) = data_all(ss,find(time_TS==p_fix_s_l(pp)+0.5*(ii-1)));
            p_fix_l(ss,ii,pp) = p_fix_s_l(pp)+0.5*(ii-1);




        end
    end
end

% outcomes:
% 1) p_mean over latencies
% 2) p_max + p_max_l -> max/min amplitude per tone
% 3) p_fix + p_fix_l -> max/min based on mean tone response (fixed latency)

%% plot timeseries with triangles like JH
% LP for plotting
fs = 1024;
filt_coef = [.5 30]; %4 Hz
filt_def = designfilt('bandpassfir', 'FilterOrder', 40, ...
    'CutoffFrequency1', filt_coef(1), 'CutoffFrequency2', filt_coef(2),...
    'SampleRate', fs);
% bypassing filtfilt fieldtrip
addpath('/appl/matlab/990/toolbox/signal/signal/')

%close all
figure('renderer','painters')
%%%%%%%%%%%%%%%%%%%%%%%
% baseline norm
time_TS = time(tidx_TS);
baseline(:) = mean(TS_sub(:,find(time_TS>-0.1 & time_TS<=0)),2);
baseline_chan=mean(TS_sub_chan(:,1:16,find(time_TS>-0.1 & time_TS<=0)),3);
TS_base = TS_sub-baseline';
%TS_base = filt_def,TS_base;
TS_base_chan = TS_sub_chan-baseline_chan;

ynh_trace = squeeze(nanmean(TS_base(YNH_idx,:)));
ynh_var = squeeze(nanstd(TS_base(YNH_idx,:)))/sqrt(length(YNH_idx));
mnh_trace = squeeze(nanmean(TS_base(MNH_idx,:)));
mnh_var = squeeze(nanmean(TS_base(MNH_idx,:)))/sqrt(length(MNH_idx));
onh_trace = squeeze(nanmean(TS_base(ONH_idx,:)));
onh_var = squeeze(nanstd(TS_base(ONH_idx,:)))/sqrt(length(ONH_idx));

% baseline
ynh_base = ynh_trace-mean(ynh_trace(1,find(time_TS>-0.1 & time_TS<=0)));
mnh_base = mnh_trace-mean(mnh_trace(1,find(time_TS>-0.1 & time_TS<=0)));
onh_base = onh_trace-mean(onh_trace(1,find(time_TS>-0.1 & time_TS<=0)));


eb1=shadedErrorBar(time(tidx_TS),filtfilt(filt_def,ynh_base),filtfilt(filt_def,ynh_var),'lineprops',{'Color',y_col})
hold on
eb2=shadedErrorBar(time(tidx_TS),filtfilt(filt_def,mnh_base),filtfilt(filt_def,mnh_var),'lineprops',{'Color',m_col})
eb3=shadedErrorBar(time(tidx_TS),filtfilt(filt_def,onh_base),filtfilt(filt_def,onh_var),'lineprops',{'Color',o_col})
p1=plot(time(tidx_TS),filtfilt(filt_def,ynh_base),'color',y_col)
p2=plot(time(tidx_TS),filtfilt(filt_def,mnh_base),'color',m_col)
p3=plot(time(tidx_TS),filtfilt(filt_def,onh_base),'color',o_col)

xlim([-.1 3.5])
ylim([-6.5 3.5])

box off
ylabel('\muV')
xlabel('Time (s)')
set(gca,'fontsize',12)

%set(gcf,'position',[101 464 668 256])
set(gcf,'position',[[101 417 668 274]])
hold on
fs = 1024;
tstim  = 0:1/fs:0.25
ramp = hann(100);win=[ramp(1:end/2-1)' ones(length(tstim)-100,1)' ramp(end/2:end)'];
stim = [sin(2*pi*256*tstim).*win zeros(size(tstim))];
stim_all = [];
for ii=1:6
    stim_all = [stim_all stim]
end
stim_all = [stim_all zeros(size([tstim tstim]))];
stime = 0:1/fs:length(stim_all)/fs-1/fs;
plot(stime,(stim_all*0.5)-5.5,'color',[0.5 0.5 0.5 0.5])


hold on
gcol =    [0.6980    0.8745    0.5412];
    %[ 0.4980    0.7882    0.4980 ]
gcoln2 = [0.2000    0.6275    0.1725];
for ii=1:6
    %p2 for old
plot(squeeze(nanmean(p_max_l(ONH_idx,ii,3))),squeeze(nanmean(p_max(ONH_idx,ii,3)))*1.2,'v','markerfacecolor',gcol,'color',gcol)
    %n2 for young
    if ii==1
        plot(squeeze(nanmean(p_max_l(YNH_idx,ii,4))),squeeze(nanmean(p_fix(YNH_idx,ii,4)))*1.27,'^','markerfacecolor',gcoln2,'color',gcoln2)
    elseif ii==5
        plot(squeeze(nanmean(p_max_l(YNH_idx,ii,4))),squeeze(nanmean(p_fix(YNH_idx,ii,4)))*2.25,'^','markerfacecolor',gcoln2,'color',gcoln2)

    else
        plot(squeeze(nanmean(p_max_l(YNH_idx,ii,4))),squeeze(nanmean(p_fix(YNH_idx,ii,4)))*1.9,'^','markerfacecolor',gcoln2,'color',gcoln2)
    end
end

hleg = legend([p1 p2 p3],'Young','Middle-aged','Older','fontsize',10);
hleg.Box = 'off'
hleg.Position = [0.2297 0.8552 0.6683 0.1037];
hleg.Orientation = 'horizontal'
fig = gcf;
saveas(fig,[savepath 'TS_all_peaks'],'svg')

%% plot P2 over tone number
figure('renderer','painters')
% gather N1, P2, p2-n1
peak_all = {}

p2 = p_fix(:,:,3);

peak_all = {p2};
tt_titles = {'P2'}
tt = 1
%p50 = P1;
idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
marks = {'^','sq','o'}
set(gcf,'renderer','painters')
for pp=1:3
    errorbar(1:6,nanmean(peak_all{tt}(idx_all{pp},:)),nanstd(peak_all{tt}(idx_all{pp},:))/sqrt(length(idx_all{pp})),'color',gcol{pp},'marker',marks{pp},'markerfacecolor',gcol{pp},'MarkerEdgecolor',gcol{pp},'MarkerSize',6.8);%,'o-','color',y_col)
    hold on
end
box off

% title(tt_titles)
xlabel('Tone number')
ylabel('\muV')
set(gca,'xtick',[1:6])
set(gca,'fontsize',12)
axis padded
xlim([0 7])
if tt==5
    hleg=legend({'Young','Mid. aged','Older'});
    hleg.Box = 'off'
    hleg.Position = [0.5728 0.6010 0.2869 0.1950];
    hleg.FontSize = 10;

end
set(gcf,'position',[440 471 362 261])

%title([tt_titles{tt}],'fontsize',12,'fontweight','normal')

fig = gcf;
saveas(fig,[savepath '/p2_tonenumber'],'svg')





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
if coff
c=colorbar;
end
end