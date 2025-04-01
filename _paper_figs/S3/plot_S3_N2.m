%% S3 figures

% plot FFR_4Hz and extract peaks

clear all
%cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
subs = dir('/work3/jonmarc/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/_derivatives/*.mat')
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat');
savepath = ['/work3/jonmarc/UHEAL_paper/_paper_figs/S3/figs/']
%% get data
for s=1:length(subs)
    
    load([subs(s).folder filesep subs(s).name])
    clc

    disp(['sub ' subs(s).name(1:5) ' loaded...'])
    sub_num(s) = str2num(subs(s).name(3:5));
    chansoi = setdiff(1:16,[5 11]);
    % get FFR
    if isfield(data,'TS')
        itpc(s,:,:) = data.itpc;
        f = data.f;
        TS_sub(s,:) = nanmean(data.TS(chansoi,:));
        TS_sub_chan(s,:,:)=data.TS;
        TS_trials = data.TS_trials;%chan x time x trial
        clear TS_trials dat

        time = data.time;
        tidx = data.tidx;
        tidx_TS = data.tidx_TS;

        age(s) =data.subinfo.age;
        gender(s) = data.subinfo.gender;

        CP(s) =  uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));
        nr_reject(s) =data.nr_reject;
        chan_labels{s} = data.chan_labels;
        chans{s} = data.channels;
        
    else

        itpc(s,:,:) =nan(16,1537);
        TS_sub(s,:,:) = nan(1,5632);
        TS_sub_chan(s,:,:) = nan(16,5632);
        age(s) = data.subinfo.age;
        gender(s) = data.subinfo.gender;
        CP(s) = uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));
        %dat_clean(s,:,:) = nan(size(dat_clean(1,:,:)));
        z_sub{s} = nan;
    end

    
end
%%
% TS_sub_chan = subject x chan x time

% reject subs
nh_idx = (~CP' & ~isnan(TS_sub_chan(:,1,1)));
data_sub = TS_sub_chan(nh_idx,:,:);

%% plot clean data
close all
plot(time(tidx_TS),squeeze(mean(TS_sub(nh_idx,:))))

% non-nan idx
nh_idx_p = find(nh_idx);

% get peaks
%% Time-series (TS) baseline
%%%%%%%%%%%%%%%%%%%%%%%

time_TS = time(tidx_TS); %-1:4s
% baseline norm from -0.1 - 0 s
baseline_idx = find(time_TS>-0.1 & time_TS<=0);

% estimate baseline
baseline(:) = mean(TS_sub(:,baseline_idx),2);
baseline_chan = mean(TS_sub_chan(:,1:16,baseline_idx),3);

% baseline corrected time-series
TS_base = TS_sub-baseline'; % mean over 14 scalp channels (chanoi) subjects x time
TS_base_chan = TS_sub_chan-baseline_chan; % subjects x chan x time


%% get age groups
% groups
YNH_idx = find(age<=25 & ~CP );
MNH_idx = find(age>25 & age<50 & ~CP )
ONH_idx = find(age>=50 & ~CP);
nh_idx = find(~CP); % all normal hearing
ages = [17 77];
% colormap
uheal_colormap;
% channels
chansoi  = setdiff(1:16,[5 11]); % all channels but T7 and T8

% plot mean over all
close all
plot(time_TS,nanmean(TS_base(nh_idx,:)))
xlim([-0.5 3.5])



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

%% plot peak latencies (max min)
close all
figure('Renderer','painters')
subplot 121
marks = {'sq','o','^','v'}
pcols = cbrewer('qual','Paired',6);
pcols = pcols([1,5,2,6],:);
for pp=1:4
    errorbar(nanmean(p_max_l(nh_idx,1:6,pp)-[0 0.5 1 1.5 2 2.5]),nanstd(p_max_l(nh_idx,1:6,pp)-[0 0.5 1 1.5 2 2.5])/sqrt(length(nh_idx)),'Marker',marks{pp},'markerfacecolor','none','color',pcols(pp,:),'markeredgecolor','none','markersize',6.8)
    hold on
end
hleg = legend('P1','N1','P2','N2')

hleg.Box = 'off'
xlabel('Tone number')
ylabel('latency (s)')
axis padded
box off
set(gca,'fontsize',12,'xtick',[1:6])
set(gcf,'position',[440 248 697 355])
hleg.Position = [0.4686 0.7146 0.1019 0.2183];
fig = gcf;
saveas(fig,[savepath ' latency_peaks'],'svg')
%%
% stats for if latency is different from t 1
for pp=1:4
    this_peak = squeeze(p_max_l(nh_idx,1:6,pp)-[0 0.5 1 1.5 2 2.5])
% reshape
r_dim = size(this_peak,1)*size(this_peak,2); % dimension
peak_tab = reshape(this_peak',r_dim,1)';
subnum_tab = reshape(repmat(sub_num(nh_idx),size(this_peak,2),1),1,r_dim);
age_tab = reshape(repmat(age(nh_idx),size(this_peak,2),1),1,r_dim);
tone_tab = reshape(repmat(1:6,size(this_peak,1),1)',1,r_dim);

a = struct;
a.peak_lat = peak_tab';
a.subnum = categorical(subnum_tab)';
a.age = age_tab';
a.tonenr = categorical(tone_tab)';
tab = sortrows(struct2table(a),'subnum');
lme{pp} = fitlme(tab,'peak_lat~tonenr*age+(1|subnum)')
end

% stats for  latency 1 vs.  rest
for pp=1:4
    this_peak = squeeze(p_max_l(nh_idx,1:6,pp)-[0 0.5 1 1.5 2 2.5])
    this_peak = [this_peak(:,1) nanmean(this_peak(:,2:6),2)];
% reshape
r_dim = size(this_peak,1)*size(this_peak,2); % dimension
peak_tab = reshape(this_peak',r_dim,1)';
subnum_tab = reshape(repmat(sub_num(nh_idx),size(this_peak,2),1),1,r_dim);
age_tab = reshape(repmat(age(nh_idx),size(this_peak,2),1),1,r_dim);
tone_tab = reshape(repmat(1:2,size(this_peak,1),1)',1,r_dim);

a = struct;
a.peak_lat = peak_tab';
a.subnum = categorical(subnum_tab)';
a.age = age_tab';
a.tonenr = categorical(tone_tab)';
tab = sortrows(struct2table(a),'subnum');
lme_12{pp} = fitlme(tab,'peak_lat~age*tonenr+(1|subnum)')
lme_12_lat{pp} = fitlme(tab,'peak_lat~tonenr+(1|subnum)')
end

% plot N2 latency for age groups
% group plots
idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
marks = {'^','sq','o'}
tt_titles = {'P1','N1','P2','N2'}
figure('Renderer','painters')
for pp=1:4
subplot(1,4,pp)
for gg=1:3
    ebgg(gg) = errorbar([1:6]+0.1*(gg-1),nanmean(p_max_l(idx_all{gg},1:6,pp)-[0 0.5 1 1.5 2 2.5]),nanstd(p_max_l(idx_all{gg},1:6,pp)-[0 0.5 1 1.5 2 2.5])/sqrt(length(idx_all{gg})),'color',gcol{gg},'marker',marks{gg},'markerfacecolor',gcol{gg},'MarkerEdgecolor','k','MarkerSize',6.8)
    hold on
end
axis padded
xlabel('Tone number')
ylabel('Latency (s)')
set(gca,'xtick',[1:6],'fontsize',12)
box off
title(tt_titles{pp},'FontWeight','normal')
end
hleg = legend([ebgg(1) ebgg(2) ebgg(3)],{'Young','Mid. aged','Older'})
hleg.Box = 'off';
hleg.Position = [0.8495 0.6874 0.0994 0.2226];
%set(gcf,'position',[18 395 1408 265])
set(gcf,'position',[18 289 1408 352])
fig = gcf;
saveas(fig,[savepath ' latency_peaks_groups'],'svg')

% latencies N2 for groups
clc
g_names = {'y','ma','old'}
for gg = 1:3
disp(['N2 Latency ' g_names{gg} ': fist tone:' num2str(nanmean(p_max_l(idx_all{gg},1,4))) '+/-' num2str(nanstd(p_max_l(idx_all{gg},1,4))) ', 2:6 tones: ' num2str(nanmean(nanmean(p_max_l(idx_all{gg},2:6,4)-[0.5 1 1.5 2 2.5],2))) '+/- ' num2str(nanstd(nanmean(p_max_l(idx_all{gg},2:6,4),2))) ])
end
%% plotting strategies on traces + peaks as a function of tones
close all


% group plots
idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
marks = {'x','<','>','^'}
% loop over strategies
strat = {p_max;p_mean;p_fix};
latency = {p_max_l;p_fix_l;p_fix_l};

titles = {'maxmin','mean','fixed'};
for ss=1:3
    figure(ss)
    set(gcf,'renderer','painters')
    for pp=1:3
        for ii=1:6
            subplot(1,6,ii)
            this_t = [0 0.5]+(0.5*(ii-1));
            this_t_idx = find(time_TS>=this_t(1) & time_TS<=this_t(2));
            plot(time_TS(this_t_idx),nanmean(data_all(idx_all{pp},this_t_idx),1),'color',gcol{pp})
            hold on
            for tt=1:4 % 4 peaks
                plot(nanmean(latency{ss}(idx_all{pp},ii,tt)),nanmean(strat{ss}(idx_all{pp},ii,tt)),'marker',marks{tt},'color',gcol{pp})

                xlim([this_t])
                ylim([-4 3])
                set(gca,'xtick',[this_t(1):0.1:this_t(2)])
                box off
            end
            set(gcf,'position',[256 128 1090 576])
            sgtitle(titles{ss})
        end
        fig = gcf;
        saveas(fig,[savepath titles{ss} '_peaks'],'svg')
    end
end

%% n100-p200
close all
figure('renderer','painters')
% gather N1, P2, p2-n1
peak_all = {}
% loop over strategies
for ss = 1%1:3 % mixed latencies
    n1p2=strat{3}(:,:,3)-strat{3}(:,:,2);
    n1 = strat{3}(:,:,2);
    p2 = strat{3}(:,:,3);
    n2 = strat{1}(:,:,4);
    p1 = strat{3}(:,:,1);
    peak_all = {p1,n1,p2,n2,n1p2};
    tt_titles = {'P1','N1','P2','N2','P2-N1'}

    for tt = 1:size(peak_all,2)
        %p50 = P1;
        idx_all = {YNH_idx,MNH_idx,ONH_idx};
        gcol = {y_col,m_col,o_col};
        marks = {'^','sq','o'}
        figure(tt)
        set(gcf,'renderer','painters')
        for pp=1:3
            errorbar(1:6,nanmean(peak_all{tt}(idx_all{pp},:)),nanstd(peak_all{tt}(idx_all{pp},:))/sqrt(length(idx_all{pp})),'color',gcol{1},'marker',marks{pp},'markerfacecolor',gcol{pp},'MarkerEdgecolor','k','MarkerSize',6.8);%,'o-','color',y_col)
            hold on
        end
        box off
        title(tt_titles)
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
        set(gcf,'position',[440 471 289 233])
        title([tt_titles{tt}],'fontsize',12,'fontweight','normal')

        fig = gcf;
        saveas(fig,[savepath tt_titles{tt} titles{ss}],'svg')


        figure(tt*10)
        set(gcf,'renderer','painters')
        peaks_all_mean = [peak_all{tt}(:,1) mean(peak_all{tt}(:,2:end),2)];
        for pp=1:3
            errorbar(1:2,nanmean(peaks_all_mean(idx_all{pp},:)),nanstd(peaks_all_mean(idx_all{pp},:))/sqrt(length(idx_all{pp})),'color',gcol{1},'marker',marks{pp},'markerfacecolor',gcol{pp},'MarkerEdgecolor','k','MarkerSize',6.8);%,'o-','color',y_col)
            hold on
        end
        box off

        xlabel('Tone number')
        ylabel('\muV')
        set(gca,'xtick',[1:6],'fontsize',12)
        axis padded
        xlim([0 3])
        %hleg=legend({'Young','Mid. aged','Older'});
        %hleg.Box = 'off'
        %hleg.Position = [0.5728 0.6010 0.2869 0.1950];
        set(gcf,'position',[440 471 289 233])
        set(gca,'xtick',[1:2],'xticklabel',{'1st','2:6'})
        title([tt_titles{tt}],'FontSize',12,'fontweight','normal')
        fig = gcf;
        saveas(fig,[savepath tt_titles{tt} ' firstvsrest' titles{ss} '_max_min'],'svg')
    end
end

%%
% LP for plotting
fs = 1024;
filt_coef = [.5 30]; %4 Hz
filt_def = designfilt('bandpassfir', 'FilterOrder', 40, ...
    'CutoffFrequency1', filt_coef(1), 'CutoffFrequency2', filt_coef(2),...
    'SampleRate', fs);
% bypassing filtfilt fieldtrip
addpath('/appl/matlab/990/toolbox/signal/signal/')

for ss=1:size(data_all,1)
    % find average peak latency for each subject
    % find mean waveform over 6 tones
    this_mean = [];
    for tt=1:6
        this_idx = [-0.1+0.5*(tt-1) 0.5+0.5*(tt-1)];
        this_mean(tt,:) = data_all(ss,find(time_TS>=this_idx(1) & time_TS<=this_idx(2)));
    end
    TS_tones(ss,:,:) =this_mean;
end
figure('Renderer','painters')
this_t = [-0.1 0.5];
this_t_idx = find(time_TS>=this_t(1) & time_TS<=this_t(2));
for gg = 1:3
% first
subplot(1,2,1)
shadedErrorBar(time_TS(this_t_idx),filtfilt(filt_def,squeeze(nanmean(TS_tones(idx_all{gg},1,:)))),...
    filtfilt(filt_def,squeeze(nanstd(TS_tones(idx_all{gg},1,:)))/sqrt(length(idx_all{gg}))),...
    'lineprops',{'-','color',gcol{gg},'linewidth',1},'transparent',1)
hold on
ylim([-4.5 2.5])
xlim([-0.1 0.5])
if gg==3
    set(gca,'fontsize',12)
    plot([0 0],[-6 5],'k--')
    plot([0.25 0.25],[-6 5],'k--')
    ylabel('\muV')
    xlabel('Time (s)')
     title(['First Tone'],'FontSize',12,'FontWeight','normal')
end
% last 5
subplot(1,2,2)
s2(gg)=shadedErrorBar(time_TS(this_t_idx),filtfilt(filt_def,squeeze(nanmean(nanmean(TS_tones(idx_all{gg},2:6,:),2),1))),...
    filtfilt(filt_def,squeeze(nanstd(nanmean(TS_tones(idx_all{gg},2:6,:),2),1))/sqrt(length(idx_all{gg}))),...
    'lineprops',{'-','color',gcol{gg},'linewidth',1},'transparent',1)
hold on
ylim([-4.5 2.5])
xlim([-0.1 0.5])
if gg==3
    set(gca,'fontsize',12)
    plot([0 0],[-6 5],'k--')
    plot([0.3 0.3],[-6 5],'k--')
    ylabel('\muV')
    xlabel('Time (s)')
    hleg=legend([s2(1).mainLine s2(2).mainLine s2(3).mainLine],{'Young','Mid. aged','Older'},'fontsize',10);
    hleg.Box = 'off'
    hleg.Position = [0.5915 0.1569 0.1866 0.1167];
    title(['Tone 2:6'],'FontSize',12,'FontWeight','normal')
end

end
set(gcf,'position',[61 215 627 480])
fig = gcf;
saveas(fig,[savepath 'ERP_firstvsrest_filt'],'svg')




%% stats full model
n1p2 = p_max(:,:,3) - p_max(:,:,2);
idx = ~isnan(n1p2(:,1));

% reshape
r_dim = size(n1p2(idx,:),1)*size(n1p2(idx,:),2); % dimension
n1p2_tab = reshape(n1p2(idx,:)',r_dim,1)';
subnum_tab = reshape(repmat(sub_num(idx),size(n1p2(idx,:),2),1),1,r_dim);
age_tab = reshape(repmat(age(idx),size(n1p2(idx,:),2),1),1,r_dim);
tone_tab = reshape(repmat(1:6,size(n1p2(idx,:),1),1)',1,r_dim);

a = struct;
a.n1p2 = n1p2_tab';
a.subnum = categorical(subnum_tab)';
a.age = age_tab';
a.tonenr = categorical(tone_tab)';
tab = sortrows(struct2table(a),'subnum');
lme = fitlme(tab,'n1p2~age*tonenr+(1|subnum)')
% lm model
%% first tone vs. rest

n1p2_mean = [n1p2(:,1) mean(n1p2(:,2:end),2)];
idx = nh_idx;
% reshape
r_dim = size(n1p2_mean(idx,:),1)*size(n1p2_mean(idx,:),2); % dimension
n1p2_tab = reshape(n1p2_mean(idx,:)',r_dim,1)';
subnum_tab = reshape(repmat(sub_num(idx),size(n1p2_mean(idx,:),2),1),1,r_dim);
age_tab = reshape(repmat(age(idx),size(n1p2_mean(idx,:),2),1),1,r_dim);
tone_tab = reshape(repmat(1:2,size(n1p2_mean(idx,:),1),1)',1,r_dim);

a = struct;
a.n1p2 = n1p2_tab';
a.subnum = categorical(subnum_tab)';
a.age = age_tab';
a.tonenr = categorical(tone_tab)';
tab = sortrows(struct2table(a),'subnum');
lme = fitlme(tab,'n1p2~age*tonenr+(1|subnum)')


%% first tone vs. rest for all peaks
clear lme
for tt=1:3%size(peak_all,2)
this_mean = [peak_all{tt}(:,1) mean(peak_all{tt}(:,2:end),2)];
idx = nh_idx;%~isnan(this_mean(:,1));
% reshape
r_dim = size(this_mean(idx,:),1)*size(this_mean(idx,:),2); % dimension
this_tab = reshape(this_mean(idx,:)',r_dim,1)';
subnum_tab = reshape(repmat(sub_num(idx),size(this_mean(idx,:),2),1),1,r_dim);
age_tab = reshape(repmat(age(idx),size(this_mean(idx,:),2),1),1,r_dim);
tone_tab = reshape(repmat(1:2,size(this_mean(idx,:),1),1)',1,r_dim);

a = struct;
a.this = this_tab';
a.subnum = categorical(subnum_tab)';
a.age = age_tab';
a.tonenr = categorical(tone_tab)';
tab = sortrows(struct2table(a),'subnum');
lme{tt} = fitlme(tab,'this~age*tonenr+(1|subnum)')

end

% display
clc
for tt=1:5
    disp(['interaction for ' tt_titles{tt} ' : t = ' num2str(lme{tt}.Coefficients.tStat(4)) ', p =' num2str(lme{tt}.Coefficients.pValue(4))])
end

%% full model 
clear lme
this_mean = [];
for tt=1:size(peak_all,2)
this_mean = peak_all{tt};
idx = nh_idx;%setdif(~isnan(this_mean(:,1)),nh_idx);

% reshape
r_dim = size(this_mean(idx,:),1)*size(this_mean(idx,:),2); % dimension
peak_tab = reshape(this_mean(idx,:)',r_dim,1)';
subnum_tab = reshape(repmat(sub_num(idx),size(this_mean(idx,:),2),1),1,r_dim);
age_tab = reshape(repmat(age(idx),size(this_mean(idx,:),2),1),1,r_dim);
tone_tab = reshape(repmat(1:6,size(this_mean(idx,:),1),1)',1,r_dim);

a = struct;
a.peak = peak_tab';
a.subnum = categorical(subnum_tab)';
a.age = age_tab';
a.tonenr = categorical(tone_tab)';
tab = sortrows(struct2table(a),'subnum');
lme{tt} = fitlme(tab,'peak~age*tonenr+(1|subnum)')

end

%% full model with peak ampitudes as factor
this_mean = [];
m_peaks = peak_all;
for tt=1:size(peak_all,2)-1
this_mean = [this_mean ;m_peaks{tt}'];

end
this_mean = this_mean';
idx = nh_idx;%~isnan(this_mean(:,1));
% reshape
r_dim = size(this_mean(idx,:),1)*size(this_mean(idx,:),2); % dimension
peak_tab = reshape(this_mean(idx,:)',r_dim,1)';
subnum_tab = reshape(repmat(sub_num(idx),size(this_mean(idx,:),2),1),1,r_dim);
age_tab = reshape(repmat(age(idx),size(this_mean(idx,:),2),1),1,r_dim);
tone_tab = reshape(repmat(repmat([1:6],1,4),size(this_mean(idx,:),1),1)',1,r_dim);
peak_name = reshape(repmat(reshape(repmat({'1_P1','2_N1','3_P2','4_N2'},6,1),1,6*4)',size(this_mean(idx,:),1),1),1,r_dim);

a = struct;
a.peak = peak_tab';
a.subnum = categorical(subnum_tab)';
a.age = age_tab';
a.tonenr = categorical(tone_tab)';
a.peak_name = categorical(peak_name)';
tab = sortrows(struct2table(a),'subnum');
lme = fitlme(tab,'peak~peak_name*age+(1|subnum)')


%% full model first vs. rest
this_mean = [];
m_peaks = peak_all;
for tt=1:size(peak_all,2)-1
this_mean = [this_mean ;m_peaks{tt}(:,1)';nanmean(m_peaks{tt}(:,2:6),2)'];

end
this_mean = this_mean';
idx = nh_idx;%~isnan(this_mean(:,1));
% reshape
r_dim = size(this_mean(idx,:),1)*size(this_mean(idx,:),2); % dimension
peak_tab = reshape(this_mean(idx,:)',r_dim,1)';
subnum_tab = reshape(repmat(sub_num(idx),size(this_mean(idx,:),2),1),1,r_dim);
age_tab = reshape(repmat(age(idx),size(this_mean(idx,:),2),1),1,r_dim);
tone_tab = reshape(repmat(repmat([1:2],1,4),size(this_mean(idx,:),1),1)',1,r_dim);
peak_name = reshape(repmat(reshape(repmat({'1_P1','2_N1','3_P2','4_N2'},2,1),1,2*4)',size(this_mean(idx,:),1),1),1,r_dim);

a = struct;
a.peak = peak_tab';
a.subnum = categorical(subnum_tab)';
a.age = age_tab';
a.tonenr = categorical(tone_tab)';
a.peak_name = categorical(peak_name)';
tab = sortrows(struct2table(a),'subnum');
lme = fitlme(tab,'peak~peak_name*age*tonenr+(1|subnum)')


%% full model first - rest :2 step model
this_mean = [];
m_peaks = peak_all;
for tt=1:size(peak_all,2)-1
this_mean = [this_mean ;m_peaks{tt}(:,1)'-nanmean(m_peaks{tt}(:,2:6),2)'];

end
this_mean = this_mean';
idx = nh_idx;%~isnan(this_mean(:,1));
% reshape
r_dim = size(this_mean(idx,:),1)*size(this_mean(idx,:),2); % dimension
peak_tab = reshape(this_mean(idx,:)',r_dim,1)';
subnum_tab = reshape(repmat(sub_num(idx),size(this_mean(idx,:),2),1),1,r_dim);
age_tab = reshape(repmat(age(idx),size(this_mean(idx,:),2),1),1,r_dim);
tone_tab = reshape(repmat(repmat([1],1,4),size(this_mean(idx,:),1),1)',1,r_dim);
peak_name = reshape(repmat(reshape({'1_P1','2_N1','3_P2','4_N2'},1,1*4)',size(this_mean(idx,:),1),1),1,r_dim);

a = struct;
a.peak = peak_tab';
a.subnum = categorical(subnum_tab)';
a.age = age_tab';
a.tonenr = categorical(tone_tab)';
a.peak_name = categorical(peak_name)';
tab = sortrows(struct2table(a),'subnum');
lme = fitlme(tab,'peak~peak_name*age+(1|subnum)')