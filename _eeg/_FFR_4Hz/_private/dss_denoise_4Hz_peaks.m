%% MCCA denoising

% plot FFR_4Hz and extract peaks

clear all
%cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/zhome/7e/f/64621/Desktop/UHEAL_paper/UHEAL_startup.m')
subs = dir('/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/_derivatives/*.mat')
load('/work3/jonmarc/UHEAL_master/UHEAL_paper/_stats/uheal_data.mat');
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
        %% dss
        clear y
        %time x chan x trials
        dat  = permute(TS_trials(chansoi,:,:),[2,1,3]);
        c0 = nt_cov(dat);
        c1 = nt_cov(mean(dat,3));
        [todss,pwr0,pwr1]=nt_dss0(c0,c1);
        z=nt_mmat(dat,todss);
        % regress out last components,keep 1:2
        tmp=nt_tsr(dat,squeeze(z(:,3:end,:)));
        dat_clean(s,:,:) = nanmean(tmp,3);
        %z_sub{s} = z;
%         close all
%         plot(mean(mean(dat_clean,2),3))
%         hold on
%         plot(mean(mean(dat,2),3))
        %TS_dss(s,:) =squeeze(mean(y(:,:,1),1));
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
        dat_clean(s,:,:) = nan(size(dat_clean(1,:,:)));
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
plot(time(tidx_TS),squeeze(mean(mean(dat_clean(nh_idx,:,:),1),3)))
hold on
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

% on DSS data
baseline_dss = mean(mean(dat_clean(:,baseline_idx,:),2),3);
baseline_chan_dss = mean(dat_clean(:,baseline_idx,:),2);

% baseline corrected time-series
TS_base = TS_sub-baseline'; % mean over 14 scalp channels (chanoi) subjects x time
TS_base_chan = TS_sub_chan-baseline_chan; % subjects x chan x time

TS_dss = mean(dat_clean,3)-baseline_dss;
TS_dss_chan = dat_clean - baseline_chan_dss;

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
hold on
plot(time_TS,nanmean(TS_dss(nh_idx,:)))


%% alternative peak find
% P1 = 0.045 -  0.065 s
% N1 = 0.085 -  0.15 s
% P2 = 0.15  -  0.25 s
% N2 = 0.2   -  0.5 s
fs = 1032;
data_all = {TS_base;TS_dss};
% peak latencies

P1_idx =[0.03 0.085];
N1_idx =[0.065  0.15];
P2_idx =[0.15 0.25];
N2_idx = [0.2 0.3];
p_idx = {P1_idx,N1_idx,P2_idx,N2_idx};
dd=2;
% loop over subjects
for ss=1:size(data_all{dd},1)
    % find average peak latency for each subject
    % find mean waveform over 6 tones
    for tt=1%:6
        this_idx = [0+0.5*(tt-1) 0.5+0.5*(tt-1)];
        this_mean(tt,:) = data_all{dd}(ss,find(time_TS>=this_idx(1) & time_TS<=this_idx(2)));
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
            p_mean(ss,ii,pp) = mean(data_all{dd}(ss,this_t));

            % max/min
            if mod(pp,2) %negative of positive peaks
                [p_max(ss,ii,pp),p_l] = max(data_all{dd}(ss,this_t));
            else
                [p_max(ss,ii,pp),p_l] = min(data_all{dd}(ss,this_t));
            end
            % convert index to time
            p_max_l(ss,ii,pp) = time_TS(this_t(p_l));

            % fixed latency based on mean over all tones
            p_fix(ss,ii,pp) = data_all{dd}(ss,find(time_TS==p_fix_s_l(pp)+0.5*(ii-1)));
            p_fix_l(ss,ii,pp) = p_fix_s_l(pp)+0.5*(ii-1);




        end
    end
end

% outcomes:
% 1) p_mean over latencies
% 2) p_max + p_max_l -> max/min amplitude per tone
% 3) p_fix + p_fix_l -> max/min based on mean tone response (fixed latency)

% compare peaks to each other
close all
titles = {'P50','N1','P2','N2'}
dss_string = {'raw','dss'};
for pp=1:4

    for ii=1:6
        figure(pp)
        subplot(1,6,ii)
        % amplitude
        scatter(p_max(nh_idx_p,ii,pp),p_mean(nh_idx_p,ii,pp))
        lsline
        hold on
        xlim([-6 6])
        ylim([-6 6])
        set(gcf,'position',[440 589 943 115])

        figure(pp+4)
        subplot(1,6,ii)
        scatter(p_max(nh_idx_p,ii,pp),p_fix(nh_idx_p,ii,pp))
        lsline
        hold on
        xlim([-6 6])
        ylim([-6 6])
        set(gcf,'position',[440 589 943 115])
        % latency
        figure(pp+8)
        subplot(1,6,ii)
        scatter(p_max_l(nh_idx_p,ii,pp),p_fix_l(nh_idx_p,ii,pp))
        lsline
        hold on
        xlim([0+0.5*(ii-1) .5+0.5*(ii-1)])
        ylim([0+0.5*(ii-1) .5+0.5*(ii-1)])
        set(gcf,'position',[440 589 943 115])

    end
    figure(pp)
    title([titles{pp} ' max vs. mean'])
    fig = gcf;
    %saveas(fig,['/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/peak_dev/peak_comp' titles{pp} '_max_vs_mean' dss_string{dd}],'epsc');
    figure(pp+4)
    title([titles{pp} ' max vs. fixed'])
    fig = gcf;
    %saveas(fig,['/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/peak_dev/peak_comp' titles{pp} '_max_vs_fixed' dss_string{dd}],'epsc');
    figure(pp+8)
    title([titles{pp} ' max vs. fixed latency'])
    fig = gcf;
    %saveas(fig,['/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/peak_dev/peak_comp' titles{pp} '_max_vs_fixed_latency' dss_string{dd}],'epsc');
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
dds_string = {'raw','dss'};
titles = {'maxmin','mean','fixed'};
for ss=1:3
    figure(ss)
    set(gcf,'renderer','painters')
    for pp=1:3
        for ii=1:6
            subplot(1,6,ii)
            this_t = [0 0.5]+(0.5*(ii-1));
            this_t_idx = find(time_TS>=this_t(1) & time_TS<=this_t(2));
            plot(time_TS(this_t_idx),nanmean(data_all{dd}(idx_all{pp},this_t_idx),1),'color',gcol{pp})
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
        %saveas(fig,['/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/peak_dev/trace_peaks_' titles{ss} '_' dds_string{dd}],'epsc')
    end
end

%% n100-p200
close all
figure('renderer','painters')
% loop over strategies
for ss = 1:3
n1p2=strat{ss}(:,:,3)-strat{ss}(:,:,2);
%p50 = P1;
idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
marks = {'o','sq','^'}
figure(ss)
set(gcf,'renderer','painters')
for pp=1:3
    errorbar(1:6,nanmean(n1p2(idx_all{pp},:)),nanstd(n1p2(idx_all{pp},:))/sqrt(length(idx_all{pp})),'color',gcol{1},'marker',marks{pp},'markerfacecolor',gcol{pp});%,'o-','color',y_col)
    hold on
end
box off
title('P200-N100')
xlabel('tone nr.')
ylabel('\muV')
set(gca,'xtick',[1:6])
xlim([0 7])
hleg=legend({'Young','Mid. aged','Older'});
hleg.Box = 'off'
hleg.Position = [0.5728 0.6010 0.2869 0.1950];
set(gcf,'position',[440 471 289 233])
title(titles{ss})

fig = gcf;
%saveas(fig,['/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/peak_dev/p2-n1_' titles{ss} '_' dds_string{dd}],'epsc')

end
%% stats full model
n1p2 = p_fix(:,:,3) - p_fix(:,:,2);
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
%n1p2=strat{ss}(:,:,3)-strat{ss}(:,:,2);
%p50 = P1;
idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
marks = {'o','sq','^'}
close all
figure(ss)
set(gcf,'renderer','painters')
for pp=1:3
    errorbar(1:2,nanmean(n1p2_mean(idx_all{pp},:)),nanstd(n1p2_mean(idx_all{pp},:))/sqrt(length(idx_all{pp})),'color',gcol{1},'marker',marks{pp},'markerfacecolor',gcol{pp});%,'o-','color',y_col)
    hold on
end
box off
title('P200-N100')
xlabel('tone nr.')
ylabel('\muV')
set(gca,'xtick',[1:6])
xlim([0 3])
hleg=legend({'Young','Mid. aged','Older'});
hleg.Box = 'off'
hleg.Position = [0.5728 0.6010 0.2869 0.1950];
set(gcf,'position',[440 471 289 233])
set(gca,'xtick',[1:2],'xticklabel',{'1st','2:6'})
title(titles{ss})
fig = gcf;
saveas(fig,['/zhome/7e/f/64621/Desktop/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/firstvsrest' titles{ss} '_' dds_string{dd}],'epsc')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% estimate peaks for averaged 2:6 tones in trace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_all = {TS_base;TS_dss};
% peak latencies

P1_idx =[0.03 0.085];
N1_idx =[0.065  0.15];
P2_idx =[0.15 0.25];
N2_idx = [0.2 0.3];
p_idx = {P1_idx,N1_idx,P2_idx,N2_idx};
dd=2;
clear p_fix p_fix_l p_mean p_mean_l p_max p_max_l
% loop over subjects
for ss=1:size(data_all{dd},1)
    % find average peak latency for each subject
    % find for first tone
    for tt=1%:6
        this_idx = [0+0.5*(tt-1) 0.5+0.5*(tt-1)];
        this_mean(tt,:) = data_all{dd}(ss,find(time_TS>=this_idx(1) & time_TS<=this_idx(2)));
    end
    % mean over tones
    this_mean = nanmean(this_mean,1);
    mean_time = time_TS(find(time_TS>=0 & time_TS<=0.5));

    % for tones 2:6
    for tt=2:6
        this_idx = [0+0.5*(tt-1) 0.5+0.5*(tt-1)];
        this_mean_t(tt,:) = data_all{dd}(ss,find(time_TS>=this_idx(1) & time_TS<=this_idx(2)));
    end
    this_mean_t = nanmean(this_mean_t,1)

    this_mean_all{dd}(ss,:,:) = [this_mean;this_mean_t];


    % loop over tones
    for ii=1:2
        for pp=1:4
            % find latency for this peak
            this_p =[p_idx{pp}(1) p_idx{pp}(2)];
            this_t = find(mean_time>=this_p(1) & mean_time<=this_p(2));

            % max/min
            if mod(pp,2) %negative of positive peaks
                [p_max(ss,ii,pp),p_l] = max(this_mean_all{dd}(ss,ii,this_t));
            else
                [p_max(ss,ii,pp),p_l] = min(this_mean_all{dd}(ss,ii,this_t));
            end
            % convert index to time
            p_max_l(ss,ii,pp) = mean_time(this_t(p_l));

            % fixed latency based on mean over all tones
            p_fix(ss,ii,pp) = this_mean_all{dd}(ss,ii,mean_time==p_max_l(ss,1,pp));
            p_fix_l(ss,ii,pp) = p_max_l(ss,ii,pp);




        end
    end
end
%% plotting strategies on traces + peaks as a function of tones
close all


% group plots
idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
marks = {'x','<','>','^'}
% loop over strategies
strat = {p_max;p_fix};
latency = {p_max_l;p_fix_l};
dds_string = {'raw','dss'};
titles = {'maxmin','fixed'};
for ss=1:2
    figure(ss)
    set(gcf,'renderer','painters')
    for pp=1:3
        for ii=1:2
            subplot(1,6,ii)
            this_t = [0 0.5];
            this_t_idx = find(mean_time>=this_t(1) & mean_time<=this_t(2));
            plot(mean_time(this_t_idx),squeeze(nanmean(this_mean_all{dd}(idx_all{pp},ii,this_t_idx),1)),'color',gcol{pp})
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
        saveas(fig,['/zhome/7e/f/64621/Desktop/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/trace_peaks_firstvsrest' titles{ss} '_' dds_string{dd}],'epsc')
    end
end
%% first tone vs. rest
n1p2 = p_fix(:,:,3) - p_fix(:,:,2);
idx = ~isnan(n1p2(:,1));
n1p2_mean = n1p2;%[n1p2(:,1) mean(n1p2(:,2:end),2)];
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
%n1p2=strat{ss}(:,:,3)-strat{ss}(:,:,2);
%p50 = P1;
idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
marks = {'o','sq','^'}
close all
figure(ss)
set(gcf,'renderer','painters')
for pp=1:3
    errorbar(1:2,nanmean(n1p2_mean(idx_all{pp},:)),nanstd(n1p2_mean(idx_all{pp},:))/sqrt(length(idx_all{pp})),'color',gcol{1},'marker',marks{pp},'markerfacecolor',gcol{pp});%,'o-','color',y_col)
    hold on
end
box off
title('P200-N100')
xlabel('tone nr.')
ylabel('\muV')
set(gca,'xtick',[1:6])
xlim([0 3])
hleg=legend({'Young','Mid. aged','Older'});
hleg.Box = 'off'
hleg.Position = [0.5728 0.6010 0.2869 0.1950];
set(gcf,'position',[440 471 289 233])
set(gca,'xtick',[1:2],'xticklabel',{'1st','2:6'})
title(titles{ss})
fig = gcf;
saveas(fig,['/zhome/7e/f/64621/Desktop/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/firstvsrest_estimate2' titles{ss} '_' dds_string{dd}],'epsc')
%% MCCA
nchans = 16
x = permute(squeeze(data_sub(:,1:16,1:end)),[3,2,1]);
xx=x(:,:); % concatenate channelwise
%xx = zscore(xx);
addpath('/work3/jonmarc/UHEAL_master/UHEAL/_scripts/_tools/NoiseTools')
C=xx'*xx;
[A,score,AA]=nt_mcca(C,nchans);
z=xx*A; % common space
%% plot first 8 components
close all
t=0:1/128:length(z)/128-1/128;
figure(2)
for ii=1:8
    subplot(8,1,ii)
plot(t-1,z(:,ii))
subtitle(['SC ' num2str(ii)])
end

%%

z_sub = squeeze(data_sub(1,:,:))'*AA{1}(:,:);
z_sub(:,100:end) = 0;
data_z = z_sub*pinv(AA{1}(:,:))
close all
subplot 121
plot(mean(data_z(:,chansoi),2))
hold on
subplot 122
plot(squeeze(mean(data_sub(1,chansoi,:),2)))



%% dds

c0 = cov(squeeze(data_sub(1,:,:)));
c1 = cov(mean(squeeze(data_sub(1,:,:))));
[todss,pwr0,pwr1]=nt_dds0(c0,c1);
nt_mat(dat,todss)

