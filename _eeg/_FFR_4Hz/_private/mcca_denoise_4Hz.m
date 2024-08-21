%% MCCA denoising

% plot FFR_4Hz and extract peaks

clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_master/UHEAL_paper/UHEAL_startup.m')
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
        % regress out last components,keep 1:3
        tmp=nt_tsr(dat,squeeze(z(:,4:end,:)));
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

% plot every subject
nh_idx_p = find(nh_idx);
for ii=1:length(nh_idx_p)
    close all
    plot(time(tidx_TS),squeeze(nanmean(dat_clean(nh_idx_p(ii),:,:),3)),'k')
    hold on
    plot(time(tidx_TS),squeeze(TS_sub(nh_idx_p(ii),:)),'color',[0.7 0.2 0.1 0.5])
    hleg = legend('clean','orig')
    hleg.Position = [ 0.7796    0.7102    0.1448    0.2432];
    set(gca,'fontsize',12)
    title([subs(nh_idx_p(ii)).name(1:5)])
    set(gcf,'position',[440 521 670 183],'renderer','painters')
    xlabel('time (s)')
    grid on
    box off
    fig = gcf;
    saveas(fig,['/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/DSS_traces/' [subs(nh_idx_p(ii)).name(1:5)]] ,'epsc')

end
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

