% plot abr and extract peaks
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
d = dir('_outputs/_derivatives/*.mat')
%%
for dd=1:length(d)
    load([d(dd).folder filesep d(dd).name])
    if isfield(data,'abr')
    sub_abr(dd,:,:) = data.abr{1};
    t_abr(dd,:) = data.time{1};
    fs(dd) = data.fs;
    sub_peaks{dd} = data.abr_peaks{1};
    AP_amp_pm(dd,:) = data.abr_peaks{1}.AP_amp-data.abr_peaks{1}.AP_neg;
    WV_amp_pm(dd,:) = data.abr_peaks{1}.WV_amp-data.abr_peaks{1}.WV_neg;
    subid{dd} = data.subid;
    subinfo{dd} = data.subinfo;
    age(dd) = data.subinfo.age;
    CP(dd) = data.subinfo.CP;
    gender(dd) = data.subinfo.gender;
    rjt_sub(dd) = 0;
     else
    sub_abr(dd,:) = nan;
    t_abr(dd,:) = nan;
    fs(dd) = nan;
    sub_peaks{dd} = [];
    AP_amp_pm(dd,:) = nan(1,16);
    WV_amp_pm(dd,:) = nan(1,16);
    subid{dd} = data.subid;
    subinfo{dd} = data.subinfo;
    age(dd) = data.subinfo.age;
    gender(dd) = data.subinfo.gender;
    CP(dd) = data.subinfo.CP;
    rjt_sub(dd) = 1;
    end
    clc
    disp([subid{dd} ' done...'])
end


%% mcca

 %% load clinical measures
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat')

% get age groups
CP = ~uheal_data.CP_new
y = find(age<=25 & CP');
m = find(age>25 & age<50 & CP');
o = find(age>=50 & CP');
nh = find(CP' );
nh_all = find(CP' & ~isnan(sub_abr(:,1,1))')

% get colormap
uheal_colormap;

%% mcca
%select only relevant time data
tidx = find(t_abr(1,:)>=-0.001 & t_abr(1,:)<=0.007);
nchans = 16
x = permute(squeeze(sub_abr(nh_all,1:16,tidx)),[3,2,1]);
xx=x(:,:); % concatenate channelwise
%xx = zscore(xx);
C=xx'*xx;
[A,score,AA]=nt_mcca(C,nchans);
z=xx*A; % common space

%% get weight maps for each subject

for ss=1:size(x,3)
    % subject specific z
    z_sub = squeeze(sub_abr(nh_all(ss),1:16,tidx))'*AA{(ss)}(:,:);
    % get mixing weights (see mcca_demo3)

    a(ss,:)=squeeze(sub_abr(nh_all(ss),1:16,tidx))*nt_normcol(z(:,1))/length(tidx);
end
%% plot first 8 components
close all
t=0:1/fs(1):length(z)/fs(1)-1/fs(1);
figure('Renderer','painters')
for ii=1:8
    subplot(8,1,ii)
plot(t_abr(1,tidx)*1e3,z(:,ii))
subtitle(['SC ' num2str(ii)])
end
set(gcf,'position',[250 122 339 589])
fig = gcf;
%saveas(fig,'figs/mcca_comp','epsc')

% sort components for ABR waveform

for ii=1:size(z,2)

    % mean
    z_w(ii) = mean(z(:,ii));
    %x-correlate with eeg?
    z_w(ii) = max(abs(xcorr(z(:,ii),squeeze(mean(sub_abr(ss,1:16,tidx),2)))));


end
[m,i] = sort(z_w,'descend'); % wave 1 components
figure;
for ii=1:9
    % wave I
    subplot(9,1,ii)
    plot(t_abr(1,tidx)'*1e3,z(:,i(ii)))
    subtitle(['SC ' num2str(i(ii))])
    
end

set(gcf,'position',[600 6 339 688],'renderer','painters')
sgtitle('xcorr')
fig = gcf;
%saveas(fig,'/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_private/figs/SC_WI','svg')
% figure(5)
% set(gcf,'position',[600 6 339 688],'renderer','painters')
% sgtitle('WV')
% fig = gcf;
%saveas(fig,'/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_private/figs/SC_WV','svg')

%% weight maps components
figure(700)
set(gcf,'renderer','painters')
for cc=1:9
    subplot(3,3,cc)
    jm_topoplot(squeeze(a(i(cc),:))',[],[],0)
    subtitle(['SC' num2str(i(cc))])
    set(gca,'Fontsize',10)

end

set(gcf,'position',[797 304 471 395])
fig = gcf;
%saveas(fig,'/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_private/figs/SC_top_WI','svg')


%%
close all
figure
data_z = nan(size(sub_abr(:,:,tidx)));

for cc = 3

        for ss=1:size(x,3)
            z_sub = squeeze(sub_abr(nh_all(ss),:,tidx))'*AA{(ss)}(:,:);
            %z_sub(:,1) = 0;
            z_sub(:,[setdiff(1:length(z),i(1:cc))]) = 0;
            %z_sub(:,[1:7 9:end]) = 0;
            data_z(nh_all(ss),:,:) = permute(z_sub*pinv(AA{ss}(:,:)),[2,1]);

        end
        % scaling
        rms_dat = squeeze(rms(sub_abr(nh_all,:,tidx),3))
        rms_z = squeeze(rms(data_z(nh_all,:,:),3));
        data_z_sc = rms_dat.*(data_z(nh_all,:,:)./rms_z);

        % mean over subjects
        figure(cc)
        subplot(1,3,2)
        plot(t_abr(1,tidx)*1e3,squeeze(nanmean(data_z_sc))')
        title(['mcca '])
        ylim([-0.3 0.3])
        xlim([-0.001 0.007]*1e3)
        box off

        %single subject
        figure(cc*10)
        subplot(1,3,2)
        plot(t_abr(1,tidx),squeeze(nanmean(data_z_sc(3,:,:)))')
        title(['mcca sub '])
       

    end
    %mean
    figure(cc)
    subplot 131
    plot(t_abr(1,:)*1e3,squeeze(nanmean(sub_abr)))
    title('raw')
    ylim([-0.3 0.3])
    xlim([-0.001 0.007]*1e3)
    box off

    set(gcf,'renderer','painters','position',[440 410 560 196])
    fig = gcf;
    %saveas(fig,['/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_private/figs/mcca_components_WIV_' num2str(cc)],'svg')

    %single subject
    figure(cc*10)
    subplot 131
    plot(t_abr(1,tidx),squeeze(nanmean(sub_abr(nh_all(3),:,tidx))))
    title('raw single subject')

    set(gcf,'renderer','painters','position',[440 410 560 196])
    fig = gcf;
    %saveas(fig,['/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_private/figs/mcca_components_WIV_single' num2str(cc)],'svg')




%% single subject

        subplot 122
        plot(t_abr(1,tidx),squeeze(nanmean(data_z_sc(3,:,:)))')
        title('mcca wI single subject')
        subplot 121
        plot(t_abr(1,tidx),squeeze(nanmean(sub_abr(nh_all(3),:,tidx))))
        title('raw single subject')


%% get peaks
for ss=1:size(data_z,1)
    [abr_peaks{ss}] =  get_abr_peaks_chan(squeeze(data_z(ss,:,:)),t_abr(1,:));
    AP_amp_pm_mcca(ss,:) = abr_peaks{ss}.AP_amp-abr_peaks{ss}.AP_neg;
    WV_amp_pm_mcca(ss,:) = abr_peaks{ss}.WV_amp-abr_peaks{ss}.WV_neg;
end

%% compare
close all
scatter(uheal_data.AP_amp_pm(nh_all,:),uheal_data.FFR_SNR(nh_all,:))
figure
scatter(AP_amp_pm_mcca(nh_all,10),uheal_data.FFR_SNR(nh_all,:))
%% gather evertyhing and save
abr_data = struct;
abr_data.subid = uheal_data.subid;
abr_data.SP_amp = nan(size(uheal_data.subid));
%uheal_data.SP_lat = nan(size(uheal_data.subid));
abr_data.AP_amp = nan(size(uheal_data.subid));
abr_data.AP_amp_pm = nan(size(uheal_data.subid));
abr_data.AP_lat =  nan(size(uheal_data.subid));
abr_data.WV_amp = nan(size(uheal_data.subid));
abr_data.WV_amp_pm = nan(size(uheal_data.subid));
abr_data.WV_lat =  nan(size(uheal_data.subid));

for s=1:length(subid)
    % get this subid
    thisID = str2double(subid{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);
    if ~isempty(sub_peaks{s})
    abr_data.SP_amp(this_idx) = sub_peaks{s}.SP_amp;
    %uheal_data.SP_lat(this_idx) = SP_lat(s);
    abr_data.AP_amp(this_idx) = sub_peaks{s}.AP_amp;
    abr_data.AP_lat(this_idx) = sub_peaks{s}.AP_latency;
    abr_data.WV_amp(this_idx) = sub_peaks{s}.WV_amp;
    abr_data.WV_lat(this_idx) = sub_peaks{s}.WV_latency;
    abr_data.AP_amp_pm(this_idx) = AP_amp_pm(s);
    abr_data.WV_amp_pm(this_idx) = WV_amp_pm(s);
    end
end
%save('/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_outputs/abr_data_table/abr_data.mat','abr_data');

%% gather traces 
abr_sub_trace = struct;
abr_sub_trace.subid = uheal_data.subid;
abr_sub_trace.sub_abr_b = sub_abr;
abr_sub_trace.t_abr = t_abr';
abr_sub_trace.CP = uheal_data.CP_new;
abr_sub_trace.rjt_sub = rjt_sub'
abr_sub_trace.age  = age';
abr_sub_trace.gender = gender';


%save('/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_outputs/abr_data_table/abr_sub_trace','-struct','abr_sub_trace')


%% plot
load('/work3/jonmarc/UHEAL_master/UHEAL/uheal_data.mat')
%plot mean 
figure
plot(t_abr(1,:),nanmean(sub_abr(:,:)))
%% 

% normal hearing
idx = ~uheal_data.CP_new
% plot age
close all
y = find(age<=25 & idx');
m = find(age>25 & age<50 & idx');
o = find(age>=50 & idx');
plot(t_abr(1,:),nanmean(sub_abr(y,:)))
hold on
plot(t_abr(1,:),nanmean(sub_abr(m,:)))
plot(t_abr(1,:),nanmean(sub_abr(o,:)))
xlim([-0.0005 0.01])
legend('Young','Mid. aged','Older')
hold off
figure

scatter(uheal_data.AP_amp_pm(idx),AP_amp_pm(idx))
title('AP')
figure
plot(age(idx),AP_amp_pm(idx),'o')
hold on
plot(uheal_data.Age(idx),uheal_data.AP_amp_pm(idx),'x')

figure

scatter(uheal_data.WV_amp_pm(idx),WV_amp_pm(idx))
title('WV')
figure
plot(age(idx),WV_amp_pm(idx),'o')
hold on
plot(uheal_data.Age(idx),uheal_data.WV_amp_pm(idx),'x')
%% manual plotting
close all
for dd=1:length(AP_amp_pm(idx))
    if ~isnan(sub_abr(dd))
    plot(t_abr(1,:),sub_abr(dd,:))
    hold on
    % wave I
    plot(sub_peaks{dd}.AP_latency,sub_peaks{dd}.AP_amp,'o')
    plot(sub_peaks{dd}.AP_neg_latency,sub_peaks{dd}.AP_neg,'x')
    % wave V
    plot(sub_peaks{dd}.WV_latency,sub_peaks{dd}.WV_amp,'o')
    plot(sub_peaks{dd}.WV_neg_latency,sub_peaks{dd}.WV_neg,'x')    
    title(subid{dd})
    pause
    hold off
    else
    end

end

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
colormap(flip(brewermap(100,'BrBG'))) % change the colormap
%colormap(brewermap(64,'YlOrRd')) % change the colormap
if coff
    c=colorbar;
else
    c=nan;
end
end