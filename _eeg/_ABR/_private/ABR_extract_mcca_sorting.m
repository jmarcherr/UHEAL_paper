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
nchans = 16
x = permute(squeeze(sub_abr(nh_all,1:16,:)),[3,2,1]);
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
plot(t_abr(1,:),z(:,ii))
subtitle(['SC ' num2str(ii)])
end
set(gcf,'position',[250 122 339 589])
fig = gcf;
%saveas(fig,'figs/mcca_comp','epsc')

% sort components for ABR wave 1 (mean between 0.8 and 1.5 ms)
t_idx_w1=find(t_abr(1,:)>=0.0008 & t_abr(1,:)<=0.0019);
for ii=1:size(z,2)
    z_w1(ii) = mean(z(t_idx_w1,ii));
end
[m,i]=sort(z_w1,'descend')
figure
for ii=1:8
    subplot(8,1,ii)
plot(t_abr(1,:),z(:,i(ii)))
subtitle(['SC ' num2str(i(ii))])
end
set(gcf,'position',[250 122 339 589])
%%
figure
data_z = nan(size(sub_abr));
for ss=1:size(x,3)
    z_sub = squeeze(sub_abr(nh_all(ss),:,:))'*AA{(ss)}(:,:);
    %z_sub(:,1) = 0;
    z_sub(:,[setdiff(1:length(z),i(1:5))]) = 0;
    %z_sub(:,[1:7 9:end]) = 0;
    data_z(nh_all(ss),:,:) = permute(z_sub*pinv(AA{ss}(:,:)),[2,1]);

end
subplot 121
plot(t_abr(1,:),squeeze(nanmean(data_z))')
subplot 122
plot(t_abr(1,:),squeeze(nanmean(sub_abr)))


%% plotting just one subject
close all
subplot 121
plot(t_abr(1,:),squeeze(nanmean(data_z(4,:,:)))')
subplot 122
plot(t_abr(1,:),squeeze(nanmean(sub_abr(4,:,:))))

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