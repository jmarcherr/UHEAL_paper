% plot abr and extract peaks
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work1/jonmarc/UHEAL_master/UHEAL_paper/UHEAL_startup.m')
d = dir('_outputs/_derivatives/*.mat')

for dd=1:length(d)
    load([d(dd).folder filesep d(dd).name])
    if isfield(data,'abr')
    sub_abr(dd,:) = data.abr{1};
    t_abr(dd,:) = data.time{1};
    fs(dd) = data.fs;
    sub_peaks{dd} = data.abr_peaks{1};
    AP_amp_pm(dd) = data.abr_peaks{1}.AP_amp-data.abr_peaks{1}.AP_neg;
    WV_amp_pm(dd) = data.abr_peaks{1}.WV_amp-data.abr_peaks{1}.WV_neg;
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
    AP_amp_pm(dd) = nan;
    WV_amp_pm(dd) = nan;
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


%% save to file
load('/work1/jonmarc/UHEAL_master/UHEAL_paper/_clin/clin_data_table/clin_data.mat')
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
save('/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_ABR/_outputs/abr_data_table/abr_data.mat','abr_data');

%% gather traces 
abr_sub_trace = struct;
abr_sub_trace.subid = uheal_data.subid;
abr_sub_trace.sub_abr_b = sub_abr;
abr_sub_trace.t_abr = t_abr';
abr_sub_trace.CP = uheal_data.CP_new;
abr_sub_trace.rjt_sub = rjt_sub'
abr_sub_trace.age  = age';
abr_sub_trace.gender = gender';


save('/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_ABR/_outputs/abr_data_table/abr_sub_trace','-struct','abr_sub_trace')


%% plot
load('/work1/jonmarc/UHEAL_master/UHEAL/uheal_data.mat')
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