
d = dir('/work3/jonmarc/UHEAL_paper/_eeg/_FFR/_outputs/_derivatives/_trials2/*.mat')
clc
for s=1:length(d)
    load([d(s).folder filesep d(s).name]);

    extract_ffr_data_trials2;

    disp([subid{s} ' done...'])
end


 %% load clinical measures
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat')

% get age groups
CP = ~uheal_data.CP_new
y = find(age<=25 & CP' & sig_avg);
m = find(age>25 & age<50 & CP' & sig_avg);
o = find(age>=50 & CP' & sig_avg);
nh = find(CP' & sig_avg);
nh_all = find(CP' & ~isnan(TS_sub(:,1,1))')
idx_all = find(~isnan(TS_sub(:,1,1))');

% get young, old and ma in idx_all
y_idx = find(age(idx_all)<=25 & CP(idx_all)');
m_idx = find(age(idx_all)>25 & age(idx_all)<50 & CP(idx_all)')
o_idx = find(age(idx_all)>=50 & CP(idx_all)')


% get colormap
uheal_colormap;