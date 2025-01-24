
d = dir('/work3/jonmarc/UHEAL_paper/_eeg/_FFR/_outputs/_derivatives/_ERP/*.mat')
clc
for s=1:length(d)
    load([d(s).folder filesep d(s).name]);

    extract_erp_data_mcca;

    disp([subid{s} ' done...'])
end


 %% load clinical measures
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat')

% get age groups
CP = ~uheal_data.CP_new
%y = find(age<=25 & CP' & sig_avg);
%m = find(age>25 & age<50 & CP' & sig_avg);
%o = find(age>=50 & CP' & sig_avg);
%nh = find(CP' & sig_avg);
nh_all = find(CP' & ~isnan(TS_sub(:,1,1))')
idx_all = find(~isnan(TS_sub(:,1,1))');

% get colormap
uheal_colormap;