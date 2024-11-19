
d = dir('/work3/jonmarc/UHEAL_paper/_eeg/_FFR/_outputs/_derivatives/_split/*.mat')
d_noise = dir('/work3/jonmarc/UHEAL_paper/_eeg/_EFR/_outputs/_derivatives/*.mat')
clc
for s=1:length(d)
    load([d(s).folder filesep d(s).name]);

    extract_ffr_data_split;
    
    load([d_noise(s).folder filesep d_noise(s).name])
    if isfield(data,'FFR_TS')
        FFR_noise(s,:,:) = data.FFR_TS;
    else
        FFR_noise(s,:,:) = nan(16,2458);
    end

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

% get colormap
uheal_colormap;