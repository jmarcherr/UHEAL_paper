% plot abr and extract peaks
function[uheal_data]= gen_ABR_mat(uheal_data,datadir)
d = dir([datadir filesep '*.mat']);
    clc
    disp(['Processing ABR data ...'])
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

end



%% gather evertyhing and save
uheal_data.SP_amp = nan(size(uheal_data.subid));
uheal_data.SP_amp_peak = nan(size(uheal_data.subid));
uheal_data.SP_lat = nan(size(uheal_data.subid));

% manual SP peaks
% check plot_fig2_SP_man.m for extracted peaks
load('/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_outputs/SP_manual_peaks/sp_man_peaks.mat')
uheal_data.SP_amp_man = sp_man_peaks(2,:)';
uheal_data.SP_lat_man = sp_man_peaks(1,:)';

uheal_data.AP_amp = nan(size(uheal_data.subid));
uheal_data.AP_amp_pm = nan(size(uheal_data.subid));
uheal_data.AP_lat =  nan(size(uheal_data.subid));
uheal_data.WV_amp = nan(size(uheal_data.subid));
uheal_data.WV_amp_pm = nan(size(uheal_data.subid));
uheal_data.WV_lat =  nan(size(uheal_data.subid));
uheal_data.abr_IV_ratio = nan(size(uheal_data.subid));
uheal_data.abr_SPAP_ratio = nan(size(uheal_data.subid));
uheal_data.abr_baseline = nan(size(uheal_data.subid));


for s=1:length(subid)
    % get this subid
    thisID = str2double(subid{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);
    if ~isempty(sub_peaks{s})


        uheal_data.SP_amp(this_idx) = sub_peaks{s}.SP_amp;%-sub_peaks{s}.baseline;
        uheal_data.SP_amp_peak(this_idx) = sub_peaks{s}.SP_amp_peak;
        uheal_data.SP_lat(this_idx) = sub_peaks{s}.SP_latency;
        uheal_data.AP_amp(this_idx) = sub_peaks{s}.AP_amp;
        uheal_data.AP_lat(this_idx) = sub_peaks{s}.AP_latency;
        uheal_data.WV_amp(this_idx) = sub_peaks{s}.WV_amp;
        uheal_data.WV_lat(this_idx) = sub_peaks{s}.WV_latency;
        uheal_data.AP_amp_pm(this_idx) = AP_amp_pm(s);
        uheal_data.WV_amp_pm(this_idx) = WV_amp_pm(s);
        uheal_data.abr_IV_ratio(this_idx) = sub_peaks{s}.AP_amp/sub_peaks{s}.WV_amp;
        %uheal_data.abr_SPAP_ratio(this_idx) = (sub_peaks{s}.SP_amp-sub_peaks{s}.baseline)/(sub_peaks{s}.AP_amp-sub_peaks{s}.baseline);
        uheal_data.abr_SPAP_ratio(this_idx) = (uheal_data.SP_amp_man(this_idx)-sub_peaks{s}.baseline)/(sub_peaks{s}.AP_amp-sub_peaks{s}.baseline);
        %uheal_data.SP_amp_man(this_idx) = uheal_data.SP_amp_man(this_idx)-sub_peaks{s}.baseline;
        uheal_data.abr_baseline(this_idx) = sub_peaks{s}.baseline;
    end
end
disp(['ABR data done'])
end