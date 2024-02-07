function [uheal_data]=gen_EFR_mat(uheal_data,datadir)

d = dir([datadir filesep '*.mat'])
clc
disp(['Processing EFR data ...'])
%% get data
for s=1:length(d)
    
    load([d(s).folder filesep d(s).name]);
    % get FFR
    if isfield(data,'FFR')

        % 3 channel average
        FFR_avg(s) = data.FFR_avg;
        SNR_avg(s) = data.FFR_SNR_avg;
        F_avg(s) = data.F_avg;
        F_crit_avg(s) = data.F_crit_avg;
        sig_idx_avg(s) = data.sig_idx_avg;
        noise_avg(s) = data.noise_avg;


    else

        % 3 channel average
        FFR_avg(s) = nan;
        SNR_avg(s) = nan;
        F_avg(s) = nan;
        F_crit_avg(s) = nan;
        sig_idx_avg(s) = 0;
        noise_avg(s) = nan;
    end
   
    subid{s} = data.subid;
    
end

uheal_data.EFR_SNR = nan(size(uheal_data.subid));
uheal_data.EFR_sig = nan(size(uheal_data.subid));
uheal_data.EFR_noise = nan(size(uheal_data.subid));

for s=1:length(SNR_avg)
    % get this subid
    thisID = str2double(subid{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);
    
    uheal_data.EFR_SNR(this_idx) = SNR_avg(s)';
    uheal_data.EFR_sig(this_idx) = sig_idx_avg(s)';
    uheal_data.EFR_noise(this_idx) = noise_avg(s)';


end
clc
disp(['EFR data done ...'])
end

