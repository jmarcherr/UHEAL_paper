
function [uheal_data] = gen_FFR_mat(uheal_data,datadir)
d = dir([datadir filesep '*.mat'])
clc
disp(['Processing FFR data ...'])
for s=1:length(d)
    load([d(s).folder filesep d(s).name]);
    extract_ffr_data;

end

uheal_data.FFR_SNR = nan(size(uheal_data.subid));
uheal_data.FFR_sig = nan(size(uheal_data.subid));
uheal_data.FFR_noise = nan(size(uheal_data.subid));

for s=1:length(SNR_avg)
    % get this subid
    thisID = str2double(subid{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);
    
    uheal_data.FFR_SNR(this_idx) = SNR_avg(s)';
    uheal_data.FFR_sig(this_idx) = sig_avg(s)';
    uheal_data.FFR_noise(this_idx) = noise_avg(s)';

end
clc
disp(['FFR data done ...'])
end

