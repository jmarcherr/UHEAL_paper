function [this_fft,FFR_avg,F_avg,SNR_avg,F_crit_avg,sig_idx_avg,this_noise]=get_fft_efr_chaoi(f,fft_sub,chaoi,fid)

% mean FFT for selected channels
this_fft = nanmean(fft_sub(chaoi,:));

% get SNR
%get noise
noisebw = 20; % noise band-width
linenoise = [100]';

nbins = find([f>=(fid-noisebw) & ...
    f<=(fid+noisebw) & f~=fid & f~=linenoise]);
%+/- 4Hz for harmonics of rep rate
nbins(find(f(nbins)>fid+0.5 & f(nbins)<fid+5)) = [];
nbins(find(f(nbins)>fid-5 & f(nbins)<fid-0.5)) = [];

% noisebw = 20;
% linenoise = [100]';
% nbins = find([f>=(fid-noisebw) & ...
%     f<=(fid+noisebw) &  f~=fid & ...
%     f~=linenoise]);

this_noise = mean(this_fft(nbins));

%F-statistic
F_avg=this_fft(find(f==fid))/this_noise;
bg_freq = nbins;
F_crit_avg = finv(0.99,2,2*length(bg_freq));
% significant?
if F_avg>=F_crit_avg
    sig_idx_avg = 1;
else
    sig_idx_avg = 0;
end

%SNR & FFR amplitude
SNR_avg = db(this_fft(find(f==fid)))-db(this_noise);
FFR_avg =this_fft(find(f==fid));

end