function [this_fft,FFR_avg,F_avg,SNR_avg,F_crit_avg,sig_idx_avg,this_noise]=get_fft_chaoi(f,fft_sub,chaoi,fid)

% mean FFT for selected channels
this_fft = nanmean(fft_sub(chaoi,:));

% get SNR
noisebw = 20; % noise band-width
nonfreqs = [fid+[2:2:10] fid-[2:2:10]]'; % harmonics of 2 Hz
nbins = find([f>=(fid-noisebw) & ...
    f<=(fid+noisebw) & f~=fid & f~=nonfreqs(1) & f~=nonfreqs(2)...
    & f~=nonfreqs(3) & f~=nonfreqs(4)...
    & f~=nonfreqs(5) & f~=nonfreqs(6)...
    & f~=nonfreqs(7) & f~=nonfreqs(8)...
    & f~=nonfreqs(9) & f~=nonfreqs(10)])

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