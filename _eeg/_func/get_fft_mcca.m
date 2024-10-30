function [f,fft_sub,f_fft_noise,FFR,F,SNR,F_crit]=get_fft_mcca(data,foi,fs)
chans = size(data,1);
for cc=1:chans  
    fid = foi;
    tmp_data = squeeze(data(cc,:));
    M = tmp_data;
    
    f_fft = [];
    %FFT
    f_fft = fft(M)/(length(M)/2);
    %Convert to power
    f_fft_pow(cc,:) = abs(f_fft.^2); %
    %Truncate negative freqencies
    fft_sub(cc,:) = squeeze((f_fft_pow(cc,1:end/2+1)));
    %Frequency vector
    f = fs/2*linspace(0,1,length(fft_sub(cc,:)));
    % get complex value at signal bin
    f_fft_cmplx = f_fft(find(f==fid));
    %get powerbin
    f_fft_sub_pow(cc) = fft_sub(cc,find(f==fid));
    
    %get noise (all other bins than the target bin
   % bw = 20;
    %nbins = find(f~=fid & f>=(fid-bw) & f<=(fid+bw));
    noisebw = 20;
    nonfreqs = [fid+[2:2:10] fid-[2:2:10]]';
    linenoise = [300 250]';
       nbins = find([f>=(fid-noisebw) & ...
            f<=(fid+noisebw) & f~=fid & f~=nonfreqs(1) & f~=nonfreqs(2)...
            & f~=nonfreqs(3) & f~=nonfreqs(4)...
            & f~=nonfreqs(5) & f~=nonfreqs(6)...
            & f~=nonfreqs(7) & f~=nonfreqs(8)...
            & f~=nonfreqs(9) & f~=nonfreqs(10)]); 

    f_fft_noise(cc) = mean(fft_sub(cc,nbins));
    
    %F-statistic
    F(cc)=f_fft_sub_pow(cc)./f_fft_noise(cc);
    bg_freq = nbins;
    F_crit(cc) = finv(0.99,2,2*length(bg_freq));
    
    %SNR
    SNR(cc) = db(f_fft_sub_pow(cc)')-db(f_fft_noise(cc)');
    FFR(cc) =f_fft_sub_pow(cc);

    
    
end