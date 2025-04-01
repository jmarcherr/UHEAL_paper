function [f,fft_sub,f_fft_noise,FFR,F,SNR,F_crit]=get_fft_efr(data,foi,fs)
chans = size(data,1);
for cc=1:chans  
    fid = foi;
    tmp_data = squeeze(data(cc,:));
    M=tmp_data;
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
    f_fft_sub_pow = fft_sub(cc,find(f==fid));
    
    %get noise
    noisebw = 20; % noise band-width
    linenoise = [100]';
    nbins = find([f>=(fid-noisebw) & ...
        f<=(fid+noisebw) & f~=fid & f~=linenoise]);
    %+/- 4Hz for harmonics of rep rate
    nbins(find(f(nbins)>fid+0.5 & f(nbins)<fid+5)) = [];
    nbins(find(f(nbins)>fid-5 & f(nbins)<fid-0.5)) = [];

    f_fft_noise(cc,:) = mean(fft_sub(cc,nbins));
    
    %F-statistic
    F(cc)=f_fft_sub_pow/f_fft_noise(cc,:);
    bg_freq = nbins;
    F_crit(cc) = finv(0.99,2,2*length(bg_freq));
    
    %SNR
    SNR(cc) = db(f_fft_sub_pow)-db(f_fft_noise(cc,:));
    FFR(cc) =f_fft_sub_pow;

    
    
end

end