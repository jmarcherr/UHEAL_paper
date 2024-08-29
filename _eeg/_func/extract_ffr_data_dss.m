if isfield(data,'FFR')
    % from cell to double
    fs(s) = data.fs;
    TS_sub(s,:,:)   = data.FFR_TS; % timeseries
    fft_sub(s,:,:)  = data.f_fft; % spectrum
    F_sub(s,:)      = data.F; % F-stat per chan
    FFR_sub(s,:)    = data.FFR; %FFR per chan
    SNR_sub(s,:)    = data.FFR_SNR; % SNR per chan
    F_critt(s,:)     = data.F_crit; %critial F
    noise_sub(s,:,:)  = data.fft_noise; %noise per chan
    noise_f_sub(s,:,:) = data.noise.f;
    time{s}         = data.time;


    % dss
    %time x chan x trials
    dat  = permute(data.trial(1:16,:,:),[2,1,3]);

    % get resonator filter at stimulus frequency
    FPEAK=326; % Hz
    Q=8; % determines width
    [b,a]=nt_filter_peak(FPEAK/(fs/2),Q);
    % covariance matrices of full band (c0) and filtered (c1)
    [c0,c1]=nt_bias_filter(dat,b,a);

    % DSS matrix
    [todss,pwr0,pwr1]=nt_dss0(c0,c1);
    p1=pwr1./pwr0; % score, proportional to power ratio of 50Hz & harmonics to full band

    % DSS components
    z=nt_mmat(meg,todss);
    % regress out last components,keep 1
    tmp=nt_tsr(dat,squeeze(z(:,2:end,:)));
    dat_clean(s,:,:) = nanmean(tmp,3);

    clear TS_trials dat
        %
    tidx{s}         = data.tidx;

    %average over Cz,Fz,FCz
    FFR_avg(s)      = data.F_avg;
    SNR_avg(s)      = data.FFR_SNR_avg;
    sig_avg(s)      = data.sig_idx_avg;
    F_avg(s)        = data.F_avg;
    F_crit_avg(s)   = data.F_crit_avg;
    noise_avg(s)    = data.noise_avg; 
    subid{s}        = data.subid;

    fft_freq(s,:) = data.fft_freq;
    stimear(s,:) = data.stimear;
    subinfo{s} = data.subinfo;
    nr_reject(s) =data.nr_reject;
    chan_labels{s} = data.chan_labels;
    chans{s} = data.channels;
    age(s) = data.subinfo.age;
    gender(s) = data.subinfo.gender;
else

    fs(s) = nan;
    TS_sub(s,:,:) = nan(1,16,2458);
    env_sub(s,:,:) = TS_sub(s,:,:);
    fft_sub(s,:,:) = nan(1,16,1025);
    F_sub(s,:) = nan(1,16);
    FFR_sub(s,:) = nan(1,16);
    SNR_sub(s,:) = nan(1,16);
    noise_sub(s,:)=nan(1,16);
    subinfo{s} = data.subinfo;
    subid{s} = data.subid;
    age(s) = data.subinfo.age;
    gender(s) = data.subinfo.gender;
    dat_clean(s,:,:) = nan(1,16,2458);

        %average over Cz,Fz,FCz
    FFR_avg(s)      = nan;
    SNR_avg(s)      = nan;
    sig_avg(s)      = 0; % non sig;
    F_avg(s)        = nan;
    F_crit_avg(s)   = nan;
    noise_avg(s)    = nan;

end

    