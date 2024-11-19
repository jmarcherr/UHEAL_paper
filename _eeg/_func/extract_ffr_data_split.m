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
    tidx{s}         = data.tidx;
    %FFR_trials = data.FFR_trials;
    %FFR_trials(s,:,:,:) = data.FFR_trials;
    FFR_train(s,:,:,:) = data.FFR_train;
    FFR_test(s,:,:,:) = data.FFR_test;
    %trial_step(s,:) = data.trial_step;
  

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
    %FFR_trials(s,:,:,:) = nan(10,16,2458);
    FFR_train(s,:,:,:) = nan(10,16,2458);
    FFR_test(s,:,:,:) = nan(10,16,2458);

        %average over Cz,Fz,FCz
    FFR_avg(s)      = nan;
    SNR_avg(s)      = nan;
    sig_avg(s)      = 0; % non sig;
    F_avg(s)        = nan;
    F_crit_avg(s)   = nan;
    noise_avg(s)    = nan;


end


%% dss

