
%% FFR analysis script
function data_ffr = EFR_analysis(data)
%% analysis script for preprocessed ffr data (EFR_preproc.mat)
% Outputs weighted average EFR as well as EFR_SNR and F statistics
try

    % check data
    if nargin < 1
        error('No iput. Please provide data from EFR_preproc');
    elseif  ~isfield(data,'hdr')
        error(['No processed FFR data for subject ' data.subid]);
    elseif isfield(data,'error')
        data.error_preproc = data.error;
        error([data.error])
    end

    
    % Get scalp EEG channels
    chaoi = [1:16];

    % Get trial ids
    trials_oi =find(data.trialinfo);

    % select data
    cfg = [];
    cfg.channel = {data.label{[chaoi]}};
    data_cond = ft_selectdata(cfg,data);
    time = data.time{1};

    % define time idx
    fs = data.fsample;
    tidx = find(time>=-0.1 & time<0.5); % pt + noise
    tidx_FFR = find(time>=0 & time<0.5);
    tidx_noisef = find(time>=0.3-1/fs & time<0.5);

    % epoch data
    epoched_data =epoch_data_3(data_cond,1,trials_oi);   %epoch x chan x time

    % artifact rejection and weighting
    [valid_trials,rjt] = threshold_rejection(epoched_data,50,tidx,trials_oi);

    % clean epoched data
    data_cc = [];
    % reject if less than half the trials remain
    if length(valid_trials)<[length(data.trial)/2]
        error(['Less than half of trials are clean. Rejecting ' data.subid]);
    else
        data_cc = epoched_data(valid_trials,:,:); % clean data
        data_cc = permute(data_cc,[2,3,1]);
    end
    data.valid_trials = valid_trials;

    clc
    nr_reject = length(rjt)/(length(trials_oi));
    fprintf('%.2f %% of trials rejected! \n',nr_reject*100)

    %% weighted average for FFR + noise
    % FFR and noise
    fs = data.fsample;
    [data_w,data_w_trial,noise_w] = weighted_average_EFR(data_cc,tidx,tidx_FFR,tidx_noisef);

    data_filt = data_w;%(:,tidx_FFR); %weighted data 
    data_trial = data_w_trial;  %weighted -0.1 to 0.5 s 
    noise_filt = noise_w;       %weighted noise

    %%
    % ------------------------ Analysis ------------------------
    foi = [120]; % AM frequency

    %get fft EFR
    [f,fft_sub,f_fft_noise,FFR,F,SNR,F_crit]=get_fft_efr(data_filt,foi,fs);

    % selected channels
    chaoi_avg = [find(strcmp(data.label,'Cz')),...
        find(strcmp(data.label,'FCz')),...
        find(strcmp(data.label,'Fz'))];

    % get channel average SNR/EFR over Cz, Fz and FCz
    [f_fft_avg,FFR_avg,F_avg,SNR_avg,F_crit_avg,sig_idx_avg,noise_avg]=get_fft_efr_chaoi(f,fft_sub,chaoi_avg,foi);
    
    %get fft noise
    [f_noise,fft_sub_noise,~,~,~,~,~]=get_fft_efr(noise_filt(:,1:length(data_filt)),foi,fs);

    %% Save FFR
    data_ffr = struct;
    data_ffr.subid = data.subid;
    data_ffr.fs         = fs;
    data_ffr.f_fft      = fft_sub;
    data_ffr.fft_freq   = f';
    data_ffr.FFR        = FFR';
    data_ffr.FFR_SNR    = SNR';
    data_ffr.FFR_TS     = data_trial;
    data_ffr.time       = data.time{1}';
    data_ffr.F          = F';
    data_ffr.F_crit     = F_crit';
    data_ffr.tidx       = tidx';
    data_ffr.fft_noise  = f_fft_noise;
    data_ffr.noise.f    = f_noise;
    data_ffr.noise.fft_sub_noise = fft_sub_noise;

    % channel average
    data_ffr.f_fft_avg  = f_fft_avg';
    data_ffr.FFR_avg    = FFR_avg';
    data_ffr.FFR_SNR_avg= SNR_avg';
    data_ffr.F_avg      = F_avg';
    data_ffr.F_crit_avg = F_crit_avg;
    data_ffr.sig_idx_avg= sig_idx_avg;
    data_ffr.noise_avg  = noise_avg;

    % parameters
    data_ffr.nr_reject = nr_reject*100;
    data_ffr.channels = chaoi';
    data_ffr.chan_labels = data.label;
    data_ffr.subinfo = data.subinfo;

catch ME
    warning(ME.message)
    data_ffr = struct;
    data_ffr.error = ['Error in FFR_analysis: ' ME.message];
    data_ffr.subid = data.subid;
    data_ffr.subinfo = data.subinfo;

end
clc
disp(['EFR data processed for subject ' data.subid])
end

