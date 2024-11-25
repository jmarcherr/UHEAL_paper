
%% FFR analysis script
function data_ffr = FFR_analysis_trials2(data)
%% analysis script for preprocessed ffr data (FFR_preproc.mat)
% Outputs weighted average FFR as well as FFR_SNR and F statistics
try

    % check data
    if nargin < 1
        error('No iput. Please provide data from FFR_preproc');
    elseif  ~isfield(data,'hdr')
        error(['No processed FFR data for subject ' data.subid]);
    elseif isfield(data,'error')
        data.error_preproc = data.error;
        error([data.error])
    end

    %conditions:
    %1) 326 @ 4Hz
    %2 or 3) 326 @ 4Hz inv

    % Get scalp EEG channels
    chaoi = [1:16];
    % Get trial ids
    cco = unique(data.trialinfo);
    % find positive and negative polarity trials
    trials_oi =find(data.trialinfo==cco(1));
    trials_inv = find(data.trialinfo==cco(2)); %inverted
    % are there differences in trial nr?
    if length(trials_oi)~=length(trials_inv)
        [a,~]=min([length(trials_oi) length(trials_inv)]);
        trials_oi = trials_oi(1:a);
        trials_inv = trials_inv(1:a);
    end
    % select data
    cfg = [];
    cfg.channel = {data.label{[chaoi]}};
    data_cond = ft_selectdata(cfg,data);
    time = data.time{1};

    % define time idx
    tidx = find(time>=0 & time<0.5); % pure tone (250 ms) + noise (250 ms)
    tidx_FFR = find(time>=-0.1 & time<0.5);
    tidx_noisef = find(time>=0.25 & time<0.5);

    % epoch data
    epoched_data =epoch_data_3(data_cond,1,trials_oi);   %epoch x chan x time
    epoched_inv = -epoch_data_3(data_cond,1,trials_inv); %inverted

    %epoched_data_p = [epoched_data-epoched_inv]/2;       %phase independent
    epoched_data = [epoched_data+epoched_inv]/2;         %phase dependent

    % artifact rejection and weighting
    [valid_trials,rjt] = threshold_rejection(epoched_data,50,tidx,trials_oi);

    % clean epoched data
    data_cc = [];
    % reject if less than half the trials remain
    %if length(valid_trials)<[486/4]
    %    error(['Less than half of trials are clean. Rejecting ' data.subid]);
    %else
        data_cc = epoched_data(valid_trials,:,:); % clean data
        data_cc = permute(data_cc,[2,3,1]);
    %end
    data.valid_trials = valid_trials;

    clc
    nr_reject = length(rjt)/(length(trials_oi));
    fprintf('%.2f %% of trials rejected! \n',nr_reject*100)

    %% weighted average for FFR + noise
    % FFR and noise
    fs = data.fsample;
    [data_w,data_w_trial,noise_w] = weighted_average_FFR(data_cc,tidx_FFR,tidx_noisef);

    data_filt = data_w(:,tidx); %weighted data 
    data_trial = data_w_trial;  %weighted -0.1 to 0.5 s 
    noise_filt = noise_w;       %weighted noise

    %%
    % ------------------------ Analysis ------------------------
    foi = [326]; % pure tone frequency

    %get fft FFR
    [f,fft_sub,f_fft_noise,FFR,F,SNR,F_crit]=get_fft(data_filt,foi,fs);

    % selected channels
    chaoi_avg = [find(strcmp(data.label,'Cz')),...
        find(strcmp(data.label,'FCz')),...
        find(strcmp(data.label,'Fz'))];

    % get channel average SNR/FFR over Cz, Fz and FCz
    [f_fft_avg,FFR_avg,F_avg,SNR_avg,F_crit_avg,sig_idx_avg,noise_avg]=get_fft_chaoi(f,fft_sub,chaoi_avg,foi);
    %get fft noise
    [f_noise,fft_sub_noise,~,~,~,~,~]=get_fft(noise_filt,foi,fs);


    %% trials 1 through 6 epoch, weighted, reject 
    % Get trial ids
    tmp = data.cfg.previous.previous.trl;
    [t_info,t_shift]=check_trial(tmp);
    if t_shift
       error(['missing trials!! checkup ' data.subid]);
    end
    cco = unique(data.trialinfo);
    for tt=1:6
    % find positive and negative polarity trials
    trials_oi =find(t_info(:,1)==cco(1) & t_info(:,2)==tt);
    trials_inv = find(t_info(:,1)==cco(2) & t_info(:,2)==tt); %inverted
    % are there differences in trial nr?
    if length(trials_oi)~=length(trials_inv)
        [a,~]=min([length(trials_oi) length(trials_inv)]);
        trials_oi = trials_oi(1:a);
        trials_inv = trials_inv(1:a);
    end
    % select data
    cfg = [];
    cfg.channel = {data.label{[chaoi]}};
    data_cond = ft_selectdata(cfg,data);
    time = data.time{1};

    % define time idx
    tidx = find(time>=0 & time<0.5); % pure tone (250 ms) + noise (250 ms)
    tidx_FFR = find(time>=-0.1 & time<0.5);
    tidx_noisef = find(time>=0.25 & time<0.5);

    % epoch data
    epoched_data =epoch_data_3(data_cond,1,trials_oi);   %epoch x chan x time
    epoched_inv = -epoch_data_3(data_cond,1,trials_inv); %inverted

    %epoched_data_p = [epoched_data-epoched_inv]/2;       %phase independent
    epoched_data = [epoched_data+epoched_inv]/2;         %phase dependent

    % artifact rejection and weighting
    [valid_trials,rjt] = threshold_rejection(epoched_data,50,tidx,trials_oi);

    % clean epoched data
    data_cc = [];

    data_cc = epoched_data(valid_trials,:,:); % clean data
    data_cc = permute(data_cc,[2,3,1]);

    %data.valid_trials = valid_trials;
    %% weighted average for FFR + noise
    % FFR and noise
    fs = data.fsample;
    [data_w,data_w_trial,noise_w] = weighted_average_FFR(data_cc,tidx_FFR,tidx_noisef);

    data_filt2{tt} = data_w(:,tidx); %weighted data
    data_trial2{tt} = data_w_trial;  %weighted -0.1 to 0.5 s
    noise_filt2{tt} = noise_w;       %weighted noise
    valid_trial2{tt} = valid_trials;
    end
    %% Save FFR
    data_ffr = struct;
    data_ffr.subid = data.subid;
    data_ffr.fs         = fs;
    data_ffr.f_fft      = fft_sub;
    data_ffr.fft_freq   = f';
    data_ffr.FFR        = FFR';
    data_ffr.FFR_SNR    = SNR';
    data_ffr.FFR_TS     = data_trial;
    data_ffr.FFR_trials = data_trial2;
    %data_ffr.FFR_trials = data_trials_tt;  % uncomment for faster execution
    %data_ffr.trial_step = steps;
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
    data_ffr.stimear = data.stimear;
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
disp(['FFR data processed for subject ' data.subid])
end

