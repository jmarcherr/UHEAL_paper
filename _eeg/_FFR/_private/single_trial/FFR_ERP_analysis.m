function data_ffr = FFR_ERP_analysis(data)
% extracts ERP for single trials of the FFR paradigm
try

    % check data
    if nargin < 1
        error('No iput. Please provide data from FFR_4Hz_preproc');
    elseif  ~isfield(data,'hdr')
        error(['No processed FFR_4Hz data for subject ' data.subid]);
    elseif isfield(data,'error')
        data.error_preproc = data.error;
        error([data.error])
    end

% Get scalp EEG channels
    chaoi = [1:16];
    chan_labels = data.label(chaoi);
    % Get trial ids
    trials_oi = find(data.trialinfo);

    % select data
    cfg = [];
    cfg.channel = {data.label{[chaoi]}};
    data_cond = ft_selectdata(cfg,data);
    time = data.time{1};

    % define time idx
    tidx = find(time>=-.1 & time<0.5); % 6 consequtive 250 ms pure-tones

    % epoch data
    epoched_data =epoch_data_3(data_cond,1,trials_oi);   %epoch x chan x time


    % artifact rejection and weighting
    [valid_trials,rjt] = threshold_rejection(epoched_data,150,tidx,trials_oi);

    % clean epoched data
    data_cc = [];
    % reject if less than half the trials remain
    if length(valid_trials)<[length(trials_oi)/2]
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
    tidx_w = time>=0 & time<0.5; % +/- 1s 
    fs = data.fsample;
    [data_w] = weighted_average_FFR_4Hz(data_cc,tidx_w);

        %% Save FFR
    data_ffr = struct;

    data_ffr.TS        = data_w;
    data_ffr.TS_trials = data_cc(:,tidx_w,:);
    data_ffr.time      = time;
    data_ffr.tidx      = tidx;
    data_ffr.tidx_TS   = tidx_w;
    data_ffr.fs        = fs;
    
    data_ffr.subid = data.subid; 
    data_ffr.nr_reject = nr_reject*100;
    data_ffr.valid_trials = valid_trials;
    data_ffr.channels = chaoi;
    data_ffr.chan_labels = chan_labels;
    data_ffr.subinfo = data.subinfo;

catch ME
    warning(ME.message)
    data_ffr = struct;
    data_ffr.error = ['Error in FFR_ERP_analysis: ' ME.message];
    data_ffr.subid = data.subid;
    data_ffr.subinfo = data.subinfo;

end
clc
disp(['FFR_ERP data processed for subject ' data.subid])
end