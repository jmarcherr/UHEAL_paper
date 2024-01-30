
%% FFR analysis script
function data_ffr = FFR_4Hz_analysis(data)
%% analysis script for preprocessed ffr_4Hz data (FFR_4Hz_preproc.mat)
% Outputs weighted average FFR_4Hz and ITPC spectrum
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
    tidx = find(time>=0 & time<3.5); % 6 consequtive 250 ms pure-tones

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
    fs = data.fsample;
    [data_w] = weighted_average_FFR_4Hz(data_cc,tidx);

    %%
    % ------------------------ Analysis ------------------------
    % itpc spectrum
    N=size(data_cc,3); % number of valid trials
    tidx_itpc = time>=0 & time<3.5;
    for cc=1:size(data_cc,1) % channels
        for it = 1:N % trials

            M=data_cc(cc,tidx_itpc,it);
            %FFT
            f_fft = fft(M)/(length(M)/2);
            %Convert to power
            pow(it,:) = abs(f_fft.^2); %
            %Truncate negative freqencies
            ft_sub(it,:) = (pow(it,1:end/2+1));
            %Frequency vector
            f = fs/2*linspace(0,1,length(ft_sub(it,:)));
            % get complex value at signal bin
            itpc(it,:) = f_fft(1:end/2+1);
        end

        itpc = itpc./abs(itpc);
        itpc      = sum(itpc);   % sum angles
        itpc      = abs(itpc)/N;   % take the absolute value and normalize
        itpc_sub(cc,:)      = squeeze(itpc);
    end

%     % itpc spectogram
%     foihz = 1:0.1:50; % frequencies to process
%     for ff=1:length(foihz)
%         cfg = [];
%         cfg.channel    = 1:16;
%         cfg.method     = 'wavelet';
%         cfg.width      = 12;
%         cfg.gwidth     = 3;
%         cfg.pad        = 'maxperlen';
%         cfg.output     = 'fourier';
%         cfg.keeptrials = 'yes'
%         cfg.foi        = foihz(ff);
%         cfg.toi        = -5:2/fs:7;
%         TFR = ft_freqanalysis(cfg, data);    %
%         %
%         itc           = [];
%         itc.label     = TFR.label;
%         itc.freq      = TFR.freq;
%         itc.time      = TFR.time;
%         itc.dimord    = 'chan_freq_time';
% 
% 
%         F = squeeze(TFR.fourierspctrm);   % copy the Fourier spectrum
%         N = size(F,1);           % number of trials
% 
%         % compute inter-trial phase coherence (itpc)
%         itc.itpc      = F./abs(F);         % divide by amplitude
%         itc.itpc      = nansum(itc.itpc,1);   % sum angles
%         itc.itpc      = abs(itc.itpc)/N;   % take the absolute value and normalize
%         itc.itpc      = squeeze(itc.itpc); % remove the first singleton dimension
% 
%         itpc_spec{ff} = itc;
%     end
    %% Save FFR
    data_ffr = struct;

    data_ffr.TS        = data_w;
    data_ffr.time      = time;
    data_ffr.tidx      = tidx;
    data_ffr.fs        = fs;
    %data_ffr.itpc_spec = itpc_spec; % itpc spectrogram
    data_ffr.itpc      = itpc_sub;  % itpc over frequency
    data_ffr.tidx_itpc = time(tidx_itpc);
    data_ffr.f          =f;
    
    data_ffr.subid = data.subid; 
    data_ffr.nr_reject = nr_reject*100;
    data_ffr.valid_trials = valid_trials;
    data_ffr.channels = chaoi;
    data_ffr.chan_labels = chan_labels;
    data_ffr.subinfo = data.subinfo;

catch ME
    warning(ME.message)
    data_ffr = struct;
    data_ffr.error = ['Error in FFR_4Hz_analysis: ' ME.message];
    data_ffr.subid = data.subid;
    data_ffr.subinfo = data.subinfo;

end
clc
disp(['FFR_4Hz data processed for subject ' data.subid])
end

