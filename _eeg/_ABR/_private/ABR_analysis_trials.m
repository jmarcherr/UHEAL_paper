
%% ABR analysis script
function data_abr = ABR_analysis_chans(data)
%% analysis script for preprocessed abr data (ABR_preproc.mat)
% Outputs baseline corrected, weighted average ABR as well as weighted
% trials
try

    % reject noisy subjects (19 out of 117 subjects in total)
    rjt_sub = {}%{'UH003','UH014','UH027','UH031','UH042','UH053','UH056','UH069',...
        %'UH070','UH077','UH086','UH088','UH089','UH090','UH108','UH109',...
       % 'UH110','UH111','UH112'};

    % check data
    if nargin < 1
        error('No iput. Please provide data from ABR_preproc');
    elseif  ~isfield(data,'hdr')
        error(['No processed ABR data for subject ' data.subid]);
    elseif sum(strcmp(data.subid,rjt_sub))
        error(['Noisy data. No processed ABR data for subject ' data.subid])
    end

    %find reference chan
    chans = 1:16;
    %% ------------Comments --------------------------------------
    %Trial ids:
    %1 = click, 80 nHL, 9.1Hz
    %2 or 3 = click, 80 nHL, 40Hz
    trial_info = unique(data.trialinfo);

    % loop over rates (9/s and 40/s)
    rate_ids = [9,40];
    rate = 1:2;
    for kk=[rate] % condition loop

        trials_oi =find(data.trialinfo==trial_info(kk)); % condition trials

        % select data
        cfg = [];
        cfg.channel = {data.label{[chans]}};
        data_cond = ft_selectdata(cfg,data); % select data
        time = data.time{1};
        % define epoch window
        tidx = find(time>=-5e-3 & time<15e-3);

        % epoch data
        epoched_data =epoch_data_3(data_cond,1,trials_oi); %epoch x chan x time

        % artifact rejection and weighting
        [valid_trials,rjt] = threshold_rejection(epoched_data,20,tidx,trials_oi);

        % clean epoched data
        data_cc = [];
        % reject if less than half the trials remain
        if length(valid_trials)<3000
            warning(['Less than half of trials are clean. Rejecting ' rateids(kk) '/s data from subject' data.subid])
            data_cc = nan;
        else
            data_cc = epoched_data(valid_trials,:,:);
            data_cc = permute(data_cc,[2,3,1]);
        end
        data.valid_trials = valid_trials;

        clc
        nr_reject(kk) = length(rjt)/(length(trials_oi));
        fprintf('%.2f %% of trials rejected! \n',nr_reject*100)
        %% weighted average
        [data_weighted{kk}] = weighted_ABR_chans(data,data_cc,tidx);
        %% delay, baseline correction
        [data_weighted{kk},time_corrected{kk},tidx_corrected{kk}] = ABR_correction_chans(data_weighted{kk},time,tidx,data.fsample,data.subid);

        %% find peaks
        [abr_peaks{kk}] = get_abr_peaks_chan(data_weighted{kk},time_corrected{kk});

        %% mcca trials
        %% loop over % trials

        % from 100% to 10% in steps of
        steps = [10:10:100];
        perm_trial = randperm(size(data_cc,3));
        %tidx_FFR = find(time>=-0.1 & data.time<0.5);
        for tt = 1:length(steps)
            nr_trials = round(steps(tt)/100*size(data_cc,3));
            %select random trials
            this_trial = perm_trial(1:nr_trials);
            this_dat = data_cc(:,:,this_trial);
            [data_w_trials] = weighted_ABR_chans_mcca(data,this_dat,tidx);
            [data_w_trials,~,~] = ABR_correction_chans(data_w_trials,time,tidx,data.fsample,data.subid);
            data_trials_tt(tt,:,:) = data_w_trials;
        end
        data_trials_weighted{kk} = data_trials_tt;
    end
    %% save processed ABR
    data_abr = struct;
    data_abr.subid = data.subid;
    data_abr.subinfo = data.subinfo;
    data_abr.stimear = data.stimear;
    data_abr.abr = data_weighted;
    data_abr.trials = data_trials_weighted;
    data_abr.time = time_corrected;
    data_abr.tidx = tidx_corrected;
    data_abr.fs = data.fsample;
    data_abr.abr_peaks = abr_peaks;
    data_abr.nr_reject = nr_reject*100;

catch ME
    warning(ME.message)
    data_abr = struct;
    data_abr.error = ['Error in ABR_analysis: ' ME.message];
    data_abr.subid = data.subid;
    data_abr.subinfo = data.subinfo;

end
clc
disp(['ABR data processed for subject ' data.subid])
end

