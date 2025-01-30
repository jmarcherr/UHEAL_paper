
%% ABR analysis script
function data_aep = AEP_analysis_NF(data)
%% analysis script for preprocessed abr data (AEP_preproc.mat)
% Outputs baseline corrected, weighted average AEP as well as weighted
% trials for each ISI condition
try

    % check data
    if nargin < 1
        error('No iput. Please provide data from AEP_preproc');
    elseif  ~isfield(data,'hdr')
        error(['No processed AEP data for subject ' data.subid]);
    end

    %Get channels of interest
        chaoi = 1:16;%setdiff(1:16,[5 11]);


        trials_oi =find(data.trialinfo); % all trials

        % select data
        cfg = [];
        cfg.channel = {data.label{[chaoi]}};
        data_cond = ft_selectdata(cfg,data);
        time = data.time{1};

        % define analysis interval
        tidx = find(time>=0 & time<=3);

        % epoch data
        epoched_data =epoch_data_3(data_cond,1,trials_oi); %epoch x chan x time

        % artifact rejection and weighting
        [valid_trials,rjt] = threshold_rejection(epoched_data,120,tidx,trials_oi);

        data_cc = epoched_data(valid_trials,:,:);
        data_cc = permute(data_cc,[2,3,1]);

        clc
        nr_reject = length(rjt)/(length(trials_oi));
        fprintf('%.2f %% of trials rejected! \n',nr_reject*100)
        % valid trials
        vts = trials_oi(valid_trials);

        % weighted average
        tidx_w = time>=0 & time<3; % +/- 1s
        fs = data.fsample;
        [data_w] = weighted_average_FFR_4Hz(data_cc,tidx_w);
        

    %% save processed AEP
    data_aep = struct;
    data_aep.fs = data.fsample;
    data_aep.subid = data.subid;
    data_aep.nr_reject = nr_reject*100;
    data_aep.time = time;
    data_aep.data_w = data_w;
    data_aep.chan_labels = data.label;
    data_aep.chanoi = chaoi;
    data_aep.subinfo = data.subinfo;
    data_aep.subid = data.subid;

catch ME
    warning(ME.message)
    data_aep = struct;
    data_aep.error = ['Error in AEP_analysis_NF: ' ME.message];
    data_aep.subid = data.subid;
    data_aep.subinfo = data.subinfo;

end
clc
disp(['AEP_NF data processed for subject ' data.subid])
end

