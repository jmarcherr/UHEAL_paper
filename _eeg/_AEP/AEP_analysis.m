
%% ABR analysis script
function data_aep = AEP_analysis(data)
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

    % get ISI conditions
    if ~isempty(data.missing_trials)
        ids = data.stimfile.id(setdiff(1:length(data.stimfile.id),data.missing_trials));
    else
        ids = data.stimfile.id;
    end
    isi = data.stimfile.isi+0.3; % adjusted ISI

    %Get channels of interest
    if strcmp(data.subid,'UH020') % noisy FC4(16 channel)
        chaoi = setdiff(1:15,[5 11]);
    else
        chaoi = setdiff(1:16,[5 11]);
    end

    % Loop over isis
    for kk=[1:4] % condition loop

        trials_oi =find(ids==kk);

        % select data
        cfg = [];
        cfg.channel = {data.label{[chaoi]}};
        data_cond = ft_selectdata(cfg,data);
        time = data.time{1};

        % define analysis interval
        tidx = find(time>=0 & time<=.5);

        % epoch data
        epoched_data =epoch_data_3(data_cond,1,trials_oi); %epoch x chan x time

        % artifact rejection and weighting
        [valid_trials,rjt] = threshold_rejection(epoched_data,80,tidx,trials_oi);

        data_cc = epoched_data(valid_trials,:,:);
        data_cc = permute(data_cc,[2,3,1]);

        clc
        nr_reject(kk) = length(rjt)/(length(trials_oi));
        fprintf('%.2f %% of trials rejected! \n',nr_reject(kk)*100)
        % valid trials per ISI
        vts{kk} = trials_oi(valid_trials);
    end

    % timelock analysis
    [time,aep_avg_filt,aep_avg,n1,n1_mean,n1_lat,p2,p2_mean,p2_lat,p1,p1_mean,p1_lat] = AEP_timelock(vts,data_cond,ids);

    %% save processed AEP
    data_aep = struct;
    data_aep.fs = data.fsample;
    data_aep.subid = data.subid;
    data_aep.nr_reject = nr_reject*100;
    data_aep.time = time;
    data_aep.aep_avg = aep_avg;
    data_aep.aep_avg_filt = aep_avg_filt;
    data_aep.p1 = p1;
    data_aep.p1_mean = p1_mean;
    data_aep.n1 = n1;
    data_aep.n1_mean = n1_mean;
    data_aep.p2 = p2;
    data_aep.p2_mean = p2_mean;
    data_aep.p1_lat = p1_lat;
    data_aep.n1_lat = n1_lat;
    data_aep.p2_lat = p2_lat;
    data_aep.subinfo = data.subinfo;
    data_aep.subid = data.subid;

catch ME
    warning(ME.message)
    data_aep = struct;
    data_aep.error = ['Error in AEP_analysis: ' ME.message];
    data_aep.subid = data.subid;
    data_aep.subinfo = data.subinfo;

end
clc
disp(['AEP data processed for subject ' data.subid])
end

