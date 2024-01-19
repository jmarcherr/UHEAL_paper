
%% ABR analysis script
function data_abr = ABR_analysis(data)
%% analysis script for preprocessed abr data (ABR_preproc.mat)
% Outputs baseline corrected, weighted average ABR as well as weighted
% trials
%cd(fileparts(matlab.desktop.editor.getActiveFilename))
try
    % check data
    if nargin < 1 ; error('No iput. Please provide data from ABR_preproc'); end
    
    if ~isfield(data,'hdr')
        error(['No processed ABR data for subject ' data.subid])
    end
    
    %find reference chan
    if data.stimear == 1 % left ear stimulation
        chans = find(strcmp(data.label,'EXG1')); % left tiptrode
    elseif data.stimear ==2 % right ear stimulation
        chans = find(strcmp(data.label,'EXG2')); % right tiptrode
    end

    % not tiptrode channels. Reject subject
    if isempty(chans)
       error('catch');
    else

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
            [data_weighted{kk},data_trials{kk}] = weighted_ABR(data,data_cc,tidx);
            %% delay, baseline correction
            [data_weighted{kk},data_trials{kk},time_corrected{kk},tidx_corrected{kk}] = ABR_correction(data_weighted{kk},data_trials{kk},time,tidx,data.fsample,data.subid);

        end
        %% save processed ABR
        data_abr = struct;
        data_abr.subid = data.subid;
        data_abr.subinfo = data.subinfo;
        data_abr.stimear = data.stimear;
        data_abr.abr = data_weighted;
        data_abr.trials = data_trials;
        data_abr.time = time_corrected;
        data_abr.tidx = tidx_corrected;
        data_abr.fs = data.fsample;
        data_abr.nr_reject = nr_reject*100;


    end   
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

