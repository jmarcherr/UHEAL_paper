function [time,aep_avg_filt,aep_avg,n1,n1_mean,n1_lat,p2,p2_mean,p2_lat] = AEP_timelock(vts,data_cond,ids)
%extract average AEP per ISI and get N100 and P200 peak amplitudes and
%latencies
    % timelock analysis
    subtrials = vts;
    % ids
    [fid,~] = sort(cell2mat(subtrials));
    cfg = [];
    cfg.trials = fid;
    cfg.keeptrials  = 'yes';

    data_cond_nofilt = ft_preprocessing(cfg,data_cond);
    cfg = [];
    cfg.keeptrials  = 'yes';
    timelock = ft_timelockanalysis(cfg,data_cond_nofilt)
    timelock.trialids = ids(fid);
    % N1 interval
    timeidx_N1 = find(timelock.time>=.08 & timelock.time<=.15);
    timeidx_N1_mean = find(timelock.time>=.1 & timelock.time<.12);
    % P2 interval
    timeidx_P2 = find(timelock.time>=0.19 & timelock.time<=0.3);
    timeidx_P2_mean = find(timelock.time>=0.19 & timelock.time<0.21);

    % get average AEP and N100
    for kk=1:4
        % Get conditional average (isi)
        aep_avg(kk,:) = squeeze(mean(mean(timelock.trial(find(timelock.trialids==kk),:,:),1),2));
        % baseline correction
        baseline = mean(aep_avg(kk,find(timelock.time>-0.1 & timelock.time<=0.01)));%
        aep_avg(kk,:) = aep_avg(kk,:)-baseline;
        % get n100
        [n1(kk),ni] = min(aep_avg(kk,timeidx_N1));
        % mean n100
        [n1_mean(kk)] = mean(aep_avg(kk,timeidx_N1_mean));
        % get P200
        [p2(kk),p2i] = max(aep_avg(kk,timeidx_P2))
        % P200 mean
        [p2_mean(kk)] = mean(aep_avg(kk,timeidx_P2_mean))
        % latencies
        n1_lat(kk) = timelock.time(timeidx_N1(ni));
        p2_lat(kk) = timelock.time(timeidx_P2(p2i));
    end

    % filtered version for plotting
    % ids
    [fid,~] = sort(cell2mat(subtrials));
    cfg = [];
    cfg.trials = fid;
    cfg.keeptrials  = 'yes';
    cfg.lpfilttype  = 'firws';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 30;

    data_cond_filt = ft_preprocessing(cfg,data_cond);
    cfg = [];
    cfg.keeptrials  = 'yes';
    timelock_filt = ft_timelockanalysis(cfg,data_cond_filt);
    timelock_filt.trialids = ids(fid);
    for kk=1:4
        % Get conditional average (isi)
        aep_avg_filt(kk,:) = squeeze(mean(mean(timelock_filt.trial(find(timelock_filt.trialids==kk),:,:),1),2));
        % baseline correction
        %baseline = mean(aep_avg_filt(kk,find(timelock_filt.time>-0.1 & timelock_filt.time<=0.01)));%
        %aep_avg_filt(kk,:) = aep_avg_filt(kk,:)-baseline;
    end
    time = timelock.time;
end