function [timelock,aep_avg,n1,n1_lat,p2,p2_lat] = AEP_timelock(vts,data_cond,ids)
%extract average AEP per ISI and get N100 and P200 peak amplitudes and
%latencies
    % timelock analysis
    subtrials = vts;
    % ids
    [fid,~] = sort(cell2mat(subtrials));
    cfg = [];
    cfg.trials = fid;
    cfg.keeptrials  = 'yes';
    cfg.lpfilttype  = 'firws';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 30;

    data_cond = ft_preprocessing(cfg,data_cond);
    cfg = [];
    cfg.keeptrials  = 'yes';
    timelock = ft_timelockanalysis(cfg,data_cond)
    timelock.trialids = ids(fid);
    % N1 interval
    timeidx_N1 = find(timelock.time>=.08 & timelock.time<.15);
    % P2 interval
    timeidx_P2 = find(timelock.time>=0.16 & timelock.time<0.3)

    % get average AEP and N100
    for kk=1:4
        % Get conditional average (isi)
        aep_avg(kk,:) = squeeze(mean(mean(timelock.trial(find(timelock.trialids==kk),:,:),1),2));
        % get n100
        [n1(kk),ni] = min(aep_avg(kk,timeidx_N1));
        % get P200
        [p2(kk),p2i] = max(aep_avg(kk,timeidx_P2))
        n1_lat(kk) = time(timeidx_N1(ni));
        p2_lat(kk) = time(timeidx_P2(p2i));
    end
end