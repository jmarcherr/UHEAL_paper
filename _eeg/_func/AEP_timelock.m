function [time,aep_avg_filt,aep_avg,n1,n1_mean,n1_lat,p2,p2_mean,p2_lat,p1,p1_mean,p1_lat,n2,n2_mean,n2_lat,p2_fix,p2_lat_fix,n1_fix,n1_lat_fix] = AEP_timelock(vts,data_cond,ids)
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
    % P1 = 0.045 -  0.065 s
    % N1 = 0.085 -  0.15 s
    % P2 = 0.15  -  0.25 s
    % N2 = 0.2   -  0.5 s

    % P1 interval
    timeidx_P1 = find(timelock.time>=0.045 & timelock.time<=0.085);
    timeidx_P1_mean = find(timelock.time>=0.045 & timelock.time<=.065);
    % N1 interval
    timeidx_N1 = find(timelock.time>=.085 & timelock.time<=.15);
    timeidx_N1_mean = find(timelock.time>=.1 & timelock.time<.12);
    % P2 interval
    timeidx_P2 = find(timelock.time>=0.19 & timelock.time<=0.3);
    timeidx_P2_mean = find(timelock.time>=0.19 & timelock.time<0.21);
    % N2 interval
    timeidx_N2 = find(timelock.time>=0.2 & timelock.time<=0.5);
    timeidx_N2_mean = find(timelock.time>=0.25 & timelock.time<0.35);
    % get average AEP and N100
    for kk=4:-1:1
        % Get conditional average (isi)
        aep_avg(kk,:) = squeeze(mean(mean(timelock.trial(find(timelock.trialids==kk),:,:),1),2));
        % baseline correction
        baseline = mean(aep_avg(kk,find(timelock.time>-0.1 & timelock.time<=0.01)));%
        aep_avg(kk,:) = aep_avg(kk,:)-baseline;
        % get p50
        [p1(kk),p1i] = max(aep_avg(kk,timeidx_P1_mean))
        % mean p50
        [p1_mean(kk)]=mean(aep_avg(kk,timeidx_P1));
        % get n100
        [n1(kk),ni(kk)] = min(aep_avg(kk,timeidx_N1));
        % mean n100
        [n1_mean(kk)] = mean(aep_avg(kk,timeidx_N1_mean));
        % get P200
        [p2(kk),p2i(kk)] = max(aep_avg(kk,timeidx_P2));
        % P200 mean
        [p2_mean(kk)] = mean(aep_avg(kk,timeidx_P2_mean));
        % get N200
        [n2(kk),n2i] = min(aep_avg(kk,timeidx_N2));
        n2_mean(kk) = mean(aep_avg(kk,timeidx_N2_mean));

        % latencies
        p1_lat(kk) = timelock.time(timeidx_P1(p1i));
        n1_lat(kk) = timelock.time(timeidx_N1(ni(kk)));
        p2_lat(kk) = timelock.time(timeidx_P2(p2i(kk)));
        n2_lat(kk) = timelock.time(timeidx_N2(n2i));


 
    end

    % get fixed latency p2/n1 for isi 1 based on ERPs from isi 2:4
    % kk4 = isi 2.3s
    for kk=1:4
        p2_fix(kk) = aep_avg(kk,timeidx_P2(round(mean(p2i(2:4)))));
        p2_lat_fix(kk) = timelock.time(timeidx_P2(round(mean(p2i(2:4)))));
        n1_fix(kk) = aep_avg(kk,timeidx_N1(round(mean(ni(2:4)))));
        n1_lat_fix(kk) = timelock.time(timeidx_N1(round(mean(ni(2:4)))));
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