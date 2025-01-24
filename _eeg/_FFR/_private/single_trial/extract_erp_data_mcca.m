if isfield(data,'TS')
    % from cell to double
    fs(s) = data.fs;
    TS_sub(s,:,:)   = data.TS; % timeseries
    TS_trials{s} = data.TS_trials;
    time{s}         = data.time;
    tidx{s}         = data.tidx;
    subid{s}        = data.subid;

    subinfo{s} = data.subinfo;
    nr_reject(s) =data.nr_reject;
    chan_labels{s} = data.chan_labels;
    chans{s} = data.channels;
    age(s) = data.subinfo.age;
    gender(s) = data.subinfo.gender;
else

    fs(s) = nan;
    TS_sub(s,:,:) = nan(1,16,512);
    TS_trials{s} = nan;
    subinfo{s} = data.subinfo;
    subid{s} = data.subid;
    age(s) = data.subinfo.age;
    gender(s) = data.subinfo.gender;




end


%% dss

