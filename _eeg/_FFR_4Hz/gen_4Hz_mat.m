function [uheal_data]=gen_4Hz_mat(uheal_data,datadir)

d = dir([datadir filesep '*.mat'])
clc
disp(['Processing 4Hz data ...'])
%% get data
for s=1:length(d)

    load([d(s).folder filesep d(s).name])

    chansoi = setdiff(1:16,[5 11]);
    % get FFR
    if isfield(data,'itpc')
        itpc(s,:,:) = data.itpc;
        f = data.f;
        TS_sub(s,:) = nanmean(data.TS(chansoi,:));
        TS_sub_chan(s,:,:)=data.TS;

        time = data.time;
        tidx = data.tidx;
        tidx_TS = data.tidx_TS;

        subinfo{s} = data.subinfo;
        if isempty(data.subinfo.age)|isempty(data.subinfo.gender)
            age(s) = nan;
            gender(s) = nan;
        else
            age(s) =data.subinfo.age;
            gender(s) = data.subinfo.gender;
        end

        nr_reject(s) =data.nr_reject;
        chan_labels{s} = data.chan_labels;
        chans{s} = data.channels;

    else

        itpc(s,:,:) =nan(16,1537);
        TS_sub(s,:,:) = nan(1,5632);
        TS_sub_chan(s,:,:) = nan(16,5632);
        subinfo{s} = data.subinfo;
        age(s) = data.subinfo.age;
        gender(s) = data.subinfo.gender;
    end
    subid{s} = data.subid;

end
fid = [2:2:20];
% mean spectrum
chansoi = setdiff(1:16,[5, 11]);
itpc_spec_mean=squeeze((nanmean(itpc(:,chansoi,:),2)));
% find harmonic frequencies
for ii=1:length(fid)
    fid_idx(ii) = find(f==fid(ii));
    itpc_freq(ii,:) = itpc_spec_mean(:,find(f==fid(ii)))';
end

% Get ratio
itpc_F0 = itpc_freq(1,:);
itpc_F19 = mean(itpc_freq(2:10,:));
itpc_ratio = log10(itpc_freq(1,:)./mean(itpc_freq(2:10,:)));


fs = 1024;
filt_coef = [.5 30]; %4 Hz
filt_def = designfilt('bandpassfir', 'FilterOrder', 40, ...
    'CutoffFrequency1', filt_coef(1), 'CutoffFrequency2', filt_coef(2),...
    'SampleRate', fs);
% bypassing filtfilt fieldtrip
addpath('/appl/matlab/990/toolbox/signal/signal/')


%%%%%%%%%%%%%%%%%%%%%%%
% baseline norm
time_TS = time(tidx_TS);
baseline(:) = mean(TS_sub(:,find(time_TS>-0.1 & time_TS<=0)),2);
TS_base = TS_sub-baseline';

%%% mean amplitude (negativity)
mean_amp = nanmean(TS_base(:,find(time_TS>=0 & time_TS<=3)),2);


% extract
uheal_data.Neg_4Hz = nan(size(uheal_data.subid));
uheal_data.ITPC_ratio = nan(size(uheal_data.subid));


for s=1:length(mean_amp)
    % get this subid
    thisID = str2double(subid{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);
    
    uheal_data.Neg_4Hz(this_idx) = mean_amp(s)';
    uheal_data.ITPC_ratio(this_idx) = itpc_ratio(s)';


end
clc
disp(['4Hz data done ...'])
end
 
