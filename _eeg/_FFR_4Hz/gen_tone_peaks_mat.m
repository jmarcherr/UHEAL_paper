function [uheal_data]=gen_tone_peaks_mat(uheal_data,datadir)

d = dir([datadir filesep '*.mat'])
clc
disp(['Processing tone_peak data ...'])
%% get data
for s=1:length(d)
    
    load([d(s).folder filesep d(s).name])
    %clc

    %disp(['sub ' d(s).name(1:5) ' loaded...'])
    sub_num(s) = str2num(d(s).name(3:5));
    chansoi = setdiff(1:16,[5 11]);
    % get FFR
    if isfield(data,'TS')
        itpc(s,:,:) = data.itpc;
        f = data.f;
        TS_sub(s,:) = nanmean(data.TS(chansoi,:));
        TS_sub_chan(s,:,:)=data.TS;
        TS_trials = data.TS_trials;%chan x time x trial
        clear TS_trials dat

        time = data.time;
        tidx = data.tidx;
        tidx_TS = data.tidx_TS;

        age(s) =data.subinfo.age;
        gender(s) = data.subinfo.gender;

        CP(s) =  uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));
        nr_reject(s) =data.nr_reject;
        chan_labels{s} = data.chan_labels;
        chans{s} = data.channels;
        
    else

        itpc(s,:,:) =nan(16,1537);
        TS_sub(s,:,:) = nan(1,5632);
        TS_sub_chan(s,:,:) = nan(16,5632);
        age(s) = data.subinfo.age;
        gender(s) = data.subinfo.gender;
        CP(s) = uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));
        %dat_clean(s,:,:) = nan(size(dat_clean(1,:,:)));
        z_sub{s} = nan;
    end
    subid{s} = data.subid;
    
end
%%
% TS_sub_chan = subject x chan x time

%% Time-series (TS) baseline
%%%%%%%%%%%%%%%%%%%%%%%

time_TS = time(tidx_TS); %-1:4s
% baseline norm from -0.1 - 0 s
baseline_idx = find(time_TS>-0.1 & time_TS<=0);

% estimate baseline
baseline(:) = mean(TS_sub(:,baseline_idx),2);
baseline_chan = mean(TS_sub_chan(:,1:16,baseline_idx),3);

% baseline corrected time-series
TS_base = TS_sub-baseline'; % mean over 14 scalp channels (chanoi) subjects x time
TS_base_chan = TS_sub_chan-baseline_chan; % subjects x chan x time

%% alternative peak find
% P1 = 0.045 -  0.065 s
% N1 = 0.065 -  0.15 s
% P2 = 0.15  -  0.25 s
% N2 = 0.2   -  0.4 s
fs = 1032;
data_all = TS_base;% {TS_base;TS_dss};
% peak latencies

P1_idx =[0.045 0.065];
N1_idx =[0.065  0.15];
P2_idx =[0.15 0.25];
N2_idx = [0.2 0.4];
p_idx = {P1_idx,N1_idx,P2_idx,N2_idx};

% loop over subjects
for ss=1:size(data_all,1)
    % find average peak latency for each subject
    % find mean waveform over 6 tones
    for tt=1%:6
        this_idx = [0+0.5*(tt-1) 0.5+0.5*(tt-1)];
        this_mean(tt,:) = data_all(ss,find(time_TS>=this_idx(1) & time_TS<=this_idx(2)));
    end
    % mean over tones
    this_mean = nanmean(this_mean,1);

    mean_time = time_TS(find(time_TS>=0 & time_TS<=0.5));
    % find mean peaks
    for pp=1:4
        this_p =[p_idx{pp}(1) p_idx{pp}(2)];
        this_t = find(mean_time>=this_p(1) & mean_time<=this_p(2));
        if mod(pp,2)
            [p_fix_s(pp),p_l] = max(this_mean(this_t));
        else
            [p_fix_s(pp),p_l] = min(this_mean(this_t));
        end
        % convert to time
        p_fix_s_l(pp) = mean_time(this_t(p_l));
    end

    % loop over tones
    for ii=1:6
        for pp=1:4
            % find latency for this peak
            this_p =[p_idx{pp}(1)+0.5*(ii-1) p_idx{pp}(2)+0.5*(ii-1)];
            this_t = find(time_TS>=this_p(1) & time_TS<=this_p(2));

            % mean over latencies
            p_mean(ss,ii,pp) = mean(data_all(ss,this_t));

            % max/min
            if mod(pp,2) %negative of positive peaks
                [p_max(ss,ii,pp),p_l] = max(data_all(ss,this_t));
            else
                [p_max(ss,ii,pp),p_l] = min(data_all(ss,this_t));
            end
            % convert index to time
            p_max_l(ss,ii,pp) = time_TS(this_t(p_l));

            % fixed latency based on mean over all tones
            p_fix(ss,ii,pp) = data_all(ss,find(time_TS==p_fix_s_l(pp)+0.5*(ii-1)));
            p_fix_l(ss,ii,pp) = p_fix_s_l(pp)+0.5*(ii-1);




        end
    end
end

% outcomes:
% 1) p_mean over latencies
% 2) p_max + p_max_l -> max/min amplitude per tone
% 3) p_fix + p_fix_l -> max/min based on mean tone response (fixed latency)

% loop over strategies
strat = {p_max;p_mean;p_fix};
latency = {p_max_l;p_fix_l;p_fix_l};

%p1
p1 = p_fix(:,:,1); % fixed
p1_lat = p_max_l(:,:,1);
p1_1vsrest = p1(:,1)-nanmean(p1(:,2:6),2);
%n1
n1 = p_fix(:,:,2); % fixed
n1_lat = p_max_l(:,:,2)
n1_1vsrest = n1(:,1)-nanmean(n1(:,2:6),2);
%p2
p2 = p_fix(:,:,3); % fixed
p2_lat = p_max_l(:,:,3)
p2_1vsrest = p2(:,1)-nanmean(p2(:,2:6),2);
%n2
n2 = p_max(:,:,4); % max/min
n2_lat = p_max_l(:,:,4)
n2_1vsrest = n2(:,1)-nanmean(n2(:,2:6),2);



% extract
uheal_data.tone_p1 = nan(size(uheal_data.subid,1),6);
uheal_data.tone_p1_lat = nan(size(uheal_data.subid,1),6);
uheal_data.tone_p1_1vsrest = nan(size(uheal_data.subid));
%n1
uheal_data.tone_n1 = nan(size(uheal_data.subid,1),6);
uheal_data.tone_n1_lat = nan(size(uheal_data.subid,1),6);
uheal_data.tone_n1_1vsrest = nan(size(uheal_data.subid));
%p2
uheal_data.tone_p2 = nan(size(uheal_data.subid,1),6);
uheal_data.tone_p2_lat = nan(size(uheal_data.subid,1),6);
uheal_data.tone_p2_1vsrest = nan(size(uheal_data.subid));
%n2
uheal_data.tone_n2 = nan(size(uheal_data.subid,1),6);
uheal_data.tone_n2_lat = nan(size(uheal_data.subid,1),6);
uheal_data.tone_n2_1vsrest = nan(size(uheal_data.subid));


for s=1:length(p1_1vsrest)
    % get this subid
    thisID = str2double(subid{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);

    %p1
    uheal_data.tone_p1(this_idx,:) = p1(s,:);
    uheal_data.tone_p1_lat(this_idx,:) = p1_lat(s,:);
    uheal_data.tone_p1_1vsrest(this_idx) = p1_1vsrest(s);
    %n1
    uheal_data.tone_n1(this_idx,:) = n1(s,:);
    uheal_data.tone_n1_lat(this_idx,:) = n1_lat(s,:);
    uheal_data.tone_n1_1vsrest(this_idx) = n1_1vsrest(s);
    %p2
    uheal_data.tone_p2(this_idx,:) = p2(s,:);
    uheal_data.tone_p2_lat(this_idx,:) = p2_lat(s,:);
    uheal_data.tone_p2_1vsrest(this_idx) = p2_1vsrest(s);
    %n2
    uheal_data.tone_n2(this_idx,:) = n2(s,:);
    uheal_data.tone_n2_lat(this_idx,:) = n2_lat(s,:);
    uheal_data.tone_n2_1vsrest(this_idx) = n2_1vsrest(s);



end
clc
disp(['tone peak data done ...'])

