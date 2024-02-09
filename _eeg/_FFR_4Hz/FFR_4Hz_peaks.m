%% plot FFR_4Hz results
% plot FFR_4Hz and extract peaks
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work1/jonmarc/UHEAL_master/UHEAL_paper/UHEAL_startup.m')
subs = dir('_outputs/_derivatives/*.mat')
load('/work1/jonmarc/UHEAL_master/UHEAL_paper/_stats/uheal_data.mat');
%% get data
for s=1:length(subs)
    
    load([subs(s).folder filesep subs(s).name])
    clc
    disp(['sub ' subs(s).name(1:5) ' loaded...'])
    sub_num(s) = str2num(subs(s).name(3:5));
    chansoi = setdiff(1:16,[5 11]);
    % get FFR
    if isfield(data,'TS')
        itpc(s,:,:) = data.itpc;
        f = data.f;
        TS_sub(s,:) = nanmean(data.TS(chansoi,:));
        TS_sub_chan(s,:,:)=data.TS;

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
    end

    
end
 
%% 
mean(nr_reject)
std(nr_reject)


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

%% get age groups
% groups
YNH_idx = find(age<=25 & ~CP );
MNH_idx = find(age>25 & age<50 & ~CP )
ONH_idx = find(age>=50 & ~CP);
nh_idx = find(~CP); % all normal hearing
ages = [17 77];
% channels
chansoi  = setdiff(1:16,[5 11]); % all channels but T7 and T8

% plot mean over all
close all
plot(time_TS,nanmean(TS_base(nh_idx,:)))
xlim([-0.5 3.5])
%% find P1,N1,P2,N2
% P1 = 0.045 -  0.065 s
% N1 = 0.085 -  0.15 s
% P2 = 0.15  -  0.25 s
% N2 = 0.2   -  0.5 s

% onset peak
for ii=1:6 % 6 tones
    for ss=1:size(TS_base,1)
        P1_idx(ii,:) =[0+0.5*(ii-1) 0.085+0.5*(ii-1)];
        P1(ss,ii) = max(TS_base(ss,find(time_TS>=P1_idx(ii,1) & time_TS<=P1_idx(ii,2))));

        N1_idx(ii,:) =[0.065+0.5*(ii-1)  0.15+0.5*(ii-1)];
        N1(ss,ii) = min(TS_base(ss,find(time_TS>=N1_idx(ii,1) & time_TS<=N1_idx(ii,2))));

        P2_idx(ii,:) =[0.15+0.5*(ii-1) 0.25+0.5*(ii-1)];
        P2(ss,ii) = max(TS_base(ss,find(time_TS>=P2_idx(ii,1) & time_TS<=P2_idx(ii,2))));

        N2_idx(ii,:) = [0.2+0.5*(ii-1) 0.5+0.5*(ii-1)];
        N2(ss,ii) = min(TS_base(ss,find(time_TS>=N2_idx(ii,1) & time_TS<= N2_idx(ii,2))));
    end
end
close all
subplot(1,2,1)
plot(nanmean(P1,1),'linewidth',2)
hold on
plot(nanmean(N1,1),'linewidth',2)
plot(nanmean(P2,1),'linewidth',2)
plot(nanmean(N2,1),'linewidth',2)
xlim([0 7])
hleg = legend('P1','N1','P2','N2');
hleg.Box = 'off';
hleg.Position = [0.1709 0.4237 0.1416 0.2393];
xlabel('tone nr.')
ylabel('\muV')
set(gca,'xtick',[1:6])
box off

