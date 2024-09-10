function [data_out,t_abr,tnidx] = ABR_correction_chans(data_in,time,tidx,fs,subid)
    
% init
delay = 1.1e-3; %Tube delay
fs = fs;
blIDX = round(delay/(1/fs));

% delay correction
tn = time(tidx)-delay;
%tnIDX = find(tn>-5e-3);

% AD correction ids
rjt_ids= setdiff(1:20,2);

%% correction for AD break
if isempty(setdiff(str2num(subid(3:end)),rjt_ids))
data_out = data_in(:,1:end-13);
%data_out_trials = data_in_trials(1:end-13,:);
else
data_out = data_in(:,14:end);
%data_out_trials = data_in_trials(14:end,:);
end
t_abr = tn(14:end);
tnidx = find(t_abr>-5e-3);
%t_abr = tn_corr(tnIDX);

%% Baseline correction
baseline = mean(data_out(:,find(t_abr>-1e-3 & t_abr<0)),2)'; % baseline from -1 ms to 0
%baseline_trials =mean(data_out_trials(find(t_abr>-1e-3 & t_abr<0),:),1);
data_out = data_out'-baseline;data_out = data_out';
%data_out_trials = data_out_trials-baseline_trials;


end