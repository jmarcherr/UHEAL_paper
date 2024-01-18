function [valid_trials,rjt] = threshold_rejection(epoched_data,thr,tidx,trials_oi)
%Threshold rejection of epoched data per channel
% !!!!!make so that if there are no inputs on 3 and 4, just take all trials 
% and times !!!!

        rjt_trials = [];
        for ii=1:size(epoched_data,1) % loop over trials
            for cc= 1:size(epoched_data,2) %loop over channels
                data_art = epoched_data(ii,cc,tidx); %single epoch
                if max(abs(data_art)) > thr
                    rjt_trials = [rjt_trials ii];

                end
            end
        end
        rjt = unique(rjt_trials);
        trials_oi = [trials_oi];
        % clean trials
        valid_trials = setxor(rjt,1:length(trials_oi)); 
       
end