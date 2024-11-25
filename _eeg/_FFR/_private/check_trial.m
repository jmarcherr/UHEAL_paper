
function [trial_info,trial_shift] = check_trial(tmp)
trial_shift = 0;
% check length
if length(tmp)<972
    % wrong length. find first break trial
    tmp2 = tmp(2:end,:)-tmp(1:end-1,:);
    idx = find(tmp2(:,1)==18021 | tmp2(:,1)==18022);
    %first break trial
    idx = idx(1);
    tmp = tmp(idx+1:end,:)
end

for ii=1:length(tmp)
    tmp2 = tmp(2:end,:)-tmp(1:end-1,:);
    % get trial number
    tt = mod(ii+trial_shift,6)
    if tt==0
        tt=6;
    end
    % check if the time between triggers makes sense
    if tt<6 && tmp2(ii)<=8193
        trial_n(ii) = tt;
    elseif tt==6 || tmp2(ii)>=18021 
        trial_n(ii) = 6;
        %trial_shift = 6;
    else
        trial_shift = trial_shift+1;
        trial_n(ii) = tt+trial_shift;
    end
end
trial_info = [tmp(:,4) trial_n'];
end