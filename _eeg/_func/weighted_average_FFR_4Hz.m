function [data_w] = weighted_average_FFR_4Hz(data_cc,tidx)
% outputs weighted average of FFR_4Hz per channel

% get length of channels
for cc=1:size(data_cc,1)
    for ii=1:size(data_cc,3) %chan x time x trials
        % FFR_4Hz
        var_tmp(cc,ii) =var(data_cc(cc,tidx,ii));
        data_weighted(cc,:,ii) = data_cc(cc,tidx,ii)./var_tmp(cc,ii);
        epoch_var_chan(cc,ii) = var_tmp(cc,ii).^(-1);

    end
    % FFR_4Hz single channel
    epoch_var_trial = epoch_var_chan(cc,:);
    summed_trial =sum(data_weighted(cc,:,:),3);
    summed_weights = sum(epoch_var_trial);
    data_w_trial(cc,:) = summed_trial./summed_weights;

    % for concatenated trials
    % Group epochs of FFR and noise together to make 1.5s trials
    catdata = [];
    catnoise = [];
    % How many epochs to group
    epochs_x_trial = 1; 

    % Create a Trial linking epochs_x_trial
    num_full_trials = floor(size(squeeze(data_weighted(cc,:,:)),2)/epochs_x_trial);

    sample_x_epoch = size(data_weighted(cc,:,:),2);

    % Reshape
    tmp_data = squeeze(data_weighted(cc,:,:));
    % FFR
    catdata = reshape(tmp_data(:,1:epochs_x_trial*num_full_trials),...
        sample_x_epoch*epochs_x_trial,...
        num_full_trials);

    % Vector of the summed weights for each trial
    epoch_var= reshape(epoch_var_chan(cc,1:num_full_trials*epochs_x_trial),...
        epochs_x_trial, num_full_trials);

    % FFR
    epoch_var = epoch_var(:);
    summed_trials =sum(catdata(:,:),2);
    summed_weights = sum(epoch_var);
    data_w(cc,:) = summed_trials./summed_weights;


end
end