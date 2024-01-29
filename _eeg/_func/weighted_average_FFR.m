function [data_w,data_w_trial,noise_w] = weighted_average_FFR(data_cc,tidx_FFR,tidx_noisef)
% outputs weighted average of FFR and noise per channel

% get length of channels
for cc=1:size(data_cc,1)
    for ii=1:size(data_cc,3) %chan x time x trials
        % FFR
        var_tmp(cc,ii) =var(data_cc(cc,tidx_FFR,ii));
        data_weighted(cc,:,ii) = data_cc(cc,tidx_FFR,ii)./var_tmp(cc,ii);
        epoch_var_chan(cc,ii) = var_tmp(cc,ii).^(-1);
        % noise
        var_tmp_noise(cc,ii) =var(data_cc(cc,tidx_noisef,ii));
        noise_weighted(cc,:,ii) = data_cc(cc,tidx_noisef,ii)./var_tmp(cc,ii);
        noise_var_chan(cc,ii) = var_tmp(cc,ii).^(-1);
    end
    % FFR single channel
    epoch_var_trial = epoch_var_chan(cc,:);
    summed_trial =sum(data_weighted(cc,:,:),3);
    summed_weights = sum(epoch_var_trial);
    data_w_trial(cc,:) = summed_trial./summed_weights;

    % for concatenated trials
    % Group epochs of FFR and noise together to make 1.5s trials
    catdata = [];
    catnoise = [];
    % How many epochs to group
    epochs_x_trial = 1; % 500 ms
    epochs_x_trial_noise = 2; % 2x250 ms

    % Create a Trial linking epochs_x_trial
    num_full_trials = floor(size(squeeze(data_weighted(cc,:,:)),2)/epochs_x_trial);
    num_full_trials_noise = floor(size(squeeze(data_weighted(cc,:,:)),2)/epochs_x_trial_noise);

    sample_x_epoch = size(data_weighted(cc,:,:),2);
    sample_x_epoch_noise = size(noise_weighted(cc,:,:),2);

    % Reshape
    tmp_data = squeeze(data_weighted(cc,:,:));
    tmp_noise = squeeze(noise_weighted(cc,:,:));
    % FFR
    catdata = reshape(tmp_data(:,1:epochs_x_trial*num_full_trials),...
        sample_x_epoch*epochs_x_trial,...
        num_full_trials);
    % Noise floor
    catnoise = reshape(tmp_noise(:,1:epochs_x_trial_noise*num_full_trials_noise),...
        sample_x_epoch_noise*epochs_x_trial_noise,...
        num_full_trials_noise);

    % Vector of the summed weights for each trial
    epoch_var= reshape(epoch_var_chan(cc,1:num_full_trials*epochs_x_trial),...
        epochs_x_trial, num_full_trials);
    noise_var = reshape(noise_var_chan(cc,1:num_full_trials_noise*epochs_x_trial_noise),...
        epochs_x_trial_noise, num_full_trials_noise);

    % FFR
    epoch_var = epoch_var(:);
    summed_trials =sum(catdata(:,:),2);
    summed_weights = sum(epoch_var);
    data_w(cc,:) = summed_trials./summed_weights;
    %noise
    noise_var = noise_var(:);
    trials_noise =sum(catnoise(:,:),2);
    weights_noise = sum(noise_var);
    noise_w(cc,:) = trials_noise./weights_noise;

end
end