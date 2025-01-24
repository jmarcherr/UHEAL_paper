
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
run('FFR_mcca_preproc.m')

%% Extract mcca mixing matrix


% mcca with kept channels
nchans = 16;
clc
disp(['running mcca for ' num2str(nchans) ' channels...'])

clear x
for ii=1:length(idx_all)
    x(:,:,ii) = permute(squeeze(TS_sub(idx_all(ii),1:nchans,tidx{1})),[3,2,1]);
end
xx=x(:,:); % concatenate channelwise
C=xx'*xx;
[A,score,AA]=nt_mcca(C,nchans);
z=xx*A; % common space

%% get fft of each component for exclusion
% fft
foi = 326;
for ii=1:size(z,2)
    this_comp = z(:,ii);
    [f_comp,fft_sub_comp(ii,:),f_fft_noise_comp,FFR_comp,F_comp(ii),SNR_comp,F_crit_comp]=get_fft_mcca(this_comp',foi,fs(1));
end
% find significant components (F-test)
c_idx=find(F_comp>=F_crit_comp);
% normalize
norm_F = F_comp(c_idx)./max(F_comp(c_idx));
% find components in top 50%
c_idx = find(norm_F>=0.5);


% select only relevant components
z_sub(:,setdiff(1:size(z,2),c_idx)) = 0;
%transform back to electrode space
% data_z(idx_all(ss),:,:) = permute(z_sub*pinv(AA{ss}(:,:)),[2,1]);

% save matrix
mcca_filter.subid_mcca = subid(idx_all);
mcca_filter.z = z;
mcca_filter.z_sub = z_sub';
mcca_filter.A = A;
mcca_filter.AA = AA;
mcca_filter.score =score;
mcca_filter.tidx_mcca = tidx{1};
mcca_filter.nchans = nchans;
mcca_filter
save('/work3/jonmarc/UHEAL_paper/_eeg/_FFR/_private/single_trial/mcca_filters/mcca_FFR.mat', '-struct','mcca_filter');
