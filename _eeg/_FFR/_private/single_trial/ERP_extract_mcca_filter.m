
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
run('ERP_mcca_preproc.m')

%% Extract mcca mixing matrix


% mcca with kept channels
nchans = 16;
clc
disp(['running mcca for ' num2str(nchans) ' channels...'])

clear x
for ii=1:length(idx_all)
    x(:,:,ii) = permute(squeeze(TS_sub(idx_all(ii),1:nchans,:)),[3,2,1]);
end
xx=x(:,:); % concatenate channelwise
C=xx'*xx;
[A,score,AA]=nt_mcca(C,nchans);
z=xx*A; % common space

%% find a way to extract ERP components


%% plot first 8 components
close all
t=0:1/fs(1):size(z,1)/fs(1)-1/fs(1);
figure('Renderer','painters')
for ii=1:8
    subplot(8,1,ii)
plot(t*1e3,z(:,ii))
subtitle(['SC ' num2str(ii)])
end
set(gcf,'position',[250 122 339 589])
fig = gcf;
%saveas(fig,'figs/mcca_comp','epsc')

% sort components for ABR waveform

for ii=1:size(z,2)

    % mean
    z_w(ii) = mean(z(:,ii));
    %x-correlate with eeg?
    z_w(ii) = max(abs(xcorr(z(:,ii),squeeze(nanmean(mean(TS_sub(:,1:16,:),2))))));


end
[m,i] = sort(z_w,'descend'); % wave 1 components
figure;
for ii=1:9
    % wave I
    subplot(9,1,ii)
    plot(t'*1e3,z(:,i(ii)))
    subtitle(['SC ' num2str(i(ii))])
    
end

set(gcf,'position',[600 6 339 688],'renderer','painters')
sgtitle('xcorr')
fig = gcf;
% %% get fft of each component for exclusion
% % fft
% foi = 326;
% for ii=1:size(z,2)
%     this_comp = z(:,ii);
%     [f_comp,fft_sub_comp(ii,:),f_fft_noise_comp,FFR_comp,F_comp(ii),SNR_comp,F_crit_comp]=get_fft_mcca(this_comp',foi,fs(1));
% end
% % find significant components (F-test)
% c_idx=find(F_comp>=F_crit_comp);
% % normalize
% norm_F = F_comp(c_idx)./max(F_comp(c_idx));
% % find components in top 50%
% c_idx = find(norm_F>=0.5);

%%
clear data_z
%transform back to electrode space
for ss=1:size(x,3)
    % subject specific z
    z_sub = squeeze(TS_sub(idx_all(ss),1:16,:))'*AA{(ss)}(:,:);
    % select only relevant components
    z_sub(:,setdiff(1:size(z,2),i(1:10))) = 0;
    sub_filt{ss} = z_sub;
   
    %transform back to electrode space
    data_z(ss,:,:) = permute(z_sub*pinv(AA{ss}(:,:)),[2,1]);

end
% scale
rms_dat = squeeze(rms(TS_sub(idx_all,:,:),3))
rms_z = squeeze(rms(data_z(:,:,:),3));
data_z_sc = rms_dat.*(data_z(:,:,:)./rms_z);
close all
plot(squeeze(nanmean(data_z_sc(:,10,:))))
hold on
plot(squeeze(nanmean(TS_sub(:,10,:))))


%%
% save matrix
mcca_filter.subid_mcca = subid(idx_all);
mcca_filter.z = z;
mcca_filter.z_sub = sub_filt;
mcca_filter.comp_keep = i(1:3);
mcca_filter.A = A;
mcca_filter.AA = AA;
mcca_filter.score =score;
mcca_filter.tidx_mcca = tidx{1}';
mcca_filter.nchans = nchans;
mcca_filter
%save('/work3/jonmarc/UHEAL_paper/_eeg/_FFR/_private/single_trial/mcca_filters/mcca_ERP.mat', '-struct','mcca_filter');


%% try out single trial

% get n1,p2 peaks

n1_idx = find(t>=0.08 & t<=0.15);
p2_idx = find(t>=0.1 & t<=0.3);
for ss=1:size(mcca_filter.subid_mcca,2)
    %find sub id
    this_sub = find(strcmp(subid{ss},mcca_filter.subid_mcca));
    if this_sub
        for tt=1:size(TS_trials{this_sub},3)
            %raw amp
            n1p2(ss,tt) = max(nanmean(TS_trials{this_sub}(:,n1_idx,tt),1))-min(nanmean(TS_trials{this_sub}(:,n1_idx,tt),1));
            % zsub
            z_trial = squeeze(TS_trials{this_sub}(1:16,:,tt))'*AA{(ss)}(:,:);
            z_trial(:,setdiff(1:size(z,2),mcca_filter.comp_keep)) = 0;
            %transform back to electrode space
            data_z = permute(z_trial*pinv(AA{ss}(:,:)),[2,1]);
            % scale
            rms_dat = squeeze(rms(TS_trials{this_sub}(:,:,tt),2));
            rms_z = squeeze(rms(data_z(:,:,:),2));
            % scaling does not work
            data_z_sc_trial(tt,:,:) = rms_dat.*(data_z(:,:,:)./rms_z);
            %data_z_sc_trial(tt,:,:) = data_z;
            n1p2_mcs{ss}(tt) = max(nanmean(data_z_sc_trial(tt,:,:),2))-min(nanmean(data_z_sc_trial(tt,:,:),2));
        end
    else
    end
    disp(['subject ' num2str(ss) ' done'])
end

%%
close all
tt=100
tidx_erp = find(t>=0 & t<=0.3)
 subplot 121
plot(t(tidx_erp),mean(squeeze(TS_trials{1}(10,tidx_erp,:))'))
subplot 122
plot(t(tidx_erp),squeeze(data_z_sc_trial(tt,10,tidx_erp)))
hold on
plot(t(tidx_erp),squeeze(TS_trials{1}(10,tidx_erp,tt)))

figure
plot(t(tidx_erp),-squeeze(mean(mean(data_z_sc_trial(:,:,tidx_erp),1),2))')

figure
scatter(1:length(n1p2_mcs{1}),n1p2_mcs{1})
