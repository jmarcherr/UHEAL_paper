
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
run('FFR_mcca_preproc_trials2.m')

%%
% plot average data for each trial config
close all
for tt = 1:6
    subplot(1,6,tt)
    plot(time{1},squeeze(mean(nanmean(FFR_trials(:,tt,10,:),1),1)),'Color',[0.6 0.6 0.6])
    ylim([-0.3 0.3])
    box off
end
    
set(gcf,'Position',[440 459 849 143],'renderer','painters')
fig = gcf;
saveas(fig,'figs/trials2/data_TS_16','svg')


%% loop over trial nrs - do MCCA
% mcca with kept channels

    nchans = 16;
    clc
    disp(['running mcca ...'])

    clear x
    for ii=1:length(idx_all)
        x(:,:,ii) = permute(squeeze(TS_sub(idx_all(ii),1:16,tidx{1})),[3,2,1]);
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
    norm_F = F_comp./max(F_comp);
    % find components in top 50%
    c_idx = find(norm_F>=0.5);

    %% SNR calculation
    % init variables
    data_z = nan(size(TS_sub(:,1:16,:)));
    FFR_mcca = nan(size(TS_sub,1),1);
    SNR_mcca = nan(size(TS_sub,1),1);
    sig_mcca = nan(size(TS_sub,1),1);
    SNR_chan_mcca = nan(size(x,3),16);
    a = nan(size(x,3),16);
    for tt = 1:6
    for ss=1:size(x,3)
        % subject specific z for tt tone in sequence of 6
        z_sub = squeeze(FFR_trials(idx_all(ss),tt,1:16,:))'*AA{(ss)}(:,:);

        % get mixing weights for first 9 components (only for 16 channels)
            for cc=1:9
                a_all(ss,cc,:) = squeeze(FFR_trials(idx_all(ss),tt,1:16,tidx{1}))*nt_normcol(z(:,cc))/size(data_z,3);
            end

        % select only relevant components
        z_sub(:,setdiff(1:size(z,2),c_idx)) = 0;
        %transform back to electrode space
        data_z(idx_all(ss),:,:) = permute(z_sub*pinv(AA{ss}(:,:)),[2,1]);

        % get mixing weights (see mcca_demo3)
        a(ss,:)=squeeze(TS_sub(idx_all(ss),1:16,tidx{1}))*nt_normcol(z(:,1))/size(data_z,3);



        % get fft
        foi = [326]; % pure tone frequency

        %get fft FFR mcca
        [f_tmp,fft_sub_tmp,~,~,~,SNR_tmp,~]=get_fft(squeeze(data_z(idx_all(ss),:,tidx{1})),foi,fs(s));
        SNR_chan_mcca(idx_all(ss),:) = SNR_tmp;
        % get channel average SNR/FFR over remaining channels
        [~,FFR_avg_tmp,~,SNR_avg_tmp,~,sig_idx_avg_tmp,~]=get_fft_chaoi(f_tmp,fft_sub_tmp,1:16,foi);
        FFR_mcca(idx_all(ss)) = FFR_avg_tmp;
        SNR_mcca(idx_all(ss)) = SNR_avg_tmp;
        sig_mcca(idx_all(ss)) = sig_idx_avg_tmp;

        %get fft FFR raw
        [f_tmp,fft_sub_tmp,f_fft_noise_tmp,FFR_tmp,F_tmp,SNR_tmp,F_crit_tmp]=get_fft(squeeze(FFR_trials(idx_all(ss),tt,:,tidx{1})),foi,fs(s));
        SNR_chan_raw(idx_all(ss),:) = SNR_tmp;
        % get channel average SNR/FFR over remaining channels
        [f_fft_avg,FFR_avg_tmp,F_avg_tmp,SNR_avg_tmp,F_crit_avg_tmp,sig_idx_avg_tmp,noise_avg_tmp]=get_fft_chaoi(f_tmp,fft_sub_tmp,1:16,foi);
        FFR_raw(idx_all(ss)) = FFR_avg_tmp;
        SNR_raw(idx_all(ss)) = SNR_avg_tmp;
        sig_raw(idx_all(ss)) = sig_idx_avg_tmp;

    end
    %% get variables for plotting
    weight_maps{tt} = a;
    SNR_mcca_it{tt} = SNR_mcca;
    SNR_avg_it{tt} = SNR_raw;
    c_idx_it{tt} = c_idx;
    sig_mcca_it{tt} = sig_mcca(idx_all);
    sig_raw_it{tt}  =
    n_sig_mcca(tt) = length(find(sig_mcca(idx_all)));
    n_sig_raw(tt) = length(find(sig_raw(idx_all)));
    n_subs = length(find(idx_all));
    norm_F_it{tt} = norm_F;
    data_z_it{tt} = data_z;
    z_it{tt} = z;
    end
%%
%% plotting
%% plot SNR of MCCA vs. raw
close all
fsize = 11
figure(1000)
set(gcf,'Renderer','painters')
cmap_mcca=cbrewer('div','BrBG',15);
cmap_trials = cbrewer('div','BrBG',10);
subplot 222
pm=plot((1:6),n_sig_mcca(1:end)./length(idx_all),'Color',cmap_mcca(12,:),'LineWidth',2)
hold on
pr=plot((1:6),n_sig_raw(1:end)./length(idx_all),'Color',cmap_mcca(2,:),'LineWidth',2)
%xlim([0 100])
ylim([.25 1])
%set(gca,'ytick',[.8 .9 1],'fontsize',fsize)
set(gca,'fontsize',fsize)
xlabel('% trials')
ylabel('% sig. subjects')
hleg = legend([pr pm],{'raw','mcca'})
hleg.Box = 'off'
hleg.Position = [0.6886 0.6671 0.1670 0.0850];
box off
fig = gcf;

%% plot SNR per trial
close all

for ii=1:6
    for ss=1:length(idx_all)
        if sig_mcca(ss) == 1
        this_mcca(ss,ii) = SNR_mcca_it{ii}(idx_all(ss));
        this_raw(ss,ii) = SNR_avg_it{ii}(idx_all(ss))
        else
            this_mcca(ss,ii) = nan;
            this_raw(ss,ii) = nan;
        end
    end
   
end
plot(1:6,nanmean(this_mcca))
%hold on
%plot(1:6,nanmean(this_raw))
%plot(1:6,nanmean(this_mcca(y,:)))
%plot(1:6,nanmean(this_mcca(o,:)))




%%
function c=jm_topoplot(var1,zlim,tit_string,coff)
load('/work3/jonmarc/UHEAL_master/UHEAL/_EEG/_func/topo_default.mat');
freq.powspctrm = var1;%nanmean(F_sub(YNH_idx,:))';
cfg = [];
cfg.comment = 'no';
cfg.marker = 'on';
cfg.maarkersymbol = '.';
cfg.layout = 'biosemi64.lay';
cfg.channel = freq.cfg.channel;
cfg.parameter = 'powspctrm';
cfg.style = 'straight';
cfg.zlim = zlim;
ft_topoplotER(cfg,freq);
title(tit_string)
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flip(brewermap(100,'BrBG'))) % change the colormap
%colormap(brewermap(64,'YlOrRd')) % change the colormap
if coff
    c=colorbar;
else
    c=nan;
end
end



