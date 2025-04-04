
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
run('FFR_mcca_preproc.m')

%% subject downsampling
% loop over subject iterations from 10:100 %
split_ps = [0.1:0.1:1];
% create datasets from idx
for it = 1:length(split_ps)


    % over 10 repititions, pick random subjects
    for itt = 1:10
        % get 10 splits of the data per percent
        datasets{it,itt} = randperm(length(idx_all),round(split_ps(it)*length(idx_all)));
    end
end

%%
close all
%loop over split iterations
for itt = 1:10
    for it = 1:10
        this_split = idx_all(datasets{it,itt});
        this_idx = 1:length(this_split);
        % mcca with kept subjects
        nsubs = length(this_split);
        clc
        disp(['iteration ' num2str(itt) '/' num2str(10)])
        disp(['running mcca for ' num2str(nsubs) ' subjects...'])

        clear x
        for ii=1:length(this_split)
            x(:,:,ii) = permute(squeeze(TS_sub(this_split(ii),1:16,tidx{1})),[3,2,1]);
        end
        xx=x(:,:); % concatenate channelwise
        C=xx'*xx;
        [A,score,AA]=nt_mcca(C,16);
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

        %% SNR calculation
        % init variables
        data_z = nan(size(TS_sub(nsubs,:,:)));
        FFR_mcca = nan(nsubs,1);
        SNR_mcca = nan(nsubs,1);
        FFR_raw = nan(nsubs,1);
        SNR_raw = nan(nsubs,1);
        sig_mcca = zeros(nsubs,1);
        sig_raw = zeros(nsubs,1);
        SNR_chan_mcca = nan(size(x,3),16);
        a = nan(size(x,3),16);
        for ss=1:size(x,3)
            % subject specific z
            z_sub = squeeze(TS_sub(this_split(ss),1:16,:))'*AA{(ss)}(:,:);

            % get mixing weights for first 9 components (only for 16 channels)
                for cc=1:9
                    a_all(ss,cc,:) = squeeze(TS_sub(this_split(ss),1:16,tidx{1}))*nt_normcol(z(:,cc))/size(data_z,3);
                end

            % select only relevant components
            z_sub(:,setdiff(1:size(z,2),c_idx)) = 0;
            %transform back to electrode space
            data_z(ss,:,:) = permute(z_sub*pinv(AA{ss}(:,:)),[2,1]);

            % get mixing weights (see mcca_demo3)
            a(ss,:)=squeeze(TS_sub(this_split(ss),1:16,tidx{1}))*nt_normcol(z(:,1))/size(data_z,3);



            % get fft
            foi = [326]; % pure tone frequency

            %get fft FFR mcca
            [f_tmp,fft_sub_tmp,~,~,~,SNR_tmp,~]=get_fft(squeeze(data_z(ss,:,tidx{1})),foi,fs(s));
            SNR_chan_mcca(ss,:) = SNR_tmp;
            % get channel average SNR/FFR over remaining channels
            [~,FFR_avg_tmp,~,SNR_avg_tmp,~,sig_idx_avg_tmp,~]=get_fft_chaoi(f_tmp,fft_sub_tmp,1:16,foi);
            FFR_mcca(ss) = FFR_avg_tmp;
            SNR_mcca(ss) = SNR_avg_tmp;
            sig_mcca(ss) = sig_idx_avg_tmp;

            %get fft FFR raw
            [f_tmp,fft_sub_tmp,f_fft_noise_tmp,FFR_tmp,F_tmp,SNR_tmp,F_crit_tmp]=get_fft(squeeze(TS_sub(this_split(ss),:,tidx{1})),foi,fs(s));
            SNR_chan_raw(ss,:) = SNR_tmp;
            % get channel average SNR/FFR over remaining channels
            [f_fft_avg,FFR_avg_tmp,F_avg_tmp,SNR_avg_tmp,F_crit_avg_tmp,sig_idx_avg_tmp,noise_avg_tmp]=get_fft_chaoi(f_tmp,fft_sub_tmp,1:16,foi);
            FFR_raw(ss) = FFR_avg_tmp;
            SNR_raw(ss) = SNR_avg_tmp;
            sig_raw(ss) = sig_idx_avg_tmp;

        end



        %% get variables for plotting
        weight_maps{it,itt} = a;
        SNR_mcca_it{it,itt} = SNR_mcca;
        SNR_avg_it{it,itt} = SNR_raw;
        c_idx_it{it,itt} = c_idx;
        sig_mcca_it{it,itt} = sig_mcca;
        sig_raw_it{it,itt} = sig_raw;
        n_sig_mcca(it,itt) = length(find(sig_mcca))/nsubs;
        n_sig_raw(it,itt) = length(find(sig_raw))/nsubs;
        n_subs(it) = length(this_idx);
        norm_F_it{it,itt} = norm_F;
        data_z_it{it,itt} = data_z;
        z_it{it,itt} = z;

    end
end

%% plotting
%% plot SNR of MCCA vs. raw
close all
fsize = 11
figure(1000)
set(gcf,'Renderer','painters')
cmap_mcca=cbrewer('div','BrBG',15);
subplot 222
pm=errorbar(11:11:110,mean(n_sig_mcca,2)',min(n_sig_mcca')-mean(n_sig_mcca,2)',max(n_sig_mcca')-mean(n_sig_mcca,2)','Color',cmap_mcca(12,:),'LineWidth',2)
hold on
pr=errorbar(11:11:110,mean(n_sig_raw,2)',min(n_sig_raw')-mean(n_sig_raw,2)',max(n_sig_raw')-mean(n_sig_raw,2)','Color',cmap_mcca(2,:),'LineWidth',2)
%xlim([0 11])
%ylim([.5 1.1])
set(gca,'fontsize',fsize)
xlabel('nr. subjects')
ylabel('% sig. subjects')
hleg = legend([pr pm],{'raw','mcca'})
hleg.Box = 'off'
hleg.Position = [0.7269 0.6416 0.1554 0.0917];
box off
fig = gcf;
axis padded
%
% SNR improvement
subplot(2,2,[1 3])
% iteration
for itt = 1:10
for ii=10:-1:1
    this_diff = [];
    for ss = 1:n_subs(ii)
        this_mcca = SNR_mcca_it{ii,itt}(ss); this_raw = SNR_avg_it{ii,itt}(ss);
        if  this_mcca>=this_raw
            this_diff(ss) = dist(SNR_mcca_it{ii,itt}(ss),SNR_avg_it{ii,itt}(ss));
        elseif this_mcca<=this_raw
            this_diff(ss) = -dist(SNR_mcca_it{ii,itt}(ss),SNR_avg_it{ii,itt}(ss));
        end
    end

    ps(it) = plot(1:n_subs(ii),sort(this_diff,'descend'),'color',cmap_mcca(ii,:),'LineWidth',1)
    hold on
    plot([1 length(n_subs(ii))],[0 0],'--k')
    xlim([-5 115])
end
end
box off
axis padded
xlabel('subject')
ylabel('\Delta SNR (dB)')
set(gca,'fontsize',fsize,'xtick',[1 50 100])
set(gcf,'position',[176 244 521 453])
%
% average snr improvement and variance across subjects
subplot 224
for ii=1:10
    this_y = [];
    for it = 1:10
    this_y(it,:) = mean(SNR_mcca_it{ii,it}-SNR_avg_it{ii,it},2);
    end
    this_y = mean(this_y);
    errorbar(ii*11,mean(this_y),std(this_y'),'o','color',cmap_mcca(13,:),'markerfacecolor',cmap_mcca(13,:))

    hold on
    %plot([1 17],[0 0],'--k')
    %ylim([-5 15])
    %xlim([0 11])
    box off
    axis padded
end
xlabel('nr. subjects')
ylabel('\Delta SNR (dB)')

set(gca,'fontsize',fsize)
fig = gcf;
saveas(fig,'figs/mcca_it/mcca_it_summary_subjects','svg')

%%
% colorbar
figure('Renderer','painters')
colormap(cmap_mcca)
c=colorbar;
c.Ticks = [[2 4 6 8 10 12 14 16]/16];
c.Limits = [2/16 1]
c.TickLabels = {'2','4','6','8','10','12','14','16'}
c.Orientation = 'Horizontal'
c.Location = 'west'
set(gcf,'position',[767 193 125 151])

fig = gcf;
%saveas(fig,'figs/mcca_it/mcca_it_summary_cb','epsc')
%% as distributions
% iteration
for itt = 1:10
for ii=10:-1:1
    this_diff = [];
    for ss = 1:n_subs(ii)
        this_mcca = SNR_mcca_it{ii,itt}(ss); this_raw = SNR_avg_it{ii,itt}(ss);
        if  this_mcca>=this_raw
            this_diff(ss) = dist(SNR_mcca_it{ii,itt}(ss),SNR_avg_it{ii,itt}(ss));
        elseif this_mcca<=this_raw
            this_diff(ss) = -dist(SNR_mcca_it{ii,itt}(ss),SNR_avg_it{ii,itt}(ss));
        end
    end

    hist(this_diff)
    hold on
    %plot([1 length(n_subs(ii))],[0 0],'--k')
    %xlim([-5 115])
end
end
box off
axis padded
xlabel('subject')
ylabel('\Delta SNR (dB)')
%set(gca,'fontsize',fsize,'xtick',[1 50 100])
set(gcf,'position',[176 244 521 453])
%% weight maps
    figure(500)
    set(gcf,'Renderer','painters')
    subplot 211
    plot(mean(weight_maps{16}),'LineWidth',1.5,'Color',cmap_mcca(12,:))
    hold on
    ylim([0.010 0.035])
    box off
    xlabel('channel')
    ylabel('weight')
    set(gca,'xtick',1:16,'xticklabels',data.chan_labels(1:16),'fontsize',10)
    subplot 212

    jm_topoplot(mean(weight_maps{16})',[],[],0)
    set(gcf,'position',[440 332 361 290])
fig = gcf;
saveas(fig,'figs/mcca_it/weight_maps','svg')
    %% scatter raw vs. mcca 16 chan

figure(300)
set(gcf,'Renderer','painters')

subplot 122
scatter(SNR_mcca_it{16}(idx_all),SNR_avg_it{16}(idx_all)','.k')
xlabel('mcca SNR (dB)');ylabel('raw SNR (dB)')
hold on
plot([-40 60],[-40 60],'color',cmap_mcca(15,:))
ylim([-40 60])
xlim([-40 60])
text(-30,40,['mcca: ' num2str(n_sig_mcca(16)) '/' num2str(n_subs) '\newline' ,...
    'raw: ' num2str(n_sig_raw(16)) '/' num2str(n_subs)],'fontsize',8)
set(gca,'Fontsize',10)
set(gcf,'position',[440 420 466 202])

% SNR improvement
subplot 121
plot(1:length(idx_all),sort(SNR_mcca_it{16}(idx_all)-SNR_avg_it{16}(idx_all)','descend'),'b')
hold on
plot([1 length(idx_all)],[0 0],'--k')
box off
xlabel('Subjects')
ylabel('\Delta SNR (dB)')
set(gca,'FontSize',10)
fig = gcf;
saveas(fig,'figs/mcca_it/scatter','svg')
%%
close all
% SNR over subjects raw
figure('Renderer','painters')
plot(1:length(idx_all),sort(SNR_avg_it{16}(idx_all),'descend'),'color',cmap_mcca(13,:),'linewidth',2)
hold on
p1=plot([0 111],[db(F_critt(1)) db(F_critt(1))],'--k')
xlabel('subjects')
ylabel('FFR SNR (dB)')
hleg = legend(p1,'SNR crit.')
hleg.Box = 'off'
box off
set(gca,'Fontsize',10)
set(gcf,'position',[440 456 220 166])
fig = gcf;
saveas(fig,'figs/mcca_it/raw_SNR_all','svg')
%% time domain
% scaling
rms_dat = squeeze(rms(TS_sub(idx_all,:,:),3))
rms_z = squeeze(rms(data_z_it{16}(idx_all,:,:),3));
data_z_sc = rms_dat.*(data_z_it{16}(idx_all,:,:)./rms_z);

figure(600)
set(gcf,'renderer','painters')
% mean over channels and scaling
subplot 121

plot(time{1},mean(squeeze(nanmean(data_z_sc))),'Color',[0.6 0.6 0.6])
ylim([-0.3 0.3])
xlim([-0.05 0.49])
set(gca,'Fontsize',10)
title('mcca denoised')
box off
subplot 122

plot(time{1},mean(squeeze(nanmean(TS_sub(idx_all,1:16,:)))),'color',[0.6 0.6 0.6])
ylim([-0.3 0.3])
xlim([-0.05 0.49])
set(gca,'Fontsize',10)

title('raw')
box off
set(gcf,'Position',[440 493 385 103],'renderer','painters')
fig = gcf;
saveas(fig,'figs/mcca_it/TS','svg')

%% weight maps components
figure(700)
set(gcf,'renderer','painters')
for cc=1:9
    subplot(3,3,cc)
        jm_topoplot(squeeze(nanmean(a_all(:,cc,:),1)),[],[],0)
        subtitle(['SC' num2str(cc)])
        set(gca,'Fontsize',10)
end
set(gcf,'position',[797 304 471 395])
fig = gcf;
saveas(fig,'figs/mcca_it/weight_maps_topo','svg')
%% plot first 8 components
figure(800)
set(gcf,'renderer','painters')
t=0:1/fs(1):length(z_it{16})/fs(1)-1/fs(1);

for ii=1:8
    subplot(8,1,ii)
plot(t,z_it{16}(:,ii),'color',[0.6 0.6 0.6])
st = subtitle(['SC ' num2str(ii)])
%st.HorizontalAlignment = 'right'
%st.Position = [0.25 0.2503 0]

set(gca,'Fontsize',10)
box off
axis off
end
set(gcf,'position',[797 304 471 395])
fig = gcf;
saveas(fig,'figs/mcca_it/SCs_TS','svg')
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




