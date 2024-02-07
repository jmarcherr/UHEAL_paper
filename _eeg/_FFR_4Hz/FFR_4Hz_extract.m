%% plot FFR_4Hz results
% plot FFR_4Hz and extract peaks
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work1/jonmarc/UHEAL_master/UHEAL_paper/UHEAL_startup.m')
subs = dir('_outputs/_derivatives/*.mat')
load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data_current_EEG/uheal_data.mat');
%% get data
for s=1:length(subs)
    
    load([subs(s).folder filesep subs(s).name])
    clc
    disp(['sub ' subs(s).name(1:5) ' loaded...'])
    sub_num(s) = str2num(subs(s).name(3:5));
    chansoi = setdiff(1:16,[5 11]);
    % get FFR
    if isfield(data,'itpc')
        itpc(s,:,:) = data.itpc;
        f = data.f;
        TS_sub(s,:) = nanmean(data.TS(chansoi,:));
        TS_sub_chan(s,:,:)=data.TS;

        time = data.time;
        tidx = data.tidx;
        tidx_TS = data.tidx_TS;
      
        subinfo{s} = data.subinfo;
        if isempty(data.subinfo.age)|isempty(data.subinfo.gender)
            age(s) = nan;
            gender(s) = nan;
        else
        age(s) =data.subinfo.age;
        gender(s) = data.subinfo.gender;
        end

        CP(s) =  uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));%data.subinfo.CP;
        nr_reject(s) =data.nr_reject;
        chan_labels{s} = data.chan_labels;
        chans{s} = data.channels;
        
    else

        itpc(s,:,:) =nan(16,1537);
        TS_sub(s,:,:) = nan(1,5632);
        TS_sub_chan(s,:,:) = nan(16,5632);
        subinfo{s} = data.subinfo;
        age(s) = data.subinfo.age;
        gender(s) = data.subinfo.gender;
    end

    
end
 
%% 
mean(nr_reject)
std(nr_reject)

%% get age groups
% mean spectrum
itpc_spec_mean=squeeze((nanmean(itpc(:,chansoi,:),2)));
% groups
YNH_idx = find(age<=25 & ~CP & ~isnan(itpc_spec_mean(:,1))');
MNH_idx = find(age>25 & age<50 & ~CP & ~isnan(itpc(:,1))')
ONH_idx = find(age>=50 & ~CP & ~isnan(itpc(:,1))');
ages = [17 77];
% colormap
uheal_colormap;
% channels
chansoi  = setdiff(1:16,[5 11]); % all channels but T7 and T8

%% age vs. ITPC frequencies
fid = [2:2:20];
close all

% plot idx
idx = find(~isnan(itpc_spec_mean(:,1)) & ~CP');
% find harmonic frequencies
for ii=1:length(fid)
    fid_idx(ii) = find(f==fid(ii));
    itpc_freq(ii,:) = itpc_spec_mean(:,find(f==fid(ii)))';
end

% Get ratio
itpc_F0 = itpc_freq(1,:);
itpc_F19 = mean(itpc_freq(2:10,:));
itpc_ratio = log10(itpc_freq(1,:)./mean(itpc_freq(2:10,:)));
[rho_ratio,pval_ratio]=corr(age(idx)',itpc_ratio(idx)','type','spearman')


%% ITPC harmonics plot
%close all
figure('renderer','painters')
subplot(1,3,1)
harm = [0:9];
noise_freqs = find(f<=20);
noise_freqs = setdiff(noise_freqs,fid_idx);
itpc_noise= nanmean(itpc_spec_mean(:,noise_freqs));
% plotting init
marks_ag = {'^','sq','o'}
col_ag = {'k','b','r'}
col_ag = {y_col,m_col,o_col}
ag_idx = {YNH_idx,MNH_idx,ONH_idx};
        for ag=1:3
            p_y = plot(harm+1,nanmean(itpc_freq(1:10,ag_idx{ag}),2),'-','color',[col_ag{ag} 0.8])%,%'marker',marks_ag{ag},'markeredgecolor',col_ag{ag},'markerfacecolor','w','linewidth',0.5);
            hold on
            for ll=1:length(harm)

                errorbar(harm(ll)+1,nanmean(itpc_freq(ll,ag_idx{ag}),2),nanstd(itpc_freq(ll,ag_idx{ag}),1)/sqrt(length(ag_idx{ag})),'color',col_ag{ag})
                eby(ag) = plot(harm(ll)+1,nanmean(itpc_freq(ll,ag_idx{ag}),2),'marker',marks_ag{ag},'color',[col_ag{ag} 0.8],'MarkerFaceColor',[col_ag{ag}],'markersize',4);
            end

        end

plot([1 10],[mean(itpc_noise) mean(itpc_noise)],'--','color',[0.75 0.75 0.75],'linewidth',1)
xlabel('Harmonic nr.')
ylabel('ITPC')
set(gca,'xtick',[1:2:10],'xticklabels',[0:2:9])

xlim([0 11])
ylim([0.05 0.45])
box off
hleg=legend([eby(1),eby(2),eby(3)],'Young','Middle-aged','Older');hleg.Box = 'off';
hleg.Position = [0.1567 0.6280 0.1875 0.2683];

set(gcf,'position',[548 464 651 164])

fig = gcf;
saveas(fig,'/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/itpc_harmonics_all','epsc')
%% mean spectrum
figure('renderer','painters')
subplot(1,3,1)
po=plot(f,mean(itpc_spec_mean(idx,:)),'color',y_col,'linewidth',1)
xlim([1 21])
box off
xlabel('Frequency (Hz)')
ylabel('ITPC')
set(gcf,'position',[441 582 443 120])
fig = gcf;

%% Time series alone
% LP for plotting
fs = 1024;
filt_coef = [.5 30]; %4 Hz
filt_def = designfilt('bandpassfir', 'FilterOrder', 40, ...
    'CutoffFrequency1', filt_coef(1), 'CutoffFrequency2', filt_coef(2),...
    'SampleRate', fs);
% bypassing filtfilt fieldtrip
addpath('/appl/matlab/990/toolbox/signal/signal/')

%close all
figure('renderer','painters')
%%%%%%%%%%%%%%%%%%%%%%%
% baseline norm
time_TS = time(tidx_TS);
baseline(:) = mean(TS_sub(:,find(time_TS>-0.1 & time_TS<=0)),2);
baseline_chan=mean(TS_sub_chan(:,1:16,find(time_TS>-0.1 & time_TS<=0)),3);
TS_base = TS_sub-baseline';
%TS_base = filt_def,TS_base;
TS_base_chan = TS_sub_chan-baseline_chan;

ynh_trace = squeeze(nanmean(TS_base(YNH_idx,:)));
ynh_var = squeeze(nanstd(TS_base(YNH_idx,:)))/sqrt(length(YNH_idx));
mnh_trace = squeeze(nanmean(TS_base(MNH_idx,:)));
mnh_var = squeeze(nanmean(TS_base(MNH_idx,:)))/sqrt(length(MNH_idx));
onh_trace = squeeze(nanmean(TS_base(ONH_idx,:)));
onh_var = squeeze(nanstd(TS_base(ONH_idx,:)))/sqrt(length(ONH_idx));

% baseline
ynh_base = ynh_trace-mean(ynh_trace(1,find(time_TS>-0.1 & time_TS<=0)));
mnh_base = mnh_trace-mean(mnh_trace(1,find(time_TS>-0.1 & time_TS<=0)));
onh_base = onh_trace-mean(onh_trace(1,find(time_TS>-0.1 & time_TS<=0)));


eb1=shadedErrorBar(time(tidx_TS),filtfilt(filt_def,ynh_base),filtfilt(filt_def,ynh_var),'lineprops',{'Color',y_col})
hold on
eb2=shadedErrorBar(time(tidx_TS),filtfilt(filt_def,mnh_base),filtfilt(filt_def,mnh_var),'lineprops',{'Color',m_col})
eb3=shadedErrorBar(time(tidx_TS),filtfilt(filt_def,onh_base),filtfilt(filt_def,onh_var),'lineprops',{'Color',o_col})
p1=plot(time(tidx_TS),filtfilt(filt_def,ynh_base),'color',y_col)
p2=plot(time(tidx_TS),filtfilt(filt_def,mnh_base),'color',m_col)
p3=plot(time(tidx_TS),filtfilt(filt_def,onh_base),'color',o_col)

xlim([-.1 3.5])
ylim([-6 3])

box off
ylabel('\muV')
xlabel('Time (s)')

set(gcf,'position',[441 538 512 164])
hold on
fs = 1024;
tstim  = 0:1/fs:0.25
ramp = hann(100);win=[ramp(1:end/2-1)' ones(length(tstim)-100,1)' ramp(end/2:end)'];
stim = [sin(2*pi*256*tstim).*win zeros(size(tstim))];
stim_all = [];
for ii=1:6
    stim_all = [stim_all stim]
end
stim_all = [stim_all zeros(size([tstim tstim]))];
stime = 0:1/fs:length(stim_all)/fs-1/fs;
plot(stime,(stim_all*0.5)-5,'color',[0.5 0.5 0.5 0.5])
hleg = legend([p1 p2 p3],'Young','Middle-aged','Older');
hleg.Box = 'off'
hleg.Position = [0.2297 0.8552 0.6683 0.1037];
hleg.Orientation = 'horizontal'
fig = gcf;
saveas(fig,'/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/TS_all','epsc')
%%
%%% mean amplitude (negativity)
mean_amp = nanmean(TS_base(:,find(time_TS>=0 & time_TS<=3)),2);
figure('renderer','painters')
subplot(1,3,2)
scatter(age(idx),mean_amp(idx),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
ll=lsline
set(ll,'linewidth',2,'color','k')
hold on
ylabel('\muV')
xlabel('Age')
xlim([15 80])

subplot(1,3,3)
errorbar([1 2],[mean(mean_amp(YNH_idx)) mean(mean_amp(ONH_idx))],...
    [std(mean_amp(YNH_idx))./sqrt(length(YNH_idx)) std(itpc_ratio(ONH_idx))/sqrt(length(ONH_idx))],...
    'k.','markersize',1)
set(gca,'xtick',[1 2],'xticklabel',{'Young' 'Older'})
xtickangle(45)
ylabel(['\muV'])
hold on
bar([1 2],[mean(mean_amp(YNH_idx)) mean(mean_amp(ONH_idx))],'FaceColor',[0.5 0.5 0.5])
box off
set(gcf,'position',[441 538 382 164])


lsline
[rho,pmean] = corr(age(idx)',mean_amp(idx))
%%%%%%%%%%%%%%%
%save fig
fig = gcf;



%% mean of all
figure('renderer','painters')
p1=plot(time_TS,filtfilt(filt_def,ynh_base),'color',y_col)
xlim([-.1 3.5])
ylim([-6 3])
hold on
plot(stime,(stim_all*0.5)-5,'color',[0.5 0.5 0.5 0.5])
axis off
set(gcf,'position',[680 923 222 83])
fig =gcf;


%% Topographies mean neg + ITPC ratio
close all
%figure('renderer','painter','position',[680 381 396 322])

top_idx ={YNH_idx,MNH_idx,ONH_idx};
spidx = [1 4 7 10;2 5 8 11];
%chanoi = setdiff(1:16,[5,6,11,13]);
chanoi = 1:16
%subplot(1,4,1)
itpc_ratio = log10(squeeze(itpc_sub(:,:,fid_idx(1)))./squeeze(nanmean(itpc_sub(:,:,fid_idx(2:end)),3)))
for ii=1:3


    zlim1 = [0.1 0.4]; % fundamental    
    zlim2 = [0.05 0.2]; % harmonics
    zlim3 = [-.8 .8]; % negativity
    zlim4 = [-0.1 0.5]; % itpc ratio
    con=0
    figure(1);set(gcf,'renderer','painters')
    subplot(3,3,ii+(ii-1)*2)
    jm_topoplot(squeeze(nanmean(itpc_ratio(top_idx{ii},chanoi),1))',[zlim4],[''],con);
    figure(2);set(gcf,'renderer','painters')
    subplot(3,3,ii+(ii-1)*2)
    jm_topoplot(nanmean(squeeze(nanmean(TS_base_chan(top_idx{ii},chanoi,find(time<=0 & time<=3)),3)))',zlim3,[''],con);
    

end
    %scale
    figure(2)
    subplot(3,3,3)
    c=jm_topoplot(zeros(16,1),[zlim3],'scale',1)
    c.FontSize =10;
    c.Label.String = 'Ampl. mV'
    c.Location = 'southoutside'
    figure(1)
    subplot(3,3,3)
    c=jm_topoplot(zeros(16,1),[zlim4],'scale',1)
    c.FontSize =10;
    c.Label.String = 'ITPC ratio'
    c.Location = 'southoutside'

figure(2)
set(gcf,'position',[680 283 520 420])
fig = gcf;
%saveas(fig,'figs/mean_neg_top','epsc')
figure(1)
set(gcf,'position',[680 283 520 420])
fig = gcf;
%saveas(fig,'/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/itpc_ratio_top','epsc')
%% find P1 N1 like Irsik
% P1 = 0.045-0.065 s

% N1 = 0.085-0.115 s
% onset peak
for ii=1:6 % 6 tones
    for ss=1:size(TS_base,1)
        P1_idx(ii,:) =[0+0.5*(ii-1) 0.085+0.5*(ii-1)];
        P1(ss,ii) = max(TS_base(ss,find(time_TS>=(0+0.5*(ii-1)) & time_TS<=(0.085+0.5*(ii-1)))));

        N1_idx(ii,:) =[0.065+0.5*(ii-1)  0.15+0.5*(ii-1)];
        N1(ss,ii) = min(TS_base(ss,find(time_TS>=(0.065+0.5*(ii-1)) & time_TS<=(0.15+0.5*(ii-1)))));

        P2_idx(ii,:) =[0.15+0.5*(ii-1) 0.25+0.5*(ii-1)];
        P2(ss,ii) = max(TS_base(ss,find(time_TS>=(0.15+0.5*(ii-1)) & time_TS<=(0.25+0.5*(ii-1)))));

        N2_idx(ii,:) = [0.2+0.5*(ii-1) 0.5+0.5*(ii-1)];
        N2(ss,ii) = min(TS_base(ss,find(time_TS>=N2_idx(ii,1) & time_TS<= N2_idx(ii,end))));
        P1P2(ss,ii) = P1(ss,ii)-P2(ss,ii);
    end
end
close all
figure('renderer','painters')
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
ylabel('\mu V')
set(gca,'xtick',[1:6])
set(gcf,'position',[440 449 466 257])
box off
fig = gcf;
saveas(fig,'/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/Peaks','epsc')
%%
% figure('renderer','painters')
% subplot(1,2,2)
% plot(time_TS,squeeze(nanmean(TS_base)))
% hold on
% for ii=1:6
%     plot(P1_idx(ii,:),ones(1,2),'k')
%     plot(N1_idx(ii,:),-1*ones(1,2),'r')
%     plot(P2_idx(ii,:),zeros(1,2),'b')
%     plot(N2_idx(ii,:),-1.5*ones(1,2),'g')
% end
% xlim([-0.1 3.5])
% 
% figure('renderer','painters')
% plot(nanmean(P2(idx,:)-N1(idx,:)))
% hold on
% plot(nanmean(P1(idx,:)-N2(idx,:)))
% plot(nanmean(P1(idx,:)-P2(idx,:)))
% %plot(nanmean(P1P2))
%close all
P1P2 = P1-P2;
P2N1 = P2-N1;
figure('renderer','painters')
subplot(1,2,1)
e1=shadedErrorBar(1:6,nanmean(P2N1(YNH_idx,:)),nanstd(P2N1(YNH_idx,:))/sqrt(length(YNH_idx)),'lineprops',{'Color',y_col})
hold on
e2=shadedErrorBar(1:6,nanmean(P2N1(MNH_idx,:)),nanstd(P2N1(MNH_idx,:))/sqrt(length(MNH_idx)),'lineprops',{'Color',m_col})
e3=shadedErrorBar(1:6,nanmean(P2N1(ONH_idx,:)),nanstd(P2N1(ONH_idx,:))/sqrt(length(ONH_idx)),'lineprops',{'Color',o_col})
title('P2-N1')
legend([e1.mainLine,e2.mainLine,e3.mainLine],{'y','m','o'},'location','best')
set(gca,'xtick',[1:6])
xlabel('tone nr.')
xlim([0 7])
set(gcf,'position',[440 449 466 257])

fig = gcf;
saveas(fig,'/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/P1P2','epsc')
%figure
%scatter(age,nanmean(P1P2,2))
%lsline
%%
figure
plot(nanmean(P2N1(YNH_idx,:)),'k')
hold on
plot(nanmean(P2N1(MNH_idx,:)),'b')
plot(nanmean(P2N1(ONH_idx,:)),'r')
title('P2-N1')
legend({'y','m','o'},'location','best')
xlabel('tone nr.')
xlim([0 7])
%set(gcf,'position',[115 723 331 283])
fig = gcf;
%saveas(fig,'figs/mean_neg_top','epsc')
%figure
%scatter(age,nanmean(P2N1,2))
%lsline

%%
close all

figure
suplot 131
scatter(age,P1)
lsline
subplot 132
scatter(age,N1)
lsline
subplot 133
scatter(age,P2)
lsline


%% extract measures and save to UHEAL_data
% save to UHEAL_data
load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data_current_EEG/uheal_data.mat');
uheal_data.neg_4Hz = nan(size(uheal_data.subid));
uheal_data.itpc_F0_4Hz = nan(size(uheal_data.subid));
uheal_data.itpc_ratio_4Hz = nan(size(uheal_data.subid));


%%
for s=1:length(mean_amp)
    % get this subid
    thisID = sub_num(s);%str2double(sub_id{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);
    
    %OBS from new analysis
    uheal_data.neg_4Hz(this_idx) = mean_amp(s)';
    uheal_data.itpc_F0_4Hz(this_idx) = itpc_F0(s)';
    uheal_data.itpc_ratio_4Hz(this_idx) = itpc_ratio(s);
end

%%
thisdir = cd;
cd('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data_current_EEG')
save('uheal_data.mat','uheal_data')

cd(thisdir)


%% extract p1,n1,p2 in singles script

%%%%%%%%%%%%%%%%%%%% end of script %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% unused
% %% topographies and mean traces_
% 
% %% individual traces
% baseline(:) = mean(TS_sub(:,find(time>-1 & time<0)),2);
%     
% close all
% TS_base = TS_sub;%-baseline';
% subplot(1,2,1)
% plot(time,TS_base(YNH_idx,:),'color',[0 0 0 .5])
% ylim([-10 10])
% subplot(1,2,2)
% hold on
% plot(time, TS_base(ONH_idx,:),'color',[1 0 0 0.5])
% ylim([-10 10])
% 
% %% drift fig
% close all
% figure('renderer','painters')
% plot(repmat([0:0.5:2.5],2,2),repmat([-6;-5],1,12),'-','color',[0.5 0.5 1])
% hold on
% plot(repmat([0.25:0.5:2.75],2,2),repmat([-6;-5],1,12),'-.','color',[1 0.7 0.5])
% 
% eb1=plot(time,ynh_base,'k')
% hold on
% eb2=plot(time,onh_base,'r')
% 
% xlim([-.1 3.5])
% ylim([-6 1.5])
% hleg = legend([eb1 eb2],'YNH','ONH');
% hleg.Box = 'off'
% hleg.Position = [0.5769 0.8964 0.2093 0.1236];
% hleg.Orientation = 'horizontal'
% box off
% ylabel('mV')
% xlabel('Time (s)')
% set(gca,'fontsize',12)
% %set(gcf,'position',[441 426 560 276])
% hA=gca;
% hA.XRuler.MinorTick ='On';
% %set(gcf,'position',[441 566 560 136])
% set(gcf,'position',[441 178 307 524])
% fig = gcf;
% box off
% %saveas(fig,'figs/TS_drift_upright_pt7filt2','epsc')
% %% each harmonic seperately
% close all
% figure('renderer','painters')
% for fid = 2:2:8
%     subplot(1,4,fid/2)
%     eb1 = shadedErrorBar(t_itpc, squeeze(nanmean(itpc_sub_spec(YNH_idx,find(f_itpc==fid),:),1)),...
%         squeeze(nanstd(itpc_sub_spec(YNH_idx,find(f_itpc==fid),:),1))/sqrt(length(YNH_idx)),'lineprops','k');
%     hold on
%     eb2=shadedErrorBar(t_itpc, squeeze(nanmean(itpc_sub_spec(ONH_idx,find(f_itpc==fid),:),1)),...
%        squeeze(nanstd(itpc_sub_spec(ONH_idx,find(f_itpc==fid),:),1))/sqrt(length(ONH_idx)),'lineprops','r');
%     ylim([0.05 .4])
%     xlim([-.5 3.5])
%     title([num2str(fid) ' Hz'])
%     
%     box off
%     
% end
% hleg = legend([eb1.mainLine eb2.mainLine], {'Young','Older'})
% hleg.Box = 'off'
% set(gcf,'position',[655 777 717 147])
% fig = gcf;
% %saveas(fig,'figs/itpc_single_freq','epsc')
% 
% %% individuals
% close all
% figure('renderer','painters')
% for fid = 2:2:8
%     subplot(1,4,fid/2)
%     eb1 = plot(t_itpc, squeeze(itpc_sub_spec(YNH_idx,find(f_itpc==fid),:)),'color',[0 0 0 .5]);
%     hold on
%     eb2=plot(t_itpc, squeeze(itpc_sub_spec(ONH_idx,find(f_itpc==fid),:)),'color',[1 0 0 .5]);
%     ylim([0 .8])
%     xlim([-.5 3.5])
%     title([num2str(fid) ' Hz'])
%     
%     box off
%     %view(90,-90)
% end
% hleg = legend([eb1(1) eb2(2)], {'Young','Older'})
% hleg.Box = 'off'
% set(gcf,'position',[441 524 793 178])
% fig = gcf;
% %saveas(fig,'figs/itpc_single_freq_ind','epsc')
% 
% 
% 
% 
% 
% %% binned phase
% close all;clear avgData
% fs = 1024;
% f_idx = [2,4,8];
% for ff=1:3
% 
%     for ss=1:size(TS_base,1)
%         numBins =round(2*pi*10);% 63;%100;
%         binWidth = 2*pi/numBins;
% 
%         % Define the window width in radians
%         windowWidth = 2*pi%0.063/180;
% 
%         % Load the EEG amplitude data
%         eegData = TS_base(ss,find(time>1 & time<=3));
%         %sin(2*pi*2*t)+sin(2*pi*4*t)+sin(2*pi*6*t)+sin(2*pi*8*t);%TS_base(ss,find(time>0.5 & time<=3.5));
% 
%         % Calculate the phase values based on the sinusoidal signal
%         time_p = (0:length(eegData)-1)/fs;
%         phaseData = 2*pi*f_idx(ff)*time_p;
% 
%         % Initialize a matrix to store the binned EEG data
%         binnedData = nan(numBins, length(eegData));
% 
%         % Loop over all the time points in the EEG data
%         for i = 1:length(eegData)
%             % Calculate the bin index for the current phase value
% 
%             binIndex = mod(floor(phaseData(i)/binWidth),numBins)+1;
% 
% 
%             binnedData(binIndex, i) = eegData(i);
% 
%         end
% 
%         % Calculate the average EEG amplitude within each bin
%         avgData(ss,:) = nanmean(binnedData, 2);
%     end
%     % Plot the results
%     figure(1)
%     subplot(3,3,ff)
%     sy=shadedErrorBar(-pi:binWidth:pi-binWidth, mean(squeeze(avgData(YNH_idx,:))),std(squeeze(avgData(YNH_idx,:)),1)/sqrt(length(YNH_idx)));
%     hold on
%     so=shadedErrorBar(-pi:binWidth:pi-binWidth, mean(squeeze(avgData(ONH_idx,:))),std(squeeze(avgData(ONH_idx,:)),1)/sqrt(length(ONH_idx)),'lineprops','r');
%     xlabel('Phase (radians)');
%     ylabel('Amplitude \muV');
%     set(gca,'xtick',[-pi,0,pi],'xticklabel',{'-\pi','0','\pi'})
%     title([num2str(f_idx(ff)) ' Hz'])
%     %ylim([-2 1.5])
%     xlim([-pi pi])
%     box off
%     set(gcf,'renderer','painters','position',[440 186 781 517])
%     if ff==3
%     hleg =legend([sy.mainLine,so.mainLine],{'Young','Older'})
%     hleg.Box = 'off';
%     hleg.Position = [0.8336 0.8593 0.1063 0.0590];
%     end
%     subplot(3,3,ff+3)
%     %bar([mean(mean(squeeze(avgData(YNH_idx,ff,:)),1)); mean(mean(squeeze(avgData(ONH_idx,ff,:)),1))])
%     %scatter(ones(length(ONH_idx),1)*2,[mean(squeeze(avgData(ONH_idx,ff,:)),2)])
%     %hold on
%     %scatter(ones(length(YNH_idx),1),mean(squeeze(avgData(YNH_idx,ff,:)),2))
%     %xlim([0 3])
% 
%     % age
%     %figure(3)
%     subplot(3,3,ff+3)
%     scatter(age,mean(squeeze(avgData(:,:)),2),'markerfacecolor','w','linewidth',2)
%     h = lsline; h.Color = 'k';h.LineWidth = 2;
%     xlabel('Age')
%     ylabel('Mean amplitude \muV')
%     xlim([15 80])
% 
% 
%     % bars
%     subplot(3,3,ff+6)
%         errorbar([1 2],[mean(mean(squeeze(avgData(YNH_idx,:)),2)) mean(mean(squeeze(avgData(ONH_idx,:)),2))],...
%         [std(mean(squeeze(avgData(YNH_idx,:)),2))/sqrt(length(YNH_idx)) std(mean(squeeze(avgData(ONH_idx,:)),2))/sqrt(length(ONH_idx))],...
%         'k.','markersize',1 )
%         hold on    
%     bar([1 2],[mean(mean(squeeze(avgData(YNH_idx,:)),2)) mean(mean(squeeze(avgData(ONH_idx,:)),2))],'FaceColor',[0.5 0.5 0.5])
%     box off
% 
%     xlim([0 3])
%     set(gca,'xtick',[1 2],'xticklabel',{'Y' 'O'})
%     ylabel('Mean amplitude \muV')
%     % INDIVIDUAL
%     figure(2)
%     subplot(1,3,ff)
%     plot(-pi:binWidth:pi-binWidth, squeeze(avgData(YNH_idx,:)),'color',[0 0 0 0.2]);
%     hold on
%     plot(-pi:binWidth:pi-binWidth, squeeze(avgData(ONH_idx,:)),'color',[1 0 0 0.2])
%     xlabel('Phase (radians)');
%     ylabel('Amplitude \muV');
%     set(gca,'xtick',[-pi,0,pi],'xticklabel',{'-\pi','0','\pi'})
%     title([num2str(f_idx(ff)) ' Hz'])
%     %ylim([-5 3])
%     xlim([-pi pi])
%     box off
%     set(gcf,'renderer','painters','position',[440 525 781 178])
% 
%     avgData_ff(:,:,ff) = avgData;
%     figure(1)
%     fig = gcf;
%     figure(2)
%     fig_ind = gcf;
% 
% end
%     %saveas(fig,'figs/binned_phase','epsc')
%     %saveas(fig,'figs/binned_phase_ind','epsc')
% 
%     %% P1, N1 and P2 from binned phase like irsik
% 
% 
% 
% %% fit function
% close all
% for ff=1:3
% 
% 
% for s=1:size(avgData_ff,1)
% tmp =squeeze(avgData_ff(s,:,ff));
% p=-pi:binWidth:pi-binWidth;
% 
% fun = 'a*((1+cos(x+p))/2)^e+b';
% %fo = fitoptions('Method','NonlinearLeastSquares',...
% %               'Lower',[-pi,1,0],...
% %               'Upper',[pi,inf,max(p)],...
% %               'StartPoint',[0 0 0]);
% %ft = fittype(fun,'problem','a','options',fo);
% %[curve2,gof2] = fit(p',tmp',ft,'problem',2)
% [curve2,gof2]=fit(p',tmp',fun,'StartPoint',[2,0,1,0],'Lower',[-inf,-2*pi,1,-inf]);
% e(s)= curve2.e;
% ar(s) = gof2.adjrsquare;
% c_sub{s} = curve2;
% end
% %%
% % best fit old and young
% [~,idxO]=max(ar(ONH_idx))
% [~,idxY]=max(ar(YNH_idx))
% %close all
% figure(ff)
% %subplot(1,2,1)
% set(gcf,'renderer','painters')
% %subplot(1,2,1)
% y_fit = plot(c_sub{YNH_idx(idxY)},'k')
% hold on
% y_dat=plot(p,squeeze(avgData_ff(YNH_idx(idxY),:,ff)),'k--')
% %figure('renderer','painters')
% %subplot(1,2,2)
% o_fit=plot(c_sub{ONH_idx(idxO)},'r')
% hold on
% o_dat=plot(p,squeeze(avgData_ff(ONH_idx(idxO),:,ff)),'r--')
% set(gca,'xtick', [-pi 0 pi],'xticklabel',{'-\pi' '0' '\pi'})
% xlabel('phase')
% ylabel('amplitude \muV')
% set(gcf,'position',[[141 482 286 221]])
% box off
% hleg = legend([y_dat,y_fit,o_dat,o_fit],{'Y data','Y fit','O data','O fit'})
% hleg.Box = 'off'
% hleg.Location = 'Best'
% fig = gcf;
% saveas(fig,['figs/fit_' num2str(f_idx(ff)) 'hz'],'epsc')
% 
% figure(ff+3)
% set(gcf,'renderer','painters')
% subplot(1,2,1)
% scatter(age,e,'markerfacecolor','w','linewidth',.5)
%     h = lsline; h.Color = 'k';h.LineWidth = 2;
%     xlim([15 80])
%     xlabel('age');ylabel('exponent')
% subplot(1,2,2)
% errorbar([1 2],[mean(e(YNH_idx)) mean(e(ONH_idx))],[std(e(YNH_idx))/sqrt(length(YNH_idx)) std(e(ONH_idx))/sqrt(length(ONH_idx))],'k.','markersize',1)
% hold on
% b=bar([1 2]',[mean(e(YNH_idx)) mean(e(ONH_idx))]','FaceColor',[0.5 0.5 0.5],'edgecolor', 'none')
% ylabel('exponent')
% box off
% hold on
% set(gca,'xtick',[1 2],'xticklabel',{'Y','O'})
% title([num2str(f_idx(ff)) ' hz'])
% set(gcf,'position',[159 242 286 150])
% fig = gcf;
% saveas(fig,['figs/fit_stats_' num2str(f_idx(ff)) 'hz'],'epsc')
% 
% 
% figure(ff+10)
% set(gcf,'renderer','painters')
% subplot(1,2,1)
% scatter(age,ar,'markerfacecolor','w','linewidth',.5)
%     h = lsline; h.Color = 'k';h.LineWidth = 2;
%     xlim([15 80])
%     xlabel('age');ylabel('adj. R^2')
% subplot(1,2,2)
% errorbar([1 2],[mean(ar(YNH_idx)) mean(ar(ONH_idx))],[std(ar(YNH_idx))/sqrt(length(YNH_idx)) std(ar(ONH_idx))/sqrt(length(ONH_idx))],'k.','markersize',1)
% hold on
% b=bar([1 2]',[mean(ar(YNH_idx)) mean(ar(ONH_idx))]','FaceColor',[0.5 0.5 0.5],'edgecolor', 'none')
% box off
% ylabel('adj. R^2')
% set(gca,'xtick',[1 2],'xticklabel',{'Y','O'})
% title([num2str(f_idx(ff)) ' hz'])
% set(gcf,'position',[159 242 286 150])
% 
% fig = gcf;
% saveas(fig,['figs/adjR2_fit_' num2str(ff) 'Hz'],'epsc')
% end
% 



%%
function c=jm_topoplot(var1,zlim,tit_string,coff)
load('/work1/jonmarc/UHEAL_master/UHEAL/_EEG/_func/topo_default.mat');
freq.powspctrm = var1;%nanmean(F_sub(YNH_idx,:))';
cfg = [];
cfg.comment = 'no';
cfg.marker = 'on';
cfg.maarkersymbol = '.';
cfg.maskparameter = 0.5;
cfg.layout = 'biosemi64.lay';
cfg.channel = freq.cfg.channel;
cfg.parameter = 'powspctrm';
cfg.style = 'straight';
cfg.zlim = zlim;
ft_topoplotER(cfg,freq);
title(tit_string)
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path
colormap(flip(brewermap(100,'RdBu'))) % change the colormap
if coff
c=colorbar;
end
end