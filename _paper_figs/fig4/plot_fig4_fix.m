% plot abr and extract peaks
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
d = dir([rootdir '/_eeg/_AEP/_outputs/_derivatives/*.mat'])
savepath = [rootdir '/_paper_figs/fig4/figs/'];

for s=1:length(d)
   
    load([d(s).folder filesep d(s).name])
    extract_aep_data
    clc
    disp([sub_id{s} ' done...'])
end


%% plot
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat')

% get age groups
idx = ~uheal_data.CP_new
y = find(age<=25 & idx');
m = find(age>25 & age<50 & idx');
o = find(age>=50 & idx');

% get colormap
uheal_colormap;
%% plot average ERPs for each ISI
close all
figure('Renderer','painter')

rates = [0.5 1 1.5 2]+0.3
YNH_idx = y;
MNH_idx = m;
ONH_idx = o;

% baseline correct average trace
for bb=1:4
    baseline(:,bb,:) = nanmean(aep_sub_filt(:,bb,find(time>-0.1 & time<=0.01)),3)
end
TS_baseline = aep_sub_filt-baseline;

% plot each isi
for kk=1:4
    % young
    subplot(1,4,kk)
    pmean_y(kk) = shadedErrorBar(time,squeeze(nanmean(TS_baseline(YNH_idx,kk,:),1)),squeeze(nanstd(TS_baseline(YNH_idx,kk,:),1))/sqrt(length(find(YNH_idx))),'lineprops',{'-','color',y_col,'linewidth',1},'transparent',1)
    hold on
    % mid.aged
    pmean_m(kk) = shadedErrorBar(time,squeeze(nanmean(TS_baseline(MNH_idx,kk,:),1)),squeeze(nanstd(TS_baseline(MNH_idx,kk,:),1))/sqrt(length(find(MNH_idx))),'lineprops',{'-','color',m_col,'linewidth',1},'transparent',1)
    % older
    pmean_o(kk) = shadedErrorBar(time,squeeze(nanmean(TS_baseline(ONH_idx,kk,:),1)),squeeze(nanstd(TS_baseline(ONH_idx,kk,:),1))/sqrt(length(find(ONH_idx))),'lineprops',{'-','color',o_col,'linewidth',1},'transparent',1)
    plot([0 0],[-6 5],'k--')
    plot([0.3 0.3],[-6 5],'k--')

    title(['ISI = ' num2str(rates(kk)) ' s' ],'FontSize',12,'FontWeight','normal')
    if kk==1
        ylabel('Amplitude \muV')
    end
    xlabel('Time (s)')
    box off
    set(gca,'Fontsize',12)
    xlim([-.1 .5])
    ylim([-5.5 3])
end

%set(gcf,'position',[144 267 1315 368])
%set(gcf,'position',[85 364 958 272]);
set(gcf,'position',[85 331 1093 305])
hleg = legend([pmean_y(1).mainLine pmean_m(1).mainLine pmean_o(1).mainLine],{'Young','Middle-aged','Older'});
hleg.Box = 'off'
hleg.Position = [0.1474 0.2191 0.1226 0.1836];hleg.FontSize = 10;%hleg.FontName = 'Arial'

fig = gcf;
print([savepath 'AEP_all'],'-dsvg')
% plot N1,P2 and P2-N1
com_sub =  p2_fix_sub-n1_fix_sub;
%%
%close all
% N100
figure('Renderer','painter')
% plot mean n1
plot(rates,nanmean(n1_fix_sub(YNH_idx,:)),'color',y_col)
hold on
plot(rates,nanmean(n1_fix_sub(MNH_idx,:)),'color',m_col)
plot(rates,nanmean(n1_fix_sub(ONH_idx,:)),'color',o_col)
% plot n1 for age groups
for ii=1:4
    %young
    eb(2,ii)=errorbar(rates(ii)-.01,nanmean(n1_fix_sub(YNH_idx,ii)),nanstd(n1_fix_sub(YNH_idx,ii))/sqrt(length(YNH_idx)),'k^','color',y_col,'markerfacecolor',y_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    %mid aged
    eb(3,ii)=errorbar(rates(ii),nanmean(n1_fix_sub(MNH_idx,ii)),nanstd(n1_fix_sub(MNH_idx,ii))/sqrt(length(MNH_idx)),'ksq','color',m_col,'markerfacecolor',m_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    %older
    eb(1,ii)=errorbar(rates(ii)+.01,nanmean(n1_fix_sub(ONH_idx,ii)),nanstd(n1_fix_sub(ONH_idx,ii))/sqrt(length(ONH_idx)),'ko','color',o_col,'markerfacecolor',o_col,'MarkerEdgecolor','k')
end

ylabel('N1 \muV')
xlabel('ISI (s)')
set(gca,'fontsize',12)
set(gcf,'position',[163 12 270 266])
%hleg = legend([eb(2) eb(3) eb(1)],'Young','Middle-aged','Older')
%hleg.Position = [0.1623 0.2148 0.5185 0.2218];
xlim([0.5 2.5])
ylim([-5 .5])
set(gca,'xtick',rates)
%hleg.Box = 'off'
%hleg.FontSize = 10;
%hleg.FontName = 'Arial'
box off
fig=gcf
saveas(fig,[savepath '/n100_age_group_fix'],'svg')

% P2
figure('Renderer','painter')
% Young
plot(rates,nanmean(p2_fix_sub(YNH_idx,:)),'color',y_col)
hold on
plot(rates,nanmean(p2_fix_sub(MNH_idx,:)),'color',m_col)
plot(rates,nanmean(p2_fix_sub(ONH_idx,:)),'-','color',o_col)
for ii=1:4
    eb(2,ii)=errorbar(rates(ii)-.01,nanmean(p2_fix_sub(YNH_idx,ii)),nanstd(p2_fix_sub(YNH_idx,ii))/sqrt(length(YNH_idx)),'^','color',y_col,'markerfacecolor',y_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    eb(2,ii)=errorbar(rates(ii),nanmean(p2_fix_sub(MNH_idx,ii)),nanstd(p2_fix_sub(MNH_idx,ii))/sqrt(length(MNH_idx)),'sq','color',m_col,'markerfacecolor',m_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    eb(1,ii)=errorbar(rates(ii)+.01,nanmean(p2_fix_sub(ONH_idx,ii)),nanstd(p2_fix_sub(ONH_idx,ii))/sqrt(length(ONH_idx)),'o','color',o_col,'markerfacecolor',o_col,'MarkerEdgecolor','k')

end


ylabel('P2 \muV')
xlabel('ISI (s)')
set(gca,'fontsize',12)
set(gcf,'position',[438 10 270 266])
set(gca,'xtick',rates)
xlim([0.5 2.5])
ylim([-1 2.5])
fig=gcf
box off
saveas(fig,[savepath 'p200_age_group_fix'],'svg')

% P2-N1 complex
% P2
figure('Renderer','painters')
% Young
%plot(rates,nanmean(com_sub(YNH_idx,:)),'color',y_col)
%hold on
%plot(rates,nanmean(com_sub(MNH_idx,:)),'color',m_col)
%plot(rates,nanmean(com_sub(ONH_idx,:)),'color',o_col)
hold on
for ii=1:4
    eb(1,ii)=errorbar(rates(ii)-.01,nanmean(com_sub(YNH_idx,ii)),nanstd(com_sub(YNH_idx,ii))/sqrt(length(YNH_idx)),'^','color',y_col,'markerfacecolor',y_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    eb(2,ii)=errorbar(rates(ii),nanmean(com_sub(MNH_idx,ii)),nanstd(com_sub(MNH_idx,ii))/sqrt(length(MNH_idx)),'sq','color',m_col,'markerfacecolor',m_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    eb(3,ii)=errorbar(rates(ii)+.01,nanmean(com_sub(ONH_idx,ii)),nanstd(com_sub(ONH_idx,ii))/sqrt(length(ONH_idx)),'o','color',o_col,'markerfacecolor',o_col,'MarkerEdgecolor','k')

end

ylabel('P2-N1 \muV')
xlabel('ISI (s)')
set(gca,'fontsize',12)
set(gcf,'position',[714 8 270 266])
set(gca,'xtick',rates)
xlim([0.5 2.5])
ylim([-1.5 8.5])
fig=gcf

box off

% get slope and intercept from p2-n1]
for ii=1:length(com_sub)
    pfit(ii,:) = polyfit(rates',com_sub(ii,:)',1);
    ffit(ii,:) = polyval(pfit(ii,:),rates);
end

plot(rates,nanmean(ffit(YNH_idx,:)),'--','color',y_col)
hold on
plot(rates,nanmean(ffit(MNH_idx,:)),'--','color',m_col)
plot(rates,nanmean(ffit(ONH_idx,:)),'--','color',o_col)
set(gca,'xtick',rates)
hleg = legend([eb(1,1) eb(2,2) eb(3,3)],{'Young','Middle-aged','Older'})
hleg.Box = 'off'
hleg.FontSize = 10;
%hleg.FontName = 'Arial'
hleg.Position = [0.1327 0.7072 0.4963 0.2105]

fig = gcf;
saveas(fig,[savepath 'P2N1_all_fix'],'svg')
    % for stats
%  for ii=1:length(com_sub)
%     pfit(ii,:) = polyfit([0:3],com_mean_sub(ii,:)',1);
%     ffit(ii,:) = polyval(pfit(ii,:),[0:3]);
%  end

 %figure
 %scatter(age,pfit(:,2))

% P50
figure('Renderer','painter')
% Young
plot(rates,nanmean(p1_sub(YNH_idx,:)),'color',y_col)
hold on
plot(rates,nanmean(p1_sub(MNH_idx,:)),'color',m_col)
plot(rates,nanmean(p1_sub(ONH_idx,:)),'-','color',o_col)
for ii=1:4
    eb(2,ii)=errorbar(rates(ii)-.01,nanmean(p1_sub(YNH_idx,ii)),nanstd(p1_sub(YNH_idx,ii))/sqrt(length(YNH_idx)),'^','color',y_col,'markerfacecolor',y_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    eb(2,ii)=errorbar(rates(ii),nanmean(p1_sub(MNH_idx,ii)),nanstd(p1_sub(MNH_idx,ii))/sqrt(length(MNH_idx)),'sq','color',m_col,'markerfacecolor',m_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    eb(1,ii)=errorbar(rates(ii)+.01,nanmean(p1_sub(ONH_idx,ii)),nanstd(p1_sub(ONH_idx,ii))/sqrt(length(ONH_idx)),'o','color',o_col,'markerfacecolor',o_col,'MarkerEdgecolor','k')

end


ylabel('P1 \muV')
xlabel('ISI (s)')
set(gca,'fontsize',12)
set(gcf,'position',[438 10 270 266])
set(gca,'xtick',rates)
xlim([0.5 2.5])
ylim([-0.5 3])
fig=gcf
box off
saveas(fig,[savepath 'p50_age_group_fix'],'svg')
%% plot stim
%close all
figure('renderer','painters')

fs_model = 10000;
subplot(4,2,[2 4])
im = linspace(0.8,2.3,4);
im = flip(kron(im,ones(1,12)))
im = [im flip(im)]
im = repmat(im,1,6)
plot(im,'-k')
xlabel('Trial nr.')
ylabel('ISI(s)')
ylim([0.25 2.5])
xlim([0 48*2])
box off
set(gca,'fontsize',12,'ytick',[1 2])%,'xtick',[1 48*2 48*4])
%set(gcf,'position',[680 855 250 168])
subplot(4,2,8)
%figure('renderer','painters')
stim = 1;
stim_all = [];
for ii=1:length(im)
    stim_all = [stim_all stim zeros(1,im(ii)*fs_model)];
end
t=0:1/fs_model:length(stim_all)/fs_model-1/fs_model;
plot(t/60,stim_all,'k')
set(gca,'fontsize',12)
ax=gca;
ax.YAxis.Visible = 'Off'
xlabel('Time (min)')
box off
xlim([0 149/60])
ylim([0 .7])
%set(gcf,'position',[680 855 250 207])
%set(gcf,'position',[0 0 297 258]);
%set(gcf,'position',[381 231 186 234])
%set(gcf,'position',[[440 433 409 273]])
%set(gcf,'position',[381 231 186 234])
set(gcf,'position',[440 329 544 304])
fig = gcf;
saveas(fig,[savepath 'stim_paradigm_new'],'svg')
 %% age corr P2-N1
 %% AEP P2-N1
var = 'AEP_p2n1_int';
labels = 'P2-N1 intercept'
lims = [-7 9];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
set(gcf,'position',[971 23 396 297])
fig = gcf;
saveas(fig,[savepath 'P2N1_int_fix'],'svg')

var = 'AEP_p2n1_slope';
labels = 'P2-N1 slope'
lims = [-4 12];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
set(gcf,'position',[971 23 396 297])
fig = gcf;
saveas(fig,[savepath 'P2N1_slope_fix'],'svg')
%% post hoc analysis
clc
disp('ISI=2.3 latencies')
disp(['mean P1 latency: ' num2str(round(nanmean(p1_sub_lat(idx,4)),3)) '+/- ' num2str(round(nanstd(p1_sub_lat(idx,4)),3)) ' s'])
disp(['mean N1 latency: ' num2str(round(nanmean(n1_sub_lat(idx,4)),3)) '+/- ' num2str(round(nanstd(n1_sub_lat(idx,4)),3)) ' s'])
disp(['mean P2 latency: ' num2str(round(nanmean(p2_sub_lat(idx,4)),3)) '+/- ' num2str(round(nanstd(p2_sub_lat(idx,4)),3)) ' s'])
disp(['mean N2 latency: ' num2str(round(nanmean(n2_sub_lat(idx,4)),3)) '+/- ' num2str(round(nanstd(n2_sub_lat(idx,4)),3)) ' s'])

% N2 for groups
idx_all = {y,m,o};
rates = {'0.8','1.3','1.8','2.3'};
gg_t = {'y','ma','o'};
disp('N2 latency for groups:')
for gg = 1:3
    for ii = 1:4
        disp(['mean N2 latency ISI ' rates{ii} 's: ' gg_t{gg} ': ' num2str(round(nanmean(n2_sub_lat(idx_all{gg},ii)),3)) '+/- ' num2str(round(nanstd(n2_sub_lat(idx_all{gg},ii)),3)) ' s'])
    end
end
% post hoc corr N1 @2.3s
this_idx = ~isnan(age') & ~isnan(n1_fix_sub(:,4));
this_idx = find(this_idx & idx);
[rhon1,pn1] =corr(zscore(age(this_idx))',zscore(n1_fix_sub(this_idx,4)),'type','spearman');

% post hoc corr P2 @0.8s
this_idx = ~isnan(age') & ~isnan(p2_fix_sub(:,1));
this_idx = find(this_idx & idx);
[rhop2,pp2] =corr(zscore(age(this_idx))',zscore(p2_fix_sub(this_idx,1)),'type','spearman');

% post hoc corr P50 mean
this_idx = ~isnan(age') & ~isnan(mean(p2_sub(:,:),2));
this_idx = find(this_idx & idx);
[rhop1m,pp1m] =corr(zscore(age(this_idx))',zscore(mean(p1_sub(this_idx,:),2)),'type','spearman');

% post hoc corr P50 mean over ISI

for ii=1:4
this_idx = ~isnan(age') & ~isnan(p1_sub(:,ii));
this_idx = find(this_idx & idx);
[rhop1(ii),pp1(ii)] =corr(zscore(age(this_idx))',zscore(p1_sub(this_idx,ii)),'type','spearman');
end
disp(['Age corr uncorrected N1@2.3 s : r = ' num2str(round(rhon1,3)) ' , p = ' num2str(round(pn1,4))])
disp(['Age corr uncorrected P2@0.8 s : r = ' num2str(round(rhop2,3)) ' , p = ' num2str(round(pp2,4))])
for ii=1:4
disp(['Age corr uncorrected P50 ' rates{ii} ' s : r = ' num2str(round(rhop1(ii),3)) ' , p = ' num2str(round(pp1(ii),4))])
end
disp(['Age corr uncorrected P50 mean : r = ' num2str(round(rhop1m,3)) ' , p = ' num2str(round(pp1m,4))])


% post hoc N2 latency
% average latencies for N2 

% average latencies for N2 for groups + corr with age