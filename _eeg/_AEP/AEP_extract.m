% plot abr and extract peaks
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
d = dir('_outputs/_derivatives/*.mat')

for s=1:length(d)
    load([d(s).folder filesep d(s).name])
    extract_aep_data
    clc
    disp([sub_id{s} ' done...'])
end


%% plot
load('/work1/jonmarc/UHEAL_master/UHEAL/uheal_data.mat')

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
    xlabel('Time (s)')
    box off
    set(gca,'Fontsize',12)
    xlim([-.1 .5])
    ylim([-5.5 3])
end

%set(gcf,'position',[144 267 1315 368])
set(gcf,'position',[85 364 958 272]);
hleg = legend([pmean_y(1).mainLine pmean_m(1).mainLine pmean_o(1).mainLine],{'Young','Middle-aged','Older'});
hleg.Box = 'off'
hleg.Position = [0.1310 0.2386 0.1491 0.2059];hleg.FontSize = 10;hleg.FontName = 'Arial'

%% plot N1,P2 and P2-N1
com_mean_sub =  p2_mean_sub-n1_mean_sub;

%close all
% N100
figure('Renderer','painter')
% plot mean n1
plot(rates,nanmean(n1_mean_sub(YNH_idx,:)),'color',y_col)
hold on
plot(rates,nanmean(n1_mean_sub(MNH_idx,:)),'color',m_col)
plot(rates,nanmean(n1_mean_sub(ONH_idx,:)),'color',o_col)
% plot n1 for age groups
for ii=1:4
    %young
    eb(2,ii)=errorbar(rates(ii)-.01,nanmean(n1_mean_sub(YNH_idx,ii)),nanstd(n1_mean_sub(YNH_idx,ii))/sqrt(length(YNH_idx)),'ksq','color',y_col,'markerfacecolor',y_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    %mid aged
    eb(3,ii)=errorbar(rates(ii),nanmean(n1_mean_sub(MNH_idx,ii)),nanstd(n1_mean_sub(MNH_idx,ii))/sqrt(length(MNH_idx)),'k^','color',m_col,'markerfacecolor',m_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    %older
    eb(1,ii)=errorbar(rates(ii)+.01,nanmean(n1_mean_sub(ONH_idx,ii)),nanstd(n1_mean_sub(ONH_idx,ii))/sqrt(length(ONH_idx)),'ko','color',o_col,'markerfacecolor',o_col,'MarkerEdgecolor','k')
end

ylabel('N100 (\muV')
xlabel('ISI (s)')
set(gca,'fontsize',12)
set(gcf,'position',[163 12 270 266])
hleg = legend([eb(2) eb(3) eb(1)],'Young','Middle-aged','Older')
hleg.Position = [0.1623 0.2148 0.5185 0.2218];
xlim([0.5 2.5])
ylim([-6 2])
set(gca,'xtick',rates)
hleg.Box = 'off'
hleg.FontSize = 10;
hleg.FontName = 'Arial'
box off
fig=gcf
%saveas(fig,'figs_paper/n100_age_group','epsc')

% P2
figure('Renderer','painter')
% Young
plot(rates,nanmean(p2_mean_sub(YNH_idx,:)),'color',y_col)
hold on
plot(rates,nanmean(p2_mean_sub(MNH_idx,:)),'color',m_col)
plot(rates,nanmean(p2_mean_sub(ONH_idx,:)),'-','color',o_col)
for ii=1:4
    eb(2,ii)=errorbar(rates(ii)-.01,nanmean(p2_mean_sub(YNH_idx,ii)),nanstd(p2_mean_sub(YNH_idx,ii))/sqrt(length(YNH_idx)),'^','color',y_col,'markerfacecolor',y_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    eb(2,ii)=errorbar(rates(ii),nanmean(p2_mean_sub(MNH_idx,ii)),nanstd(p2_mean_sub(MNH_idx,ii))/sqrt(length(MNH_idx)),'sq','color',m_col,'markerfacecolor',m_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    eb(1,ii)=errorbar(rates(ii)+.01,nanmean(p2_mean_sub(ONH_idx,ii)),nanstd(p2_mean_sub(ONH_idx,ii))/sqrt(length(ONH_idx)),'o','color',o_col,'markerfacecolor',o_col,'MarkerEdgecolor','k')

end


ylabel('P200 (\muV)')
xlabel('ISI (s)')
set(gca,'fontsize',12)
set(gcf,'position',[438 10 270 266])
set(gca,'xtick',rates)
xlim([0.5 2.5])
ylim([-1 3])
fig=gcf
box off
% saveas(fig,'figs_paper/p200_age_group','epsc')

% P2-N1 complex
% P2
figure('Renderer','painters')
% Young
plot(rates,nanmean(com_mean_sub(YNH_idx,:)),'color',y_col)
hold on
plot(rates,nanmean(com_mean_sub(MNH_idx,:)),'color',m_col)
plot(rates,nanmean(com_mean_sub(ONH_idx,:)),'color',o_col)

for ii=1:4
    eb(1,ii)=errorbar(rates(ii)-.01,nanmean(com_mean_sub(YNH_idx,ii)),nanstd(com_mean_sub(YNH_idx,ii))/sqrt(length(YNH_idx)),'^','color',y_col,'markerfacecolor',y_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    eb(2,ii)=errorbar(rates(ii),nanmean(com_mean_sub(MNH_idx,ii)),nanstd(com_mean_sub(MNH_idx,ii))/sqrt(length(MNH_idx)),'sq','color',m_col,'markerfacecolor',m_col,'MarkerEdgecolor','k','MarkerSize',6.8)
    eb(3,ii)=errorbar(rates(ii)+.01,nanmean(com_mean_sub(ONH_idx,ii)),nanstd(com_mean_sub(ONH_idx,ii))/sqrt(length(ONH_idx)),'o','color',o_col,'markerfacecolor',o_col,'MarkerEdgecolor','k')

end

ylabel('P200-N100 (\muV)')
xlabel('ISI (s)')
set(gca,'fontsize',12)
set(gcf,'position',[714 8 270 266])
set(gca,'xtick',rates)
xlim([0.5 2.5])
ylim([-2 7])
fig=gcf

box off

% get slope and intercept from p2-n1]
for ii=1:length(com_mean_sub)
    pfit(ii,:) = polyfit(rates',com_mean_sub(ii,:)',1);
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
hleg.FontName = 'Arial'
hleg.Position = [0.1327 0.7072 0.4963 0.2105]

fig = gcf;

    % for stats
 for ii=1:length(com_mean_sub)
    pfit(ii,:) = polyfit([0:3],com_mean_sub(ii,:)',1);
    ffit(ii,:) = polyval(pfit(ii,:),[0:3]);
 end

 %figure
 %scatter(age,pfit(:,2))

 %% extract measures
