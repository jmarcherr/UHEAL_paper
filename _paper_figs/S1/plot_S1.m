% script for plotting clincial measures 
% audiogram
% MEMR
% TEOAE
% Questionaires and clincial outcomes
close all
clear all
cd('/work3/jonmarc/UHEAL_paper')
UHEAL_startup
addpath('/work3/jonmarc/UHEAL_paper/_stats/_corr')

mac_screen =  [1           1        1440         900];
%% init colormap
ages = [18 77];
        % new color map
        cmap_all = [];
        % young
        cmap = cbrewer('seq','Greys',(25-17)+5);
        cmap(1:5,:) = [];
        cmap_all = [cmap_all; cmap];
        % middle aged
        cmap = cbrewer('seq','Blues',(50-25)+5);
        cmap(1:5,:) = [];
        cmap_all = [cmap_all ; cmap];
        % older
        cmap = cbrewer('seq','Reds',(max(ages)-(50))+5);
        cmap(1:5,:) = [];
        cmap_all = [cmap_all ; cmap];

        % combine
        cmap = cmap_all;
        cmap(find(cmap>1))=1;cmap(find(cmap<0))=0;
        ageidx = linspace(min(ages),max(ages),size(cmap,1));
        
        % group colors
        y_col = cmap(13-5,:);
        m_col = cmap(30,:);
        o_col = cmap(end-7,:);
%% load table
load('/work3/jonmarc/UHEAL_paper/_clin/clin_data_table/clin_data.mat')
% old groups
%y_idx = find(uheal_data.Age<=25 & ~uheal_data.CP_new); 
%m_idx = find(uheal_data.Age>26 & uheal_data.Age<50 & ~uheal_data.CP_new)
%o_idx = find(uheal_data.Age>=50 & ~uheal_data.CP_new);
% new CP subgroups
CP_X = (mean(uheal_data.aud(:,1:6),2)<=20 & all(uheal_data.aud(:,1:6)<=15,2) & all(uheal_data.aud(:,1:6)>=-5,2));
uheal_data.CP_X =CP_X;
% new groups
y_idx = find(uheal_data.Age<=25 & uheal_data.CP_X); 
m_idx = find(uheal_data.Age>25 & uheal_data.Age<50 & uheal_data.CP_X)
o_idx = find(uheal_data.Age>=50 & uheal_data.CP_X);

uheal_data.CP_new = ~CP_X;
clc
y = find(uheal_data.Age<=25 & ~uheal_data.CP_new);
disp(['YOUNG, n=' num2str(length(y)) ': mean age:' num2str(round(mean(uheal_data.Age(y)),1)) ' +/- ' num2str(round(std(uheal_data.Age(y)),1)) ...
    ', ' num2str(length(find(uheal_data.gender(y)==2))) ' Males'])
m = find(uheal_data.Age>25 & uheal_data.Age<50 & ~uheal_data.CP_new);
disp(['M AGED, n=' num2str(length(m)) ': mean age:' num2str(round(mean(uheal_data.Age(m)),1)) ' +/- ' num2str(round(std(uheal_data.Age(m)),1))...
    ', ' num2str(length(find(uheal_data.gender(m)==2))) ' Males'])
o = find(uheal_data.Age>=50 & ~uheal_data.CP_new);
disp(['OLDER, n=' num2str(length(o)) ': mean age:' num2str(round(mean(uheal_data.Age(o)),1)) ' +/- ' num2str(round(std(uheal_data.Age(o)),1))...
    ', ' num2str(length(find(uheal_data.gender(o)==2))) ' Males'])
%% Audiogram plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots audiograms vs. age as well as age-distribution in various ways. Can
% be left out if only extraction of data is wanted (plot_on=0, save_figs=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for plotting function
age_sub = uheal_data.Age;
aud = uheal_data.aud;
freq_aud = uheal_data.audfreq;
gender_sub = uheal_data.gender;
savepath = '/work3/jonmarc/UHEAL_paper/_paper_figs/S1/figs/'
save_figs = 1;

%% Audiogram NH
close all
NH_idx = find(uheal_data.CP_X==1);
fig = plot_audiogram_groups(NH_idx,uheal_data.Age,uheal_data.aud,uheal_data.audfreq,cmap)
% save figure
%title(['NH, n=' num2str(length(NH_idx))])
if save_figs
    saveas(fig,[savepath 'audiograms_NH'],'svg')
end

hold on
py=plot(uheal_data.audfreq(1,:),nanmean(uheal_data.aud(y_idx,:)),'^','color',y_col,'linewidth',2)%,'markerfacecolor',[cmap(13-5,:)])
pm=plot(uheal_data.audfreq(1,:),nanmean(uheal_data.aud(m_idx,:)),'sq','color',m_col,'linewidth',2)%,'markerfacecolor',[cmap(30,:) 0.5])
po=plot(uheal_data.audfreq(1,:),nanmean(uheal_data.aud(o_idx,:)),'o','color',o_col,'linewidth',2)%,'markerfacecolor',[cmap(end-7,:) 0.5])
hleg=legend([py,pm,po],{'Young','Middle-aged','Older'});
hleg.Box='off';
hleg.Position = [0.2300 0.2887 0.3935 0.2475];
xtickangle(0)
set(gca,'xtick',[250,500,1000,2000,4000,8000,16000],'xticklabel',{'.25','.5','1','2','4','8','16'},'FontSize',12)
%set(gcf,'position',[[305 412 478 299]])
set(gcf,'position',[305 473 422 225])
if save_figs
    saveas(fig,[savepath 'audiograms_NH_groups'],'svg')
end

%% PTA hf
figure
var = 'PTA_hf';
labels = 'PTA_{hf} (dB HL)'
lims = [-40 120];
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
fig = gcf;
if save_figs
    saveas(fig,[savepath 'PTAHF_NH_groups'],'svg')
end


%% PTA lf
var = 'PTA_lf';
labels = 'PTA_{lf} (dB HL)'
lims = [-20 25];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
ax1.YTick = [-10 0 10 20 30];ax2.YTick = [-10:10:30];
%set(hleg,'Visible','off')
%set(gcf,'position',[552 636 325 169])
fig = gcf;
if save_figs
    saveas(fig,[savepath 'PTALF_NH_groups'],'svg')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% CORTICAL
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat')
% old groups
%y_idx = find(uheal_data.Age<=25 & ~uheal_data.CP_new); 
%m_idx = find(uheal_data.Age>26 & uheal_data.Age<50 & ~uheal_data.CP_new)
%o_idx = find(uheal_data.Age>=50 & ~uheal_data.CP_new);
% new CP subgroups
CP_X = (mean(uheal_data.aud(:,1:6),2)<=20 & all(uheal_data.aud(:,1:6)<=15,2) & all(uheal_data.aud(:,1:6)>=-5,2));
uheal_data.CP_X =CP_X;
% new groups
y_idx = find(uheal_data.Age<=25 & uheal_data.CP_X); 
m_idx = find(uheal_data.Age>25 & uheal_data.Age<50 & uheal_data.CP_X)
o_idx = find(uheal_data.Age>=50 & uheal_data.CP_X);

uheal_data.CP_new = ~CP_X;

% 4Hz neg
var = 'Neg_4Hz';
labels = '\muV'
lims = [-3 2];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
%set(gcf,'position',[462 556 297 197])
%set(ax1,'FontSize',11);set(ax2,'FontSize',11)
%set(gcf,'position',[229 226 396 318])
fig = gcf;
saveas(fig,[savepath 'Neg_age'],'svg')

%4Hz itpc_ratio
var = 'ITPC_ratio';
labels = 'ITPC ratio'
lims = [-1 1];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
%set(gcf,'position',[462 556 297 197])
%set(ax1,'FontSize',11);set(ax2,'FontSize',11)
%set(gcf,'position',[229 226 396 318])
fig = gcf;
saveas(fig,[savepath 'itpc_age'],'svg')

% AEP int
var = 'AEP_p2n1_int';
labels = 'P2-N1 intercept'
lims = [-8 8];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
%set(gcf,'position',[971 23 396 297])
fig = gcf;
saveas(fig,[savepath 'P2N1_int'],'svg')
%% TEOAE
figure('Renderer','painters')
        marks_ag = {'^','sq','o'}
        col_ag = {'k','b','r'}
        col_ag = {y_col,m_col,o_col}
        ag_idx = {y_idx,m_idx,o_idx};
        levels = [1000:1000:5000];
        for ag=1:3
            p_y = semilogx(levels,nanmean(uheal_data.teoae_SNR(ag_idx{ag},:),1),'-','color',col_ag{ag},'marker',marks_ag{ag},'markeredgecolor',col_ag{ag},'markerfacecolor','w');
            hold on
            for ll=1:5

                errorbar(levels(ll),nanmean(uheal_data.teoae_SNR(ag_idx{ag},ll),1),nanstd(uheal_data.teoae_SNR(ag_idx{ag},ll),1)./sqrt(length(ag_idx{ag})),'color',col_ag{ag})
                
                eby(ag) = semilogx(levels(ll),nanmean(uheal_data.teoae_SNR(ag_idx{ag},ll),1),'marker',marks_ag{ag},'color',col_ag{ag},'MarkerFaceColor',col_ag{ag});
            end

        end
        set(gca,'fontsize',12)
        set(gcf,'position',[275 161 353 205])
        plot([0 6000],[6000 6000],'--','color',[.75 .75 .75])
        set(gca,'xtick',[levels],'xticklabel',{'1','2','3','4','5'},'FontSize',12,'ytick',[5:5:20])
        
xlim([800 5500])
ylim([3 25])
box off

tmp = gca;
xlabel('Frequency (kHz)','FontName','Arial')
ylabel('SNR (dB)')
set(gca,'fontsize',12)

fig = gcf;
if save_figs
    saveas(fig,[savepath 'TEOAE_NH_groups'],'svg')
end

% alternative size
set(gcf,'position',[275 182 239 184])
fig = gcf;
if save_figs
    saveas(fig,[savepath 'TEOAE_NH_groups_small'],'svg')
end
%% test fig
figure('Renderer','painters')
set(gcf,'position',[275 161 353 205])
plot(1:100,1:100)
set(gca,'fontsize',12)
xlabel('Frequency (kHz)')
ylabel('Hearing Level (dB)')
fig = gcf;
if save_figs
    saveas(fig,[savepath 'test_fig'],'svg')
end
%% MEMR plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots average MEMR data and linear coefficient as function of age. Can
% be left out if only extraction of data is wanted (plot_on=0, save_figs=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for plotting
reflex_sub = uheal_data.memr_reflex;
growth_sub_alt = uheal_data.memr_reflex_growth;
levels = uheal_data.memr_levels;
f_center = uheal_data.memr_fcenter;
age_sub = uheal_data.Age;
MEM_slope = uheal_data.memr_slope;



%close all
figure('Renderer','painter')
subplot(1,3,2)
cmap = cbrewer('seq','YlGnBu',7+2);
cmap = cmap(2:end,:);

Yidx = y_idx;%find(age_sub<=25 & ~uheal_data.CP_new');
Midx = m_idx;%find(age_sub>25 & age_sub<60 & ~uheal_data.CP_new');
Oidx = o_idx;%find(age_sub>=60 & ~uheal_data.CP_new');
NH_idx = find(~uheal_data.CP_new');
% mean
MEMR_mean = squeeze(nanmean(reflex_sub,1));
% growth curve

marks_ag = {'^','sq','o'}
col_ag = {'k','b','r'}
col_ag = {y_col,m_col,o_col}
ag_idx = {Yidx,Midx,Oidx};
for ag=1:3
    p_y = plot(levels(1:end-1),nanmean(growth_sub_alt(ag_idx{ag},1:end-1),1),'-','color',col_ag{ag},'marker',marks_ag{ag},'markeredgecolor',col_ag{ag},'markerfacecolor','w');
    hold on
    for ll=1:length(levels)-1

        errorbar(levels(ll),nanmean(growth_sub_alt(ag_idx{ag},ll),1),nanstd(growth_sub_alt(ag_idx{ag},ll),1)/sqrt(length(ag_idx{ag})),'color',col_ag{ag})
        eby(ag) = plot(levels(ll),nanmean(growth_sub_alt(ag_idx{ag},ll),1),'marker',marks_ag{ag},'color',col_ag{ag},'MarkerFaceColor',col_ag{ag});
    end

end

hold on
set(gca,'fontsize',12)
box off
set(gcf,'position',[228 839/2 432/2 209])
set(gca,'Xtick',levels,'Xticklabels',{'','80','','90','','100',''},'ytick',[0,.2,.4,.6],'yticklabels',{'0','.2','.4','.6'})
xlim([72 108])
ylim([0 .8])
xlabel('Elicitor level (dB)')
ylabel('\Sigma |\Delta Absorbance|')
hlegg=legend([eby(1) eby(2) eby(3)],{'Young','Middle-aged','Older'})
hlegg.Box = 'off';
hlegg.Position = [0.4307 0.6877 0.2265 0.2823];
xtickangle(0);
% mean results
%figure(2)
subplot(1,3,1)
for i=1:length(levels-1)
    % Plot
    semilogx(f_center, nanmean(reflex_sub(NH_idx,:,i),1),'color',cmap(i,:));
    axis([250 8000 -0.08 0.099])
    a=gca;
    a.XTick = [250,500,1000,2000,4000,8000];
    a.XTickLabel = [{'.25'},{'.5'},{'1'},{'2'},{'4'},{'8'}];
    hold on
    set(gca,'fontsize',12);
    box off
    for x = 1:length(levels)-1
        txt(x) = {[num2str(levels(x)),' dB SPL']};
    end
    set(gcf,'position',[228 839/2 432 209]);
    hleg = legend(txt);
    hleg.Box = 'off';
    hleg.Position = [0.5547 0.1231 0.3704 0.8038];
    xlabel('Frequency (kHz)');
    ylabel('\Delta Absorbance');
end
set(gcf,'Position',[228 420 618 225]);
hleg.Position= [0.6509 0.1470 0.2589 0.8038];
xtickangle(0)
fig = gcf;
if save_figs
    saveas(fig,[savepath 'MEMR_NH_groups'],'svg')
end

% MEMR slope

var = 'memr_slope';
labels = 'MEMR growth'
lims = [-0.015 0.06];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
%hleg.Position = [0.6 0.8 0.01 0.1];hleg.Location = 'best'
ax1.YTick = [0 0.02 0.04 0.06];ax1.YTickLabels = {'0','.02','.04','.06'}
ax2.YTick = [0 0.02 0.04 0.06];ax2.YTickLabels = {'0','.02','.04','.06'}%
%hleg.Box = 'on'
%set(hleg,'Visible','off')
fig = gcf;
if save_figs
    saveas(fig,[savepath 'MEMR_slope_NH_groups'],'svg')
end


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% plot audiogram function
function fig = plot_audiogram_groups(idx,age_sub,aud,aud_freq,cmap)
    figure('renderer','painter')
for s=1:length(idx)

            p1 = semilogx(aud_freq(idx(s),:),aud(idx(s),:)','-','color',[cmap(age_sub(idx(s))-17,:) 0.8]);
            plot_aud_param(p1,aud_freq(idx(s),:));
            hold on
            
            cb=colorbar;
            
            cb.FontSize = 10;
            cb.Limits = [0 1];
            age_ticks = [20:10:70];
            cb.Ticks = (age_ticks-17)/(max(age_sub)-min(age_sub));%[linspace(0,1,6)];
            cb.TickLabels = age_ticks;%{linspace(18,max(age_sub),5)};
            cb.Label.String = 'Age';
            cb.Label.Rotation = 90;
            cb.Label.FontSize = 12;
            cb.Label.FontName = 'Arial';
            colormap(cmap)
end
set(gcf,'position',[305 412 432 299]);
fig = gcf
fig=gcf;
end