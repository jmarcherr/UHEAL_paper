close all
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat')

data = uheal_data;
%% age distribution
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


mac_screen =  [1           1        1440         900]
%% PTA hf
var = 'PTA_hf';
labels = 'PTA_{hf} (dB HL)'
lims = [-40 120];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')



%% PTA lf
var = 'PTA_lf';
labels = 'PTA_{lf} (dB HL)'
lims = [-20 25];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
ax1.YTick = [-10 0 10 20 30];ax2.YTick = [-10:10:30];
%set(hleg,'Visible','off')
%set(gcf,'position',[552 636 325 169])

%% MEMR
var = 'memr_slope';
labels = 'MEMR growth'
lims = [-0.015 0.06];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.6 0.8 0.01 0.1];hleg.Location = 'best'
ax1.YTick = [0 0.02 0.04 0.06];ax1.YTickLabels = {'0','.02','.04','.06'}
ax2.YTick = [0 0.02 0.04 0.06];ax2.YTickLabels = {'0','.02','.04','.06'}%
%hleg.Box = 'on'
%set(hleg,'Visible','off')

%% FFR
close all
data_ffr = uheal_data;
data_ffr.FFR_SNR(find(data_ffr.FFR_sig==0))=nan;
var = 'FFR_SNR'
labels = 'FFR_{SNR} dB'
lims = [-10 75]
figure

[h,ax1,ax2]=stat_plots_uh_age_pt(data_ffr,var,labels,lims)

scatter(uheal_data.Age(find(uheal_data.FFR_sig==0 & uheal_data.gender==1)),uheal_data.FFR_SNR(find(uheal_data.FFR_sig==0 & uheal_data.gender==1)),'o','markeredgecolor',[0.5 0.5 0.5],'SizeData',10,'linewidth',.05)
scatter(uheal_data.Age(find(uheal_data.FFR_sig==0 & uheal_data.gender==2)),uheal_data.FFR_SNR(find(uheal_data.FFR_sig==0 & uheal_data.gender==2)),'o','markeredgecolor',[0.5 0.5 0.5],'SizeData',10,'linewidth',.05)
set(gcf,'position',[100 100 297 199]*(96/72))
%set(hleg,'Visible','off')

%% ABR measures
close all
% AP
var = 'AP_amp_pm';
labels = '\muV'
lims = [-0.2 1.1];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')

% WV
var = 'WV_amp_pm';
labels = '\muV'
lims = [-0.2 1.1];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')

% SP
var = 'SP_amp';
labels = '\muV'
lims = [-0.19 0.29];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')

%% clinical measures
close all
% digit span
var = 'rds';
labels = 'R. digit-span score'
lims = [0 30];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
% SSQ
var = 'ssq12_mean';
labels = 'SSQ-12 score'
lims = [0 12];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
% NESI log 10
uheal_data.nesi(isinf(uheal_data.nesi))=nan; data=uheal_data;
data.nesi(find(~isnan(data.nesi))) = log10(data.nesi(find(~isnan(data.nesi))));
var = 'nesi';
labels = 'log_{10}(NESI)'
lims = [-2 3];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
% TTS
var = 'tts';
labels = 'TTS score'
lims = [-1 8];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')

%% cortical EEG
close all
% 4Hz neg
var = 'Neg_4Hz';
labels = '\muV'
lims = [-3 2];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
%set(gcf,'position',[462 556 297 197])
set(ax1,'FontSize',11);set(ax2,'FontSize',11)

%4Hz itpc_ratio
var = 'ITPC_ratio';
labels = 'ITPC ratio'
lims = [-1 1];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
%set(gcf,'position',[462 556 297 197])
set(ax1,'FontSize',11);set(ax2,'FontSize',11)
%% cortical AEP
close all
% n100
var = 'AEP_n1';
labels = 'N100 \muV'
lims = [-16 6];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
%set(gcf,'position',[462 556 297 197])
%set(ax1,'FontSize',11);set(ax2,'FontSize',11)

%p200
var = 'AEP_p2';
labels = 'P200 \muV'
lims = [-4 7];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
%set(gcf,'position',[462 556 297 197])
%set(ax1,'FontSize',11);set(ax2,'FontSize',11)

%% AEP P2-N1
var = 'AEP_p2n1_int';
labels = 'P2-N1 intercept'
lims = [-8 8];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
%set(gcf,'position',[100 100 297 197]*(96/72))

%%
clear tmp_data
close all
clear p
nanidx = ~isnan(data.AEP_p2(:,1));
tmp_data.AG = repmat([0],size(data.Age))
tmp_data.AG(y,:) = ones(size(y));
tmp_data.AG(m,:) = ones(size(m))*2;
tmp_data.AG(o,:) = ones(size(o))*3;
tmp_data.AG = repmat(tmp_data.AG(nanidx,:),4,1);
tmp_data.Age = repmat(data.Age(nanidx),4,1)
tmp_data.gender = repmat(data.gender(nanidx),4,1)
tmp_data.subid = repmat(data.subid(nanidx),4,1);
tmp_data.PTAlf = repmat(data.PTA_lf(nanidx),4,1);
n1p2_cmp = ([data.AEP_p2(nanidx,:)-data.AEP_n1(nanidx,:)])

tmp_data.n1p2_cmp = n1p2_cmp(:);
isize = size(find(nanidx));

tmp_data.isi = [zeros(isize);ones(isize);ones(isize)*2;ones(isize)*3];   
tmp_table = struct2table(tmp_data')
%writetable(tmp_table,'/Users/jmarcher/Documents/_postdoc/UHEAL_paper/_statsR/isi_AEP/isi_table.csv')
%%
for ii=1:length(unique(tmp_data.subid))
p(ii,:) = polyfit(tmp_data.isi(tmp_data.subid==tmp_data.subid(ii)),tmp_data.n1p2_cmp(tmp_data.subid==tmp_data.subid(ii)),1);
f(ii,:) = polyval(p(ii,:),[0:3])
end
% p(:,2) = intercept

%aoctool(x,y,group,alpha,xname,yname,gname)
[h,atab,ctab,stats] = aoctool(tmp_data.isi,tmp_data.n1p2_cmp,tmp_data.subid,0.05,'isi','n1p2_cmp','subid','off');
data.n1p2_int = nan(size(data.n1p2_cmp5))
data.n1p2_int(nanidx,:) = p(:,2)
data.n1p2_slope(nanidx,:) = p(:,1);
data.isifit = nan(size(data.AEP_n1))
data.isifit(nanidx,:) = f;


var = ['n1p2_int'];
labels = 'P2-N1 intercept';
lims = [-4 10];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
set(gcf,'position',[462 556 297 197])
%ylim([-4 20])

var = ['n1p2_slope'];
labels = 'P2-N1 slope';
lims = [-4 10];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
set(gcf,'position',[462 556 297 197])
%ylim([-4 20])

figure
plot([1:4],nanmean(data.AEP_p2(y,:)-data.AEP_n1(y,:),1),'ko-')
hold on
plot([1:4],nanmean(data.isifit(y,:),1),'--k')
plot([1:4],nanmean(data.AEP_p2(m,:)-data.AEP_n1(m,:),1),'bo-')
plot([1:4],nanmean(data.isifit(m,:),1),'--b')
plot([1:4],nanmean(data.AEP_p2(o,:)-data.AEP_n1(o,:),1),'ro-')
plot([1:4],nanmean(data.isifit(o,:),1),'--r')

%% EFR
data_efr = uheal_data;
var = 'EFR_SNR';
data_efr.EFR_SNR(find(data_efr.EFR_sig==0))=nan;
labels = 'EFR_{SNR} dB'
lims = [-20 50]
figure

[h,ax1,ax2]=stat_plots_uh_age_pt(data_efr,var,labels,lims)

scatter(uheal_data.Age(find(uheal_data.EFR_sig==0 & uheal_data.gender==1)),uheal_data.EFR_SNR(find(uheal_data.EFR_sig==0 & uheal_data.gender==1)),'o','markeredgecolor',[0.5 0.5 0.5],'SizeData',10,'linewidth',.1)
scatter(uheal_data.Age(find(uheal_data.EFR_sig==0 & uheal_data.gender==2)),uheal_data.EFR_SNR(find(uheal_data.EFR_sig==0 & uheal_data.gender==2)),'o','markeredgecolor',[0.5 0.5 0.5],'SizeData',10,'linewidth',.1)
set(gcf,'position',[100 100 297 199]*(96/72))
%set(hleg,'Visible','off')



