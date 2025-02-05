% plot abr and extract peaks
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
d = dir([rootdir '/_eeg/_ABR/_outputs/_derivatives/*.mat'])
load([rootdir '/_stats/uheal_data.mat'])
savepath = [rootdir '/_paper_figs/fig2/figs/'];
for dd=1:length(d)
    load([d(dd).folder filesep d(dd).name])
    if isfield(data,'abr')
    sub_abr(dd,:) = data.abr{1};
    t_abr(dd,:) = data.time{1};
    fs(dd) = data.fs;
    sub_peaks{dd} = data.abr_peaks{1};
    AP_amp_pm(dd) = data.abr_peaks{1}.AP_amp-data.abr_peaks{1}.AP_neg;
    WV_amp_pm(dd) = data.abr_peaks{1}.WV_amp-data.abr_peaks{1}.WV_neg;
    subid{dd} = data.subid;
    subinfo{dd} = data.subinfo;
    age(dd) = data.subinfo.age;
    CP(dd) = data.subinfo.CP;
    gender(dd) = data.subinfo.gender;
    rjt_sub(dd) = 0;
     else
    sub_abr(dd,:) = nan;
    t_abr(dd,:) = nan;
    fs(dd) = nan;
    sub_peaks{dd} = [];
    AP_amp_pm(dd) = nan;
    WV_amp_pm(dd) = nan;
    subid{dd} = data.subid;
    subinfo{dd} = data.subinfo;
    age(dd) = data.subinfo.age;
    gender(dd) = data.subinfo.gender;
    CP(dd) = data.subinfo.CP;
    rjt_sub(dd) = 1;
    end
    clc
    disp([subid{dd} ' done...'])
end




%% gather evertyhing and save
abr_data = struct;
abr_data.subid = uheal_data.subid;
abr_data.SP_amp = nan(size(uheal_data.subid));
%uheal_data.SP_lat = nan(size(uheal_data.subid));
abr_data.AP_amp = nan(size(uheal_data.subid));
abr_data.AP_amp_pm = nan(size(uheal_data.subid));
abr_data.AP_lat =  nan(size(uheal_data.subid));
abr_data.WV_amp = nan(size(uheal_data.subid));
abr_data.WV_amp_pm = nan(size(uheal_data.subid));
abr_data.WV_lat =  nan(size(uheal_data.subid));

for s=1:length(subid)
    % get this subid
    thisID = str2double(subid{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);
    if ~isempty(sub_peaks{s})
    abr_data.SP_amp(this_idx) = sub_peaks{s}.SP_amp;
    %uheal_data.SP_lat(this_idx) = SP_lat(s);
    abr_data.AP_amp(this_idx) = sub_peaks{s}.AP_amp;
    abr_data.AP_lat(this_idx) = sub_peaks{s}.AP_latency;
    abr_data.WV_amp(this_idx) = sub_peaks{s}.WV_amp;
    abr_data.WV_lat(this_idx) = sub_peaks{s}.WV_latency;
    abr_data.AP_amp_pm(this_idx) = AP_amp_pm(s);
    abr_data.WV_amp_pm(this_idx) = WV_amp_pm(s);
    end
end
%save('/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_outputs/abr_data_table/abr_data.mat','abr_data');

%% gather traces 
abr_sub_trace = struct;
abr_sub_trace.subid = uheal_data.subid;
abr_sub_trace.sub_abr_b = sub_abr;
abr_sub_trace.t_abr = t_abr';
abr_sub_trace.CP = uheal_data.CP_new;
abr_sub_trace.rjt_sub = rjt_sub'
abr_sub_trace.age  = age';
abr_sub_trace.gender = gender';


%save('/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_outputs/abr_data_table/abr_sub_trace','-struct','abr_sub_trace')

%% get groups

close all
figure('renderer','painters')
YNH_idx = find(abr_sub_trace.age<=25 & ~abr_sub_trace.CP & ~abr_sub_trace.rjt_sub); 
MANH_idx = find(abr_sub_trace.age>25 & abr_sub_trace.age<50 & ~abr_sub_trace.CP & ~abr_sub_trace.rjt_sub);
ONH_idx = find(abr_sub_trace.age>=50 & ~abr_sub_trace.CP & ~abr_sub_trace.rjt_sub);
HI_idx = find(abr_sub_trace.CP & ~abr_sub_trace.rjt_sub & abr_sub_trace.age>30);
%subplot(1,3,[1 2])
sub_abr_YNH = abr_sub_trace.sub_abr_b(YNH_idx,:);
sub_abr_MNH = abr_sub_trace.sub_abr_b(MANH_idx,:);
sub_abr_ONH = abr_sub_trace.sub_abr_b(ONH_idx,:);


sub_abr_group{1} = sub_abr_YNH;
sub_abr_group{2} = sub_abr_MNH;
sub_abr_group{3} = sub_abr_ONH;

abr_data_nh = IndexedStructCopy(abr_data,find(~abr_sub_trace.CP & ~abr_sub_trace.rjt_sub));

[fig5a]=plot_abr_group(abr_sub_trace.t_abr(:,1)',sub_abr_group,abr_data)
saveas(fig5a,[savepath 'abr_grouped'],'svg')


%% Stats
% AP
var = 'AP_amp_pm';
labels = '\muV'
lims = [-0.2 1.1];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
fig = gcf;
saveas(fig,[savepath 'AP_age'],'svg')

% WV
var = 'WV_amp_pm';
labels = '\muV'
lims = [-0.2 1.1];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
fig = gcf;
saveas(fig,[savepath 'WV_age'],'svg')
% SP
var = 'SP_amp';
labels = '\muV'
lims = [-0.19 0.29];
figure
[h,ax1,ax2]=stat_plots_uh_age_pt(uheal_data,var,labels,lims)
%hleg.Position = [0.5320 0.7793 0.2155 0.1272]
%set(hleg,'Visible','off')
fig = gcf;
saveas(fig,[savepath 'SP_age'],'svg')


%% plotting functions
function [fig1,fig2]=plot_abr(t_abr,sub_abr,age)
cm = cbrewer('qual','Set1',10)
cmap = cm([1 2 10],:);
rate_colors = {'k',[0.5 0.5 0.5]};
%rate_colors = {cmap(1,:),cmap(2,:),cmap(3,:)};
%subplot(1,3,[1 2])
% loop over rates
for kk=1:2
    
    abr_var=squeeze(nanstd(sub_abr(:,kk,:)))/sqrt(length(find(~isnan(sub_abr(:,1,1)))));
    abr_mean = squeeze(nanmean(sub_abr(:,kk,:)))';
    %sb=shadedErrorBar(t_abr,abr_mean,abr_var,...
    %'lineprops',['-' rate_colors{kk}],'transparent',0);
%     if kk==2
%     sb.mainLine.Color = [0.5 0.5 0.5]
%     sb.FaceColor = [0.4 0.4 0.4];
%     sb.edge(1).Color = [0.3 0.3 0.3];
%     sb.edge(2).Color = [0.3 0.3 0.3]
%     end

    hold on
    p(kk)=plot(t_abr,abr_mean,'color',rate_colors{kk},'linewidth',2);

    hold on
end
        xlim([-.5e-3 8e-3])%changed from 6e-3 to 8e-3
        ylim([-.2 .4])
        hold on
        box off
        %grid on
        
        
 plot(t_abr,zeros(size(t_abr)),'k--');
   

hleg = legend([p(1) p(2)],'9 Hz','40Hz');
hleg.Box = 'off';
%hleg.Position = [[0.1976 0.7609 0.4444 0.0964]];

set(gca,'fontsize',12,'xtick',[0:1:8]*1e-3)
%set(gca,'TickLabelInterpreter','latex');
ylabel('Amplitude\muV');
xlabel('Time (s)')  
set(gcf,'position',[293 510 315 215])
hleg.Position = [0.4472 0.7352 0.2794 0.1884];
fig1=gcf;
%axis tight
figure
%subplot(1,3,3)
hist(age,100)
xlim([0 99])
ylim([0 8])
xlabel('Age')
ylabel('n')
set(gca,'YAxisLocation','left','fontsize',12,'xtick',[25,50,75])
set(gcf,'position',[612 510 174 215])
 box off
 %set(gca,'TickLabelInterpreter','latex');

fig2=gcf;
end

function [fig1,fig2]=plot_abr_group(t_abr,sub_abr,abr_data)

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
rate_colors = {y_col,m_col,o_col};

%subplot(1,3,[1 2])
% loop over groups
for kk=1:3
    
    abr_var=squeeze(nanstd(sub_abr{kk}))/sqrt(length(find(~isnan(sub_abr{kk}(:,1)))));
    abr_mean = squeeze(nanmean(sub_abr{kk}(:,:)));

     shadedErrorBar(t_abr,abr_mean,abr_var,...
     'lineprops',{'color',rate_colors{kk}},'transparent',0);
    hold on
    p(kk)=plot(t_abr,abr_mean,'color',rate_colors{kk},'linewidth',1);

end
        xlim([-.5e-3 8e-3])%changed from 6e-3 to 8e-3
        ylim([-.25 .4])
        hold on
        box off
        %grid on
        
        
 plot(t_abr,zeros(size(t_abr)),'k--');

 %plot peak labels
 t=text(nanmean(abr_data.AP_lat),0.365,'AP','HorizontalAlignment','center','FontName','Arial','Color',[0.3 0.3 0.3])
 text(0.0005,0.15,'SP','HorizontalAlignment','center','FontName','Arial','Color',[0.3 0.3 0.3])
 text(nanmean(abr_data.WV_lat),-0.05,'WV','HorizontalAlignment','center','FontName','Arial','Color',[0.3 0.3 0.3])



set(gca,'fontsize',12,'xtick',[0:1:8]*1e-3)
ylabel('Amplitude \muV');
xlabel('Time (s)')  
set(gcf,'position',[293 510 315 215])

hleg = legend([p(1) p(2) p(3)],'Young','Middle-aged','Older');
hleg.Box = 'off';
hleg.Position = [0.3901 0.7189 0.4444 0.2744];
hleg.FontSize = 10;
fig1=gcf;

end
function [fig1,fig2]=plot_abr_group_tiny(t_abr,sub_abr)

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
rate_colors = {y_col,m_col,o_col};

%subplot(1,3,[1 2])
% loop over groups
for kk=1:3
    
    abr_var=squeeze(nanstd(sub_abr{kk}))/sqrt(length(find(~isnan(sub_abr{kk}(:,1,1)))));
    abr_mean = squeeze(nanmean(sub_abr{kk}(:,1,:)))';

 %    shadedErrorBar(t_abr,abr_mean,abr_var,...
 %    'lineprops',{'color',rate_colors{kk}},'transparent',0);
    hold on
    p(kk)=plot(t_abr,abr_mean,'color',rate_colors{kk},'linewidth',2);

end
        xlim([-.5e-3 8e-3])%changed from 6e-3 to 8e-3
        ylim([-.2 .4])
        hold on
        box off
        %grid on
        
        
 plot(t_abr,zeros(size(t_abr)),'k--');
   

%hleg = legend([p(1) p(2) p(3)],'Young','Middle-aged','Older');
%hleg.Box = 'off';
%hleg.Position = [0.3869 0.6910 0.4444 0.2744];

set(gca,'fontsize',20,'xtick',[0:2:8]*1e-3,'Xticklabel',[0:2:8]);%[{'0'},{''},{'2'},{''},{'4'},{''},{'6'},{''},{'8'},])
ylabel('\muV');
xlabel('Time (ms)')  
set(gcf,'position',[293 510 315 215])

fig1=gcf;

end

function T = IndexedStructCopy(S, Condition, FieldList)
if nargin == 2
   FieldList = fieldnames(S);
end 
for iField = 1:numel(FieldList)
   Field    = FieldList{iField};
   T.(Field) = S.(Field)(Condition);
end
end
