%load('uheal_data.mat')
%var = 'PTA_hf';
%labels = 'PTA_{hf} (dB HL)'
%lims = [-40 120];


% Function for plotting uheal data stats

% plots box plots over groups and scatter plots as a function of age

function [h,ax1,ax2] = stat_plots_uh_age_pt(data,var,labels,lims)
%% get data and var

if strcmp(var,'AEP_n1')
    % sig points
    this_data = data.(var); % gets relevant field
    this_data = this_data(:,4);
    %this_data = ;
    %this_data(find(data.FFR_sig==0))=nan;
    %data.Age(~isnan(this_data)) = nan;
    %data.CP_new(~isnan(this_data)) = nan;
    %data.gender(~isnan(this_data)) = nan;
    % non sig points
    %non_sig = data.(var);
    %non_sig(find(data.FFR_sig==0)) = nan;
elseif strcmp(var,'AEP_p2')
    this_data =data.(var);
    this_data = this_data(:,1);
else
    this_data = data.(var);
end

uheal_data = data;
%% define color map

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

% find groups
% female
ynhf = find(uheal_data.Age<=25 & ~uheal_data.CP_new & uheal_data.gender==1);
mnhf = find(uheal_data.Age>25 & uheal_data.Age<50 & ~uheal_data.CP_new & uheal_data.gender==1);
onhf = find(uheal_data.Age>=50 & ~uheal_data.CP_new & uheal_data.gender==1);

% male
ynhm = find(uheal_data.Age<=25 & ~uheal_data.CP_new & uheal_data.gender==2);
mnhm = find(uheal_data.Age>25 & uheal_data.Age<50 & ~uheal_data.CP_new & uheal_data.gender==2);
onhm = find(uheal_data.Age>=50 & ~uheal_data.CP_new & uheal_data.gender==2);

% all
y = find(uheal_data.Age<=25 & ~uheal_data.CP_new);
m = find(uheal_data.Age>25 & uheal_data.Age<50 & ~uheal_data.CP_new);
o = find(uheal_data.Age>=50 & ~uheal_data.CP_new);

% find outliers
outym = isoutlier(this_data(ynhm));outyf = isoutlier(this_data(ynhf));oy =isoutlier(this_data(y));
outmm = isoutlier(this_data(mnhm));outmf = isoutlier(this_data(mnhf));om = isoutlier(this_data(m));
outom = isoutlier(this_data(onhm));outof = isoutlier(this_data(onhf));oo = isoutlier(this_data(o));

%% get distirubtions
figure
% male
h=distributionPlot({this_data(ynhm(~outym)),this_data(mnhm(~outmm)),this_data(onhm(~outom))},'widthDiv',[2 1],'histOri','left','color',[0.9 0.9 0.9],'showMM',0)
tmp = h(3);
for ii=1:3
    curvemx(ii,:) = tmp{1}.Children(ii).XData(1,:)';
    curvemy(ii,:) = tmp{1}.Children(ii).YData(1,:)';
end
close gcf
h=distributionPlot(gca,{this_data(ynhf(~outyf)),this_data(mnhf(~outmf)),this_data(onhf(~outof))},'widthDiv',[2 2],'histOri','right','color',[0.7 0.7 0.7],'showMM',0)
tmp = h(3);
for ii=1:3
    
    curvefx(ii,:) = tmp{1}.Children(ii).XData(2,:)';
    curvefy(ii,:) = tmp{1}.Children(ii).YData(2,:)';
end

close gcf
figure
h=distributionPlot(gca,{this_data(y(~oy)),this_data(m(~om)),this_data(o(~oo))},'widthDiv',[2 2],'histOri','right','color',[0.7 0.7 0.7],'showMM',0)
tmp = h(3);
for ii=1:3
    
    curvexa(ii,:) = tmp{1}.Children(ii).XData(2,:)';
    curveya(ii,:) = tmp{1}.Children(ii).YData(2,:)';
end

fig =gcf;
close gcf

%% plotting

figure
fig1 = gcf;
ax1=subplot(1,2,1)
hold on
% plot distribution curves for all
plot(curvexa',curveya','color','k')

% unused for data points
%plotSpread({uheal_data.PTA_hf(ynh),uheal_data.PTA_hf(mnh),uheal_data.PTA_hf(onh)},'showMM',0,'distributionMarkers','o','XValues',[1 2 3]-0.25,'distributionColors','k')
%pc1=plotSpread({uheal_data.PTA_hf(ynh(~outy)),uheal_data.PTA_hf(mnh(~outm)),uheal_data.PTA_hf(onh(~outo))},'showMM',0,'distributionMarkers','o','XValues',[1 2 3]-0.25,'distributionColors','k','categoryIdx',[ones(size(ynh(~outy))); zeros(size(mnh(~outm))); ones(size(onh(~outo)))],'categoryMarkers',{'o','^'})
%pc=plotSpread({uheal_data.PTA_hf(y),uheal_data.PTA_hf(m(~om)),uheal_data.PTA_hf(o(~oo))},'showMM',0,'distributionMarkers','o','XValues',[1 2 3]-0.25,'distributionColors','k','categoryIdx',[uheal_data.gender(y)==1;uheal_data.gender(m(~om))==1 ; uheal_data.gender(o(~oo))==1],'categoryMarkers',{'o','^'},'categoryColors',[0 0 0;0.5 0.5 0.5])
%plotSpread({uheal_data.PTA_hf(ynhf(~outyf)),uheal_data.PTA_hf(mnhf(~outmf)),uheal_data.PTA_hf(onhf(~outof))},'showMM',0,'distributionMarkers','^','XValues',[1 2 3]-0.25,'distributionColors',[0.5 0.5 0.5])

% plot boxplot for 3 groups
boxplot(this_data(y),'positions',[1],'boxstyle','filled','color',y_col,'symbol','.','widths',.3)
boxplot(this_data(m),'positions',[2],'boxstyle','filled','color',m_col,'symbol','.','widths',.3)
boxplot(this_data(o),'positions',[3],'boxstyle','filled','color',o_col,'symbol','.','widths',.3)
set(gca,'fontsize',12,'xtick',[1 2 3],'xticklabels',{'\leq 25';'26-49';'\geq 50'},'TickLabelInterpreter','tex')%,'fontname','Arial')
ylabel(labels)
xtickangle(45)
ylim(lims)
xlim([0.5 3.5])
box off
set(gcf,'position',[680 852 197 246])

% for patch
%patch('XData',curvefx','YData',curvefy')


% scatter plots
% get lsline
figure(1000)
s=scatter(uheal_data.Age(~uheal_data.CP_new),this_data(~uheal_data.CP_new))
ll=lsline;xls = ll.XData;yls=ll.YData;
fig = gcf;
close gcf

%return to original plot
fig1;
ax2=subplot(1,2,2)

% scatter of all points
hold on
% female
scatter(uheal_data.Age(~uheal_data.CP_new & uheal_data.gender==1),this_data(~uheal_data.CP_new & uheal_data.gender==1),'k.','markeredgecolor','k','markerfacecolor','w','markerfacealpha',1,'SizeData',30,'linewidth',.5)
% male
s=scatter(uheal_data.Age(~uheal_data.CP_new & uheal_data.gender==2),this_data(~uheal_data.CP_new & uheal_data.gender==2),'k.','markerfacecolor','w','markerfacealpha',1,'SizeData',30,'linewidth',.5)


% plot ls line
plot(xls,yls,'color',[0.6 0 0.1 0.5],'linewidth',2)

% get legend
%mm=plot(-100,-100,'ko','markerfacecolor','w','MarkerSize',2)
%fm=plot(-100,-100,'k^','markeredgecolor','k','markerfacecolor','w','MarkerSize',2)

%hleg = legend([mm,fm],{'Male','Female'})
%hleg.Box = 'off';hleg.FontSize = 7;hleg.Location = 'Best';%;hleg.Orientation = 'horizontal';
%hleg.Position = [0.5110 0.2823 0.2155 0.1272]
set(gca,'fontsize',12)

% set labels and dimensions
set(gcf,'position',[100 100 297 169]*(96/72))
xlabel('Age')
xtickangle(45)
xlim([15 87])


% link and scale
linkaxes([ax1,ax2],'y');
set(ax2,'Position',[0.1892+0.36 0.2686 0.2754 0.6564])
set(ax1,'position',[0.1892 0.2686 0.2754 0.6564])
set(gcf,'renderer','painters')

h = gcf;
end
