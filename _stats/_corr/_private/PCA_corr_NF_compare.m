close all
clear all
%cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat')
%%

% preprocessing
%nha
nh_idx = find(uheal_data.CP_new==0);
Field_list = fieldnames(uheal_data);
for ii=1:numel(Field_list)
    if (size(uheal_data.(Field_list{ii}),2)==1)% && size(uheal_data.(Field_list{ii}),2)==117)
        uheal_tmp.(Field_list{ii}) = uheal_data.(Field_list{ii})
    end
end
uheal_nh = rmfield(uheal_tmp,'memr_fcenter');
uheal_nh.AEP_p2n1_comp = uheal_data.AEP_p2n1_comp(:,4);
uheal_nh.AEP_n1 = uheal_data.AEP_n1(:,4);
uheal_nh.AEP_p2 = uheal_data.AEP_p2(:,1);
%%
%uheal_nh.AP_amp_pm = uheal_data.AP_amp_pm
%uheal_nh.FFR_SNR = uheal_data.FFR_SNR
%uheal_nh.FFR_sig = uheal_data.FFR_sig
%uheal_nh.AEP_p2n1_int = uheal_data.AEP_p2n1_int;
%uheal_tmp
    %%

uheal_nh = IndexedStructCopy(uheal_nh,nh_idx);
% groups
ynh = find(uheal_nh.Age<25);
mnh = find(uheal_nh.Age>=25 & uheal_nh.Age<50);
onh = find(uheal_nh.Age>=50);

% FFR SNR
sig_idx = find(uheal_nh.FFR_sig==0);
uheal_nh.FFR_SNR(sig_idx) = nan;
%%
% select variables
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
Varnames = {'AEP_NF','ABR_NF','FFR_NF'}
X = {zscor_xnan(uheal_nh.AEP_NF_int);
    zscor_xnan(uheal_nh.ABR_NF_int)
    zscor_xnan(uheal_nh.FFR_NF_int);
    };
X = reshape(cell2mat(X),105,length(Varnames));
X_imp = knnimpute(X);
[coeff,score,latent,tsquared,explained]=pca(X_imp);

% skree plot
close all
bar(1:length(Varnames),explained,'FaceColor',[0.5 0.5 0.5])
hold on
plot(1:length(Varnames),explained,'k.-')
box off
xlabel('Dim.')
ylabel('Var. explained (%)')
set(gca,'fontsize',12)
set(gcf,'Position',[100 100 202 140]/(96/72))
xlim([0.1 6.9])


%% Correlation plots
close all
% AP, FFR, ITPC, 4 Hz neg
fnames = {'AEP_NF_int','ABR_NF_int','FFR_NF_int'}
f_idx = uheal_nh.gender ==1;
m_idx = uheal_nh.gender ==2;
clear rho p
for ii=1:length(fnames)
for jj=1:length(fnames)
    thisidx = ~isnan(uheal_nh.(fnames{ii})) & ~isnan(uheal_nh.(fnames{jj}));
    thisidx_m = ~isnan(uheal_nh.(fnames{ii})) & ~isnan(uheal_nh.(fnames{jj})) & uheal_nh.gender ==1;
    thisidx_f = ~isnan(uheal_nh.(fnames{ii})) & ~isnan(uheal_nh.(fnames{jj})) & uheal_nh.gender ==2;
    [rho(jj+(ii-1)*length(fnames)),p(jj+(ii-1)*length(fnames))]=corr(zscore(uheal_nh.(fnames{ii})(thisidx)),zscore(uheal_nh.(fnames{jj})(thisidx)),'type','spearman');
    subplot(length(fnames),length(fnames),jj+(ii-1)*length(fnames))
    scatter(zscore(uheal_nh.(fnames{ii})(thisidx)),zscore(uheal_nh.(fnames{jj})(thisidx)),'k.')
    hold on
    xlim([-3 3])
    ylim([-3 3])
    ll = lsline
    xl = ll.XData;
    yl = ll.YData;
    hold off
    %ss = scatter(zscore(uheal_nh.(fnames{ii})(thisidx)),zscore(uheal_nh.(fnames{jj})(thisidx)),'.','markerEdgecolor',[1 0.3 0.1])
    hold on
    this_pt = find(thisidx);
    this_x = zscore(uheal_nh.(fnames{ii})(thisidx));
    this_y = zscore(uheal_nh.(fnames{jj})(thisidx));
    for ss=1:length(this_pt)
        pp=plot(this_x(ss),this_y(ss),'.','MarkerEdgeColor',cmap(uheal_nh.Age(this_pt(ss))-17,:),'MarkerSize',5);
    end
    %scatter(zscore(uheal_nh.(fnames{ii})(thisidx)),zscore(uheal_nh.(fnames{jj})(thisidx)),'k.')
    plot([xl],[yl],'color',[0.75 0.75 0.75])
    xlim([-3 3])
    ylim([-3 3])
    if jj+(ii-1)*length(fnames) == length(fnames)*length(fnames)-(length(fnames)-1)
        set(gca,'fontsize',12)
    else
        set(gca,'xticklabels',{},'yticklabels',{})
    end
end
end

%

cidx = [2 3 6]
%cidx = [2:5,8:10,14:15,20]
for tt=1:length(cidx)
subplot(length(fnames),length(fnames),cidx(tt))
hold off
plot(0,1)
axis off
%text(0,1,[num2str(round(rho(cidx(tt)),4)) '\newline' num2str(round(p_corr(cidx(tt)),4))])
end
%
diag_idx = find(eye(length(fnames),length(fnames)));%[1,8,15,22,29,36];
for tt=1:length(diag_idx)
    subplot(length(fnames),length(fnames),diag_idx(tt))
    hold off
    thisidx = ~isnan(uheal_nh.(fnames{tt}));
    h=histfit(zscore(uheal_nh.(fnames{tt})(thisidx)),10,'kernel')
    h(1).FaceColor = [1 1 1];
    h(1).EdgeColor = [1 1 1];
    h(2).Color = [0 0 0]
    xlim([-3 3])
    ylim([0 35])
    set(gca,'xticklabels',{},'yticklabels',{})
    box off
end
%set(gcf,'position',[823 828 494 501])
%
pidx = setdiff(cidx,diag_idx)
p_corr = p*((length(p)-length(fnames))/2);
c_y = [239 210 84]/255;
c_b = [158 204 239]/255;
for pp=1:length(pidx)
    subplot(length(fnames),length(fnames),pidx(pp))
    %text(1,1,num2str(p_corr(pidx(pp))))
    %hold on
    %text(1,2,num2str(rho(pidx(pp))))
    %hold off
    if p_corr(pidx(pp))<0.05
        circ = nsidedpoly(1000, 'Center', [0 0], 'Radius', abs(rho(pidx(pp))/max(rho)));
        ccol = {'r','b'}
        
        if rho(pidx(pp))<0
        circp=plot(circ,'EdgeColor',c_b,'FaceColor','none','Linewidth',2);%circp.EdgeAlpha = 0.2;
        else
            circp=plot(circ,'EdgeColor',c_y,'FaceColor','none','Linewidth',2);%circp.EdgeAlpha = 0.2;
        end
        
        hold on
        text(0,0,['p = ' num2str(round(p_corr(pidx(pp)),3)) '\newline' 'r = ' num2str(round(rho(pidx(pp)),3))])
        xlim([-1 1])
        ylim([-1 1])
        axis off
    end
end
%
pnames = {'AEP_{NF}','ABR_{NF}','FFR_{NF}'}
%yidx = [1,6,11,16,21]
yidx = [1,4,7]%,17,22,31]
for ii=1:length(fnames)
    subplot(length(fnames),length(fnames),ii)
    t = title(pnames{ii},'FontSize',14,'FontName','Arial')
    subplot(length(fnames),length(fnames),yidx(ii))
    % if ii==4
    %     ylab = ylabel('Neg. \newline 4Hz')
    % else
    %     ylab = ylabel(pnames{ii})
    % end
    % ylab.FontSize =14;
    % ylab.FontName = 'Arial';
    % ylab.FontWeight = 'bold'
    % ylab.Rotation = 90;
    % ylab.HorizontalAlignment = 'Right'
    % ylab.VerticalAlignment = 'middle'
    % ylab.Position = ylab.Position-[1 0 0]
    % ylab.Color = [0 0 0];
    %tmp.YLabel.Position = tmp.YLabel.Position-[0.8 0 1];
    
end

fig=gcf;
set(fig,'renderer','painters')
%cd(root
set(gcf,'position',[133.3333 133.3333 621.6667 533.3333])
saveas(fig,['figs/corr_eeg_NF_compare'],'svg')
%% correlation with questionaires
uheal_nh.nesi_log = log10(uheal_nh.nesi); uheal_nh.nesi_log(isinf(uheal_nh.nesi_log))=nan;
clinnames = {'nesi_log','tts','ssq12_mean','PTA_lf','PTA_hf','memr_slope'}
eegnames = {'AP_amp_pm','FFR_SNR','ITPC_ratio','Neg_4Hz'};
peegnames = {'AP','FFR','ITPC','Neg. 4Hz','N1','P2'};
pclinnames = {'NESI','TTS','SSQ12','PTA_{lf}','PTA_{hf}','MEMR'}
clear p rho
close all
figure('renderer','painters')
%set(gcf,'position',[119 784 356 311])
cc=0;
for ii=1:length(clinnames)
    for jj=1:length(eegnames)
        cc = cc+1
        thisidx = ~isnan(uheal_nh.(clinnames{ii})) & ~isnan(uheal_nh.(eegnames{jj}))
        [rho(ii,jj),p(ii,jj)]=corr(zscore(uheal_nh.(clinnames{ii})(thisidx)),zscore(uheal_nh.(eegnames{jj})(thisidx)));
        p_corr = p*size(p,1)*size(p,2);

        subplot(6,6,cc)
        scatter(zscore(uheal_nh.(clinnames{ii})(thisidx)),zscore(uheal_nh.(eegnames{jj})(thisidx)),'.k')
        lsline
        xlim([-3 3])
        ylim([-3 3])
        box on
        if cc==36
            set(gca,'yaxislocation','right')
        else
        set(gca,'xticklabel',{},'yticklabel',{})
        end
    end
end
for ii=1:6
    subplot(6,6,ii)
    title(peegnames{ii})
end
cc = 0;
for ii=[1,7,13,19,25,31]
    cc = cc+1
    subplot(6,6,ii)
    ylab = ylabel(pclinnames{cc})
    ylab.Rotation = 0;ylab.HorizontalAlignment = 'right'
    set(gca,'FontWeight','bold')
    set(get(gca,'YLabel'),'Visible','on')
end

% corr plot

p_plot = p';p_plot = p_plot(:);p_plot = p_plot*length(p_plot);
rho_plot = rho';rho_plot = rho_plot(:)
figure('renderer','painters')
%set(gcf,'position',[119 784 356 311])
for pp = 1:length(p_plot)

    subplot(6,6,pp)
    
    if p_plot(pp)<0.05
            circ = nsidedpoly(1000, 'Center', [0 0], 'Radius', abs(rho_plot(pp)/max(abs(rho_plot))));
        ccol = {'r','b'}
        if rho_plot(pp)<0
        circp=plot(circ,'EdgeColor',cmap(1,:),'FaceColor','none','Linewidth',2);%circp.EdgeAlpha = 0.2;
        else
            circp=plot(circ,'EdgeColor',cmap(2,:) ...
                ,'FaceColor','none','Linewidth',2);%circp.EdgeAlpha = 0.2;
        end
        xlim([-1.5 1.5])
        ylim([-1.5 1.5])
        
        %axis off
        box on
    else
        %axis off
        box on
      set(gca,'xticklabel',{},'yticklabel',{})
    end    
    set(gca,'xticklabel',{},'yticklabel',{},'xtick',[],'ytick',[])

end
for ii=1:6
    subplot(6,6,ii)
    tit = title(peegnames{ii})
    %tit.Rotation = 20
    %tit.HorizontalAlignment = 'left'
end
cc = 0;
for ii=[1,7,13,19,25,31]
    cc = cc+1
    subplot(6,6,ii)
    ylab = ylabel(pclinnames{cc})
    ylab.Rotation = 0;ylab.HorizontalAlignment = 'right'
    set(gca,'FontWeight','bold')
    set(get(gca,'YLabel'),'Visible','on')
    
end


%%
% functions
function T = IndexedStructCopy(S, Condition, FieldList)
if nargin == 2
   FieldList = fieldnames(S);
end 
for iField = 1:numel(FieldList)
   Field    = FieldList{iField};
   T.(Field) = S.(Field)(Condition);
end
end

function [rate_colors,cmap] = colormapuheal
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
end

function [rho,pval] = corr_plot_UH(x,y,labels,fidx,midx)
figure('renderer','painter')
xlab = labels.xlab;
ylab = labels.ylab;
subplot(2,2,1)


s=scatter(x,y,'o','markerfacecolor','none','MarkerEdgecolor','none')
% get lsline data
ll=lsline
ls_x = ll.XData;
ls_y = ll.YData;

% female
scatter(x(fidx),y(fidx),'o','markerfacecolor',[0.8 0.5 0.5],'markeredgecolor','k')
hold on
% male
scatter(x(midx),y(midx),'o','markerfacecolor',[0.5 0.5 0.5],'markeredgecolor','k')
ll=plot(ls_x,ls_y,'-');
set(ll,'linewidth',2,'color','k');
set(gca,'fontsize',12)
ylabel(ylab)
xlabel(xlab)
%set(gcf,'Position',[114 210 280 418]);
set(gcf,'Position',[114 210 177 418]);
% get non-nan values
idx = ~isnan(x) & ~isnan(y)
[rho,pval]=corr(x(idx),y(idx))
fig = gcf;
hold on


if pval<0.05
    annotation('textbox',[.2 .35 .1 .1], ...
    'String',['p = ' num2str(pval) ',\newline \rho = ' num2str(round(rho,3))],'EdgeColor','none')
else
    annotation('textbox',[.2 .35 .1 .1], ...
    'String',['p = n.s. ',num2str(pval), '\rho = ' num2str(round(rho,3))],'EdgeColor','none')    
end

set(gcf,'renderer','painter')

end

function [s,ll] = scatteruh(x,y,data)
% get lsline
figure(1000)
s=scatter(x,y)
ll=lsline;xls = ll.XData;yls=ll.YData;
fig = gcf;
close gcf

%return to original plot
figure

% scatter of all points
hold on
% female
scatter(x(data.gender==1),y(data.gender==1),'k.','markeredgecolor','k','markerfacecolor','w','markerfacealpha',1,'SizeData',20,'linewidth',.5)


% male
s=scatter(x(data.gender==2),y(data.gender==2),'k.','markerfacecolor','w','markerfacealpha',1,'SizeData',20,'linewidth',.5)
ll=lsline
ll(1).Color = 'r'
ll(2).Color = 'b'
% plot ls line
%plot(xls,yls,'color',[0.6 0 0.1 0.5],'linewidth',2)

% get legend
mm=plot(-100,-100,'k.','markerfacecolor','k','MarkerSize',2)
fm=plot(-100,-100,'k.','markeredgecolor','k','markerfacecolor','w','MarkerSize',2)

%hleg = legend([mm,fm],{'Male','Female'})
%hleg.Box = 'off';hleg.FontSize = 7;hleg.Location = 'Best';%;hleg.Orientation = 'horizontal';
%hleg.Position = [0.5110 0.2823 0.2155 0.1272]
%set(gca,'fontsize',12)

% set labels and dimensions
set(gcf,'position',[938 371 104 141])%%[552 636 297 169])
xlabel('Age')
xtickangle(45)
xlim([15 87])
end
