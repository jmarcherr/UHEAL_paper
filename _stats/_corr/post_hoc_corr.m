% post hoc correlations

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

%EFR SNR
sig_idx = find(uheal_nh.EFR_sig==0);
uheal_nh.EFR_SNR(sig_idx) = nan;

%% post hoc correlations
% AP vs. FFR, EFR, MEMR
% FFR vs. EFR, MEMR
[cmap_g,cmap] = colormapuheal;
fnames = {'AP_amp_pm','FFR_SNR','EFR_SNR','memr_slope'}

close all
clear rho p
clear rho p
for ii=1:4
for jj=1:4
    thisidx = ~isnan(uheal_nh.(fnames{ii})) & ~isnan(uheal_nh.(fnames{jj}));
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
    if jj+(ii-1)*length(fnames) == 21
        set(gca,'fontsize',12)
    else
        set(gca,'xticklabels',{},'yticklabels',{})
    end
end
end

%
cidx = [2:4,7,8,12]
for tt=1:length(cidx)
subplot(length(fnames),length(fnames),cidx(tt))
hold off
plot(0,1)
axis off
%text(0,1,[num2str(round(rho(cidx(tt)),4)) '\newline' num2str(round(p_corr(cidx(tt)),4))])
end
%
diag_idx = find(eye(length(fnames),length(fnames)));%,29,36];
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
set(gcf,'position',[823 828 494 501])
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
        text(0,0,['p = ' num2str(round(p(pidx(pp)),3)) '\newline' 'r = ' num2str(round(rho(pidx(pp)),3))])
        xlim([-1 1])
        ylim([-1 1])
        axis off
    else
        text(0,0,['p = ' num2str(round(p(pidx(pp)),3)) '\newline' 'r = ' num2str(round(rho(pidx(pp)),3))])
        xlim([-1 1])
        ylim([-1 1])
        axis off
    end
end
%
pnames = {'AP','FFR','EFR','MEMR'}
yidx = [1,6,11,16,21]
for ii=1:length(fnames)
    subplot(length(fnames),length(fnames),ii)
    t = title(pnames{ii},'FontSize',14,'FontName','Arial')
    subplot(length(fnames),length(fnames),yidx(ii))
%     if ii==4
%         ylab = ylabel('Neg. \newline 4Hz')
%     else
%         ylab = ylabel(pnames{ii})
%     end
%     ylab.FontSize =14;
%     ylab.FontName = 'Arial';
%     ylab.FontWeight = 'bold'
%     ylab.Rotation = 0;
%     ylab.HorizontalAlignment = 'Right'
%     ylab.VerticalAlignment = 'middle'
%     ylab.Position = ylab.Position-[0.6 0 0]
%     ylab.Color = [0 0 0];
    %tmp.YLabel.Position = tmp.YLabel.Position-[0 0 1];
    
end

fig=gcf;
set(fig,'renderer','painters')
%cd(root
set(gcf,'position',[100 100 397 400]*(96/72))
saveas(gcf,'/work3/jonmarc/UHEAL_paper/_stats/_corr/figs/subcort_corr','svg')



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