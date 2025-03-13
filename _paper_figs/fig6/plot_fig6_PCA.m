close all
clear all
%cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat')
savepath = ['/work3/jonmarc/UHEAL_paper/_paper_figs/fig6/figs/']
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
ynh = find(uheal_nh.Age<=25);
mnh = find(uheal_nh.Age>25 & uheal_nh.Age<50);
onh = find(uheal_nh.Age>=50);

% FFR SNR
sig_idx = find(uheal_nh.FFR_sig==0);
uheal_nh.FFR_SNR(sig_idx) = nan;
%%
% select variables
zscor_xnan = @(x) bsxfun(@rdivide, bsxfun(@minus, x, mean(x,'omitnan')), std(x, 'omitnan'));
Varnames = {'AP','FFR','ITPC_{Q}','AEP_{\mu}','P2_{adapt}','P2-N1_{int}'}
X = {zscor_xnan(uheal_nh.AP_amp_pm);
    zscor_xnan(uheal_nh.FFR_SNR);
    zscor_xnan(uheal_nh.ITPC_ratio);
    zscor_xnan(uheal_nh.Neg_4Hz);
    zscor_xnan(uheal_nh.tone_p2_1vsrest);
    zscor_xnan(uheal_nh.AEP_p2n1_int)};
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

%% contributions plot
% Variance in variable I explained by principal component J
for i = 1:length(Varnames)
    for j = 1:2
        varI(i,j) = coeff(i,:)*(latent.*coeff(i,:)')
        varIfromJ(i,j) = coeff(i,j)*latent(j)*coeff(i,j)
        percVarIfromJ(i,j) = varIfromJ(i,j)/varI(i,j)
    end
end
b_col = cbrewer('qual','Set1',4);
b_col = b_col(2,:);
cols = [repmat([0.5 0.5 0.5],2,1);repmat(b_col,4,1)]
% contributions dim 1
close all
for pp=1:2
    figure('renderer','painters')
    %if pp==2
        [s,i]=sort(varIfromJ(:,pp),'ascend')

    %else
    %    [s,i]=sort(varIfromJ(:,pp),'descend')
    %end

b = bar(s*100,'FaceColor','flat')
b.CData = cols(i,:);
for ii=1:length(Varnames)
sortnames{ii} = Varnames{i(ii)}
end

set(gca,'fontsize',12)
%set(gcf,'Position',[680 958 202 140])

%xlim([0.1 6.9])
%ylim([0 0.7])
box off
%if pp==2
    f = gcf
%f.Position = [939 586 103 141]
%f.Position = [100 100 103 111]*(96/72)
f.Position = [519 215 264 231];
xlab = xlabel(['PC' num2str(pp) ' (%)'])
xlab.FontSize = 12
xlab.FontName = 'none'
xlab.HorizontalAlignment = 'center'

%xlab.Position = xlab.Position% + [0 0 0];
b.Horizontal = 'on'
set(gca,'yticklabels',sortnames,'Yaxislocation','right','xaxislocation','top','FontName','Arial','Fontsize',12)
    ylim([0.5 length(Varnames)+0.8])
    xlim([0 0.8]*100)
%[939 586 145 141]
%ytickangle(45)
%else
%     ylabel('Contribution (%)')
%     set(gca,'xticklabels',sortnames)
%     xtickangle(45)
%     xlim([0.1 6.9])
%     ylim([0 0.7])
% f = gcf
% f.Position = [680 824 300 274]*1.1
%end
%title(['Dimension ' num2str(pp)])
saveas(f,[savepath 'pca_contri_PC' num2str(pp)],'svg')
end

%% PCA plot
close all
[cmap_g,cmap] = colormapuheal
clear s
figure('renderer','painters')
rscore = rescale(score,-1,1);
%s(1) = scatter(rscore(ynh,1),rscore(ynh,2),'k^','markerfacecolor',[cmap{1}],'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none')
hold on
%s(2) = scatter(rscore(mnh,1),rscore(mnh,2),'ksq','markerfacecolor',[cmap{2}],'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none')
%s(3) = scatter(rscore(onh,1),rscore(onh,2),'ko','markerfacecolor',[cmap{3}],'MarkerFaceAlpha',0.5,'MarkerEdgeColor','none')
QX = coeff(:,1);
QY = coeff(:,2);
%{'AP','FFR','TE-ITPC_{Q}','TE-AEP_{mean}','TE-AEP_{P2 adapt.}','NE-AEP_{P2-N1}'}
b_col = cbrewer('qual','Set1',4);
b_col = b_col(2,:);
for ii = 1:length(QX)
    plot([0 QX(ii)],[0 QY(ii)],'color',[0.5 0.5 0.5]-0.5,'linewidth',1)
    % AP
    if ii==1
        t=text(QX(ii)-0.15,QY(ii)+0.1,Varnames{ii},'Color','k')
        
    elseif ii==2 % FFR
                t=text(QX(ii)-0.1,QY(ii)+0.1,Varnames{ii},'Color','k')
    elseif ii==3 % ITPC
                t=text(QX(ii)-0.3,QY(ii)-0,Varnames{ii},'Color',b_col)
    elseif ii==4 % TE mean
                t=text(QX(ii)+0.05,QY(ii)-0,Varnames{ii},'Color',b_col)
    elseif ii==5 % P2
                t=text(QX(ii)-0.1,QY(ii)+0.18,Varnames{ii},'Color',b_col)
    elseif ii==6 % p2-n1   
                t=text(QX(ii)+0.05,QY(ii)-0.0,Varnames{ii},'Color',b_col)
    else
        t=text(QX(ii)+0.05,QY(ii),Varnames{ii})
    end
    t.FontName = 'Arial';
    t.FontSize = 11;
    t.FontWeight = 'Bold';
end
p = nsidedpoly(100000, 'Center', [0 0], 'Radius', 1);
p=plot(p,'FaceColor','none');p.EdgeAlpha = 0.2;
xlim([-1 1])
ylim([-1 1])
box off
grid on
%plot([-1 1],[0 0],'--','color',[0.5 0.5 0.5])
%plot([0 0],[-1 1],'--k')
f = gcf
%f.Position = [680 824 300 274]*1.1
%f.Position = [100 100 282 243]*(96/72)
f.Position = [[387 464 460 410]];
xlabel(['PC1 (' num2str(round(explained(1))) ' %)'])
ylabel(['PC2 (' num2str(round(explained(2))) ' %)'])
set(gca,'fontsize',12)
%hleg = legend(s,'Young','Mid. aged','Older')
%hleg.Box = 'on'
%hleg.Position =[0.1471 0.1423 0.3000 0.1715];


fig=gcf;
%cd(root)
saveas(fig,[savepath 'pca_eeg'],'svg')


%% PCA C1 vs. age

[s,ll]=scatteruh(uheal_nh.Age,score(:,1),uheal_nh)
    f = gcf
    set(f,'renderer','painters')
f.Position = [938 371 104 141];%[939 586 192 141]
ylabel('Dim.1')
[rho_age,p_age]=corr(uheal_nh.Age,score(:,1))
[rho_age2,p_age2]=corr(uheal_nh.Age,score(:,2))
ll(1).Color = cmap_g{2}
ll(2).Color = cmap_g{3}


[s,ll]=scatteruh(uheal_nh.Age,score(:,2),uheal_nh)
    f = gcf
    set(f,'renderer','painters')
f.Position = [938 371 104 141];%[939 586 192 141]
ylabel('Dim.2')
[rho_sex1,p_sex1] = corr(uheal_nh.gender,score(:,1))
[rho_sex,p_sex]=corr(uheal_nh.gender,score(:,2))
ll(1).Color = cmap_g{2}
ll(2).Color = cmap_g{3}

%% PC1 and PC2 vs clincal measures
fnames = {'Age','gender','nesi','tts','ssq12_mean','memr_slope','PTA_hf'}
pnames = {'Age','Sex','NESI','TTS','SSQ12','MEMR','PTA_{HF}'}
clear p rho p_corr
c_y = [239 210 84]/255;
c_b = [158 204 239]/255;
close all
figure('renderer','painters')
for jj=1:2 % PC1 and 2
    for ii=1:length(fnames) % NESI, TTS, SSQ12, DS, MEMR
        thisidx = ~isnan(uheal_nh.(fnames{ii})) & ~isnan(score(:,jj));
        if ii==2
            [rho(jj,ii),p(jj,ii)]=corr(uheal_nh.(fnames{ii})(thisidx),...
                score(thisidx,jj),'type','Spearman');
        else
            [rho(jj,ii),p(jj,ii)]=corr(zscore(uheal_nh.(fnames{ii})(thisidx)),...
                score(thisidx,jj),'type','Spearman');
        end
        p_corr(jj,ii) = p(jj,ii)*length(fnames)*2;
        
        subplot(2,length(fnames),ii+(jj-1)*length(fnames))
        box on
        set(gca,'xticklabel',{},'yticklabel',{},'xtick',[],'ytick',[])
        %titles
        if jj==2
            t_box = xlabel(pnames{ii})
            t_box.FontSize = 12;
            t_box.FontName = 'Arial'
            t_box.Rotation = 45;
            t_box.HorizontalAlignment = 'right'
            t_box.FontWeight = 'Bold'
        else
        end
        %ylab
        if ii==1
            ylab = ylabel(['PC' num2str(jj)])
            ylab.FontSize = 12;
            ylab.FontName = 'Arial'
            ylab.Rotation = 0;
            ylab.HorizontalAlignment = 'right'
            ylab.FontWeight = 'Bold'
            ylab.Position = ylab.Position - [0 .15 0]
        end
        if p_corr(jj,ii)<0.05
            circ = nsidedpoly(1000, 'Center', [0 0], 'Radius', abs(rho(jj,ii)/max(rho)));
            ccol = {'r','b'}
            if rho(jj,ii)<0
                circp=plot(circ,'EdgeColor',c_b,'FaceColor','none','Linewidth',2);%circp.EdgeAlpha = 0.2;
            else
                circp=plot(circ,'EdgeColor',c_y,'FaceColor','none','Linewidth',2);%circp.EdgeAlpha = 0.2;
            end
            %axis off
            box  on
            
            %ylabel(num2str(jj))
            %xlim([-1 1])
            %ylim([-1 1])
            set(gca,'xticklabel',{},'yticklabel',{},'xtick',[],'ytick',[])
            %t_box = xlabel(pnames{ii})
            %                    t_box.FontSize = 12;
            %t_box.FontName = 'Arial'
            %t_box.Rotation = 45;
            % t_box.HorizontalAlignment = 'left'
            if ii==1
                ylab = ylabel(['PC' num2str(jj)])
                ylab.FontSize = 12;
                ylab.FontName = 'none'
                ylab.Rotation = 0;
                ylab.HorizontalAlignment = 'right'
                ylab.FontWeight = 'Bold'
                ylab.Position = ylab.Position - [0 .4 0]
            end
        end
        
    end
    
end
rho
p_corr

set(gcf,'Position',[[128 193 779 286]])
fig = gcf;
saveas(fig,[savepath 'pca_clin'],'svg')
%% Correlation plots
close all
% AP, FFR, ITPC, 4 Hz neg
%{'AP','FFR','ITPC','Mean Neg.','P2_{adapt.}','N1','P2'}
fnames = {'AP_amp_pm','FFR_SNR','ITPC_ratio','Neg_4Hz','tone_p2_1vsrest','AEP_p2n1_int'}
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
        pp=plot(this_x(ss),this_y(ss),'.','MarkerEdgeColor',[0.2 0.2 0.2],'MarkerSize',2);
        %cmap(uheal_nh.Age(this_pt(ss))-17,:)
    end
    %scatter(zscore(uheal_nh.(fnames{ii})(thisidx)),zscore(uheal_nh.(fnames{jj})(thisidx)),'k.')
    plot([xl],[yl],'color',[0.8 0.8 0.8])
    xlim([-3 3])
    ylim([-3 3])
    if jj+(ii-1)*length(fnames) == 35
        set(gca,'fontsize',12)
    else
        set(gca,'xticklabels',{},'yticklabels',{})
    end
end
end

%
diag_idx = find(eye(length(fnames),length(fnames)));
cidx = [];
for ii=1:length(diag_idx)
    cidx = [cidx diag_idx(ii)+1:length(fnames)*ii]
end

for tt=1:length(cidx)
subplot(length(fnames),length(fnames),cidx(tt))
hold off
plot(0,1)
axis off

end
%
diag_idx = find(eye(length(fnames),length(fnames)));
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
        %text(0,0,['p = ' num2str(round(p_corr(pidx(pp)),3)) '\newline' 'r = ' num2str(round(rho(pidx(pp)),3))])
        xlim([-1 1])
        ylim([-1 1])
        axis off
    end
end
%
pnames = {'AP_{ }','FFR_{ }','ITPC_{Q}','AEP_{\mu}','P2_{adapt}','P2-N1_{int}'}
%yidx = [1,6,11,16,21]
yidx = []
for ii=1:length(pnames)
yidx =[yidx (ii-1)*length(pnames)+1];
end
for ii=1:length(fnames)
    subplot(length(fnames),length(fnames),ii)
    t = title(pnames{ii},'FontSize',11,'FontName','Arial','fontweight','Bold','Rotation',0,'VerticalAlignment','bottom')
    subplot(length(fnames),length(fnames),yidx(ii))
%     if ii==4
%         ylab = ylabel('Neg. \newline 4Hz')
%     else
         ylab = ylabel(pnames{ii})
%     end
    ylab.FontSize =11;
    ylab.FontName = 'Arial';
    ylab.FontWeight = 'Bold'
    ylab.Rotation = 0;
    ylab.HorizontalAlignment = 'Right'
    ylab.VerticalAlignment = 'middle'
    ylab.Position = ylab.Position-[0.6 0 0]
    ylab.Color = [0 0 0];
    %tmp.YLabel.Position = tmp.YLabel.Position-[0 0 1];
    
end

fig=gcf;
set(fig,'renderer','painters')
%cd(root

set(gcf,'position',[133.3333 133.3333 614.6667 533.6667])
set(gcf, 'papersize',[10 20]);
%                PaperSize: [8.5000 11]
              %  PaperType: 'usletter'
              % PaperUnits: 'inches'
saveas(fig,[savepath 'corr_eeg'],'svg')
%%
clc
cnames = {'AP:FFR','AP:ITPC','AP:AEPmean','AP:P2adapt','AP:P2N1',...
          'FFR:ITPC','FFR:AEPmean','FFR:P2adapt','FFR:P2N1', ...
          'ITPC:AEPmean','ITPC:P2adapt','ITPC:P2N1',...
          'AEPmean:P2adapt','AEPmean:P2N1',...
          'P2adapt:P2N1'}
for ii=1:length(pidx)
    disp([cnames{ii} '     :r = ' num2str(rho(pidx(ii))) ', p = ' num2str(p_corr(pidx(ii)))])
end
%% correlation with questionaires
uheal_nh.nesi_log = log10(uheal_nh.nesi); uheal_nh.nesi_log(isinf(uheal_nh.nesi_log))=nan;
clinnames = {'Age','nesi_log','tts','ssq12_mean','PTA_lf','PTA_hf','memr_slope'}
eegnames = {'AP_amp_pm','FFR_SNR','ITPC_ratio','Neg_4Hz','tone_p2_1vsrest','AEP_p2n1_int'};
peegnames = {'AP_{ }','FFR_{ }','ITPC_{Q}','Mean_{ }','P2_{adapt}','P2-N1_{ }'};
pclinnames = {'Age','NESI','TTS','SSQ12','PTA_{lf}','PTA_{hf}','MEMR'}
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

        subplot(length(clinnames),length(eegnames),cc)
        scatter(zscore(uheal_nh.(clinnames{ii})(thisidx)),zscore(uheal_nh.(eegnames{jj})(thisidx)),'.k')
        lsline
        xlim([-3 3])
        ylim([-3 3])
        box on
        if cc==length(clinnames)*length(eegnames)
            set(gca,'yaxislocation','right')
        else
        set(gca,'xticklabel',{},'yticklabel',{})
        end
    end
end
for ii=1:length(peegnames)
    subplot(length(clinnames),length(eegnames),ii)
    title(peegnames{ii})
end
cc = 0;
yidx = []
for ii=1:length(clinnames)
    yidx = [yidx (ii-1)*length(eegnames)+1]
end
for ii=1:length(yidx)
    cc = cc+1
    subplot(length(clinnames),length(eegnames),yidx(ii))
    ylab = ylabel(pclinnames{cc})
    ylab.Rotation = 0;ylab.HorizontalAlignment = 'right'
    set(gca,'FontWeight','bold')
    set(get(gca,'YLabel'),'Visible','on')
end

fig = gcf;
saveas(fig,[savepath 'clin_corr_scatter'],'svg')
% corr plot

p_plot = p';p_plot = p_plot(:);p_plot = p_plot*length(p_plot);
rho_plot = rho';rho_plot = rho_plot(:)
figure('renderer','painters')
%set(gcf,'position',[119 784 356 311])
for pp = 1:length(p_plot)

    subplot(length(clinnames),length(eegnames),pp)
    
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
for ii=1:length(eegnames)
    subplot(length(clinnames),length(eegnames),ii)
    t = title(peegnames{ii})
end
cc = 0;
idx_i = [];
for ii=1:length(clinnames)
    idx_i = [idx_i 1+(ii-1)*length(eegnames)]
end
for ii=idx_i
    cc = cc+1
    subplot(length(clinnames),length(eegnames),ii)
    ylab = ylabel(pclinnames{cc})
    ylab.Rotation = 0;ylab.HorizontalAlignment = 'right'
    set(gca,'FontWeight','bold')
    set(get(gca,'YLabel'),'Visible','on')
    
end
set(gcf,'position',[603 206 456 420])
fig = gcf;
saveas(fig,[savepath 'clin_corr_pr'],'svg')
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
