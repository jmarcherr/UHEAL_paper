% script for plotting clincial measures 
% audiogram
% MEMR
% TEOAE
% Questionaires and clincial outcomes
close all
clear all
cd('/work3/jonmarc/UHEAL_paper')
UHEAL_startup


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
%m_idx = find(uheal_data.Age>28 & uheal_data.Age<50 & ~uheal_data.CP_new)
%o_idx = find(uheal_data.Age>=50 & ~uheal_data.CP_new);
% new CP subgroups
CP_X = (mean(uheal_data.aud(:,1:6),2)<=20 & all(uheal_data.aud(:,1:6)<=15,2) & all(uheal_data.aud(:,1:6)>=-5,2));
uheal_data.CP_X =CP_X;
% new groups
y_idx = find(uheal_data.Age<=25 & uheal_data.CP_X); 
m_idx = find(uheal_data.Age>26 & uheal_data.Age<50 & uheal_data.CP_X)
o_idx = find(uheal_data.Age>=50 & uheal_data.CP_X);

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
savepath = '/work3/jonmarc/UHEAL_paper/_clin/figs/'
save_figs = 1;
%%
for plot_on=1
    if plot_on
        close all
        %% All audiograms with age-colorbar
        idx = find(~isnan(uheal_data.Age));
        fig = plot_audiogram_groups(idx,uheal_data.Age,uheal_data.aud,uheal_data.audfreq,cmap)
        % save figure
        title(['all, n=' num2str(length(idx))])
        if save_figs
            saveas(fig,[savepath 'audiogram_all'],'svg')
        end
        %% only NH
        close all
        NH_idx = find(uheal_data.CP_X==1);
        fig = plot_audiogram_groups(NH_idx,uheal_data.Age,uheal_data.aud,uheal_data.audfreq,cmap)
        % save figure
        %title(['NH, n=' num2str(length(NH_idx))])
        if save_figs
            saveas(fig,[savepath 'audiograms_NH_subg'],'svg')
        end

        hold on
        py=plot(uheal_data.audfreq(1,:),nanmean(uheal_data.aud(y_idx,:)),'^','color',y_col,'linewidth',2)%,'markerfacecolor',[cmap(13-5,:)])
        pm=plot(uheal_data.audfreq(1,:),nanmean(uheal_data.aud(m_idx,:)),'sq','color',m_col,'linewidth',2)%,'markerfacecolor',[cmap(30,:) 0.5])
        po=plot(uheal_data.audfreq(1,:),nanmean(uheal_data.aud(o_idx,:)),'o','color',o_col,'linewidth',2)%,'markerfacecolor',[cmap(end-7,:) 0.5])
        hleg=legend([py,pm,po],{'Young','Middle-aged','Older'});
        hleg.Box='off';
        hleg.Position = [0.2300 0.2887 0.3935 0.2475];
        xtickangle(0)
        set(gca,'xtick',[250,500,1000,2000,4000,8000,16000],'xticklabel',{'.25','.5','1','2','4','8','16'})
        %set(gcf,'position',[[305 412 478 299]])
        set(gcf,'position',[305 473 422 225])
        if save_figs
            saveas(fig,[savepath 'audiograms_NH_groups_subg'],'svg')
        end

        %% NH seperate
        
        %% only HI
        close all
        HI_idx = find(uheal_data.CP_new==1);
        fig = plot_audiogram_groups(HI_idx,age_sub,aud,freq_aud,cmap)

        % save figure
        title(['HI, n=' num2str(length(HI_idx))])
        if save_figs
            saveas(fig,['figs/figs_NH/audiograms_HI_subg_' date],'svg')
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %% age plots (only NH)
        % all
        close all
        figure('renderer','painter')
        [N,X]=hist(age_sub(NH_idx),max(age_sub(NH_idx))-min(age_sub(NH_idx)))
        Bh = bar(X,N,'facecolor',[0.5 0.5 0.5],'EdgeColor', [0.5 0.5 0.5],'BarWidth',0.5);
        xlabel('Age')
        ylabel('n')
        ax = ancestor(gca, 'axes');
        xrule = ax.XAxis;

        hold on

        plot([25 25],[0 8],'k--')
        plot([56 56],[0 8],'k--')
        set(gca,'fontsize',16)
        xlim([16 70])
        xrule.FontSize = 14;
        set(gcf,'position',[305 403 191 299])
        box off
        fig= gcf
        if save_figs
            saveas(fig,['figs/figs_NH/age_hist_all_subg_' date],'svg')
        end

        %% age plots sex divded (only NH)
        % all
        close all
        figure('renderer','painter')
        %subplot(1,2,1)
        maleidx = find(gender_sub==2 & uheal_data.CP_X==1);
        femaleidx =find(gender_sub==1 & uheal_data.CP_X==1);

        [N,X]=hist(age_sub(maleidx),5);
        [Nf,Xf] = hist(age_sub(femaleidx),5);

        Bh = bar([Xf],[N;Nf],'stacked','BarWidth',.8);
        Bh(1).FaceColor = [0.5 0.5 0.5];Bh(1).EdgeColor = [0 0 0];
        Bh(2).FaceColor = [1 1 1];Bh(2).EdgeColor = [0 0 0];
        xlabel('Age')
        ylabel('Participants')
        ax = ancestor(gca, 'axes');
        xrule = ax.XAxis;
        hold on

        set(gca,'fontsize',12)
        xlim([16 70])
        xrule.FontSize = 12;
        %set(gcf,'position',[305 403 191 299])
        set(gcf,'position',[305 403 227 225])
        hleg = legend('Male','Female');
        hleg.FontSize = 10;
        hleg.Box = 'off';
        hleg.Position = [0.4634    0.7620    0.5236    0.1288];
        box off
        ylim([0 45])

        
        fig= gcf
        if save_figs
            saveas(fig,['figs/figs_NH/age_hist_male_female_subg' date],'svg')
        end

        %% new dist


        %% grouped male female
        figure(3)
        Y   = gender_sub(find(age_sub<=25 & ~uheal_data.CP_new));
        O1  = gender_sub(find(age_sub>28 & age_sub<=50 & ~uheal_data.CP_new));
        O2  = gender_sub(find(age_sub>50 & ~uheal_data.CP_new));

        b1=bar([1:3],[length(find(Y==1)) length(find(O1==1)) length(find(O2==1));...
            length(find(Y==2)) length(find(O1==2)) length(find(O2==2))]','stacked')

        hold on
        xlabel('Group')
        xtickangle(45);
        ylim([0 70])
        ylabel('Nr. of participants')
        set(gca,'xtick',[1:3],'xticklabels',{'Young','Middle-aged','Older'},'fontsize',12)
        hleg = legend(b1,'Female','Male')
        b1(1).FaceColor = [0.2 0.2 0.5];
        b1(2).FaceColor = [0.5 0.7 0.5];
        hleg.Position = [0.6188 0.7419 0.2824 0.1689];
        set(gcf,'position',[305 412 432 299])
        fig= gcf
        if save_figs
            saveas(fig,['figs/age_hist_gender_subg' date],'svg')
        end
        %% grouped male female and HI

        HI_idx = (uheal_data.CP_new==1)

        figure(4)
        Y   = gender_sub(find(age_sub<=25));
        O1  = gender_sub(find(age_sub(~HI_idx)>25 & age_sub(~HI_idx)<=60));
        O1HI = gender_sub(find(age_sub(HI_idx)>25 & age_sub(HI_idx)<=60));
        O2  = gender_sub(find(age_sub(~HI_idx)>60));
        O2HI = gender_sub(find(age_sub(HI_idx)>60 ));

        b1=bar([1:5],[length(find(Y==1)) length(find(O1==1)) length(find(O2==1)) length(find(O1HI==1)) length(find(O2HI==1));...
            length(find(Y==2)) length(find(O1==2)) length(find(O2==2)) length(find(O1HI==2)) length(find(O2HI==2))]','stacked')

        hold on
        xlabel('Group')
        ylabel('Nr. of participants')
        set(gca,'xtick',[1:5],'xticklabels',{'Y','O1','O2','O1HI','O2HI'},'fontsize',16)
        legend(b1,'Female','Male')
        xtickangle(45)
        b1(1).FaceColor = [0.2 0.2 0.5]
        b1(2).FaceColor = [0.5 0.7 0.5]
        set(gcf,'position',[305 412 432 299])
        ylim([0 max(b1(1).YData)+max(b1(2).YData)+5]);
        ymax = max(b1(1).YData)+max(b1(2).YData)+5;
        fig= gcf
        if save_figs
            saveas(fig,['figs/age_hist_gender_HI' date],'svg')
        end

        %% grouped male female and HI from hvidovre
        for i=1

            % find who has been to hvidovre (from who-what-when 10-05-21)
            DRCMR_idx = zeros(length(gender_sub),1);
            DRCMR_idx = find(uheal_data.DRCMR),

            %
            O_idx = intersect(find(~HI_idx),find(DRCMR_idx));
            OHI_idx =intersect(find(HI_idx),find(DRCMR_idx));


            figure(5)
            Y   = gender_sub(find(age_sub(DRCMR_idx)<=25));
            O1  = gender_sub(find(age_sub(O_idx)>25 & age_sub(O_idx)<=60));
            O1HI = gender_sub(find(age_sub(OHI_idx)>25 & age_sub(OHI_idx)<=60));
            O2  = gender_sub(find(age_sub(O_idx)>60));
            O2HI = gender_sub(find(age_sub(OHI_idx)>60 ));


            b1=bar([1:5],[length(find(Y==1)) length(find(O1==1)) length(find(O2==1)) length(find(O1HI==1)) length(find(O2HI==1));...
                length(find(Y==2)) length(find(O1==2)) length(find(O2==2)) length(find(O1HI==2)) length(find(O2HI==2))]','stacked')

            hold on
            xlabel('Group')
            ylabel('Nr. of participants')
            set(gca,'xtick',[1:5],'xticklabels',{'Y','O1','O2','O1HI','O2HI'},'fontsize',16)
            legend(b1,'Female','Male')
            xtickangle(45)
            b1(1).FaceColor = [0.2 0.2 0.5]
            b1(2).FaceColor = [0.5 0.7 0.5]
            %b1(3).FaceColor = [0.2 0.2 0.6];
            %b1(4).Facecolor = [0.5 0.7 0.4];
            set(gcf,'position',[305 412 432 299])
            ylim([0 ymax]);
            fig= gcf
            if save_figs
                saveas(fig,['figs/age_hist_gender_DRCMR' date],'svg')
            end
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% MEMR plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots average MEMR data and linear coefficient as function of age. Can
% be left out if only extraction of data is wanted (plot_on=0, save_figs=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_figs = 0;

% for plotting
reflex_sub = uheal_data.memr_reflex;
growth_sub_alt = uheal_data.memr_reflex_growth;
levels = uheal_data.memr_levels;
f_center = uheal_data.memr_fcenter;
age_sub = uheal_data.Age;
MEM_slope = uheal_data.memr_slope;

for plot_on=1
    if plot_on

        close all
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
        xlabel('Elicitor level [dB]')
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
            axis([250 8000 -0.08 0.08])
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
            xlabel('Frequency [Hz]');
            ylabel('\Delta Absorbance');
        end
        set(gcf,'Position',[228 420 618 209]);
        hleg.Position= [0.6509 0.1470 0.2589 0.8038];
        xtickangle(0)
        fig= gcf
        if save_figs
            saveas(fig,['figs/figs_NH/memr_summary' date],'svg')
        end

        %% MEMR slope

        figure('Renderer','painter')
        MEM_slope(42,:) = nan;
        jit = randn(size(MEM_slope,1)).*0.2;
        for s=1:length(MEM_slope)
            subplot(1,3,1)
            %plot(1+jit(s),db(FFR(s))','o','markerfacecolor',[cmap(age(s)-17,:)],'color',[cmap(age(s)-17,:)])
            plot(1+jit(s),MEM_slope(s,1)','o','markerfacecolor',[0.6 0.6 0.6],'color','k')
            xlim([0 2])
            hold on
            xlabel('')
            ylabel('Linear coefficient')
            set(gca,'xtick',[],'fontsize',12,'Xcolor','none')
            box off
            ylim([-0.02 0.05])

        end
        %close all
        figure('Renderer','painter')
        %subplot(1,3,[2 3])
        scatter(age_sub(NH_idx),MEM_slope(NH_idx,1),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
        set(gca,'fontsize',12)
        ylabel('Linear coefficient')

        xlabel('Age')
        ylim([-0.02 0.05])
        %set(gcf,'position',[228 839 432 209])
        %set(gcf,'position',[228 416 441 609])
        set(gcf,'Position',[228 420 618/3 209]);
        ll=lsline
        set(ll,'linewidth',2,'color','k')

        fig= gcf
        if save_figs
            saveas(fig,'figs/figs_NH/memr_age','svg')
        end
    end
end



%% TEOAE plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots average TEOAE data. Can be left out if only extraction of data 
% is wanted (plot_on=0, save_figs=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

teoae_sub_amp = uheal_data.teoae_amp;
teoae_SNR = uheal_data.teoae_SNR;
teoae_sig = uheal_data.teoae_sig;
age_sub =uheal_data.Age;

for plot_on=0
    if plot_on
        % colormap
        cmap_all = [];
        % young
        cmap = cbrewer('seq','Greys',(25-17)+5);
        cmap(1:5,:) = [];
        cmap_all = [cmap_all; cmap];
        % middle aged
        cmap = cbrewer('seq','Blues',(60-25)+5);
        cmap(1:5,:) = [];
        cmap_all = [cmap_all ; cmap];
        % older
        cmap = cbrewer('seq','Reds',(max(ages)-(60))+5);
        cmap(1:5,:) = [];
        cmap_all = [cmap_all ; cmap];

        % combine
        cmap = cmap_all;
        cmap(find(cmap>1))=1;cmap(find(cmap<0))=0;
        ages = [18 77];
        ageidx = linspace(min(ages),max(ages),size(cmap,1));
        %% mean

        close all
        % valid measurements
        figure('Renderer','painter')
        idx = find(~isnan(teoae_sub_amp(:,1,1,1)) & ~uheal_data.CP_new);
        bb=bar(squeeze(nanmean(teoae_sub_amp(idx,:,1),1)),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
        hold on
        sig = zeros(size(teoae_sub_amp(idx,:,:,:),1),5);
        jit_age = linspace(-0.15,0.15,60);
        for s=1:size(teoae_sub_amp(idx,:,:,:),1)
            %jit_s = randn*0.1;
            for kk=1:5

                if  teoae_SNR(s,kk)>=6
                    p_sig=plot(kk+jit_age(age_sub(s)-17),squeeze(teoae_sub_amp(idx(s),kk,1)),'o',...
                        'color',cmap(age_sub(s)-17,:),'markerfacecolor',cmap(age_sub(s)-17,:),...
                        'markersize',4)
                    sig(s,kk)=1;

                else
                    p_nonsig=plot(kk+jit_age(age_sub(s)-17),squeeze(teoae_sub_amp(idx(s),kk,1)),'+',...
                        'color',cmap(age_sub(s)-17,:),'markersize',4)

                end
            end
        end

        % area under the curve instead
        %nf= area(1:5,nanmean(teoae_sub_amp(:,:,2),1),-40,'EdgeColor','none','FaceColor',[0.5 0.5 0.5],'ShowBaseLine','off')
        %nf=errorbar(1:5,nanmean(teoae_sub_amp(:,:,2),1),nanstd(teoae_sub_amp(:,:,2),1)/sqrt((length(idx))),'color',[1 0.5 0.5]);

        %ylim([-35 30])
        %plot(0:6,ones(1,7)*-10,'--k')
        % %significant
        for kk=1:5
            this_sig = length(find(sig(:,kk)))/size(sig,1);
            text(kk-0.1,27,[num2str(round(this_sig,2)*100) '%'])
        end
        set(gca,'fontsize',16)
        ylim([-30 24])
        p_sig=plot(1,-60,'ko');
        p_nonsig = plot(1,-60,'k+');
        hleg = legend([p_sig,p_nonsig],{'Significant','N.S','Noise Floor'})
        hleg.Position = [0.2044    0.2280    0.3013    0.1829];
        hleg.Box = 'off';
        hleg.FontSize = 12;

        ylabel('TEOAE (dB SPL)')
        xlabel('Frequency bands (+/- .5 kHz)')
        set(gcf,'Position',[441 386 468 339])
        box off
        fig= gcf
        if save_figs
            saveas(fig,['figs/figs_NH/teoae_overview_amp' date],'svg')
        end


        %uheal_data.teoae = teoae_amp;
        %% alternative (SNR)
        close all

        % valid measurements
        figure('Renderer','painter')
        idx = find(~isnan(teoae_sub_amp(:,1,1,1)) & ~uheal_data.CP_new);
        bar(squeeze(nanmean(teoae_SNR(idx,:),1)),'FaceColor',[0.5 0.5 0.5],'EdgeColor',[0.5 0.5 0.5])
        hold on
        sig = zeros(size(teoae_sub_amp(idx,:,:,:),1),5);
        for s=1:size(teoae_SNR(idx,:),1)
            jit_s = randn*0.1;
            for kk=1:5
                if teoae_SNR(idx(s),kk)>=6
                    p_sig=plot(kk+jit_age(age_sub(s)-17),squeeze(teoae_SNR(idx(s),kk)),'o',...
                        'color',cmap(age_sub(s)-17,:),'markerfacecolor',cmap(age_sub(s)-17,:),...
                        'markersize',4)
                    sig(s,kk)=1;
                else
                    p_sig=plot(kk+jit_age(age_sub(s)-17),squeeze(teoae_SNR(idx(s),kk)),'+',...
                        'color',cmap(age_sub(s)-17,:),'markerfacecolor',cmap(age_sub(s)-17,:),...
                        'markersize',4)
                    sig(s,kk)=0;
                end
            end
        end

        %nf=errorbar(1:5,nanmean(teoae_sub_amp(:,:,2),1),nanstd(teoae_sub_amp(:,:,2),1)/sqrt((length(idx))),'r')
        ylim([-15 42])
        plot(0:6,ones(1,7)*6,'--k')
        % %significant
        for kk=1:5
            this_sig = length(find(sig(:,kk)))/size(sig,1);
            text(kk-0.1,39,[num2str(round(this_sig,2)*100) '%'])
        end
        set(gca,'fontsize',16)
        p_sig=plot(1,-60,'ko');
        p_nonsig = plot(1,-60,'k+');
        hleg = legend([p_sig,p_nonsig],{'Significant','N.S'})

        hleg.Position = [0.1830 0.1956 0.3013 0.1829];
        hleg.Box = 'off';
        hleg.FontSize = 12;

        ylabel('TEOAE SNR (dB)')
        xlabel('Frequency bands (+/- .5 kHz)')
        box off
        set(gcf,'Position',[441 386 468 339])

        fig= gcf
        if save_figs
            saveas(fig,'figs/figs_NH/teoae_overview_SNR','svg')
        end
        %% mean wave form
        close all
        figure('Renderer','painter')
        %plot(squeeze(teoae_sub_resp(1,:,1)),teoae_sub_resp(:,:,3),'color',[0.5 0.5 0.5])
        hold on
        plot(squeeze(teoae_sub_resp(1,:,1)),nanmean(teoae_sub_resp(:,:,3),1),'k')
        shadedErrorBar(squeeze(teoae_sub_resp(1,:,1)),nanmean(teoae_sub_resp(:,:,3),1),nanstd(teoae_sub_resp(:,:,3),1)/sqrt(size(teoae_sub_resp,1)))
        set(gca,'fontsize',16)
        set(gcf,'Position',[441 552 468 173])
        xlabel('Time [ms]')
        ylabel('mPA')
        xlim([4 12.5])
        box off
        %plot(squeeze(teoae_sub_resp(1,:,1)),teoae_sub_resp(:,:,2))
        fig= gcf
        if save_figs
            saveas(fig,['figs/teoae_waveform' date],'svg')
        end
    end
end

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