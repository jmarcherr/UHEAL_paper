%% plot clinical measures

close all
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/zhome/7e/f/64621/Desktop/UHEAL_paper/UHEAL_startup.m')
clindir = '/work3/jonmarc/UHEAL_master/UHEAL_paper/_clin/_clindata';
freq_aud = [250 500 1000 2000 4000 8000 9000 10000 11200 12500 14000 16000];
d=dir(fullfile([clindir filesep 'UH*.mat']));

% plot?
%params
plot_on=0;      % 1= plot
save_figs = 0;  % 1= save figs
%% get clin data

for s=1:length(d)
   results{s} =  get_data(d,s);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% generate UHEAL_DATA mat %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%init
% get subject parameters 
for s=1:length(d)
    age_sub(s) = results{s}.age_sub;
    aud(s,:,:) = results{s}.aud;
    aud_L(s,:,:) = results{s}.aud_L;
    aud_R(s,:,:) = results{s}.aud_R;
    aud_freq(s,:) = results{s}.aud_freq;
    PTA_HF(s) = nanmean(results{s}.aud(7:end)); % PTA 9-16k
    %PTA_LF(s) = nanmean(results{s}.aud(1:6));   % PTA 250-8k
    PTA_LF(s) = nanmean(results{s}.aud(1:4)); % PTA 250-2k
    CP_sub(s,:) = results{s}.CP_sub;
    gender_sub(s,:) = results{s}.gender_sub;
    HV_sub(s) = results{s}.HV_sub;
    subid(s) = str2double(results{s}.sub_id(3:5));
end
%
%load('/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/scraped/uheal_data_table/uheal_data.mat')
%this_dir = cd;
uheal_data = struct;
uheal_data.subid = subid';
uheal_data.Age = age_sub';
uheal_data.CP = CP_sub; % CP from excel 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% new HI criteria

%from carcagno 2020
%Selection criteria included audiometric thresholds for both ears below 20
%dB HL at octave frequencies from 0.125 to 2 kHz, and below 40 dB HL at 4
%kHz. No selection criteria were imposed for frequencies above 4 kHz.

for ii=1:4 % .125-2k
    aud_tmp = aud(:,ii);
    crit_l(ii,:) = aud_tmp<=20;
end
crit_4 = aud(:,5)<=40;
crit_all = [crit_l' crit_4];
crit_sum = sum(crit_all,2);


find(crit_sum==5);
uheal_data.CP_new = ones(1,117);
uheal_data.CP_new(find(crit_sum==5))=0;
uheal_data.CP_new = uheal_data.CP_new';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subgroup CP for more strict analysis
% all thresholds from .25-8kHz -5 db HL=><=15 dB HL 
CP_SG = double(mean(aud(:,1:6),2)<=20 & ...
    all(aud(:,1:6)<=15,2) & ...
    all(aud(:,1:6)>=-5,2));
uheal_data.CP_SG = CP_SG;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uheal_data.gender = gender_sub;
uheal_data.DRCMR = HV_sub';
uheal_data.aud = aud;
uheal_data.audL = aud_L;
uheal_data.audR = aud_R;
uheal_data.audfreq = aud_freq;
uheal_data.PTA_lf = PTA_LF';
uheal_data.PTA_hf = PTA_HF';

% MEMR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract MEMR values to uheal_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for s=1:length(d)
    reflex_sub(s,:,:) = results{s}.memr.reflex_sub;
    growth_sub_alt(s,:) = results{s}.memr.growth_sub;
    MEM_slope(s,:) =results{s}.memr.MEM_slope;
    disp(num2str(s))
end
    levels = results{1}.memr.levels;
    f_center = results{1}.memr.f_center;
    freq = results{1}.memr.freq;
    
 uheal_data.memr_slope = MEM_slope(:,1);  
 uheal_data.memr_slope(42) = nan; % bad data for sub 42
 uheal_data.memr_reflex_growth = growth_sub_alt;
 uheal_data.memr_reflex = reflex_sub;
 uheal_data.memr_fcenter = f_center;
 uheal_data.memr_levels = levels;

% TEOAE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract TEOAE values to uheal_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
for s=1:length(d)
    % Load scaped data
    teoae_sub_amp(s,:,:) = results{s}.teoae.amp;
    teoae_sub_resp(s,:,:,:) = results{s}.teoae.resp;
    for kk=1:5
        teoae_amp(s,kk) = squeeze(teoae_sub_amp(s,kk,1));
        teoae_SNR(s,kk) = teoae_amp(s,kk)-squeeze(teoae_sub_amp(s,kk,2));
        if teoae_SNR(s,kk)>=6 % 6 dB SNR criteria
            teoae_sig(s,kk)=1;
        else
            teoae_sig(s,kk)=0;
        end
    end
end


uheal_data.teoae_amp = teoae_amp;
uheal_data.teoae_sig = teoae_sig;
uheal_data.teoae_SNR = teoae_SNR;

% Reverse digit span, nesi, ssq
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract clincal values to uheal_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for s=1:length(d)
    % rds
    if ~isfield(results{s}.rds,'Fresp')
        fds(s) = nan;
    else
        fds(s) = sum([results{s}.rds.Fresp.Total_score{:}]);
    end
    if ~isfield(results{s}.rds,'Bresp')
        rds(s) = nan;
    else
        rds(s) = sum([results{s}.rds.Bresp.Total_score{:}])
    end
    %nesi
    if isempty(results{s}.nesi)
        nesi(s) = nan;
    else
    nesi(s) = results{s}.nesi;
    end
    if isempty(results{s}.tts);
        tts(s) = nan;
    else
        tts(s)= results{s}.tts;
    end
    % ssq
    if isempty(results{s}.ssq)
        ssq(s,:) = nan(1,12);
    else
        ssq(s,:) = results{s}.ssq(:,2);
        %ssq(s,:)= normalize(results{s}.ssq(:,2));
    end
end
% save to uheal struct
uheal_data.rds = rds';
uheal_data.fds = fds';
uheal_data.nesi =nesi';
uheal_data.tts = tts';
uheal_data.ssq12_mean = nanmean(ssq,2);
uheal_data.ssq12 = ssq;

% Extract ACALOS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract ACALOS values to uheal_data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fs = [500,1000,2000,4000];
x= [-10:2.5:120];
%addpath('/work1/jonmarc/UHEAL_master/UHEAL/_scripts/_tools/ACALOS')
for s=1:length(d)
    if ~isempty(results{s}.acalos);
        tmp_data = results{s}.acalos.raw.RawData;
        freq_idx = tmp_data(:,1);
        
        for ff=1:4
        this_f = find(freq_idx==fs(ff));
        acalos_f{ff} = tmp_data(this_f,:,:);
        dataLoudFit = reshape(acalos_f{ff}(:,2:3)', 1, []);
        
        select_fitting_function = 'BTUX';
        
        [fitparams,rawData] = fit_loudness_function(dataLoudFit,select_fitting_function);
        data_tmp(ff,2)=fs(ff);
        data_tmp(ff,[3,4,5])=fitparams;
        %% Estimated parameters
        
        [y, failed] = loudness_function_bh2002(x, fitparams);
        iZero = (max(find(y==0)));
        iUCL = (min(find(y==50)));
        %%%%%%%%%%%%%%%%%%%%% HTL
        if  ~isempty(iZero)
            HTL(ff)=x(iZero+1);
        else
            if ~isempty(x(min(find(y<=1))+1))
                HTL(ff) =  x(min(find(y<=1))+1);
            else
                HTL(ff) = x(find(min(y)));
            end
        end
        %%%%%%%%%%%%%%%%%%%%% UCL
        if  ~isempty(iUCL)
            UCL(ff)=x(iUCL-1);
        else
            UCL(ff) = x(y==max(y));
        end
        
        Lcut =fitparams(1);
        m_lo = fitparams(2);
        %        HTL_A(ii) = Lcut - 22.5/m_lo;
        HTL_A(ff)= fitparams(1)-22.5/m_lo;
        MCL(ff) = x(max(find(y<=25)));
        
        y_plot(:,ff) = y;
        x_plot(:,ff) = x;
        
    end
    % save results
    acalos_yplot(s,:,:) = y_plot;
    HTL_A_sub(s,:) = HTL_A;
    MCL_sub(s,:) = MCL;
    AC_slope(s,:) = data_tmp(:,4)';
    L_cut(s,:) = data_tmp(:,3)';
    m_high(s,:) = data_tmp(:,5)';
    %acalos_ = results{s}.acalos.proc;
    else
        acalos_yplot(s,:,:) = nan(size(acalos_yplot(1,:,:)));
        HTL_A_sub(s,:) = nan(size(HTL_A_sub(1,:)));
        MCL_sub(s,:) = nan(size(MCL_sub(1,:)));
        AC_slope(s,:) = nan(size(AC_slope(1,:)));
        L_cut(s,:) = nan(size(L_cut(1,:)));
        m_high(s,:) = nan(size(m_high(1,:)));
    end
    clc
    disp(['subject ' num2str(s) ' processed'])
end

uheal_data.acalos_AC_slope = AC_slope;
uheal_data.acalos_MCL = MCL_sub;
uheal_data.acalos_L_cut = L_cut;
uheal_data.acalos_m_high = m_high;
uheal_data.acalos_HTL_A_sub = HTL_A_sub;

% remove fields
uheal_data = rmfield(uheal_data,'DRCMR');
uheal_data = rmfield(uheal_data,'CP');
%% save uheal data
savedir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_clin/clin_data_table/'

save([savedir 'clin_data.mat'],'uheal_data')

 
%% % plotting %%%%


%% Audiogram plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots audiograms vs. age as well as age-distribution in various ways. Can
% be left out if only extraction of data is wanted (plot_on=0, save_figs=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load('uheal_data_table_current/uheal_data.mat')
cd('/work1/jonmarc/UHEAL_master/UHEAL_paper/_clin/_clindata')
figdir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_clin/figs'
ages = [18 77];
% groups
y_idx = find(uheal_data.Age<=25 & ~uheal_data.CP_new); 
m_idx = find(uheal_data.Age>25 & uheal_data.Age<50 & ~uheal_data.CP_new)
o_idx = find(uheal_data.Age>=50 & ~uheal_data.CP_new);

save_figs = 1
%%
for plot_on=1
    if plot_on
        close all
        uheal_colormap
        
        %% All audiograms with age-colorbar
        idx = find(~isnan(age_sub));
        fig = plot_audiogram_groups(idx,age_sub,aud,freq_aud,cmap)
        % save figure
        title(['all, n=' num2str(length(idx))])
        if save_figs
            saveas(fig,[figdir '/audiogram_all_' date],'epsc')
        end
        %% only NH
        close all
        NH_idx = find(uheal_data.CP_new==0);
        fig = plot_audiogram_groups(NH_idx,age_sub,aud,freq_aud,cmap)
        % save figure
        %title(['NH, n=' num2str(length(NH_idx))])
        if save_figs
            saveas(fig,[figdir '/audiograms_NH_' date],'epsc')
        end

        hold on
        py=plot(freq_aud,nanmean(aud(y_idx,:)),'o','color',y_col,'linewidth',2)%,'markerfacecolor',[cmap(13-5,:)])
        pm=plot(freq_aud,nanmean(aud(m_idx,:)),'^','color',m_col,'linewidth',2)%,'markerfacecolor',[cmap(30,:) 0.5])
        po=plot(freq_aud,nanmean(aud(o_idx,:)),'sq','color',o_col,'linewidth',2)%,'markerfacecolor',[cmap(end-7,:) 0.5])
        hleg=legend([py,pm,po],{'Young','Middle-aged','Older'});
        hleg.Box='off';
        hleg.Position = [0.2323 0.3188 0.3935 0.2475];
        xtickangle(0)
        set(gca,'xtick',[250,500,1000,2000,4000,8000,16000],'xticklabel',{'.25','.5','1','2','4','8','16'})
        
        if save_figs
            saveas(fig,[figdir '/audiograms_NH_groups' date],'epsc')
        end
        
        %% only HI
        close all
        HI_idx = find(uheal_data.CP_new==1);
        fig = plot_audiogram_groups(HI_idx,age_sub,aud,freq_aud,cmap)

        % save figure
        title(['HI, n=' num2str(length(HI_idx))])
        if save_figs
            saveas(fig,[figdir '/audiograms_HI' date],'epsc')
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
            saveas(fig,[figdir '/age_hist_all_' date],'epsc')
        end

        %% age plots sex divded (only NH)
        % all
        close all
        figure('renderer','painter')
        maleidx = find(gender_sub==2 & uheal_data.CP_new==0);
        femaleidx =find(gender_sub==1 & uheal_data.CP_new==0);

        [N,X]=hist(age_sub(maleidx),6);
        [Nf,Xf] = hist(age_sub(femaleidx),6);

        Bh = bar([Xf],[N;Nf],'stacked','BarWidth',.8);
        Bh(1).FaceColor = [0.5 0.5 0.5];Bh(1).EdgeColor = [0.5 0.5 0.5];
        Bh(2).FaceColor = [0.8 0.5 0.5];Bh(2).EdgeColor = [0.8 0.5 0.5];
        xlabel('Age')
        ylabel('Nr. participants')
        ax = ancestor(gca, 'axes');
        xrule = ax.XAxis;
        hold on

        set(gca,'fontsize',16)
        xlim([16 77])
        xrule.FontSize = 14;
        %set(gcf,'position',[305 403 191 299])
        set(gcf,'position',[305 403 227 299])
        hleg = legend('Male','Female');
        hleg.FontSize = 10;
        hleg.Box = 'off';
        hleg.Position = [0.4634    0.7620    0.5236    0.1288];
        box off
        fig= gcf
        if save_figs
            saveas(fig,[figdir '/age_hist_male_female_' date],'epsc')
        end

        %% grouped male female
        figure(3)
        Y   = gender_sub(find(age_sub<=25 & ~uheal_data.CP_new'));
        O1  = gender_sub(find(age_sub>25 & age_sub<=50 & ~uheal_data.CP_new'));
        O2  = gender_sub(find(age_sub>50 & ~uheal_data.CP_new'));

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
            saveas(fig,[figdir '/age_hist_gender' date],'epsc')
        end

    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%% MEMR plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots average MEMR data and linear coefficient as function of age. Can
% be left out if only extraction of data is wanted (plot_on=0, save_figs=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for plot_on=0
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
            axis([freq(1) freq(end) -0.08 0.08])
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
            %saveas(fig,['figs/figs_NH/memr_summary' date],'epsc')
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
            saveas(fig,'figs/figs_NH/memr_age','epsc')
        end
    end
end

%% TEOAE plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots average TEOAE data. Can be left out if only extraction of data 
% is wanted (plot_on=0, save_figs=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for plot_on=0
    if plot_on
        % colormap
                % new color map
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
        nf=errorbar(1:5,nanmean(teoae_sub_amp(:,:,2),1),nanstd(teoae_sub_amp(:,:,2),1)/sqrt((length(idx))),'color',[1 0.5 0.5])

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
        hleg = legend([p_sig,p_nonsig,nf],{'Significant','N.S','Noise Floor'})
        hleg.Position = [0.2044    0.2280    0.3013    0.1829];
        hleg.Box = 'off';
        hleg.FontSize = 12;

        ylabel('TEOAE (dB SPL)')
        xlabel('Frequency bands (+/- .5 kHz)')
        set(gcf,'Position',[441 386 468 339])
        box off
        fig= gcf
        if save_figs
            saveas(fig,['figs/figs_NH/teoae_overview_amp' date],'epsc')
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
            saveas(fig,'figs/figs_NH/teoae_overview_SNR','epsc')
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
            saveas(fig,['figs/teoae_waveform' date],'epsc')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Digit-span, NESI, TTS, SSQ plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots average clincial data. Can be left out if only extraction of data
% is wanted (plot_on=0, save_figs=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for plot_on=0
    if plot_on
        figure('Renderer','painter')
        scatter(age_sub,rds,'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
        ll=lsline
        title('Digit-span')
        xlabel('Age')
        ylabel('Digit-span score')
        set(ll,'linewidth',2,'color','k')
        set(gca,'fontsize',12)
        set(gcf,'position',[441 459 270 266])

        fig= gcf
        if save_figs
            saveas(fig,['figs/rds_age' date],'epsc')
        end

        figure('Renderer','painter')
        scatter(age_sub,nesi,'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
        ll=lsline
        title('NESI')
        xlabel('Age')
        ylabel('NESI score')
        set(ll,'linewidth',2,'color','k')
        set(gca,'fontsize',12)
        set(gcf,'position',[441 459 270 266])

        fig= gcf
        if save_figs
            saveas(fig,['figs/NESI_age' date],'epsc')
        end

        figure('Renderer','painter')
        for ii=1:12
            subplot(3,4,ii)
            scatter(age_sub,ssq(:,ii)','o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
            ll=lsline
            title(['SSQ-12 - ' num2str(ii)])
            xlabel('Age')
            ylabel('SSQ score')
            set(ll,'linewidth',2,'color','k')
            set(gca,'fontsize',12)
            set(gcf,'position',[184 106 1049 596])
        end
        fig= gcf
        if save_figs
            saveas(fig,['figs/ssq_age' date],'epsc')
        end

        % tts
        figure('Renderer','painter')
        scatter(age_sub,tts,'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
        ll=lsline
        title('TTS')
        xlabel('Age')
        ylabel('TTS score')
        set(ll,'linewidth',2,'color','k')
        set(gca,'fontsize',12)
        set(gcf,'position',[441 459 270 266])

        fig= gcf
        if save_figs
            saveas(fig,'figs/tts_age','epsc')
        end
    end
end

%% ACALOS plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plots average clincial data. Can be left out if only extraction of data
% is wanted (plot_on=0, save_figs=0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for plot_on =1
    if plot_on
        close all
        YNH = find(age_sub<25);
        ONH = find(age_sub>45);
        figure('Renderer','painter')
        for ii=1:4
            plot(x,squeeze(nanmean(acalos_yplot(NH_idx,:,ii),1)),'linewidth',2,'color',cmap(ii*10,:))
            hold on
        end
        xlim([0 120])
        xlabel('dB SPL')
        ylabel('CU')
        box off
        hleg = legend(num2str(fs(:)),'location','best');
        hleg.Box = 'off'
        set(gca,'fontsize',12)
        %set(gcf,'position',[150 641 216 359])
        set(gcf,'position',[150 781 232 219])
        fig= gcf
        if save_figs
            saveas(fig,['figs/ACALOS_all_2' date],'epsc')
        end
        %%
        figure('Renderer','painter')
        for ii=1:4
            subplot(1,4,ii)
            plot(x,squeeze(nanmean(acalos_yplot(YNH,:,ii))),'k')
            hold on
            plot(x,squeeze(nanmean(acalos_yplot(ONH,:,ii))),'r')
            xlim([0 105])
            xlabel('dB SPL')
            ylabel('CU')
            box off
            title([num2str(fs(ii)) ' Hz'])
            set(gcf,'position',[681 838 560 164])
        end
        hleg = legend('Young','Older');
        hleg.Box = 'off'
        fig= gcf
        if save_figs
            saveas(fig,['figs/ACALOS_yvo' date],'epsc')
        end

        figure('Renderer','painter')
        for ii=1:4

            subplot(1,4,ii)
            scatter(age_sub,AC_slope(:,ii),'o','markerfacecolor',[0.6 0.6 0.6],'MarkerEdgecolor','k')
            ll=lsline
            title([num2str(fs(ii)) ' Hz'])
            xlabel('age')
            ylabel('Slope (CU/dB)')
            ylim([0.2 0.8])
            box off
            set(ll,'linewidth',2,'color','k')
            set(gcf,'position',[681 838 560 164])
            hold on
            scatter(age_sub(find(CP_sub)),AC_slope(find(CP_sub),ii),'rx')
        end
        fig= gcf
        if save_figs
            saveas(fig,['figs/ACALOS_slope' date],'epsc')
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions

% plot audiogram function
function fig = plot_audiogram_groups(idx,age_sub,aud,aud_freq,cmap)
    figure('renderer','painter')
for s=1:length(idx)

            p1 = semilogx(aud_freq,aud(idx(s),:)','-','color',[cmap(age_sub(idx(s))-17,:) 0.8]);
            plot_aud_param(p1,aud_freq);
            hold on
            
            cb=colorbar;
            
            cb.FontSize = 12;
            cb.Limits = [0 1];
            age_ticks = [20:10:70];
            cb.Ticks = (age_ticks-17)/(max(age_sub)-min(age_sub));%[linspace(0,1,6)];
            cb.TickLabels = age_ticks;%{linspace(18,max(age_sub),5)};
            cb.Label.String = 'Age';
            cb.Label.Rotation = 90;
            cb.Label.FontSize = 16;
            cb.Label.FontName = 'Arial';
            colormap(cmap)
end
set(gcf,'position',[305 412 432 299]);
fig = gcf
fig=gcf;
end


