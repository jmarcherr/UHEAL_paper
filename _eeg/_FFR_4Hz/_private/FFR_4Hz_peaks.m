%% plot FFR_4Hz results
% plot FFR_4Hz and extract peaks
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/zhome/7e/f/64621/Desktop/UHEAL_paper/UHEAL_startup.m')
subs = dir('/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/_derivatives/*.mat')
load('/work3/jonmarc/UHEAL_master/UHEAL_paper/_stats/uheal_data.mat');
%% get data
for s=1:length(subs)
    
    load([subs(s).folder filesep subs(s).name])
    clc
    disp(['sub ' subs(s).name(1:5) ' loaded...'])
    sub_num(s) = str2num(subs(s).name(3:5));
    chansoi = setdiff(1:16,[5 11]);
    % get FFR
    if isfield(data,'TS')
        itpc(s,:,:) = data.itpc;
        f = data.f;
        TS_sub(s,:) = nanmean(data.TS(chansoi,:));
        TS_sub_chan(s,:,:)=data.TS;

        time = data.time;
        tidx = data.tidx;
        tidx_TS = data.tidx_TS;

        age(s) =data.subinfo.age;
        gender(s) = data.subinfo.gender;

        CP(s) =  uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));
        nr_reject(s) =data.nr_reject;
        chan_labels{s} = data.chan_labels;
        chans{s} = data.channels;
        
    else

        itpc(s,:,:) =nan(16,1537);
        TS_sub(s,:,:) = nan(1,5632);
        TS_sub_chan(s,:,:) = nan(16,5632);
        age(s) = data.subinfo.age;
        gender(s) = data.subinfo.gender;
        CP(s) = uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));
    end

    
end
 
%% 
mean(nr_reject)
std(nr_reject)


%% Time-series (TS) baseline
%%%%%%%%%%%%%%%%%%%%%%%

time_TS = time(tidx_TS); %-1:4s
% baseline norm from -0.1 - 0 s
baseline_idx = find(time_TS>-0.1 & time_TS<=0);

% estimate baseline
baseline(:) = mean(TS_sub(:,baseline_idx),2);
baseline_chan = mean(TS_sub_chan(:,1:16,baseline_idx),3);

% baseline corrected time-series
TS_base = TS_sub-baseline'; % mean over 14 scalp channels (chanoi) subjects x time
TS_base_chan = TS_sub_chan-baseline_chan; % subjects x chan x time

%% get age groups
% groups
YNH_idx = find(age<=25 & ~CP );
MNH_idx = find(age>25 & age<50 & ~CP )
ONH_idx = find(age>=50 & ~CP);
nh_idx = find(~CP); % all normal hearing
ages = [17 77];
% colormap
uheal_colormap;
% channels
chansoi  = setdiff(1:16,[5 11]); % all channels but T7 and T8

% plot mean over all
close all
plot(time_TS,nanmean(TS_base(nh_idx,:)))
xlim([-0.5 3.5])
%% find P1,N1,P2,N2
% P1 = 0.045 -  0.065 s
% N1 = 0.085 -  0.15 s
% P2 = 0.15  -  0.25 s
% N2 = 0.2   -  0.5 s

% onset peak
close all
figure('renderer','painters')
  for ss=1:size(TS_base,1)
      close all
    for ii=1:6 % 6 tones
  
        P1_idx(ii,:) =[0.03+0.5*(ii-1) 0.085+0.5*(ii-1)];
        [P1(ss,ii),P1_l(ss,ii)] = max(TS_base(ss,find(time_TS>=P1_idx(ii,1) & time_TS<=P1_idx(ii,2))));

        N1_idx(ii,:) =[0.065+0.5*(ii-1)  0.15+0.5*(ii-1)];
        [N1(ss,ii),N1_l(ss,ii)] = min(TS_base(ss,find(time_TS>=N1_idx(ii,1) & time_TS<=N1_idx(ii,2))));

        P2_idx(ii,:) =[0.15+0.5*(ii-1) 0.25+0.5*(ii-1)];
        [P2(ss,ii),P2_l(ss,ii)] = max(TS_base(ss,find(time_TS>=P2_idx(ii,1) & time_TS<=P2_idx(ii,2))));

        N2_idx(ii,:) = [0.2+0.5*(ii-1) 0.5+0.5*(ii-1)];
        [N2(ss,ii),N2_l(ss,ii)] = min(TS_base(ss,find(time_TS>=N2_idx(ii,1) & time_TS<= N2_idx(ii,2))));
        
%         % plot individual traces
%         subplot(6,1,ii)
%         this_t = [0 0.5]+(0.5*(ii-1));
%         this_t_idx = find(time_TS>=this_t(1) & time_TS<=this_t(2));
%         plot(time_TS(this_t_idx),TS_base(ss,this_t_idx),'color','k')
%         hold on
%         plot(time_TS(this_t_idx(P1_l(ss,ii)))+0.03,P1(ss,ii),'rx')
%         plot(time_TS(this_t_idx(N1_l(ss,ii)))+0.065,N1(ss,ii),'r<')
%         plot(time_TS(this_t_idx(P2_l(ss,ii)))+0.15,P2(ss,ii),'r>')
%         xlim([this_t])
%         set(gca,'xtick',[this_t(1):0.1:this_t(2)])
%         grid on
%         box off
%         if ii==1
%         title([num2str(subs(ss).name(1:5))])
%         end
%         set(gcf,'position',[440 104 231 600])
%         fig = gcf;
%         saveas(fig,['/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/sub_traces/' subs(ss).name(1:5)],'epsc')
    end
    %pause
  end
  

  for pp=1:3
      subplot(3,1,pp)
  plot(time_TS,TS_base(idx_all{pp},:),'color',[0.5 0.5 0.5 0.5])
  hold on
  plot(time_TS,nanmean(TS_base(idx_all{pp},:)),'color',gcol{pp},'linewidth',2)
  xlim([-0.5 3.5])
  ylim([-6 6])
  end
%%
close all
figure('renderer','painters')
  % group plots
    for ii=1:6
    
    subplot(1,6,ii)
    idx_all = {YNH_idx,MNH_idx,ONH_idx};
    gcol = {y_col,m_col,o_col};
    marks = {'o','sq','^'}

    for pp=1:3
        this_t = [0 0.5]+(0.5*(ii-1));
        this_t_idx = find(time_TS>=this_t(1) & time_TS<=this_t(2));
        plot(time_TS(this_t_idx),nanmean(TS_base(idx_all{pp},this_t_idx),1),'color',gcol{pp})
        hold on
        plot(time_TS(this_t_idx(round(nanmean(P1_l(idx_all{pp},ii))))),nanmean(P1(idx_all{pp},ii)),'x','color',gcol{pp})
        plot(time_TS(this_t_idx(round(nanmean(N1_l(idx_all{pp},ii)))))+0.065,nanmean(N1(idx_all{pp},ii)),'<','color',gcol{pp})
        plot(time_TS(this_t_idx(round(nanmean(P2_l(idx_all{pp},ii)))))+0.15,nanmean(P2(idx_all{pp},ii)),'>','color',gcol{pp})
        xlim([this_t])
        ylim([-4 3])
        set(gca,'xtick',[this_t(1):0.1:this_t(2)])
        box off
    end
    set(gcf,'position',[256 128 1090 576])
    end
            fig = gcf;
        %saveas(fig,'/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/replot_ffr','epsc')
%%
close all
subplot(1,2,1)
plot(nanmean(P1,1),'linewidth',2)
hold on
plot(nanmean(N1,1),'linewidth',2)
plot(nanmean(P2,1),'linewidth',2)
plot(nanmean(N2,1),'linewidth',2)
xlim([0 7])
hleg = legend('P1','N1','P2','N2');
hleg.Box = 'off';
hleg.Position = [0.1709 0.4237 0.1416 0.2393];
xlabel('tone nr.')
ylabel('\muV')
set(gca,'xtick',[1:6])
box off


% plot sequences ontop

%%  plot p2 relative to p1
close all
figure('renderer','painters')
p2p1=P2-P1;
n1p1 = N1-P1;
n2p1 = N2-P1;
subplot 241
box off
plot(1:6,nanmean(p2p1,1),'linewidth',2)
hold on
plot(1:6,nanmean(n1p1,1),'linewidth',2)
plot(1:6,nanmean(n2p1,1),'linewidth',2)
xlabel('tone nr.')
ylabel('\muV')
set(gca,'xtick',[1:6])
xlim([0 7])
titles={'P200-P50','N100-P50','N200-P50'};
hleg =legend(titles)
hleg.Box = 'off'
title('Relative peak amplitude')
box off
ylim([-5 2])
p_all = {p2p1;n1p1;n2p1};
%
for ii=1:3
subplot(2,4,ii+1)
errorbar(1:6,nanmean(p_all{ii}(YNH_idx,:)),nanstd(p_all{ii}(YNH_idx,:))/sqrt(length(YNH_idx)),'color',y_col,'marker','o','markerfacecolor',y_col);%,'o-','color',y_col)
hold on
errorbar(1:6,nanmean(p_all{ii}(MNH_idx,:)),nanstd(p_all{ii}(MNH_idx,:))/sqrt(length(MNH_idx)),'color',y_col,'marker','sq','markerfacecolor',m_col);%,'o-','color',y_col)
errorbar(1:6,nanmean(p_all{ii}(ONH_idx,:)),nanstd(p_all{ii}(ONH_idx,:))/sqrt(length(ONH_idx)),'color',y_col,'marker','^','markerfacecolor',o_col);%,'o-','color',y_col)
%ylim([-2 0.5])
box off
xlabel('tone nr.')
ylabel('\muV')
set(gca,'xtick',[1:6])
xlim([0 7])
title(titles{ii})
if ii==1
hleg=legend({'Young','Mid. aged','Older'})
end
hleg.Box = 'off'
end
subplot(2,4,[5:8])
idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
for ii=1:3
plot(time_TS,nanmean(TS_base(idx_all{ii},:)),'color',gcol{ii})
hold on
end
box off
xlim([-0.5 3.5])
ylabel('\muV')
xlabel('Time (s)')

fs = 1024;
tstim  = 0:1/fs:0.25
ramp = hann(100);win=[ramp(1:end/2-1)' ones(length(tstim)-100,1)' ramp(end/2:end)'];
stim = [sin(2*pi*256*tstim).*win zeros(size(tstim))];
stim_all = [];
for ii=1:6
    stim_all = [stim_all stim]
end
stim_all = [stim_all zeros(size([tstim tstim]))];
stime = 0:1/fs:length(stim_all)/fs-1/fs;
plot(stime,(stim_all*0.2)-3,'color',[0.5 0.5 0.5 0.5])

set(gcf,'position',[214 275 1059 429])

fig = gcf;
%saveas(fig,'/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/P50_rel_peaks','epsc')

%% n100-p200
close all
figure('renderer','painters')
n1p2=P2-N1;
p50 = P1;
idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
marks = {'o','sq','^'}
for pp=1:3
    errorbar(1:6,nanmean(n1p2(idx_all{pp},:)),nanstd(n1p2(idx_all{pp},:))/sqrt(length(idx_all{pp})),'color',gcol{1},'marker',marks{pp},'markerfacecolor',gcol{pp});%,'o-','color',y_col)
    hold on
end
box off
title('P200-N100')
xlabel('tone nr.')
ylabel('\muV')
set(gca,'xtick',[1:6])
xlim([0 7])
hleg=legend({'Young','Mid. aged','Older'});
hleg.Box = 'off'
hleg.Position = [0.5728 0.6010 0.2869 0.1950];
set(gcf,'position',[440 471 289 233])

fig = gcf;
saveas(fig,'/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/P2-N1','epsc')

%% P50
close all
figure('renderer','painters')

idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
marks = {'o','sq','^'}
for pp=1:3
    errorbar(1:6,nanmean(P1(idx_all{pp},:)),nanstd(P1(idx_all{pp},:))/sqrt(length(idx_all{pp})),'color',gcol{1},'marker',marks{pp},'markerfacecolor',gcol{pp});%,'o-','color',y_col)
    hold on
end
box off
title('P50')
xlabel('tone nr.')
ylabel('\muV')
set(gca,'xtick',[1:6])
xlim([0 7])
hleg=legend({'Young','Mid. aged','Older'});
hleg.Box = 'off'
hleg.Position = [0.5432 0.1985 0.3599 0.2017];
set(gcf,'position',[440 471 289 233])

fig = gcf;
saveas(fig,'/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/P50','epsc')
%% corr with age
close all
figure('renderer','painters')
p_all = {p2p1;n1p1;n2p1};
idx_all = {YNH_idx,MNH_idx,ONH_idx};
gcol = {y_col,m_col,o_col};
titles={'P200-P50','N100-P50','N200-P50'};
for tt=1:6 % 6 tones
    
    for  ii=1:3
        subplot(3,6,tt+(ii-1)*6)


        scatter(age(nh_idx),p_all{ii}(nh_idx,tt),'markeredgecolor',gcol{2},'marker','.','sizeData',50)
        ll=lsline
        ll.LineWidth =2;
        ll.Color = 'k';
        hold on
        if ii==1
            tit=title(['Tone nr. ' num2str(tt)]);
            tit.FontName = 'Arial'
        end
        if tt==1
            ylab=ylabel(titles{ii});
            ylab.Rotation = 0;
            ylab.Position = ylab.Position-[25 0 0];
            ylab.FontWeight = 'bold'
            ylab.FontName = 'Arial'
        end


    end
end
set(gcf,'position',[303 163 829 540])
fig = gcf;
saveas(fig,'/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/P50_rel_peaks_scatter','epsc')

%% tone 2- tone 1 response
close all
figure('renderer','painters')
for pp =1:3
    subplot(1,3,pp)
    scatter(age(nh_idx),p_all{pp}(nh_idx,2)-p_all{pp}(nh_idx,1),'markeredgecolor',gcol{2},'marker','.','sizeData',50)
    ll=lsline
    ll.LineWidth =2;
    ll.Color = 'k';
    title(titles{pp})
    if pp==1
        ylabel('Tone nr.2 - Tone nr. 1')
    end
    ylim([-6 6])
    %set(gca,'fontsize',12)
end
set(gcf,'position',[440 552 494 151])
fig = gcf;
saveas(fig,'/work3/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/figs/P50_rel_peaks_t2-t1','epsc')


