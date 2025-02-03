%% plot AEP_NF results
% plot AEP_NF and extract peaks
clear all
cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
subs = dir('_outputs/*.mat')
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat');
fig_save_path = '/work3/jonmarc/UHEAL_paper/_eeg/_NFest/_ABR/figs/';

%% get data
for s=1:length(subs)
    
    load([subs(s).folder filesep subs(s).name])
    clc
    disp(['sub ' subs(s).name(1:5) ' loaded...'])
    sub_num(s) = str2num(subs(s).name(3:5));
    
    % get FFR
    if isfield(data,'data_w')
        fs = data.fs;
        TS_sub(s,:) = nanmean(data.data_w(:,:));
        if size(data.data_w,1)<14
            TS_sub_chan(s,1:size(data.data_w,1),:) = data.data_w;
            TS_sub_chan(s,size(data.data_w,1)+(14-size(data.data_w,1)),:) = nan(1,1536);
        else
            TS_sub_chan(s,:,:)=data.data_w;
        end

        time = data.time;
      
        subinfo{s} = data.subinfo;
        if isempty(data.subinfo.age)|isempty(data.subinfo.gender)
            age(s) = nan;
            gender(s) = nan;
        else
        age(s) =data.subinfo.age;
        gender(s) = data.subinfo.gender;
        end

        CP(s) =  uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));%data.subinfo.CP;
        nr_reject(s) =data.nr_reject;
        
    else

        TS_sub(s,:,:) =nan(1,1536);
        TS_sub_chan(s,:,:) = nan(16,1536);
        subinfo{s} = data.subinfo;
        age(s) = data.subinfo.age;
        gender(s) = data.subinfo.gender;
    end
end

 %% 
mean(nr_reject)
std(nr_reject)

%% get age groups

% groups
YNH_idx = find(age<=25 & ~CP & ~isnan(TS_sub(:,1))');
MNH_idx = find(age>25 & age<50 & ~CP & ~isnan(TS_sub(:,1))')
ONH_idx = find(age>=50 & ~CP & ~isnan(TS_sub(:,1))');
idx_all = find(~CP & ~isnan(TS_sub(:,1))');
ages = [17 77];
% colormap
uheal_colormap;
% channels
chansoi  =1:16;% setdiff(1:16,[5 11]); % all channels but T7 and T8   
fs = data.fs;


%% get FFT
tidx = time>=0 & time<3;
pow_sub_chan = nan(117,16,769);coeff = nan(117,16,2);
pow_sub = nan(117,769);coeff = nan(117,2);
SNR = nan(117,4);NCL = nan(117,4);
for ss = 1:length(idx_all)
    for cc=1:length(chansoi) % channels
        M=squeeze(TS_sub_chan(idx_all(ss),cc,find(tidx)));

        %FFT
        f_fft = fft(M)/(length(M)/2);

        %Convert to power
        pow = abs(f_fft.^2); %
        %Truncate negative freqencies
        ft_sub = (pow(1:end/2+1));
        pow_sub_chan(idx_all(ss),cc,:) = squeeze(ft_sub);

        %Frequency vector
        f = fs/2*linspace(0,1,length(ft_sub));
        % estimate noise floor
        foi = 2;
        nbins = [];%2:2:20];
        aband = f(find(f>7 &  f<12));
        fitF = f(find(f>0.7 & f<20));
        fitF = setdiff(fitF,[nbins aband]);
        feedback = logical(0);
        [~,coeff_chan(idx_all(ss),cc,:)] = NNfloorEstim(squeeze(ft_sub)',f,foi,fitF,feedback);
        
    end
    % mean over channels
    M=squeeze(nanmean(TS_sub_chan(idx_all(ss),chansoi,find(tidx)),2));
    %FFT
    f_fft = fft(M)/(length(M)/2);

    %Convert to power
    pow = abs(f_fft.^2); %
    %Truncate negative freqencies
    ft_sub = (pow(1:end/2+1));
    pow_sub(idx_all(ss),:) = squeeze(ft_sub);
    foi = 2:2:8;
    [~,coeff(idx_all(ss),:)] = NNfloorEstim(squeeze(ft_sub)',f,foi,fitF,feedback);

end


%% plot mean spectrum
close all
gg = {YNH_idx;MNH_idx;ONH_idx};
cols = {y_col,m_col,o_col};
% over channels
for cc = 1:14
    subplot(2,3,1)
semilogx(f,db(squeeze(nanmean(pow_sub_chan(:,cc,:),1))))
xlim([0.7 40])
xlabel('frequency Hz')
title('log')
hold on
subplot(2,3,2)
plot(f,squeeze(nanmean(pow_sub_chan(:,cc,:))))
title('linear')
xlim([0.7 20])
xlabel('frequency Hz')
hold on
%hleg = legend(chan_labels{1}(chansoi),'FontSize',8);
hleg = legend(data.chan_labels{data.chanoi},'FontSize',8);

end
set(gcf,'renderer','painters','position',[440 187 726 420])
hleg.Position = [0.6576 0.4338 0.0992 0.4917];

% mean

    subplot(2,3,4)
semilogx(f,db(squeeze(nanmean(pow_sub))),'k')
xlim([0.7 40])
xlabel('frequency Hz')
subplot(2,3,5)
plot(f,squeeze(nanmean(pow_sub)),'k')
%title('linear')
xlim([0.7 20])

xlabel('frequency Hz')
fig = gcf;
saveas(fig,[fig_save_path 'spectrum_chan'],'svg')
%
% mean over groups
figure('renderer','painters','Position',[440 187 726 420])

for g = 1:3
    subplot(2,3,1)
    semilogx(f,db(nanmean(pow_sub(gg{g},:))),'color',cols{g})
    hold on
    xlim([0.7 40])
    xlabel('frequency Hz')
    title('log')

    subplot(2,3,2)
    plot(f,nanmean(pow_sub(gg{g},:)),'Color',cols{g})
    hold on
    title('linear')
    hleg = legend('Young','Mid-aged','Older')
    hleg.Position = [0.6368 0.8207 0.1768 0.1048];
    xlim([0.7 20])
    xlabel('frequency Hz')
    


end
fig = gcf;
saveas(fig,[fig_save_path 'spectrum_groups'],'svg')

% scatter plots
figure('renderer','painters')
subplot(2,3,1)
%slope
scatter(age(idx_all),20*coeff(idx_all,1),'k.')
[r,p]=corr(age(idx_all)',20*coeff(idx_all,1),'type','spearman');
text(45,-50,['r = ' num2str(round(r,2)) newline 'p = ' num2str(round(p,3))],'fontsize',8)
ylim([-60 10])
lsline
xlabel('age')
ylabel('slope (dB/Hz)')
subplot(2,3,2)
%intercept
scatter(age(idx_all),20*coeff(idx_all,2),'k.')
[r,p]=corr(age(idx_all)',20*coeff(idx_all,2),'type','spearman');
text(50,-5,['r = ' num2str(round(r,2)) newline 'p = ' num2str(round(p,3))],'fontsize',8)
ylim([-60 10])
lsline
xlabel('age')
ylabel('intercept (dB \muV)')
%subplot(2,3,3)
%estimate
% scatter(age(idx_all),db(est(idx_all,1)),'k.')
% [r,p]=corr(age(idx_all)',db(est(idx_all,1)),'type','spearman');
% text(50,-13.5,['r = ' num2str(round(r,2)) newline 'p = ' num2str(round(p,3))],'fontsize',8)
% lsline
% xlabel('age')
% ylabel('estimate (dB \muV)')

for g = 1:3
        subplot(2,3,4)
    plot(ones(size(coeff(gg{g},1)))*g,20*coeff(gg{g},1),'o','Color',cols{g})
    title('slope')
    box off
    hold on
    xlim([0 4])
    ylim([-60 10])
    set(gca,'xtick',[1:3],'xticklabels',{'Y','Ma','O'})
    subplot(2,3,5)
    plot(ones(size(coeff(gg{g},2)))*g,20*coeff(gg{g},2),'o','Color',cols{g})
    title('intercept')
    hold on
    box off
    xlim([0 4])
    ylim([-60 10])
    set(gca,'xtick',[1:3],'xticklabels',{'Y','Ma','O'})
    % subplot(2,3,6)
    % plot(ones(size(est(gg{g})))*g,db(est(gg{g})),'o','Color',cols{g})
    % box off
    % title('est')
    % hold on
    % xlim([0 4])
    % set(gca,'xtick',[1:3],'xticklabels',{'Y','Ma','O'})
end
fig = gcf;
saveas(fig,[fig_save_path 'coeff_est_corr'],'svg')
%
% plot model for each group (mean)
%close all
figure('renderer','painters')
for g = 1:3
    subplot(2,3,1)

    [est_g,coeff_g] = NNfloorEstim(nanmean(pow_sub(gg{g},:)),f,foi,fitF,feedback);
    modelx = log10(fitF);
    modely = coeff_g(1)*modelx + coeff_g(2);
    semilogx(f,db(nanmean(pow_sub(gg{g},:))),'color',[cols{g} 0.3])
    hold on
    semilogx(fitF,20*modely,'--','color',cols{g})
    xlim([0.7 25])
    xlabel('frequency Hz')
    title('log')
    box off

    subplot(2,3,2)
    plot(f,db(nanmean(pow_sub(gg{g},:))),'Color',[cols{g} 0.3])
    hold on
    plot(fitF,20*modely,'--','color',[cols{g}])
    title('linear')
    box off
    %legend([pp{1} pp{2} pp{3}],'Young','Mid-aged','Older')
    xlim([0.7 25])
    xlabel('frequency Hz')
    set(gcf,'position',[440 187 726 420])



end
fig = gcf;
saveas(fig,[fig_save_path 'NF_est'],'svg')

