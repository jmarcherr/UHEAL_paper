function [uheal_data]=gen_aep_mat(uheal_data,datadir)

d = dir([datadir filesep '*.mat'])
clc
disp(['Processing aep data ...'])

for s=1:length(d)
    load([d(s).folder filesep d(s).name])
    extract_aep_data
    clc
    disp([sub_id{s} ' done...'])
end


com_sub = p2_sub-n1_sub; %obs

% for stats
for ii=1:length(com_sub)
    %pfit(ii,:) = polyfit([0:3],com_sub(ii,:)',1);
    %ffit(ii,:) = polyval(pfit(ii,:),[0:3]);
    pfit(ii,:) = polyfit([0.8:0.5:2.3],com_sub(ii,:)',1);
    ffit(ii,:) = polyval(pfit(ii,:),[0.8:0.5:2.3]);
end


 %% save latency and amplitudes of n100
 %p50
uheal_data.AEP_p1 = nan(size(uheal_data.subid,1),4);
 %n100
uheal_data.AEP_n1 = nan(size(uheal_data.subid,1),4);
 %p200
uheal_data.AEP_p2 = nan(size(uheal_data.subid,1),4);

uheal_data.AEP_p2n1_comp = nan(size(uheal_data.subid,1),4);
uheal_data.AEP_p2n1_int = nan(size(uheal_data.subid,1),1);
uheal_data.AEP_p2n1_slope = nan(size(uheal_data.subid,1),1);


%% save
for s=1:length(n1_mean_sub)
    % get this subid
    thisID = str2double(sub_id{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);
    uheal_data.AEP_p1(this_idx,:) = p1_sub(s,:)
    uheal_data.AEP_n1(this_idx,:) = n1_sub(s,:);
    uheal_data.AEP_p2(this_idx,:) = p2_sub(s,:);
    uheal_data.AEP_p2n1_comp(this_idx,:) = com_sub(s,:);
    uheal_data.AEP_p2n1_int(this_idx,:) = pfit(s,2);
    uheal_data.AEP_p2n1_slope(this_idx) = pfit(s,1);
end
disp(['AEP data done!'])
end