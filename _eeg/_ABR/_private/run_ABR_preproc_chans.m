% run ABR preproc
clear all
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
cd(rootdir)
ft_defaults
addpath(datadir)
% save dir
savedir = '/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_private/_outputs/_derivatives';

% get avalible subjects
d = dir([datadir filesep 'UH*']);
subids = {d.name};
numsub = length(subids);

%% HPC cluster parameters
clust=parcluster('dcc');    % load the MDCS cluster profile
clust.AdditionalProperties.MemUsage = '8GB';
clust.AdditionalProperties.WallTime = '20:00';
clust.saveProfile;
parpool(clust, 20);
% run preproc and save
%%
parfor ss=1:numsub
    % run preprocessing
    data_preproc = ABR_preproc_chans(datadir,subids{ss});
    % run analysis
    data_abr = ABR_analysis_chans(data_preproc);
    % save data
    saveABRdata(data_abr,savedir,subids{ss});
    % ...
end

%% %%%%%%


