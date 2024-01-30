% run EFR preproc
clear all
cd('/work1/jonmarc/UHEAL_master/UHEAL_paper')
UHEAL_startup
cd(rootdir)
ft_defaults
addpath(datadir)
% save dir
savedir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_EFR/_outputs/_derivatives';

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
    data_preproc = EFR_preproc(datadir,subids{ss});
    % run analysis
    data_ffr = EFR_analysis(data_preproc);
    % save data
    saveEFRdata(data_ffr,savedir,subids{ss});
    % ...
end

%% %%%%%%



