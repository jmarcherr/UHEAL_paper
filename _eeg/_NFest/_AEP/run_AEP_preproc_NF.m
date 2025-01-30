% run ABR preproc
clear all
cd(['/work3/jonmarc/UHEAL_paper'])
UHEAL_startup
cd(rootdir)
ft_defaults
addpath(datadir)
% save dir
savedir = '/work3/jonmarc/UHEAL_paper/_eeg/_NFest/_AEP/_outputs';
savedir_preproc = '/work3/jonmarc/UHEAL_paper/_eeg/_NFest/_AEP/_outputs/preproc_aep';

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
    data_preproc = AEP_preproc_NF(datadir,subids{ss});
    % run analysis
    data_aep = AEP_analysis_NF(data_preproc);
    % save data
    saveAEPdata(data_aep,savedir,subids{ss});
    % ...
end

%% %%%%%%



