% run ABR preproc
clear all
cd('/work1/jonmarc/UHEAL_master/UHEAL_paper')
UHEAL_startup
cd(rootdir)
ft_defaults
addpath(datadir)
% save dir
savedir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_ABR/_outputs/_derivatives'

% get avalible subjects
d = dir([datadir filesep 'UH*']);
numsub = length(d);

%% HPC cluster parameters
clust=parcluster('dcc');    % load the MDCS cluster profile
clust.AdditionalProperties.MemUsage = '8GB';
clust.AdditionalProperties.WallTime = '20:00';
clust.saveProfile;
parpool(clust, 20);
% run preproc and save
%%
parfor dd=1:length(d)
    % run preprocessing
    data_preproc = ABR_preproc(datadir,d,dd);
    % run analysis
    data_abr = ABR_analysis(data_preproc);
    % save data
    saveABRdata(data_abr,savedir,d,dd);
    % ...
end

%% %%%%%%



