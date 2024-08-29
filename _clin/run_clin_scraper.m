% run clin_scraper
clear all
%cd('/work1/jonmarc/UHEAL_master/UHEAL_paper')
run('/zhome/7e/f/64621/Desktop/UHEAL_paper/UHEAL_startup') % path to your folder
cd(rootdir)
ft_defaults
addpath(datadir)
% save dir
savedir = '/zhome/7e/f/64621/Desktop/UHEAL_paper/_clin/_clindata';

% get avalible subjects
d = dir([datadir filesep 'UH*']);
subids = {d.name};
numsub = length(subids);
clc
fprintf('Found %d subject/session folders.\n',numsub)

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
    [dataalm] = UHEAL_clin_scraper(subids{ss},datadir);

    % save data
    save_clindata(dataalm,savedir,subids{ss});
    % ...
end
