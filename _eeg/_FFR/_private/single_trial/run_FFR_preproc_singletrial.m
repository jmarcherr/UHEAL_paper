% run FFR preproc
clear all
cd('/work3/jonmarc/UHEAL_paper')
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
cd(rootdir)
ft_defaults
addpath(datadir)
% save dir
savedir = '/work3/jonmarc/UHEAL_paper/_eeg/_FFR/_outputs/_derivatives/_trials';

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
    data_preproc = FFR_preproc(datadir,subids{ss});
    % run analysis
    data_ffr = FFR_analysis_trials(data_preproc);
    % save data
    saveFFRdata(data_ffr,savedir,subids{ss});
    % ...
end

% Get mcca for subjects
%get_mcca...
% save mixing matrix

%run preproc again for each subject
% apply filter to single trials
% correlate


%% %%%%%%



