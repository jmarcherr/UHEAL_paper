% run ABR preproc
clear all

cd('/work1/jonmarc/UHEAL_master/UHEAL_paper')
UHEAL_startup
cd(rootdir)
ft_defaults
addpath(datadir)
% save dir
savedir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_preproc/_ABR/_outputs'

% go to data dir
cd(datadir)
d = dir('UH*');
numsub = length(d);


% %% HPC cluster parameters
% clust=parcluster('dcc');    % load the MDCS cluster profile
% clust.AdditionalProperties.MemUsage = '8GB';
% clust.AdditionalProperties.WallTime = '20:00';
% clust.saveProfile;
%%
parpool(clust, 20);
for ii=1:length(d)
    % run preprocessing
    data = ABR_preproc(datadir,d,dd);
    % save data
    log = saveABRdata;
end

function saveABRdata(savedir,d,dd)
    cd(datadir)
    cd ..
    % save to folder
    cd([savedir filesep '_preprocdata_ABR'])
    %%  Save mat
    savefile = [d(dd).name '_ABR.mat'];
    save(savefile,'data','-v7.3');
    % back to root
    cd(rootdir)
end