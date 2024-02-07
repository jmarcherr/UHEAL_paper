%% extract audiograms, and subinfo and save as mat.

clear all;close all
cd(fileparts(matlab.desktop.editor.getActiveFilename))

datafol_remote = '/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data';

%addpath(datafol);
addpath(datafol_remote);
cd(datafol_remote);
% dependent scripts
addpath('/work1/jonmarc/UHEAL_master/UHEAL/_scripts');
root = cd;
cont = dir(datafol_remote);
cont = cont(~ismember({cont.name},{'.','..','scraped'}));


sfol = cont([cont.isdir]); %isolate data folders
addpath(sfol.name)

%make the output directory
cd(datafol_remote)
outputfol = 'scraped';
if ~exist(outputfol)
mkdir(outputfol)
end
cd(outputfol)
outputfol = fullfile(cd);

fprintf('Found %d subject/session folders.\n',size(sfol,1))
rootdir = cd;
%% loop
for i = 1:size(sfol,1)%1:size(sfol,1) %loop across subject/session folders




    % extract audiogram and subinfo

    load([sfol(i).name '.mat'])
    aud = dataalm.aud;
    subinfo =dataalm.subinfo;

    % save in datafolder
    cd(datafol_remote)
    cd(sfol(i).name)%enter each folder
    fprintf('Entered folder: %s \n',sfol(i).name)
    save('aud.mat','aud')
    save('subinfo.mat','subinfo')
    cd(rootdir)
end