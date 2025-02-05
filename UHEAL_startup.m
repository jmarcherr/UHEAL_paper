% startup
rootdir = fileparts(which('UHEAL_startup.m'));
%datadir = ([rootdir, filesep 'UHEAL_data']); % relative datafolder
datadir = '/work3/jonmarc/UHEAL_master/UHEAL/UHEAL_data/'; % hardcoded datafolder

set(0,'DefaultAxesFontName','SansSerif')
try %#ok
    rng(1); 
end  


fprintf('\n project directory now added to the current path \n')

addpath([ rootdir '/_external/fieldtrip-master'])
if ~exist(fileparts(which('ft_defaults.m')))
    fprintf('remember to add fieldtrip to you path! \n')
end


% add project paths
addpath(genpath(fullfile('_eeg/')))
addpath(genpath(fullfile('_clin/')))

addpath (fullfile('_external/cbrewer/cbrewer'))
addpath (fullfile('_external/NoiseTools'))
addpath (fullfile('_external'))
addpath( genpath(fullfile('_external/distributionPlot')))

fprintf('\n directory addded to the path')
% Determine where your m-file's folder is.
%folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
%addpath(genpath(folder));


