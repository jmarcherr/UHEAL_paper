% startup
rootdir = fileparts(which('UHEAL_startup.m'));
%datadir = ([rootdir, filesep 'UHEAL_data']); % relative datafolder
datadir = '/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/'; % hardcoded datafolder


try %#ok
    rng(1); 
end  


fprintf('\n project directory now added to the current path \n')

addpath([ rootdir '/_tools/fieldtrip-master'])
if ~exist(fileparts(which('ft_defaults.m')))
    fprintf('remember to add fieldtrip to you path! \n')
end



addpath(fullfile('_eeg/'))
%addpath(fullfile('_func/'))
addpath (fullfile('_eeg/_preproc'))
addpath (fullfile('_eeg/_func'))
addpath (fullfile('_eeg/_analysis'))
addpath (fullfile('_tools/cbrewer/cbrewer'))
addpath (fullfile('_tools'))

fprintf('\n directory addded to the path')



