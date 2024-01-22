
function data = AEP_preproc(dataset_root,subid)
%% Preprocessing of noise burst AEP for ISI = 0.8:0.5:2.3s conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs:
%       dataset_root = directory with subject preprocessed data
%       subid = UHEAL data identifyer
%    Outputs:
%        Preprocessed AEP data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1 || isempty(dataset_root); dataset_root = '/work1/jonmarc/UHEAL_master/UHEAL/UHEAL_data/'; end
if nargin < 2 || isempty(subid); error('Please provide subject identifier'); end

try
    % load subject info for this participant
    filename_subjectdir = dir(fullfile([dataset_root filesep 'scraped' filesep subid '.mat']));
    if isempty(filename_subjectdir)
        error(sprintf('Subject data is not available for subject indentifier %s.',subid))
    end

    load([dataset_root filesep 'scraped' filesep subid '.mat']);

    % load stim file
    stimfile_name = dir(fullfile(dataset_root,subid,'AEP_stim_*'));
    if ~isempty(stimfile_name)
        load(fullfile(stimfile_name.folder, stimfile_name.name));
    else
        % could not get stimfile
        warning(sprintf('stimfile could not be identified for subject %s',subid));
        stim = nan;
    end
    % get BDF name of subject dd
    filename_bdf = dir(fullfile(dataset_root,subid,'*.bdf'));

    if sum(strcmp(subid,{'UH001'})) % wrong AEP data for subjects 1
        error(sprintf('AEP data is not available for subject indentifier %s',subid))
    end

    filename_bdf = fullfile(filename_bdf.folder , filename_bdf.name);

    % Triggers indicate stimulus onset
    triggers = [100,102]; % 100 or 102 depending on lab
    hdr = ft_read_header(filename_bdf);

    cfg=[];
    cfg.layout =  'biosemi64.lay';
    cfg.continuous = 'yes';
    cfg.channel     = 'all';
    cfg.dataset = filename_bdf;
    cfg.trialdef.eventtype    = 'STATUS';
    cfg.trialdef.eventvalue   = triggers;
    cfg.trialdef.prestim      = .1; % -100 ms
    cfg.trialdef.poststim     = .5; % 500 ms
    cfg = ft_definetrial(cfg);

    % Redefine trigger values
    for tt=1:length(triggers)
        cfg.trl(cfg.trl(:,4)==triggers(tt),4) = 1; % all the same trigger
    end
    % get missing trials
    tmp = (cfg.trl(2:end,:)-cfg.trl(1:end-1,:));
    missing_trials = find(tmp(:,1)>=3.8e4);

    %initial preprocessing
    data_int = ft_preprocessing(cfg);

    % Resampling before filtering
    cfgres = [];
    cfgres.resamplefs = 4096;
    cfgres.detrend    = 'no';
    data = ft_resampledata(cfgres,data_int);

    % Rereferencing (P9 and P10)
    cfg.channel     = {'eeg','EXG1','EXG2' '-Status'};
    cfg.reref       = 'yes';
    cfg.refchannel = {'P9','P10'}; 
    cfg.layout      =  'biosemi64.lay';
    cfg.continuous  = 'yes';
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50 100  150];
    cfg.lpfilttype  = 'firws';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 1000;
    cfg.hpfilter    = 'yes';
    cfg.hpfilttype  = 'firws';
    cfg.hpfreq      = .7; 
    cfg.demean      = 'yes';

    % rereferenced data struct
    data = ft_preprocessing(cfg);

    % add stimear
    data.missing_trials = missing_trials;
    data.stimfile  = stim;
    data.subid = subid;
    data.subinfo = dataalm.subinfo;

catch ME
    %No data for current subject

    warning(ME.message)
    data.error = ['error in AEP_preproc: ' ME.message];
    data.subid = subid;
    data.subinfo = dataalm.subinfo;
end
end




