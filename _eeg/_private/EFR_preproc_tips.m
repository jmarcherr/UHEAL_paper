
function data = EFR_preproc_tips(dataset_root,subid)
%% Preprocessing of 120 Hz AM noise EFR for positive and negative polarity conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs:
%       dataset_root = directory with subject preprocessed data
%       subid = UHEAL data identifyer
%    Outputs:
%        Preprocessed EFR data

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

    % load stim ear
    stimfile_name = dir(fullfile(dataset_root,subid,'ffr_SW_stim*'));
    if ~isempty(stimfile_name)
        load(fullfile(stimfile_name.folder, stimfile_name.name));
        stimear = stim.ear(1);
    else
        % could not get stimear
        warning(sprintf('stimulus ear could not be identified for subject %s',subid));
        stimear = nan;
    end

    % get BDF name of subject dd
    filename_bdf = dir(fullfile(dataset_root,subid,'*.bdf'));

    if isempty(filename_bdf) && sum(strcmp(subid,{'UH003','UH031'})) % no tiptrode channels
        error(sprintf('EFR data is not available for subject indentifier %s',subid))
    end


    filename_bdf = fullfile(filename_bdf.folder , filename_bdf.name);

    % Triggers indicate stimulus onset
    triggers = [100,102]; % noise onset.
    hdr = ft_read_header(filename_bdf);

    cfg=[];
    cfg.layout =  'biosemi64.lay';
    cfg.continuous = 'yes';
    cfg.channel     = 'all';
    cfg.dataset = filename_bdf;
    cfg.trialdef.eventtype    = 'STATUS';
    cfg.trialdef.eventvalue   = triggers;
    cfg.trialdef.prestim      = .1; %-100 ms
    cfg.trialdef.poststim     = .6; % 500 ms
    cfg = ft_definetrial(cfg);

    % Redefine trigger values
    for tt=1:length(triggers)
        cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
    end

    %initial preprocessing
    data_int = ft_preprocessing(cfg);

    % Resampling before filtering
    cfgres = [];
    cfgres.resamplefs = 4096;
    cfgres.detrend    = 'no';
    data_res = ft_resampledata(cfgres,data_int);

    % Rereferencing (ipsilateral tiptrode)
    cfg.channel     = {'eeg','EXG1','EXG2' '-Status'};
    cfg.reref       = 'yes';
    if stimear ==1
        cfg.refchannel  = {'EXG1'}; %vertex electrodes%linked mastoids
    else
        cfg.refchannel  = {'EXG2'};
    end
    cfg.layout      =  'biosemi64.lay';
    cfg.continuous  = 'yes';
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50 100  150];
    cfg.hpfilter    = 'yes';
    cfg.hpfilttype  = 'firws';
    cfg.hpfreq      = 80; % changed from 80

    % re-referenced and filtered data struct
    data = ft_preprocessing(cfg,data_res);

    % subinfo
    data.subid = subid;
    data.subinfo = dataalm.subinfo;
    data.stimear = stimear;

catch ME
    %No data for current subject

    warning(ME.message)
    data.error = ['error in EFR_preproc: ' ME.message];
    data.subid = subid;
    data.subinfo = dataalm.subinfo;
end
end




