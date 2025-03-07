
function data = FFR_preproc(dataset_root,subid)
%% Preprocessing of tonal FFR for positive and negative polarity conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs:
%       dataset_root = directory with subject preprocessed data
%       subid = UHEAL data identifyer
%    Outputs:
%        Preprocessed FFR data

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
        error(sprintf('FFR data is not available for subject indentifier %s',subid))
    end

    filename_bdf = fullfile(filename_bdf.folder , filename_bdf.name);

    % Triggers indicate stimulus onset
    triggers = [10,20,22]; %10 = positive polarity, [20,22] = negative polarity.
    hdr = ft_read_header(filename_bdf);

    cfg=[];
    cfg.layout =  'biosemi64.lay';
    cfg.continuous = 'yes';
    cfg.channel     = 'all';
    cfg.dataset = filename_bdf;
    cfg.trialdef.eventtype    = 'STATUS';
    cfg.trialdef.eventvalue   = triggers;
    cfg.trialdef.prestim      = .1; %-100 ms
    cfg.trialdef.poststim     = .5; % 500 ms
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

    % Rereferencing (Cz, Fz, FCz)
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

    % repair bad channels
    if strcmp(subid,'UH020')
        noisy_chan = 'FC4';
    elseif strcmp(subid,'UH073')
        noisy_chan = 'FC3';
    else
        noisy_chan = [];
    end

    if noisy_chan
        %frepair noisy channels
        cfg = [];
        cfg.method        = 'distance';
        cfg.neighbourdist = 1;
        cfg.layout = 'biosemi64.lay';
        cfg.channel  = data.label;
        neighbours = ft_prepare_neighbours(cfg);
        missing = noisy_chan;
        cfg.method         = 'average';
        cfg.layout = 'biosemi64.lay';
        cfg.neighbours = neighbours;
        cfg.badchannel = missing;
        cfg.feedback      = 'yes'
        data_rep = ft_channelrepair(cfg,data)
        data = rmfield(data_rep,'elec');
    end



    % add stimear
    data.stimear = stimear;
    data.subid = subid;
    data.subinfo = dataalm.subinfo;

catch ME
    %No data for current subject

    warning(ME.message)
    data.error = ['error in FFR_preproc: ' ME.message];
    data.subid = subid;
    data.subinfo = dataalm.subinfo;
end
end




