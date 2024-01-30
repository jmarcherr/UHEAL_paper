
function data = FFR_4Hz_preproc(dataset_root,subid)
%% Preprocessing of tonal sequence (3s) 2Hz presentation rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Inputs:
%       dataset_root = directory with subject preprocessed data
%       subid = UHEAL data identifyer
%    Outputs:
%        Preprocessed FFR_4Hz data

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

    if isempty(filename_bdf) || strcmp(subid,'UH001') % wrong data from UH001
        error(sprintf('FFR_4Hz data is not available for subject indentifier %s',subid))
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
    cfg.trialdef.prestim      = 5; %-100 ms
    cfg.trialdef.poststim     = 7; % 500 ms
    cfg = ft_definetrial(cfg);

    % Redefine trigger values
    for tt=1:length(triggers)
        cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
        trl_org = cfg.trl;
    end
    if size(cfg.trl,1)==972
        missing_trials = 0;
    else
        warning('Subject has missing trials')
        missing_trials = 1;
    end

    % find 6 trial sequence
    if sum(strcmp(subid,{'UH003','UH031'}))
        fs=2048;
    else
        fs=4096;
    end
    for ii=1:length(trl_org)-1
        this_stamp = trl_org(ii,1); % this zero point
        next_stamp(ii) = trl_org(ii+1,1)-this_stamp;
    end
    break_trials = find(next_stamp/fs>18000/fs);
    break_trials = [1 break_trials+1];
    trl_tmp = cfg.trl(break_trials,:);
    cfg.trl = trl_tmp;

    % initial preprocessing
    data_int = ft_preprocessing(cfg);
    fs = data_int.fsample;
    for ii=1:length(trl_tmp)-1
        these_trials = break_trials(ii):break_trials(ii+1)-1;
        this_break = break_trials(ii+1)-1;
        time_stamps{ii}= (trl_org([these_trials],1)-trl_org(these_trials(1),1))/fs;
    end
    these_trials = break_trials(ii+1):length(trl_org);
    this_break = length(trl_org);
    time_stamps{ii+1}=(trl_org([these_trials],1)-trl_org(these_trials(1),1))/fs;

    % Resampling
    cfgres = [];
    cfgres.resamplefs = 4096/4;
    cfgres.detrend    = 'no';
    data = ft_resampledata(cfgres,data_int);

    %%
    %Rereferencing (l/r mastoid)
    cfg.channel     = {'eeg','-Status'};%chaoi;;%chaoi;
    cfg.reref       = 'yes';
    cfg.refchannel  = {'P9','P10'};

    cfg.layout      =  'biosemi64.lay';
    cfg.continuous  = 'yes';
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50 100];
    cfg.lpfilttype  = 'firws';
    cfg.lpfilter    = 'yes';
    cfg.lpfreq      = 100;

    cfg.hpfilter    = 'yes';
    cfg.hpfilttype  = 'firws';
    cfg.hpfreq      = .7;%80; % changed from .1

    % rereferenced data struct
    data = ft_preprocessing(cfg,data);
    data.missing_trials = missing_trials;
    data.trl_org = trl_org;
    data.time_stamps = time_stamps;
    data.break_trials = break_trials;
    data.subid = subid;
    data.subinfo = dataalm.subinfo;

catch ME
    %No data for current subject

    warning(ME.message)
    data.error = ['error in FFR_4Hz_preproc: ' ME.message];
    data.subid = subid;
    data.subinfo = dataalm.subinfo;
end
end


