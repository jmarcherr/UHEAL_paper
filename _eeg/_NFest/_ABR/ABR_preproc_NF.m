
function data = ABR_preproc_NF(dataset_root,subid)
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

    if ~isempty(filename_bdf) && sum(strcmp(subid,{'UH003','UH031','UH067','UH091'})) %no ABR for subjects 3,31,67 and 91
        error(sprintf('ABR data is not available for subject indentifier %s',subid))
    end

    filename_bdf = fullfile(filename_bdf.folder , filename_bdf.name);

    % Triggers indicate stimulus onset
    triggers = [50,60,62]; % 100 or 102 depending on lab
    hdr = ft_read_header(filename_bdf);

    cfg=[];
    cfg.layout =  'biosemi64.lay';
    cfg.continuous = 'yes';
    cfg.channel     = 'all';
    cfg.dataset = filename_bdf;
    cfg.trialdef.eventtype    = 'STATUS';
    cfg.trialdef.eventvalue   = triggers;
    cfg.trialdef.prestim      = 0; % -100 ms
    cfg.trialdef.poststim     = 3; % 500 ms
    cfg = ft_definetrial(cfg);

    % Redefine trigger values
    for tt=1:length(triggers)
        cfg.trl(cfg.trl(:,4)==triggers(tt),4) = 1; % all the same trigger
    end
    % create new cfg with 3s blocks (disregarding triggers)
    sample_step = cfg.trl(1,2)-cfg.trl(1,1);
    first_trigger = cfg.trl(1,1);
    last_trigger = cfg.trl(end,1);
    % number of trials
    ntrig = floor((last_trigger-first_trigger)/sample_step);
    cfg_tmp = cfg.trl(1,:);
    for ii=1:ntrig-1
        this_line = cfg_tmp;
        this_line(1,1) = first_trigger+(ii-1)*sample_step+(ii-1)*1;
        this_line(1,2) = first_trigger+ii*sample_step+(ii-1)*1;
        cfg_new(ii,:) = this_line;
    end

    cfg.trl = cfg_new;

    %initial preprocessing
    data_int = ft_preprocessing(cfg);

    % Resampling before filtering
    cfgres = [];
    cfgres.resamplefs = 512;
    cfgres.detrend    = 'no';
    data = ft_resampledata(cfgres,data_int);

    % Rereferencing (P9 and P10)
    cfg.channel     = {'eeg','EXG1','EXG2' '-Status'};
    cfg.reref       = 'yes';
    cfg.refchannel = {'P9','P10'}; 
    cfg.layout      =  'biosemi64.lay';
    cfg.continuous  = 'yes';
    %cfg.dftfilter   = 'yes';
    %cfg.dftfreq     = [50 100  150];
    %cfg.lpfilttype  = 'firws';
    %cfg.lpfilter    = 'yes';
    %cfg.lpfreq      = ;
    cfg.hpfilter    = 'yes';
    cfg.hpfilttype  = 'firws';
    cfg.hpfreq      = .7; 
    %cfg.demean      = 'yes';

    % rereferenced data struct
    data = ft_preprocessing(cfg,data);

    % add stimear
    data.subid = subid;
    data.subinfo = dataalm.subinfo;

catch ME
    %No data for current subject

    warning(ME.message)
    data.error = ['error in ABR_preproc_NF: ' ME.message];
    data.subid = subid;
    data.subinfo = dataalm.subinfo;
end
end




