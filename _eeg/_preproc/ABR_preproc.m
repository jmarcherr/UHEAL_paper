
function ABR_par(datadir,d,dd)
%% Preprocessing of click ABR for 9/s and 40/s conditions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   outputs data.mat of ABR trials and saves to _preprocABR folder for
%   subsequent analysis
%   Inputs:
%   datadir = directory with subject preprocessed data (UHEAL_scaper)
%   d = structure of UHEAL data
%   dd = currect subject nr.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define root and data directories
rootdir = cd;
cd(datadir)
cd(d(dd).name)

% get BDF name of subject dd
bdf = dir('*.bdf')
if ~isempty(bdf) || strcmp(d(dd).name,'UH091') || strcmp(d(dd).name,'UH067') %no ABR for subjects 91 and 67
    dataset = bdf.name;
    % load stim file
    stim_file = [];
    stim_file = dir('click_abr_stim*');
    if ~isempty(stim_file)
        load(stim_file.name);
        % Get stim ear
        stimear = stim.ear(1);
    end

    % load dataalm for this participant
    cd ..
    cd('scraped')
    load([d(dd).name '.mat'])
    eeg_lab = dataalm.subinfo.lab;

    % back to datadir
    cd(datadir)
    % Subject dir
    cd(d(dd).name)
    %% ------------Event extraction --------------------------------------
    triggers = [50,60,62]; % 50 = 9.1/s, [60,62] = 40/s.


    hdr = ft_read_header(dataset);
    cfg=[];
    cfg.layout =  'biosemi64.lay';
    cfg.continuous = 'yes';
    cfg.channel     = 'all';
    cfg.dataset = dataset;
    cfg.trialdef.eventtype    = 'STATUS';
    cfg.trialdef.eventvalue   = triggers;
    cfg.trialdef.prestim      = 10e-3; %10 ms
    cfg.trialdef.poststim     = 20e-3; %20 ms
    cfg = ft_definetrial(cfg);

    % redefine trigger values
    for tt=1:length(triggers)
        cfg.trl(cfg.trl(:,4)==triggers(tt),4) = tt;
    end

    %% Rereferencing (Cz, Fz, FCz)

    cfg.dataset = dataset;
    cfg.channel     = {'eeg','EXG1','EXG2' '-Status'};%chaoi;;%chaoi;
    cfg.reref       = 'yes';
    cfg.refchannel = {'Cz','Fz','FCz'};
    cfg.layout      =  'biosemi64.lay';
    cfg.continuous  = 'yes';
    % line noise filter
    cfg.dftfilter   = 'yes';
    cfg.dftfreq     = [50 100  150];
    % low-pass filter
    cfg.lpfilttype  = 'but';
    cfg.lpfilter    = 'yes';
    cfg.lpfiltord   = 2;
    if dd==31 || dd==3 % recorded with 2048 fs
        cfg.lpfreq      = 1000;%3000;
    else
        cfg.lpfreq      = 3000;
    end
    % high-pass filter
    cfg.hpfilter    = 'yes';
    cfg.hpfilttype  = 'but';
    cfg.hpfreq      = 100;
    cfg.demean = 'yes';
    cfg.hpfiltord   = 4;

    % rereferenced data struct
    data = ft_preprocessing(cfg);

    % Resample to 16 kHz
    if dd==31 || dd==3 % recorded with 2048 fs
        cfgres = [];
        cfgres.resamplefs = 16384;
        cfgres.detrend    = 'no';
        data = ft_resampledata(cfgres, data);
    end

    cd(datadir)
    cd ..
    % save to folder
    cd(['_EEG' filesep '_preprocdata_ABR'])
    %%  Save mat
    savefile = [d(dd).name '_ABR.mat'];
    save(savefile,'data','-v7.3');
    % back to root
    cd(rootdir)
else
% no data for this subjet
warning(['No data for subject ' d(dd).name])
%back to root
cd(rootdir)
end
end





