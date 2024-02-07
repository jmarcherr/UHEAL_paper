%% generate UHEAL_data mat
clear all
run('/work1/jonmarc/UHEAL_master/UHEAL_paper/UHEAL_startup.m')
%load clin data
load('/work1/jonmarc/UHEAL_master/UHEAL_paper/_clin/clin_data_table/clin_data.mat')
% run ABR extract
abr_dir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_ABR/_outputs/_derivatives'
uheal_data = gen_ABR_mat(uheal_data,abr_dir);
% run FFR extract
ffr_dir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR/_outputs/_derivatives'
uheal_data = gen_FFR_mat(uheal_data,ffr_dir);
% run EFR extract
efr_dir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_EFR/_outputs/_derivatives'
uheal_data = gen_EFR_mat(uheal_data,efr_dir);
% run FFR 4Hz extract
ffr4hz_dir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/_derivatives';
uheal_data = gen_4Hz_mat(uheal_data,ffr4hz_dir);
% run AEP extract
aep_dir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_eeg/_AEP/_outputs/_derivatives';
uheal_data = gen_aep_mat(uheal_data,aep_dir);
% save uheal_data
savedir = '/work1/jonmarc/UHEAL_master/UHEAL_paper/_stats'
save([savedir filesep 'uheal_data.mat'],'uheal_data')

% convert to table
uheal_table = struct2table(rmfield(uheal_data,{'memr_reflex','memr_fcenter','memr_levels'}));
writetable(uheal_table,[savedir filesep 'uheal_table.csv'])