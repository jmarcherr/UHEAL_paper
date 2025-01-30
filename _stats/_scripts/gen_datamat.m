%% generate UHEAL_data mat
clear all
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
%load clin data
load('/work3/jonmarc/UHEAL_paper/_clin/clin_data_table/clin_data.mat')
% run ABR extract
abr_dir = '/work3/jonmarc/UHEAL_paper/_eeg/_ABR/_outputs/_derivatives'
uheal_data = gen_ABR_mat(uheal_data,abr_dir);
% run FFR extract
ffr_dir = '/work3/jonmarc/UHEAL_paper/_eeg/_FFR/_outputs/_derivatives'
uheal_data = gen_FFR_mat(uheal_data,ffr_dir);
% run EFR extract
efr_dir = '/work3/jonmarc/UHEAL_paper/_eeg/_EFR/_outputs/_derivatives'
uheal_data = gen_EFR_mat(uheal_data,efr_dir);
% run FFR 4Hz extract
ffr4hz_dir = '/work3/jonmarc/UHEAL_paper/_eeg/_FFR_4Hz/_outputs/_derivatives';
uheal_data = gen_4Hz_mat(uheal_data,ffr4hz_dir);
% run AEP extract
aep_dir = '/work3/jonmarc/UHEAL_paper/_eeg/_AEP/_outputs/_derivatives';
uheal_data = gen_aep_mat(uheal_data,aep_dir);
% run AEP NF extract
aep_NF_dir = '/work3/jonmarc/UHEAL_paper/_eeg/_NFest/_AEP/_outputs/';
uheal_data = gen_aepNF_mat(uheal_data,aep_NF_dir)
% run ABR NF extract
abr_NF_dir = '/work3/jonmarc/UHEAL_paper/_eeg/_NFest/_ABR/_outputs/';
uheal_data = gen_abrNF_mat(uheal_data,abr_NF_dir)
% save uheal_data
savedir = '/work3/jonmarc/UHEAL_paper/_stats'
save([savedir filesep 'uheal_data.mat'],'uheal_data')

% convert to table
uheal_table = struct2table(rmfield(uheal_data,{'memr_reflex','memr_fcenter','memr_levels'}));
writetable(uheal_table,[savedir filesep 'uheal_table.csv'])