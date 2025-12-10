close all
clear all
%cd(fileparts(matlab.desktop.editor.getActiveFilename))
run('/work3/jonmarc/UHEAL_paper/UHEAL_startup.m')
load('/work3/jonmarc/UHEAL_paper/_stats/uheal_data.mat')
savepath = ['/work3/jonmarc/UHEAL_paper/_paper_figs/fig6/figs/']
%%

% preprocessing
%nha
nh_idx = find(uheal_data.CP_new==0);
Field_list = fieldnames(uheal_data);
for ii=1:numel(Field_list)
    if (size(uheal_data.(Field_list{ii}),2)==1)% && size(uheal_data.(Field_list{ii}),2)==117)
        uheal_tmp.(Field_list{ii}) = uheal_data.(Field_list{ii})
    end
end
uheal_nh = rmfield(uheal_tmp,'memr_fcenter');
uheal_nh.AEP_p2n1_comp = uheal_data.AEP_p2n1_comp(:,4);
uheal_nh.AEP_n1 = uheal_data.AEP_n1(:,4);
uheal_nh.AEP_p2 = uheal_data.AEP_p2(:,1);
%%
%uheal_nh.AP_amp_pm = uheal_data.AP_amp_pm
%uheal_nh.FFR_SNR = uheal_data.FFR_SNR
%uheal_nh.FFR_sig = uheal_data.FFR_sig
%uheal_nh.AEP_p2n1_int = uheal_data.AEP_p2n1_int;
%uheal_tmp
    %%

uheal_nh = IndexedStructCopy(uheal_nh,nh_idx);
% groups
ynh = find(uheal_nh.Age<=25);
mnh = find(uheal_nh.Age>25 & uheal_nh.Age<50);
onh = find(uheal_nh.Age>=50);

% FFR SNR
sig_idx = find(uheal_nh.FFR_sig==0);
uheal_nh.FFR_SNR(sig_idx) = nan;


%% review 1: find rejected participants per group

% ABR
y_rjt_abr = length(find(isnan(uheal_nh.AP_amp) & uheal_nh.Age<=25))
m_rjt_abr = length(find(isnan(uheal_nh.AP_amp) & uheal_nh.Age>25 & uheal_nh.Age<50))
o_rjt_abr = length(find(isnan(uheal_nh.AP_amp) & uheal_nh.Age>=50))
% FFR
y_rjt_ffr = length(find(~uheal_nh.FFR_sig & uheal_nh.Age<=25))
m_rjt_ffr = length(find(~uheal_nh.FFR_sig & uheal_nh.Age>25 & uheal_nh.Age<50))
o_rjt_ffr = length(find(~uheal_nh.FFR_sig & uheal_nh.Age>=50))
% EFR
y_rjt_efr = length(find(~uheal_nh.EFR_sig & uheal_nh.Age<=25))
m_rjt_efr = length(find(~uheal_nh.EFR_sig & uheal_nh.Age>25 & uheal_nh.Age<50))
o_rjt_efr = length(find(~uheal_nh.EFR_sig & uheal_nh.Age>=50))

% ITPC / 4Hz
uheal_nh.Age(find(isnan(uheal_nh.Neg_4Hz)))
% AEP
uheal_nh.Age(find(isnan(uheal_nh.AEP_p2n1_comp)))

%% variance of N1 peak


% functions
function T = IndexedStructCopy(S, Condition, FieldList)
if nargin == 2
   FieldList = fieldnames(S);
end 
for iField = 1:numel(FieldList)
   Field    = FieldList{iField};
   T.(Field) = S.(Field)(Condition);
end
end