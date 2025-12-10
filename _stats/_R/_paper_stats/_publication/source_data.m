%% generate source data from Marcher-Rørsted et al. "Peripheral and central effects of auditory aging"


clear all
datapath = '/work3/jonmarc/UHEAL_paper/_stats/'; % insert data path here
% load data
load([datapath 'uheal_data.mat'])
data=uheal_data;
% extract only nh
nh_idx = find(data.CP_new==0);
Field_list = fieldnames(data);
for ii=1:numel(Field_list)
    %if (size(data.(Field_list{ii}),2)==1)%
        data_tmp.(Field_list{ii}) = data.(Field_list{ii})
    %end
end
rmvars = {'memr_reflex_growth','memr_reflex','memr_fcenter','memr_levels'...
            'audL','audR','CP_new',...
            'abr_IV_ratio','abr_baseline','SP_amp_peak','SP_amp','SP_lat','AP_amp','WV_amp',...
            'FFR_noise','EFR_noise','ssq12','fds','teoae_amp',...
            'acalos_MCL','acalos_L_cut','acalos_m_high','acalos_HTL_A_sub',...
            'tone_p1_1vsrest','tone_n1_1vsrest','tone_p2_1vsrest','tone_n2_1vsrest',...
            'AEP_NF_slope','AEP_NF_int','AEP_p2n1_comp','ABR_NF_slope','ABR_NF_int','FFR_NF_slope','FFR_NF_int'};
for rr = 1:length(rmvars)
data_tmp = rmfield(data_tmp,rmvars{rr});
end
data_nh = data_tmp;
data_nh = IndexedStructCopy(data_nh,nh_idx)

Oldnames ={'CP_SG','gender','aud','audfreq','SP_amp_man','SP_lat_man','Neg_4Hz','ITPC_ratio'};
Newnames = {'Subgroup_idx','Sex','Audiogram','Audiogram_Freq','SP_amp','SP_lat','tone_AEP_mean','ITPC_Q'}
% rename fields
for nn = 1:length(Oldnames)
%data_nh = table2struct(renamevars(struct2table(data_nh), Oldnames{nn}, Newnames{nn}));
[data_nh.(Newnames{nn})] = data_nh.(Oldnames{nn}); 
% This adds a new field, newName, with the same data as in the field oldName
idx = find(strcmp(fieldnames(data_nh),Oldnames{nn}))
data_nh = rmfield(data_nh,Oldnames{nn});
data_nh = orderfields(data_nh,[1:idx-1,length(fieldnames(data_nh)),idx:length(fieldnames(data_nh))-1]); 
% In this example, there are 11 fields and I changed the name of the 4th field. 
%This command reorders the new field into the place of the old field.

end

source_data = data_nh;

% save uheal_data
savedir = '/work3/jonmarc/UHEAL_paper/_stats/_R/_paper_stats/_publication'
save([savedir filesep 'source_data.mat'],'source_data')

% convert to table
source_table = struct2table(source_data);
writetable(source_table,[savedir filesep 'source_data_table.csv'])

%%
% audiogram
plot(data.audfreq(1,:),data.aud(find(data.CP_new==0),:),'k')
xlabel('Frequency (Hz)')
ylabel('Hearing Level (dB)')

% plot scatter plots
varx = 'Age';
vary = {'PTA_lf','PTA_hf','memr_slope','SP_amp','AP_amp_pm','WV_amp_pm',...
    'FFR_SNR','AEP_p2n1_int','Neg_4Hz','ITPC_ratio'};
varlabels = {'PTA_{lf}','PTA_{hf}','MEMR growth','SP amplitude','AP amplitude','WV amplitude',...
    'FFR_{SNR}','P2N1 intercept','AEP_{\mu}','ITPC_{Q}'};
varunits = {'HL (dB)','HL (dB)','MEMR growth','\muV','\muV','\muV',...
    'dB','\muV','\muV','ITPC_{Q}'};

% PTA lf
figure
for ss=1:length(vary)
subplot(4,3,ss)
scatter(getfield(data_nh,varx),getfield(data_nh,vary{ss}),'k.')
lsline
xlabel('Age')
ylabel(varunits{ss})
title(varlabels{ss},'FontWeight','normal')
end

% functions
function T = IndexedStructCopy(S, Condition, FieldList)
if nargin == 2
   FieldList = fieldnames(S);
end 
for iField = 1:numel(FieldList)
   Field    = FieldList{iField};
   T.(Field) = S.(Field)(Condition,:);
end
end
