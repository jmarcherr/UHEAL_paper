%% plot data from Marcher-Rørsted et al. "Peripheral and central effects of auditory aging"


clear all
datapath = '/work3/jonmarc/UHEAL_paper/_stats/'; % insert data path here
% load data
load([datapath 'uheal_data.mat'])
data=uheal_data;
% extract only nh
nh_idx = find(data.CP_new==0);
Field_list = fieldnames(data);
for ii=1:numel(Field_list)
    if (size(data.(Field_list{ii}),2)==1)%
        data_tmp.(Field_list{ii}) = data.(Field_list{ii})
    end
end
data_nh = rmfield(data_tmp,'memr_fcenter');
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
   T.(Field) = S.(Field)(Condition);
end
end
