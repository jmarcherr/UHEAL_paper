% run ABR preproc

for ii=1:subjects
    % run preprocessing
    data = ABR_preproc(datadir,d,dd);
    % save data
    log = saveABRdata;
end