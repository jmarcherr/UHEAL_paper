%% functions
function saveAEPdata(data,savedir,subid)


    %%  Save mat
    savefile = [savedir filesep subid '_FFR.mat'];
    save(savefile,'data','-v7.3');

end