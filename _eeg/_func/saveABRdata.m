%% functions
function saveABRdata(data,savedir,subid)


    %%  Save mat
    savefile = [savedir filesep subid '_ABR.mat'];
    save(savefile,'data','-v7.3');

end