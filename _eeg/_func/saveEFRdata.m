%% functions
function saveEFRdata(data,savedir,subid)


    %%  Save mat
    savefile = [savedir filesep subid '_EFR.mat'];
    save(savefile,'data','-v7.3');

end
