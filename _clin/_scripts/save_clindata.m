%% functions
function save_clindata(dataalm,savedir,subid)


    %%  Save mat
    savefile = [savedir filesep subid '.mat'];
    save(savefile,'dataalm','-v7.3');

end