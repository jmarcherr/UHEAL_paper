%% functions
function saveABRdata(data,savedir,d,dd)


    %%  Save mat
    savefile = [savedir filesep d(dd).name '_ABR.mat'];
    save(savefile,'data','-v7.3');

end