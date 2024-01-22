function [data_w,data_trials] = weighted_ABR(data,data_cc,tidx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

        fs = data.fsample;
        time = data.time{1};
        var_tmp =[];data_weighted=[];epoch_var=[];
        if size(data_cc,3)<3000
            var_tmp = nan;
            data_weighted=nan;
            epoch_var=nan;
        else
            for cc=1:size(data_cc,1)
                for ii=1:size(data_cc,3) %chan x time x trials
                    var_tmp(cc,ii) =var(data_cc(cc,tidx,ii));
                    data_weighted(cc,:,ii) = data_cc(cc,tidx,ii)./var_tmp(cc,ii);
                    epoch_var(cc,ii) = var_tmp(cc,ii).^(-1);
                    
                end
            end
        end

        % Correct channels
        corrdata =-squeeze(data_weighted(cc,:,:)); %(inverting channel)


        epoch_var = epoch_var(:);
        summed_trials =sum(corrdata,2);
        summed_weights = sum(epoch_var);
        data_w = summed_trials./summed_weights;
        data_trials = corrdata./epoch_var';
end