function [data_w] = weighted_ABR_chans_mcca(data,data_cc,tidx)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

        fs = data.fsample;
        time = data.time{1};
        var_tmp =[];data_weighted=[];epoch_var=[];
        % if size(data_cc,3)<3000
        %     var_tmp = nan;
        %     data_weighted=nan;
        %     epoch_var=nan;
        % else
            for cc=1:size(data_cc,1)
                for ii=1:size(data_cc,3) %chan x time x trials
                    var_tmp(cc,ii) =var(data_cc(cc,tidx,ii));
                    data_weighted(cc,:,ii) = data_cc(cc,tidx,ii)./var_tmp(cc,ii);
                    epoch_var(cc,ii) = var_tmp(cc,ii).^(-1);
                    
                end
            end
        %end

        % Correct channels
        %corrdata =-squeeze(data_weighted(cc,:,:)); %(inverting channel)


        %epoch_var = epoch_var(:);
        summed_trials =sum(data_weighted,3);
        summed_weights = sum(epoch_var,2);
        data_w = summed_trials./summed_weights;
        %data_trials = data_weighted./epoch_var';
end