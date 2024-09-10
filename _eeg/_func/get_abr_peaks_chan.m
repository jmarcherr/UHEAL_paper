function [abr_peaks] = get_abr_peaks_chan(data_in,t_abr)
%find ABR peaks automatic

ap_peak = find(t_abr>0e-3 & t_abr<2.3e-3);
wv_peak = find(t_abr>3.5e-3 & t_abr<9e-3); % 4-7 ms before
%sp_peak = find(t_abr>0e-3 & t_abr<.7e-3);
for cc=1:size(data_in,1)
if ~isempty(data_in)

        % AP
        [abr_peaks.AP_amp(cc),AP_lat] = max(squeeze(data_in(cc,ap_peak)));
        this_ap = AP_lat-1+ap_peak(1);
        abr_peaks.AP_latency(cc) = t_abr(this_ap);
        [abr_peaks.AP_neg(cc),AP_neg_lat] = min(squeeze(data_in(this_ap:this_ap+15)));
        abr_peaks.AP_neg_latency(cc) = t_abr(AP_neg_lat-1+this_ap);

        % use AP to define SP
        % SP
        [abr_peaks.SP_amp(cc)] = mean(squeeze(data_in(cc,find(t_abr>0 & t_abr<0.8e-3))));
        %abr_peaks.SP_latency = t_abr(SP_lat+sp_peak(1));

        % WV
        [abr_peaks.WV_amp(cc),WV_lat] = max(squeeze(data_in(cc,wv_peak)));
        this_wv = WV_lat-1 + wv_peak(1);
        abr_peaks.WV_latency(cc) = t_abr(this_wv);
        [abr_peaks.WV_neg(cc),WV_neg_lat] = min(squeeze(data_in(this_wv:this_wv+20)));
        abr_peaks.WV_neg_latency(cc) = t_abr(WV_neg_lat-1+this_wv);
    else
        abr_peaks.SP_amp(cc) =nan;
        abr_peaks.AP_amp(cc) = nan;
        abr_peaks.AP_latency(cc) = nan;
        abr_peaks.AP_neg(cc) = nan;
        abr_peaks.AP_neg_latency(cc) = nan;
        abr_peaks.WV_amp(cc) = nan;
        abr_peaks.WV_latency(cc) = nan;
        abr_peaks.WV_neg(cc) = nan;
        abr_peaks.WV_neg_latency(cc) = nan;
end

end
