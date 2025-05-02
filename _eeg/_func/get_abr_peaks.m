function [abr_peaks] = get_abr_peaks(data_in,t_abr)
%find ABR peaks automatic

ap_peak = find(t_abr>0e-3 & t_abr<2.3e-3);
wv_peak = find(t_abr>3.5e-3 & t_abr<9e-3); % 4-7 ms before
baseline = find(t_abr>0e-3 & t_abr<1e-3); % lowest amplitude within the firstmsec (liberman 2016)
sp_peak = find(t_abr>0e-3 & t_abr<.8e-3);

if ~isempty(data_in)

        % AP
        [abr_peaks.AP_amp,AP_lat] = max(squeeze(data_in(ap_peak)));
        this_ap = AP_lat-1+ap_peak(1);
        abr_peaks.AP_latency = t_abr(this_ap);
        [abr_peaks.AP_neg,AP_neg_lat] = min(squeeze(data_in(this_ap:this_ap+15)));
        abr_peaks.AP_neg_latency = t_abr(AP_neg_lat-1+this_ap);

        % use AP to define SP
        % SP
        % mean 0-0.8ms
        [abr_peaks.SP_amp] = mean(squeeze(data_in(find(t_abr>0 & t_abr<0.8e-3))));
        % max 0-0.8ms
        [abr_peaks.SP_amp_peak,SP_lat] = max(squeeze(data_in(sp_peak)));
        abr_peaks.SP_latency = t_abr(SP_lat+sp_peak(1));

        % WV
        [abr_peaks.WV_amp,WV_lat] = max(squeeze(data_in(wv_peak)));
        this_wv = WV_lat-1 + wv_peak(1);
        abr_peaks.WV_latency = t_abr(this_wv);
        [abr_peaks.WV_neg,WV_neg_lat] = min(squeeze(data_in(this_wv:this_wv+20)));
        abr_peaks.WV_neg_latency = t_abr(WV_neg_lat-1+this_wv);

        % Baseline
        [abr_peaks.baseline,~] = min(squeeze(data_in(baseline)));

    else
        abr_peaks.SP_amp =nan;
        abr_peaks.SP_amp_peak = nan;
        abr_peaks.SP_latency = nan;
        abr_peaks.AP_amp = nan;
        abr_peaks.AP_latency = nan;
        abr_peaks.AP_neg = nan;
        abr_peaks.AP_neg_latency = nan;
        abr_peaks.WV_amp = nan;
        abr_peaks.WV_latency = nan;
        abr_peaks.WV_neg = nan;
        abr_peaks.WV_neg_latency = nan;
        abr_peaks.baseline = nan;
end

end
