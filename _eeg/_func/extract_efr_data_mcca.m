
    % get FFR
    if isfield(data,'FFR')
        fs(s) = data.fs;
        FFR(s,:) = data.FFR;
        FFR_SNR(s,:) = data.FFR_SNR;
        F(s,:) = data.F;
        F_crit(s,:) = data.F_crit;
        sig_idx(s,:) = F(s,:)>F_crit(s,:);
        sig_idx_cz(s) = sig_idx(s,10);
        subid{s} = data.subid;
        f_fft(s,:,:) = data.f_fft;
        fft_freq(s,:) = data.fft_freq;
        TS(s,:,:) = data.FFR_TS;
        time = data.time;
        tidx = data.tidx;
        fft_noise(s,:) = data.fft_noise;
        noisef(s,:,:) = data.noise.fft_sub_noise;
        subinfo{s} = data.subinfo;

        % 3 channel average
        FFR_avg(s) = data.FFR_avg;
        SNR_avg(s) = data.FFR_SNR_avg;
        F_avg(s) = data.F_avg;
        F_crit_avg(s) = data.F_crit_avg;
        sig_idx_avg(s) = data.sig_idx_avg;
        noise_avg(s) = data.noise_avg;

        if isempty(data.subinfo.age)|isempty(data.subinfo.gender)
            age(s) = nan;
            gender(s) = nan;
        else
        age(s) =data.subinfo.age;
        gender(s) = data.subinfo.gender;
        end

        nr_reject(s) =data.nr_reject;
        chan_labels{s} = data.chan_labels;
        chans{s} = data.channels;
        
    elseif ~isfield(data.subinfo,'age')
        subinfo{s} = data.subinfo;
        age(s) = nan;
    else
        FFR(s,:) = nan(1,16);
        FFR_SNR(s,:) = nan(1,16);
        F(s,:)  = nan(1,16);
        F_crit(s,:) = nan(1,16);
        sig_idx(s,:) = zeros(1,16);
        sig_idx_cz(s) = 0;
        f_fft(s,:,:) = nan(16,4097);
        fft_freq(s,:) = nan(1,4097);
        TS(s,:,:) = nan(16,2458);

        fft_noise(s,:) = nan(1,16);
        noisef(s,:,:) = nan(16,4097);
        subinfo{s} = data.subinfo;

        % 3 channel average
        FFR_avg(s) = nan;
        SNR_avg(s) = nan;
        F_avg(s) = nan;
        F_crit_avg(s) = nan;
        sig_idx_avg(s) = 0;
        noise_avg(s) = nan;
        subinfo{s} = data.subinfo;
        age(s) = nan;
    end
    
        CP(s) = uheal_data.CP_new(find(uheal_data.subid==sub_num(s)));