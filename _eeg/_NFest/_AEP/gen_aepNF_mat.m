function [uheal_data]=gen_aepNF_mat(uheal_data,datadir)

d = dir([datadir filesep '*.mat'])
clc
disp(['Processing aep data ...'])

for s=1:length(d)
    load([d(s).folder filesep d(s).name])
    clc
    
    % get FFR
    if isfield(data,'data_w')
        fs = data.fs;
        TS_sub(s,:) = nanmean(data.data_w(:,:));

        TS_sub_chan(s,:,:)=data.data_w;


        time = data.time;
      
        subinfo{s} = data.subinfo;
        if isempty(data.subinfo.age)|isempty(data.subinfo.gender)
            age(s) = nan;
            gender(s) = nan;
        else
        age(s) =data.subinfo.age;
        gender(s) = data.subinfo.gender;
        end

        nr_reject(s) =data.nr_reject;
        sub_id{s} = data.subid;
        sub_num(s) = str2num(sub_id{s}(end-2:end));
        
    else

        TS_sub(s,:,:) =nan(1,1536);
        TS_sub_chan(s,:,:) = nan(16,1536);
        subinfo{s} = data.subinfo;
        age(s) = data.subinfo.age;
        gender(s) = data.subinfo.gender;
        sub_id{s} = data.subid;
        sub_num(s) = str2num(sub_id{s}(end-2:end));
    end
end

% channels
chansoi  =1:16;% setdiff(1:16,[5 11]); % all channels but T7 and T8   
fs = data.fs;


%% get FFT
tidx = time>=0 & time<3;
pow_sub_chan = nan(117,16,769);coeff = nan(117,16,2);
pow_sub = nan(117,769);coeff = nan(117,2);
SNR = nan(117,4);NCL = nan(117,4);
for ss = 1:length(d)
    for cc=1:length(chansoi) % channels
        M=squeeze(TS_sub_chan(ss,cc,find(tidx)));

        %FFT
        f_fft = fft(M)/(length(M)/2);

        %Convert to power
        pow = abs(f_fft.^2); %
        %Truncate negative freqencies
        ft_sub = (pow(1:end/2+1));
        pow_sub_chan(ss,cc,:) = squeeze(ft_sub);

        %Frequency vector
        f = fs/2*linspace(0,1,length(ft_sub));
        % estimate noise floor
        foi = 2;
        nbins = [5];%2:2:20];
        aband = [7 12];
        fitF = f(find(f>0.7 & f<20));
        fitF = setdiff(fitF,[nbins aband]);
        feedback = logical(0);
        [~,coeff_chan(ss,cc,:)] = NNfloorEstim(squeeze(ft_sub)',f,foi,fitF,feedback);
        
    end
    % mean over channels
    M=squeeze(nanmean(TS_sub_chan(ss,chansoi,find(tidx)),2));
    %FFT
    f_fft = fft(M)/(length(M)/2);

    %Convert to power
    pow = abs(f_fft.^2); %
    %Truncate negative freqencies
    ft_sub = (pow(1:end/2+1));
    pow_sub(ss,:) = squeeze(ft_sub);
    foi = 2:2:8;
    [~,coeff(ss,:)] = NNfloorEstim(squeeze(ft_sub)',f,foi,fitF,feedback);

end

 %% save latency and amplitudes of n100
 %p50
uheal_data.AEP_NF_slope = nan(size(uheal_data.subid,1),1);
uheal_data.AEP_NF_int = nan(size(uheal_data.subid,1),1);


%% save
for s=1:length(coeff)
    % get this subid
    thisID = str2double(sub_id{s}(3:5))
    this_idx = find(uheal_data.subid==thisID);
    uheal_data.AEP_NF_slope(this_idx,:) = coeff(s,1);
    uheal_data.AEP_NF_int(this_idx,:) = coeff(s,2);

end
disp(['AEP NF data done!'])
end