if isfield(data,'aep_avg')
    aep_sub(s,:,:) = data.aep_avg;
    aep_sub_filt(s,:,:)=data.aep_avg_filt;
    n1_sub(s,:) = data.n1;
    n1_mean_sub(s,:)=data.n1_mean;
    n1_sub_lat(s,:) = data.n1_lat;
    p2_mean_sub(s,:) = data.p2_mean;
    p2_sub(s,:) = data.p2;
    p2_sub_lat(s,:) = data.p2_lat;

    time = data.time;
    age(s) = data.subinfo.age;
    gender(s) = data.subinfo.gender;
    sub_id{s} = data.subid;
    sub_num(s) = str2num(sub_id{s}(end-2:end));
    nr_reject(s,:) = data.nr_reject;
else
    aep_sub(s,:,:) = nan(4,308);
    aep_sub_filt(s,:,:) = nan(4,308);
    n1_sub(s,:) = nan(1,4);
    n1_mean_sub(s,:)=nan(1,4);
    n1_sub_lat(s,:) = nan(1,4);
    p2_mean_sub(s,:) = nan(1,4);
    p2_sub(s,:) = nan(1,4);
    p2_sub_lat(s,:) = nan(1,4);

    %time = nan;
    age(s) = data.subinfo.age;
    gender(s) = data.subinfo.gender;
    sub_id{s} = data.subid;
    sub_num(s) = str2num(sub_id{s}(end-2:end));
    %nr_reject(s,:) = nan;

end