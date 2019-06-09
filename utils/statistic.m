function statistic(data_name)
    load(['./result/',data_name,'.mat']);
    res_std=[std(res_RLS); std(res_LapRLS); std(res_LapRLS_pcg); std(res_nystrom); std(res_nystrom_pcg)];
    tmp=[res_RLS(:,1)-res_LapRLS_pcg(:,1), res_LapRLS(:,1)-res_LapRLS_pcg(:,1), res_LapRLS_pcg(:,1)-res_LapRLS_pcg(:,1), res_nystrom(:,1)-res_LapRLS_pcg(:,1), res_nystrom_pcg(:,1)-res_LapRLS_pcg(:,1)];
    mean(tmp)./std(tmp)*sqrt(size(tmp,1))>2.015
    res_mean(:,1)'
    res_std(:,1)'
    disp([num2str(res_mean(:,2)),['|';'|';'|';'|';'|'],num2str(res_mean(:,3))])
end