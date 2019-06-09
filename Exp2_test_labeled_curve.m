function Exp2_test_labeled_curve(data_name)
    addpath('./methods');
    addpath('./utils');
    addpath('./utils/');
    rand('state', 0)
    totally=tic();
    %% Parameters ----------
    options = parameters(data_name)
    %% Load Data ----------
    dataset_path=['./datasets/', data_name];
    [y, X]=libsvmread(dataset_path);

    train_fraction=0.7;

    %% Data split ---------- 
    all_idx=randperm(length(y), length(y));
    test_idx=all_idx(1:ceil((1-train_fraction)*length(all_idx)));
    X_test=X(test_idx,:);
    y_test=y(test_idx);

    train_idx=setdiff(all_idx,test_idx);
    X_train=X(train_idx,:);
    y_train=y(train_idx);

    clear X;
    clear all_idx;
    clear train_idx;
    clear test_idx;

    %% Generating grapha Laplacian L and predefinition of kernel----------
    L=construct_laplacian_graph(data_name, X_train, options.NN);
    kernel=load_kernel(data_name, length(y), options.KernelParam);

    times=3;
    res_RLS=[];
    res_LapRLS_pcg=[];
    res_nystrom_pcg=[];
    labeled_partitions=2.^(0:6)./100;
    for labeled_fraction=labeled_partitions   
        record=zeros(5,times);
        disp(['labeled_fraction:', num2str(labeled_fraction)]);
        for t=1:times
            %% Preprocess ----------
            train_size=length(y_train);
            labeled_unsorted=randperm(train_size, ceil(labeled_fraction*train_size));
            unlabeled_idx=setdiff(1:train_size, labeled_unsorted);
            labeled_idx=setdiff(1:train_size, unlabeled_idx);
            sample_idx=1:ceil(0.1*train_size);

            vectors = struct('labeled_idx',labeled_idx, ...
                         'unlabeled_idx',unlabeled_idx, ...
                         'sample_idx',sample_idx,...
                         'y_train',y_train,...
                         'y_test',y_test); 
            if train_size < 30000
                result=zeros(3,4);            
                K_nn=kernel(X_train, 1:length(y_train), X_train, 1:length(y_train), t, 'knn');
                K_tn=kernel(X_test, 1:length(y_test), X_train, 1:length(y_train), t, 'ktn');
                K_mm=K_nn(labeled_idx, labeled_idx);
                K_tm=K_tn(:,labeled_idx);
                [result(1,1), result(1,2), result(1,3), result(1,4)] = Method_RLS(K_mm, K_tm, vectors, options);
                clear K_mm K_tm

                K_ns=K_nn(:,sample_idx);
                K_ss=K_nn(sample_idx,sample_idx)+1e-6*eye(length(sample_idx));
                K_ms=K_ns(labeled_idx,:);
                [result(2,1), result(2,2), result(2,3), result(2,4)] = Method_LapRLS_PCG(K_nn, K_tn,K_ns, K_ss, L, vectors, options);

                
                K_ts = K_tn(:, sample_idx);              
                clear K_tn
                [result(3,1), result(3,2), result(3,3), result(3,4)] = Method_LapRLS_Nystrom_PCG(K_ns, K_ms, K_ss, K_ts, L, vectors, options);

                record(1,t)=result(1,1);
                record(2,t)=result(2,1);
                record(3,t)=result(3,1); 
                clear K_nn K_tn K_mm K_tm K_ns K_ss
            else
                result=zeros(3,4);
                K_mm=kernel(X_train, labeled_idx, X_train, labeled_idx, t, ['kmm_',num2str(labeled_fraction)]);   
                K_tm=kernel(X_test, 1:length(y_test), X_train, labeled_idx, t, ['ktm_',num2str(labeled_fraction)]);  
                [result(1,1), result(1,2), result(1,3), result(1,4)] = Method_RLS(K_mm, K_tm, vectors, options);
                clear K_mm K_tm;

                K_ns=kernel(X_train, 1:length(y_train), X_train, sample_idx, t, ['kns_',num2str(labeled_fraction)]); 
                K_ms=K_ns(labeled_idx,:);
                K_ss=K_ns(sample_idx,:)+1e-6*eye(length(sample_idx));
                K_ts=kernel(X_test, 1:length(y_test), X_train, sample_idx, t, ['kts_',num2str(labeled_fraction)]);
                [result(3,1), result(5,2), result(5,3), result(5,4)] = Method_LapRLS_Nystrom_PCG(K_ns, K_ms, K_ss, K_ts, L, vectors, options);
                
                record(1,t)=result(1,1);
                record(2,t)=result(2,1);
                record(3,t)=result(3,1); 
                clear K_ns K_ms K_ss K_ts
            end
        end
        mean_std=[mean(record,2),std(record,0,2)];
        res_RLS=[res_RLS;mean_std(1,:)];
        res_LapRLS_pcg=[res_LapRLS_pcg;mean_std(2,:)];
        res_nystrom_pcg=[res_nystrom_pcg;mean_std(3,:)];
    end
    str=['./result/labeled_curve_',data_name];
    save(str, 'res_RLS', 'res_LapRLS_pcg', 'res_nystrom_pcg');
    disp(['totally cost ', num2str(toc(totally))]);
    
    draw_labeled_curve(data_name);
end

