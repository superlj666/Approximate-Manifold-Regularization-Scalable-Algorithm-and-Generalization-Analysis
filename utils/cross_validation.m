function cross_validation(data_name, folds)
    rand('state', 0);
    all_start=tic();
    options=struct('KernelParam',-12:2:-6,...[-6,-4,-2,0,2,4,6],...
        'NN',4,...[2,4],...
        'lambda_A', [0.001,0.01,0.1],...[0.0001,0.001, 0.01],...
        'lambda_I', [0.01,0.1,1]);%[0.00001, 0.0001,0.001]);
    train_fraction=0.7;
    labeled_fraction=0.1;

    cv_records=[];
    bif_records=[];
    for t=1:3    
        %% Load Data ----------
        dataset_path=['/home/lijian/datasets/', data_name];
        [y, X]=libsvmread(dataset_path);
        %% Data split ----------    
        % train test split
        all_idx=randperm(length(y), length(y));
        if length(all_idx) < 286
            test_idx=all_idx(1:length(all_idx)-200);
        else
            test_idx=all_idx(1:ceil((1-train_fraction)*length(all_idx)));
        end
        X_test=X(test_idx,:);
        train_idx=setdiff(all_idx,test_idx);
        train_idx=train_idx(1:length(train_idx)-mod(length(train_idx),200));
        if length(train_idx) > 2000
            train_idx=train_idx(1:2000);
        end
        X_train=X(train_idx,:);
        y_train=y(train_idx)/abs(max(y));
        y_test=y(test_idx)/abs(max(y));

        % label split
        train_size=length(y_train);
        labeled_unsorted=randperm(train_size, ceil(labeled_fraction*train_size));
        unlabeled_idx=setdiff(1:train_size, labeled_unsorted);    
        labeled_idx=setdiff(1:train_size, unlabeled_idx);    
        y_train(unlabeled_idx)=0;

        % fold split
        unlabeled_folds_idx=buffer(unlabeled_idx, folds);
        labeled_folds_idx=buffer(labeled_idx, folds);
        l=length(labeled_idx);
        u=length(unlabeled_idx);

        %% tune parameters with cv
        cv_best_sigma=0;
        cv_best_k=0;
        cv_best_lambda_A=0;
        cv_best_lambda_I=0;
        cv_best_error=2^32;
        cv_best_alpha=[];

        bif_best_sigma=0;
        bif_best_k=0;
        bif_best_lambda_A=0;
        bif_best_lambda_I=0;
        bif_best_error=2^32;
        bif_best_alpha=[];

        cv_counter=0;
        bif_counter=0;
        for sigma=options.KernelParam
            K_train_train=kernel(X_train, X_train, sigma, [data_name,'_train_',t]);
            K_labeled_labeled=K_train_train(labeled_idx, labeled_idx);
            sample_idx=randperm(l+u, ceil(sqrt(l+u)));
            K_XU=K_train_train(:, sample_idx);
            K_UU=K_train_train(sample_idx, sample_idx);
            for k=options.NN
                for l_A=options.lambda_A
                    for l_I=options.lambda_I                        
                        L_all=laplacian_graph(X_train, k+1);
                        % H_LS based on origin form
                        diag_v=zeros(train_size,1);
                        diag_v(labeled_idx,1)=1;
                        J=spdiags(diag_v,0,train_size,train_size);
                        H_LS_bif=(J/l+l_I/(l+u)^2*L_all)*K_train_train+l_A/l*eye(train_size);              
                        alpha_all=(H_LS_bif\y_train)./l;  
                        y_lap_train=K_train_train*alpha_all;

                        % H_LS based on bif                    
                        H_LS_start=tic();           
                        A_inv=sparse(l+u, l+u);
                        H_ll_inv=pinv(K_labeled_labeled./l+l_A/l*eye(l));
                        A_inv(labeled_idx, labeled_idx)=H_ll_inv;
                        I_uu=l/l_A;
                        for unlabel_idx_i=unlabeled_idx
                            A_inv(unlabel_idx_i,unlabel_idx_i)=I_uu;
                        end 
                        X=l_I/(l+u)^2*K_XU;
                        Z=K_XU'*L_all;
                        inv_part=pinv(K_UU+Z*A_inv*X);
    %                     H_LS_bif1=J*K_train_train*J/l+l_A/l*eye(train_size)+l_I/(l+u)^2*(K_train_train*L_all);
    %                     H_LS_bif_inv1=pinv(H_LS_bif1);
                        H_LS_time=toc(H_LS_start);
                        bif_counter=bif_counter+H_LS_time;   


                        cv_fold_errors=zeros(1,folds);
                        bif_fold_errors=zeros(1,folds);
                        for fold=1:folds                
                            % cv
                            cv_start=tic();
                            fold_one_idx=[unlabeled_folds_idx(fold,:),labeled_folds_idx(fold,:)];
                            local_labeled_idx=length(unlabeled_folds_idx(fold,:))+1:length(fold_one_idx);
                            fold_one_unlabeled_idx=unlabeled_folds_idx(fold,:);
                            fold_one_labeled_idx=labeled_folds_idx(fold,:);
                            y_fold_one=y_train(fold_one_labeled_idx);

                            %fold_other_idx=setdiff([unlabeled_idx,labeled_idx],fold_one_unlabeled_idx);
                            fold_other_idx=setdiff([unlabeled_idx,labeled_idx],fold_one_idx);
                            %fold_other_idx=sort([fold_other_idx, fold_one_labeled_idx(1:ceil(2*length(fold_one_labeled_idx)/3))]);
                            y_fold_other=y_train(fold_other_idx);

                            m=length(labeled_folds_idx(fold,:));
                            n=length(unlabeled_folds_idx(fold,:));
                            
                            J1=J(fold_other_idx, fold_other_idx);
                            H_LS=(J1/l+l_I/(l+u)^2*L_all(fold_other_idx, fold_other_idx))*K_train_train(fold_other_idx, fold_other_idx)+...
                                l_A/l*eye(length(fold_other_idx));
                            
%                             H_LS1=H_LS_bif(fold_other_idx, fold_other_idx);
%                             H_LS2=(J1/l+l_I/(l+u)^2*L_other_other)*K_train_train(fold_other_idx, fold_other_idx)+...
%                                 l_A/l*eye(length(fold_other_idx));
                            alpha=H_LS\y_fold_other./l;

                            K_one_other=K_train_train(fold_one_labeled_idx, fold_other_idx);
                            y_pre_cv=K_one_other*alpha;
                            cv_fold_errors(1,fold)=sqrt(sum((y_pre_cv-y_fold_one).^2)/length(y_fold_one));%predict(y_pre_cv, y_fold_one);   
                            cv_time=toc(cv_start);
                            cv_counter=cv_counter+cv_time;

                            % bif
                            bif_start=tic();
                            K_train_one=K_train_train(:,fold_one_idx);
                            y_lap_one=K_train_one'*alpha_all;
                            %y_lap_one=K_train_train(fold_one_idx, fold_other_idx)*alpha_all(fold_other_idx);
                            
                            mu_si=2*(y_lap_train-y_train);
                            mu_si(fold_one_unlabeled_idx)=0;
                            mu_si=mu_si(fold_one_idx);
                            L_one_one=L_all(fold_one_idx, fold_one_idx);
                            %B_LS=H_LS_bif_inv1(fold_one_labeled_idx,:)*(-K_train_one*mu_si./(2*m)-l_A/l*y_lap_train-l_I/(m+n)^2*K_train_one*L_one_one*y_lap_one);
                            v=(-K_train_one*mu_si./(2*m)-l_A/l*y_lap_train-l_I/(m+n)^2*K_train_one*L_one_one*y_lap_one);                        
                            B_LS=A_inv(fold_one_labeled_idx,:)*v-A_inv(fold_one_labeled_idx,:)*(X*(inv_part*(Z*(A_inv*v))));
                            H_mm_inv=pinv(K_train_train(fold_one_labeled_idx, fold_one_labeled_idx)+l_A*eye(m));
                            G_diag=diag(l_A*H_mm_inv,0);
                            G=ones(m, 1)-G_diag;

                            y_pre_bif=y_lap_one(local_labeled_idx)+B_LS./(1-folds)+(B_LS./(0.33*(1-folds)^2))./(1.-G);
                            bif_fold_errors(1,fold)=sqrt(sum((y_pre_bif-y_fold_one).^2)/length(y_fold_one));%predict(y_pre_bif, y_fold_one); 
                            bif_time=toc(bif_start);
                            bif_counter=bif_counter+bif_time;                        
                        end
                        cv_mean_error=mean(cv_fold_errors);                    
                        if cv_best_error > cv_mean_error
                            cv_best_sigma=sigma;
                            cv_best_k=k;
                            cv_best_lambda_A=l_A;
                            cv_best_lambda_I=l_I;
                            cv_best_error=cv_mean_error;
                            cv_best_alpha=alpha_all;
                        end

                        bif_mean_error=mean(bif_fold_errors);                    
                        if bif_best_error > bif_mean_error
                            bif_best_sigma=sigma;
                            bif_best_k=k;
                            bif_best_lambda_A=l_A;
                            bif_best_lambda_I=l_I;
                            bif_best_error=bif_mean_error;
                            bif_best_alpha=alpha_all;
                        end
                    end
                end
            end
        end

%         save(['result/',data_name,'_cv'], 'cv_best_sigma', 'cv_best_k', 'cv_best_lambda_A', 'cv_best_lambda_I', 'cv_best_error', 'cv_counter', 'cv_best_alpha');
%         save(['result/',data_name,'_bif'], 'bif_best_sigma', 'bif_best_k', 'bif_best_lambda_A', 'bif_best_lambda_I', 'bif_best_error', 'bif_counter', 'cv_best_alpha');
        disp(['cv cost ', num2str(cv_counter)]);
        disp(['bif cost ', num2str(bif_counter)]);

        disp(['best sigma: ', num2str(cv_best_sigma),' ', num2str(bif_best_sigma)]);
        disp(['best lambda_A: ', num2str(cv_best_lambda_A), ' ', num2str(bif_best_lambda_A)]);
        disp(['best lambda_I: ', num2str(cv_best_lambda_I), ' ', num2str(bif_best_lambda_I)]);
        K_test_train=kernel(X_test, X_train, cv_best_sigma, [data_name,'_test_',t]);
        cv_pre_error=predict(K_test_train*cv_best_alpha, y_test)
        K_test_train=kernel(X_test, X_train, bif_best_sigma, [data_name,'_test_',t]);
        bif_pre_error=predict(K_test_train*bif_best_alpha, y_test)
        %t_statistic(data_name, 10);
        cv_records=[cv_records;[cv_pre_error,cv_counter]];
        bif_records=[bif_records;[bif_pre_error,bif_counter]];
    end
    
    result=[mean(cv_records),std(cv_records);mean(bif_records),std(bif_records)];
    result=result(:,[1,3,2])
    save(['result/',data_name,'_cv_bif_',num2str(folds)], 'cv_records', 'bif_records', 'result');
    disp(['totally use ', num2str(toc(all_start))]);
end

function error=predict(y_pre, y_true)
    if length(unique(y_true))==2
        smaller_label=min(unique(y_true));
        bigger_label=max(unique(y_true));
        y_pre(y_pre>mean(unique(y_true)))=bigger_label;
        y_pre(y_pre<=mean(unique(y_true)))=smaller_label;
        error=sum(y_pre~=y_true)/length(y_pre);
    else
        error=sqrt(sum((y_pre-y_true).^2)/length(y_true));
    end
end

function kernel=kernel(X_row, X_col, sigma, prefix)
%     str=['result/',prefix,'_cv_',num2str(sigma),'.mat'];
%     if exist(str,'file')
%         load(str);
%     else
        norms_row = sum(X_row'.^2);
        norms_col = sum(X_col'.^2);

        kernel = exp((-norms_row'*ones(1,size(X_col,1)) - ones(size(X_row,1),1)*norms_col + 2*(X_row*X_col')).*2^sigma);
%         save(str, 'kernel');
%     end
end

function L=laplacian_graph(X_train, K)
%     str=['result/',prefix,'_cv_',num2str(K-1),'.mat'];
%     if exist(str,'file')
%         load(str);
%     else
        n_sample=size(X_train,1);
        % knn
%         start=tic();
        knn=knnsearch(X_train,X_train,'K',K+1);

        %% directed graph -------
        W=sparse(n_sample, n_sample);
        for i=1:n_sample
            for j=2:K+1
                col=knn(i,j);
                if ~W(i,col)
                    W(i,col)=exp(-sum((X_train(i,:)-X_train(col,:)).^2)/2);
                    W(col,i)=W(i,col);
                end
            end
        end
        D=sum(W)';
        L=D-W;
%         D=spdiags(sqrt(1./D), 0, n_sample, n_sample);
%         L=speye(n_sample)-D*W*D;
%         time=toc(start);
%         disp(['Generate L using: ', num2str(time)]);
%         save(str, 'L');
%     end
end