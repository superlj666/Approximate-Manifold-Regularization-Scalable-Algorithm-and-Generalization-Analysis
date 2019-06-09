function [r_mse, iter, time, all_time, r_mse_arr] = Method_LapRLS_PCG(K_nn, K_tn,K_ns, K_ss, L, vectors, options)    
    tic();
    labeled_idx=vectors.labeled_idx;
    unlabeled_idx=vectors.unlabeled_idx;
    sample_idx=vectors.sample_idx;
    y_train=vectors.y_train;
    y_relabel=y_train;
    y_relabel(unlabeled_idx)=0;
    y_test=vectors.y_test;
    
    lambda_A=options.lambda_A;
    lambda_I=options.lambda_I;
    tol=options.tol;
    maxit=options.T; 

    diag_v=zeros(length(y_train),1);
    diag_v(labeled_idx,1)=1;
    J=spdiags(diag_v,0,length(y_train),length(y_train));
    J_lambdaI_L=J+lambda_I*L;
    pinv_part=lambda_A*K_ss+K_ns'*(J_lambdaI_L*K_ns)+(1e-6)*eye(length(sample_idx));
    
    %% direct way    
    function y=afun(x) 
        y=J_lambdaI_L*(K_nn*x)+lambda_A*x;
    end
    function y=mfun(v) 
        y=(v-J_lambdaI_L*(K_ns*(pinv_part\(K_ns'*v))))/lambda_A;
    end
    
    if options.performance
        % our pcg
        start=tic();
        [x_arr, iter]=ours_cg(@afun,y_relabel, tol, maxit, []);
        time=toc(start);
        disp(['direct pcg using: ', num2str(time)]);    
        r_mse_arr=zeros(iter,1);
        for it=1:iter
            y_pre=K_tn*x_arr(:,it);
            r_mse=sqrt(sum((y_pre-y_test).^2)/length(y_test));
            r_mse_arr(it,1)=r_mse;
        end
        all_time=toc;
        disp(['direct pcg totally using: ', num2str(all_time)]);
    else
        % matlab pcg
        start=tic();
        [alpha, flag, rev, iter]=cgs(@afun,y_relabel,tol,maxit, @mfun);
        time=toc(start);
        disp(['direct pcg using: ', num2str(time)]);
        y_pre=K_tn*alpha;
        r_mse=sqrt(sum((y_pre-y_test).^2)/length(y_test));
        all_time=toc;
        disp(['direct pcg totally using: ', num2str(all_time)]);
        r_mse_arr=[];
    end
end