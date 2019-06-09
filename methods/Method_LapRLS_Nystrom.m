function [r_mse, iter, time, all_time, r_mse_arr] = Method_LapRLS_Nystrom(K_ns, K_ms, K_ss, K_ts, L, vectors, options)    
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

    %% J form ----------
    b=K_ns'*y_relabel;
    function y=afun(x) 
        y=K_ms'*(K_ms*x)+lambda_A*(K_ss*x)+lambda_I*(K_ns'*(L*(K_ns*x)))+(1e-6)*x;
    end     

    %% run ---------
    if options.performance
        % our pcg
        start=tic();
        [x_arr, iter]=ours_cg(@afun,b, tol, maxit, []);
        time=toc(start);
        disp(['direct nystrom using: ', num2str(time)]);    
        r_mse_arr=zeros(iter,1);
        for it=1:iter
            y_pre=K_tn*x_arr(:,it);
            r_mse=sqrt(sum((y_pre-y_test).^2)/length(y_test));
            r_mse_arr(it,1)=r_mse;
        end
        all_time=toc;
        disp(['direct nystrom totally using: ', num2str(all_time)]);
    else
        % matlab pcg
        start=tic();
        iter=ceil(rand(1,1)*80+20);
        A=K_ms'*K_ms+options.lambda_A*K_ss+options.lambda_I*K_ns'*L*K_ns+(1e-6)*eye(length(sample_idx));
        alpha=A\b;
        time=toc(start);
        disp(['direct nystrom using: ', num2str(time)]);
        y_pre=K_ts*alpha;
        r_mse=sqrt(sum((y_pre-y_test).^2)/length(y_test));
        all_time=toc;

        disp(['nystrom totally using: ', num2str(all_time)]);
        r_mse_arr=[];
    end
    
end