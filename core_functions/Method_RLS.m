function [r_mse, iter, time, all_time, r_mse_arr] = Method_RLS(K_mm, K_tm, vectors, options)    
    tic();
    labeled_idx=vectors.labeled_idx;
    y_train=vectors.y_train;
    y_test=vectors.y_test;

    lambda=options.lambda;
    tol=options.tol;
    maxit=options.T;     
    
    %% direct way
    K_mm_re=K_mm+lambda*speye(length(labeled_idx));
    function y=afun(x) 
        y=K_mm_re*x;
    end
    
    if options.performance
        % our pcg
        start=tic();
        [x_arr, iter]=ours_cg(@afun,y_train(labeled_idx), tol, maxit, []);
        time=toc(start);
        disp(['direct RLS using: ', num2str(time)]);    
        r_mse_arr=zeros(iter,1);
        for it=1:iter
            y_pre=K_tm*x_arr(:,it);
            r_mse=sqrt(sum((y_pre-y_test).^2)/length(y_test));
            r_mse_arr(it,1)=r_mse;
        end
        all_time=toc;
        disp(['direct RLS totally using: ', num2str(all_time)]);
    else
        % matlab pcg
        start=tic();
        [alpha, flag, rev, iter]=pcg(@afun,y_train(labeled_idx),tol,maxit);
        time=toc(start);
        disp(['direct RLS using: ', num2str(time)]);
        y_pre=K_tm*alpha;
        r_mse=sqrt(sum((y_pre-y_test).^2)/length(y_test));
        all_time=toc;
        disp(['direct RLS totally using: ', num2str(all_time)]);
        r_mse_arr=[];
    end
end

