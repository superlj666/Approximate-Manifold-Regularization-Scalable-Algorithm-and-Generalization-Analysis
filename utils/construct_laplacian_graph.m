function [L,W]=construct_laplacian_graph(data_name, X_train, K)
    tic();
    str=['data/',data_name,'/laplacian_',num2str(K),'.mat'];
    if exist(str,'file')
        load(str);
    else
        if ~exist(['data/',data_name], 'dir')
            mkdir(['data/',data_name], 'dir');
        end
        n_sample=size(X_train,1);
        % knn
        start=tic();
        knn=knnsearch(X_train,X_train,'K',K+1);
        time=toc(start);
        disp([num2str(K), '-NN using: ', num2str(time)]);

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
        D=spdiags(sqrt(1./D), 0, n_sample, n_sample);
        L=speye(n_sample)-D*W*D;
        save(str, 'L');
    end

    disp(['generate L using: ',num2str(toc)]);
end