 % For small datasets (n < 20000, 3G) -----

function f = load_kernel(data_name, kernel_size, sigma)
    if kernel_size > 1000000
        %todo
        f = @(X, row, col) res(X(row,:), X(col,:), sigma);
    else
        f = @(X_row, row, X_col, col, seq, postfix) full_kernel(data_name, sigma, X_row, row,X_col, col, seq, postfix);
    end
end

function kernel=full_kernel(data_name, sigma, X_row, row,X_col, col, seq, postfix)    
    str=['data/',data_name,'/Gaussian_',num2str(sigma),'_',num2str(seq),'_',postfix];
    if exist(str,'file')
        fid=fopen(str,'r');
        kernel=fread(fid, [length(row),length(col)], 'double');
        fclose(fid);
    else
        if ~exist(['data/',data_name], 'dir')
            mkdir(['data/',data_name], 'dir');
        end
        tic();
        X_row=X_row(row,:);
        X_col=X_col(col,:);
        norms_row = sum(X_row'.^2);
        norms_col = sum(X_col'.^2);

        kernel = exp((-norms_row'*ones(1,size(X_col,1)) - ones(size(X_row,1),1)*norms_col + 2*(X_row*X_col'))/(2*sigma^2));
        
q
        
        disp([data_name, ' generate kernel matrix using: ',num2str(toc)]);
    end
end

function D = res(X1, X2, sigma)
    sq1 = sum(X1.^2,2);
    sq2 = sum(X2.^2,2)';
    D = X1*X2';

    clear X2
    clear X1

    D = -2.0*D;
    D = bsxfun(@plus, D, sq2);
    clear sq2

    D = bsxfun(@plus, D, sq1);
    clear sq1
    D = bsxfun(@times, D, -1/(2*sigma^2));

    if isa(D, 'gpuArray')
        D = arrayfun(@exp, D);
    else
        D = exp(D);
    end
end