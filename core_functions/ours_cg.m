function [x_arr,iter] = ours_cg(afun, b, tol, maxit, mfun)
    x_arr=zeros(length(b),maxit);
    x = zeros(length(b), 1);
    r = b - afun(x);
    if isempty(mfun) % cg
        p = r;
        rsold = r' * r;

        for iter = 1:maxit
            Ap = afun(p);
            alpha = rsold / (p' * Ap);
            x = x + alpha * p;
            x_arr(:,iter)=x;
            r = r - alpha * Ap;
            rsnew = r' * r;
            if sqrt(rsnew) < tol
                disp(['convergent at ', num2str(iter)]);
                break;
            end
            if rsnew==rsold
                disp(['residual does not change at ', num2str(iter)]);
                break;
            end
            p = r + (rsnew / rsold) * p;
            rsold = rsnew;
        end
        if iter==maxit && sqrt(rsnew) >= tol
            disp(['achive maxit iteration and not convergent']);
        end
    else % pcg
        z = mfun(r);
        p = z;
        rsold = r' * z;

        for iter = 1:maxit
            Ap = afun(p);
            alpha = rsold / (p' * Ap);
            x = x + alpha * p;
            r = r - alpha * Ap;            
            x_arr(:,iter)=x;
            if sqrt(r'*r) < tol
                disp(['convergent at ', num2str(iter)]);
                break;
            end
            z=mfun(r);
            rsnew = r' * z;
            if rsnew==rsold
                disp(['residual does not change at ', num2str(iter)]);
                break;
            end
            p = r + (rsnew / rsold) * p;
            rsold = rsnew;
        end
        if iter==maxit && sqrt(rsnew) >= tol
            disp(['achive maxit iteration and not convergent'])
        end        
    end
end