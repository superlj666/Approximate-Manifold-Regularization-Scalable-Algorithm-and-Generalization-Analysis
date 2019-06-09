function f=kernel_mul_ssl_call(kernel_method)
    f=@(left_idx, right_idx, L, y) actual_f(kernel_method, L, left_idx, right_idx,y);
end
function res=actual_f(kernel, L, left_idx, right_idx, y) 
    n=length(left_idx);
    s=length(right_idx);

    step = 2^30/s;
    ms = ceil(linspace(0, n, ceil(n/step)+1));
    if isempty(y)
        if isempty(L) %K_ms^T * K_ms
            K_sls=zeros(s,s);
            for i=1:ceil(n/step)
                K_ls = kerenl(left_idx(ms(i)+1:ms(i+1)), right_idx);
                K_sls = K_sls + K_ls'*K_ls;
            end
            res=K_sls;
        else
            K_sns = zeros(s, s);
            for i=1:ceil(n/step)
                KL=zeros(s,length(ms(i)+1:ms(i+1)));
                for j=1:ceil(n/step)
                    K_step_j_s = kernel(left_idx(ms(j)+1:ms(j+1)), right_idx);
                    KL=KL+K_step_j_s'*L;
                end
                K_step_i_s=kernel(left_idx(ms(i)+1:ms(i+1)), right_idx);
                K_sns=K_sns+KL*K_step_i_s;
            end
            res=K_sns;
        end
    else
        z = zeros(s, 1);
        ms = ceil(linspace(0, n, ceil(n/step)+1));
        for i=1:ceil(n/step)
            Kr = kernel(left_idx(ms(i)+1:ms(i+1)), right_idx);
            z = z + Kr'*y;
        end
        res=z;
    end
end