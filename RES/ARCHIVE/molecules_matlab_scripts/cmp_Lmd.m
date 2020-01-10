function cmp_Lmd(T, n, lmd, lmd_th, cut_err, tit, ticks, log_scale)
    N = length(lmd_th);
    d = zeros(1,N);
    %for i = 1:N
        %d(i) = max(abs(lmd(i) / lmd_th(i) - 1), abs(lmd_th(i) / lmd(i) - 1));
    %end
    d = rel_err(lmd, lmd_th);
    d(d < eps) = eps;
    lg_d = log(d)/log(10);
    
    if(log_scale)
        y = d;
    else
        y = lg_d;
    end
    plot_lmd(T, n, y, cut_err, tit, ticks, log_scale);
end

