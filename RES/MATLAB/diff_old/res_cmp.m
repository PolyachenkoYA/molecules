function res_cmp(T, n, lmd, lmd_th, cut_err, tit, ticks, log_scale)
    N = length(lmd_th);
    d = zeros(1,N);
    d = eps_err(lmd, lmd_th);
    %d = (lmd_th - lmd) ./ max(lmd, lmd_th);
    %d(abs(d) < eps) = eps .*  sign(d);
    lg_d = log(d)/log(10);
    
    if(log_scale) % this is done the right way!
        z = d;
    else
        z = lg_d;
    end
    
    z_2d = reshape(z, [length(T) length(n)]);
    
    figure;
    ax2 = axes;
    res_plot(ax2, T, n, z_2d, log(cut_err)/log(10), tit, ticks, log_scale);
end

