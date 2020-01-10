function [D_th, D_int, D_Einst, lin_fit, r2, v_cor] = get_diff_coef(params, time, lin_i, cubo_time_cut, x, v, Nfrm)
    D_th = sqrt(params.Tmp / pi) * 2 / (3 * pi * params.n);
    [r2, v_cor] = comp_diff(Nfrm, x, v);
    D_int = sum(v_cor(time < cubo_time_cut)) * (1/params.dumpDT) / 3;
    lin_fit = polyfit(time(lin_i), r2(lin_i), 1);
    D_Einst = lin_fit(1) / 6;
end

