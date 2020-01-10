function l = Lmd_lg_t(fit_b, fit_l, T)
    l = (3/4) * sqrt(pi * T / 2) * exp((fit_b(2) - fit_l(2)) / (fit_l(1) - fit_b(1)));
end

