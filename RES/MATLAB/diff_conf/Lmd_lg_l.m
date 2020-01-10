function l = Lmd_lg_l(fit_b, fit_l)
    l = (1/4) * sqrt(3 * pi / (1 + exp(-2))) * exp((fit_b(2) * fit_l(1) - fit_b(1) * fit_l(2)) / (2 * (fit_l(1) - fit_b(1))));
end

