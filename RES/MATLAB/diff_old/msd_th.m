function r2 = msd_th(t, tau, T)
    r2 = (t / tau + exp(-t / tau) - 1) * (6*T*tau^2); 
end

