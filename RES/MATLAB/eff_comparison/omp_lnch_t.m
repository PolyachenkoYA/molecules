function t = omp_lnch_t(p)
    pl = round(log(p)/log(2)) + 1;
    N = 1e5;
    work_times = [0.0365, 0.105, 0.135, 5.15, 7.55];
    t = zeros(size(pl));
    for i1 = 1:size(pl,1)
        for i2 = 1:size(pl,2)
            t(i1, i2) = work_times(pl(i1, i2)) / N * 1e9;
        end
    end    
end