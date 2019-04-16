function d2 = get_mean_disp(x, params)
    x0 = zeros(params.Ntot, 3);
    for i1 = 1:params.Ntot
        for i2 = 1:3
            x0(i1, i2) = mean(x(:, i1, i2));
        end
    end
    d2 = zeros(1,params.framesN);
    for i = 1:params.framesN
        d2(i) = sum(sum((squeeze(x(i,:,:)) - x0).^2));
    end
    d2 = sqrt(sum(d2) / (params.Ntot * params.framesN));
end

