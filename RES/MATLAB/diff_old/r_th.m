function res = r_th(T, mode)
    switch(mode)
        case 1
            res = ones(1, length(T));
        case 2
            res = sqrt(1 + 2 ./ (3 .* T));
        case 3
            res = ((1 + sqrt(1 + 3 .* T / 2)) / 2).^(-1/6);
    end
end

