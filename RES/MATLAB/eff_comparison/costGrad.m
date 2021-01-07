function [J, gJ] = costGrad(tht, lmd, pts)
% tht = [x0, y0, R]
    dim  = length(tht);
    p = pts(:, 2);
    N = pts(:, 3);
    if(dim == 3)
        f = 1 ./ (tht(1)./p + (tht(2) + tht(3) .* p) ./ N);
    else
        f = 1 ./ ((tht(1) + tht(2)./N)./p + (tht(3) + tht(4)./N + tht(5).*p)./N);
    end
    d = f - pts(:, 1);
    f2 = f.^2;

    g1 = d .* f2;
    g11 = g1 ./ p;
    g12 = g1 ./ N;
    J = (d' * d + (tht' * tht + 2*sum(abs(tht))) * lmd)/2;
    if(dim == 3)
        gJ(1) = -sum(g11);
        gJ(2) = -sum(g12);
        gJ(3) = -sum(g12 .* p);
    else
        gJ(1) = -sum(g11);
        gJ(2) = -sum(g12 ./ p);
        gJ(3) = -sum(g12);
        gJ(4) = -sum(g12 ./ N);
        gJ(5) = -sum(g12 .* p);
    end
    gJ = gJ + (tht + sign(tht))' * lmd;
    
    gJ = gJ';         
end


