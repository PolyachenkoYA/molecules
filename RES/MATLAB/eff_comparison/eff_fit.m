function e_th = eff_fit(data, M, dim, max_iter, N_fit_ind, p_fit_ind, lmd)
    %tht_init = getInitGuess(data, p_a, N_a);
    data_fit = data(N_fit_ind, p_fit_ind) / M;
    p_grid_fit = p_grid(N_fit_ind, p_fit_ind);
    N_grid_fit = N_grid(N_fit_ind, p_fit_ind);
    Npts_fit = numel(data_fit);
    pts_fit = [reshape(data_fit, [Npts_fit, 1]),...
               reshape(p_grid_fit, [Npts_fit, 1]),...
               reshape(N_grid_fit, [Npts_fit, 1])];

    if(dim == 5)
        a = [0.1, 0.3, 1, 3];
        b = [0.1, 0.3, 1, 3];
        c = [0.1, 0.3, 1, 3];
        d = [0.1, 0.3, 1, 3];
        e = [0.1, 0.3, 1, 3];        
    else
        %a = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
        %b = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
        %c = [0.01, 0.03, 0.1, 0.3, 1, 3, 10, 30];
        a = [0.5, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5];
        b = [0.5, 0.7, 0.8, 0.9, 1, 1.1, 1.2, 1.3, 1.5];
        c = [1, 1.5, 1.7, 1.9, 2, 2.1, 2.3, 2.5, 3];        
    end
    Na = length(a);
    Nb = length(b);
    Nc = length(c);    
    if(dim == 5)
        Nd = length(d);
        Ne = length(e);
        Ntot_fit = Na * Nb * Nc * Nd * Ne;
    else
        Ntot_fit = Na * Nb * Nc;
    end
    tht_init_grid = zeros(Ntot_fit, dim);
    if(dim == 5)
        for ai = 1:Na 
            for bi = 1:Nb
                for ci = 1:Nc
                    for di = 1:Nd
                        for ei = 1:Ne
                            ind = ei + (di-1 + (ci-1 + (bi-1 + (ai-1)*Nb)*Nc)*Nd)*Ne;
                            tht_init_grid(ind, 1) = a(ai);
                            tht_init_grid(ind, 2) = b(bi);
                            tht_init_grid(ind, 3) = c(ci);
                            tht_init_grid(ind, 4) = d(di);
                            tht_init_grid(ind, 5) = e(ei);
                        end
                    end
                end
            end
        end
    else
        for ai = 1:Na 
            for bi = 1:Nb
                for ci = 1:Nc
                    ind = ci + (bi-1 + (ai-1)*Nb)*Nc;
                    tht_init_grid(ind, 1) = a(ai);
                    tht_init_grid(ind, 2) = b(bi);
                    tht_init_grid(ind, 3) = c(ci);
                end
            end
        end        
    end
    tht = grid_srch(tht_init_grid, pts_fit, max_iter, lmd);
    % tht = [1.374, 2.3083, 0.1502]';
    % tht = [30.5054, 0, 4.6454]';
    if(dim == 5)
        e_th = (1 ./ ((tht(1) + tht(2)./N_grid)./p_grid + (tht(3) + tht(4)./N_grid + tht(5).*p_grid)./N_grid)) * M;
    else
        e_th = (1 ./ (tht(1) ./ p_grid + tht(2) ./ N_grid + tht(3) .* (p_grid ./ N_grid))) * M;
    end            
end

