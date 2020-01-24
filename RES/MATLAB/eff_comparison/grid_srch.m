function tht = grid_srch(tht_init_grid, pts_fit, max_iter, lmd)
    N = size(tht_init_grid, 1);
    dim = size(tht_init_grid, 2);
    J_arr = zeros(N, 1);
    tht_ans = zeros(N, dim);
    parfor i = 1:N
        [tht_ans(i, :), J_arr(i)] = fminunc(@(tht)(costGrad(tht, lmd, pts_fit)), tht_init_grid(i, :)',...
                                 optimset('GradObj', 'on', 'MaxIter', max_iter));    
        disp(i/N);
    end
    [~, min_ind] = min(J_arr);
    tht = tht_ans(min_ind, :)';
end

