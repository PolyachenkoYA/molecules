function eff_draw(path, filename, draw_time)
    filepath = fullfile(path, [filename '.dat']);
    p_arr = 1:6;
    N_arr = 2.^(4:17);
    N_p = length(p_arr);
    N_N = length(N_arr);

    data = load(filepath);
    [p_grid, N_grid] = meshgrid(p_arr, N_arr);

    getFig('$N_p$', 'N', '$\eta = T/dt \cdot N^2 / t_{wallclock}$', 'log', 'log', 'log');
    surf(p_grid, N_grid, data, 'DisplayName', '$\eta(N, p)$',...
        'EdgeColor', 'interp', 'FaceColor', 'interp');
    colorbar;
    set(gca,'ColorScale','log');
    savefig(fullfile('DATA', [filename '_eff_full.fig']));
    
    getFig('$N_p$', 'N', '$\eta = T/dt \cdot N^2 / t_{wallclock}$', 'log', 'log', 'log');
    surf(p_grid(:, 1:4), N_grid(:, 1:4), data(:, 1:4), 'DisplayName', '$\eta(N, p)$',...
        'EdgeColor', 'interp', 'FaceColor', 'interp');
    colorbar;
    set(gca,'ColorScale','log');
    savefig(fullfile('DATA', [filename '_eff_main.fig']));    

    if(draw_time)
        getFig('$N_p$', 'N', '$t(N, N_p)$', 'log', 'log', 'log');
        T = max(5120 * (2048 ./ N_grid)^2, 4) / 512;
        surf(p_grid, N_grid, (N_grid.^2 ./ data) * T*512, 'DisplayName', '$\t_{wallclock}(N, p)$',...
            'EdgeColor', 'interp', 'FaceColor', 'interp');
        colorbar;
        set(gca,'ColorScale','log');
        savefig(fullfile(path, [filename '_time.fig']));
    end
end

