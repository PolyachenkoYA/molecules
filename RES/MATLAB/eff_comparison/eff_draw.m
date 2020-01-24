function eff_draw(data, e_th, filepath, filename, p_arr, N_arr, cut_ind, draw_time, save_pics, draw_th, draw_full, p_normalize)        
    N_p = length(p_arr);
    N_N = length(N_arr);
    if(isempty(cut_ind))
        cut_ind = 1:N_p;
    end
    [p_grid, N_grid] = meshgrid(p_arr, N_arr);
    if(p_normalize)
        data = data ./ p_grid;
        e_th = e_th ./ p_grid;    
        eff_name = '\eta_p';
        eff_th_name = '\eta^{th}_p(N, p)';
        eff_full_name = '\eta_p = (T/dt \cdot N^2 / t_{wallclock}) / N_p';
    else
        eff_name = '\eta';
        eff_th_name = '\eta_p(N, p)';
        eff_full_name = '\eta = T/dt \cdot N^2 / t_{wallclock}';
    end

    if(draw_full)
        getFig('$N_p$', 'N', ['$' filename '; \hspace{10pt} ' eff_full_name '$'], 'log', 'log', 'log');
        if(draw_th)
            plot3(p_grid, N_grid, data, 'x',...
                  'HandleVisibility', 'off', 'Color', getMyColor(1), 'LineWidth', 2);    
            surf(p_grid, N_grid, e_th,...
                'DisplayName', ['$' eff_th_name '$'], 'EdgeColor', 'interp', 'FaceColor', 'interp');
        else
            surf(p_grid, N_grid, data, 'DisplayName', ['$' eff_name '$'],...
                 'EdgeColor', 'interp', 'FaceColor', 'interp');
        end
        legend off;
        if(save_pics)
            savefig(fullfile(filepath, [filename '_eff_full.fig']));
        end
    end
    
    getFig('$N_p$', 'N', ['$' filename '; \hspace{10pt} ' eff_full_name '$'], 'log', 'log', 'log');
    if(draw_th)
        plot3(p_grid(:, cut_ind), N_grid(:, cut_ind), data(:, cut_ind), 'x',...
              'HandleVisibility', 'off', 'Color', getMyColor(1), 'LineWidth', 2);    
        surf(p_grid(:, cut_ind), N_grid(:, cut_ind), e_th(:, cut_ind),...
            'DisplayName', ['$' eff_th_name '$'], 'EdgeColor', 'interp', 'FaceColor', 'interp');    
    else
        surf(p_grid(:, cut_ind), N_grid(:, cut_ind), data(:, cut_ind), 'DisplayName', ['$' eff_name '$'],...
            'EdgeColor', 'interp', 'FaceColor', 'interp');
    end
    legend off;
    if(save_pics)
        savefig(fullfile(filepath, [filename '_eff_main.fig']));    
    end

    if(draw_time)
        getFig('$N_p$', 'N', ['$' filename '; \hspace{10pt} t(N, N_p)$'], 'log', 'log', 'log');
        T = max(5120 * (2048 ./ N_grid)^2, 4) / 512;
        surf(p_grid, N_grid, (N_grid.^2 ./ data) * T*512, 'DisplayName', '$\t_{wallclock}(N, p)$',...
            'EdgeColor', 'interp', 'FaceColor', 'interp');
        if(save_pics)
            savefig(fullfile(filepath, [filename '_time.fig']));
        end
    end
end

