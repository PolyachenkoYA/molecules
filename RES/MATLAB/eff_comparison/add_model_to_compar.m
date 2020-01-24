function [pi, Ni] = add_model_to_compar(p, N, mdl, ax_pfix, ax_Nfix, clr)
    draw_th = ~any(size(mdl.data) ~= size(mdl.e_th));
    
    if(~isempty(p))
        pi = (mdl.p_arr == p);
        if(any(pi))
            plot(ax_pfix, mdl.N_arr, mdl.data(:, pi), 'o--',...
                'DisplayName', mdl.model_name, 'LineWidth', 1.5, 'Color', clr);
            if(draw_th)
                plot(ax_pfix, mdl.N_arr, mdl.e_th(:, pi), ...
                    'DisplayName', [mdl.model_name ' th'], 'Color', clr);
            end
        end
    end
    
    if(~isempty(N))
        [m, Ni] = min(abs(log(mdl.N_arr/N)));
        if(m > 0.5)
            Ni = 0;
        end
        if(any(Ni))
            plot(ax_Nfix, mdl.p_arr, mdl.data(Ni, :), 'o--',...
                'DisplayName', mdl.model_name, 'LineWidth', 1.5, 'Color', clr);
            if(draw_th)
                plot(ax_Nfix, mdl.p_arr, mdl.e_th(Ni, :), ...
                    'DisplayName', [mdl.model_name ' th'], 'Color', clr);
            end        
        end
    end
end

