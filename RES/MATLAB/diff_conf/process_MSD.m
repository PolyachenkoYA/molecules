function [k_exp, Lmd, consist] = process_MSD(base_path, n, T, small_time, big_time, sazerlend_time, print_res, graph_mod, close_figs)
    if(graph_mod)
        line_w_was = get(0, 'DefaultLineLineWidth');
        set(0, 'DefaultLineLineWidth', 2);
    end

    consist.th_line = 1 / (6 * T);
    consist.th_log = sqrt(3*T*(1 + exp(-2)));
    Lmd.th_0 = Lmd_th(r_th(T, 1), n);
    Lmd.th_S = Lmd_th(r_th(T, 2), n);
    Lmd.th_T = Lmd_th(r_th(T, 3), n);
    % Lmd.th_0 = 1 / (pi * n * sqrt(2));
    % Lmd.th_S = Lmd.th_0 / (thS_r(T)^2);
    % Lmd.th_T = Lmd.th_0 / (thT_r(T)^2);

    model_name = fullfile(base_path, ['n' num2str(n) '_T' num2str(T)]);
    if(round(T) == T)
       model_name = [model_name '.0'];
    end
    filename = fullfile(model_name, 'data.out');
    data = dlmread(filename);
    t = data(:, 1);
    msd = data(:, 2);    
    
    t_line = t;
    msd_line = msd;
    for i = 1:5
        line_fit = polyfit(t_line, msd_line, 1);
        line_time_ind = (t > big_time * line_fit(1) * consist.th_line);
        msd_line = msd(line_time_ind);
        t_line = t(line_time_ind);            
    end    
    k_exp = line_fit(1); % slope 1st estimation
    tau = k_exp * consist.th_line;
    ballistic_time_ind = (t < small_time * tau);
    %ballistic_time_ind(1) = 0; % in case t=0 or msd=0
    sazerlend_time_ind = ((t > tau * big_time) .* (t < tau * sazerlend_time)) > 0;
    msd_ballistic = msd(ballistic_time_ind);
    t_ballistic = t(ballistic_time_ind);
    msd_sazerlend = msd(sazerlend_time_ind);
    t_sazerlend = t(sazerlend_time_ind);    
    
    log_end_fit = polyfit(log(t_line), log(msd_line), 1);
    log_beg_fit = polyfit(log(t_ballistic), log(msd_ballistic), 1);    
    r1_fit = polyfit(t_ballistic, sqrt(msd_ballistic), 1);
    sazerlend_fit = polyfit(t_sazerlend, msd_sazerlend, 1);
    
    %a = abs(polyval(line_fit, 0));
    a = abs(polyval(sazerlend_fit, 0));
    consist.exp_line = a / k_exp^2;
    consist.exp_log = exp((log_beg_fit(2)*(log_end_fit(1)-2) - log_end_fit(2)*(log_beg_fit(1)-2)) / (2*(log_end_fit(1) - log_beg_fit(1))));
    consist.d_line = eps_err(consist.exp_line, consist.th_line);
    consist.d_log = eps_err(consist.exp_log, consist.th_log);
    msd_th_my = msd_th(t, tau, T);
    Lmd.a = Lmd_a(a);
    Lmd.k = Lmd_k(k_exp, T);
    Lmd.log_t = Lmd_lg_t(log_beg_fit, log_end_fit, T);
    Lmd.log_l = Lmd_lg_l(log_beg_fit, log_end_fit);
    
    if(print_res)
        disp(['Lmd = ' num2str(Lmd.a) '; ' num2str(Lmd.k) '; ' num2str(Lmd.log_t) '; ' num2str(Lmd.log_l)]);
        disp(['Lmd_0 = ' num2str(Lmd.th_0) '; Lmd_th_S = ' num2str(Lmd.th_S) '; Lmd_th_T = ' num2str(Lmd.th_T)]);
        disp(['consist_line = ' num2str(eps_err(Lmd.a, Lmd.k)) ';  consist_log = ' num2str(eps_err(Lmd.log_t, Lmd.log_l))]);
    end
    if(graph_mod)
        
        fig_r2_name = 'r2';
        fig_r2 = figure;
        ax_r2 = axes;
        plot(ax_r2, t, msd, 'DisplayName', 'MSD');
        hold(ax_r2, 'on');
        plot(ax_r2, t, polyval(line_fit, t), ':', 'DisplayName', 'asymptote');
        plot(ax_r2, t, msd_th_my, '--', 'DisplayName', 'theory');
        plot(ax_r2, [0, max(t)], [0 0], 'Color', [0 0 0], 'DisplayName', 'zero', 'LineWidth', 1);
        plot(ax_r2, [tau tau], [-a, 1.3*a], '-', 'Color', [0 0 0], 'DisplayName', '\tau');
        x_draw = big_time * tau;
        plot(ax_r2, [x_draw x_draw], [-a/2, msd_th(x_draw, tau, T) + a/2], '--', 'Color', [0 0 0], 'DisplayName', 'line time');
        x_draw = small_time * tau;
        plot(ax_r2, [x_draw x_draw], [-a/2, msd_th(x_draw, tau, T) + a/2], '-.', 'Color', [0 0 0], 'DisplayName', 'ballistic time');
        grid(ax_r2, 'on');
        xlabel(ax_r2, 'time (L-J)');
        ylabel(ax_r2, 'MSD (L-J)');
        title(ax_r2, ['n = ' num2str(n) ';  T = ' num2str(T) ';  t_{lin} = ' num2str(big_time) '\tau;  c_{line} = ' num2str(eps_err(Lmd.a, Lmd.k)) '; c_{log} = ' num2str(eps_err(Lmd.log_t, Lmd.log_l))]);
        legend(ax_r2, 'show', 'Location', 'best');
        if(close_figs)
            savefig(fig_r2, fullfile(model_name, [fig_r2_name '.fig']));
            saveas(fig_r2, fullfile(model_name, [fig_r2_name '_whole.png']));            
        end        
        x_draw = t_line(1)*1.5;
        xlim(ax_r2, [0 x_draw]);
        ylim(ax_r2, [-a*1.1 msd_th(x_draw, tau, T)]);
        if(close_figs)
            legend(ax_r2, 'off');
            saveas(fig_r2, fullfile(model_name, [fig_r2_name '_main.png']));
            close(fig_r2);
        end

        fig_log_name = 'log';
        fig_log = figure;
        ax_log = axes;
        plot(ax_log, t, msd, 'DisplayName', 'MSD');
        hold(ax_log, 'on');
        x_draw = [t_ballistic(end) t_line(end)];
        plot(ax_log, x_draw, exp(polyval(log_end_fit, log(x_draw))), ':', 'DisplayName', 'ballistic fit');
        x_draw = [t_ballistic(1) t_line(1)];
        plot(ax_log, x_draw, exp(polyval(log_beg_fit, log(x_draw))), ':', 'DisplayName', 'linear fit');
        plot(ax_log, t, msd_th_my, '--', 'DisplayName', 'theory');
        plot(ax_log, [tau tau], [a / 100, 5*a], '-', 'Color', [0 0 0], 'DisplayName', '\tau');
        x_draw = big_time * tau;
        plot(ax_log, [x_draw x_draw], [msd_th(x_draw, tau, T)/5, msd_th(x_draw, tau, T)*5], '--', 'Color', [0 0 0], 'DisplayName', 'line time');
        x_draw = tau * 2;
        plot(ax_log, [x_draw x_draw], [msd_th(x_draw, tau, T)/5, msd_th(x_draw, tau, T)*5], ':', 'Color', [0 0 0], 'DisplayName', 'theory intercept');
        x_draw = small_time * tau;
        plot(ax_log, [x_draw x_draw], [msd_th(x_draw, tau, T)/5, msd_th(x_draw, tau, T)*5], '-.', 'Color', [0 0 0], 'DisplayName', 'ballistic time');
        grid(ax_log, 'on');
        xlabel(ax_log, 'time (L-J)');
        ylabel(ax_log, 'MSD (L-J)');
        set(ax_log, 'XScale', 'log');
        set(ax_log, 'YScale', 'log');
        title(ax_log, ['n = ' num2str(n) ';  T = ' num2str(T) ';  t_{lin} = ' num2str(big_time) '\tau;  c_{line} = ' num2str(eps_err(Lmd.a, Lmd.k)) '; c_{log} = ' num2str(eps_err(Lmd.log_t, Lmd.log_l))]);    
        legend(ax_log, 'show', 'Location', 'best');    
        if(close_figs)
            savefig(fig_log, fullfile(model_name, [fig_log_name '.fig']));
            saveas(fig_log, fullfile(model_name, [fig_log_name '_whole.png']));                   
        end        
        x_draw = [t_ballistic(end)/1.5 t_line(1)*2.5];
        xlim(ax_log, x_draw);
        ylim(ax_log, msd_th(x_draw, tau, T));
        if(close_figs)
            legend(ax_log, 'off');
            saveas(fig_log, fullfile(model_name, [fig_log_name '_main.png']));                   
            close(fig_log);
        end
                
        fig_r1_name = 'r1';
        fig_r1 = figure;
        ax_r1 = axes;
        plot(ax_r1, t, sqrt(msd), 'DisplayName', 'MSD^{1/2}');
        hold(ax_r1, 'on');
        x_draw = [0 t_line(1)];
        plot(ax_r1, x_draw, polyval(r1_fit, x_draw), ':', 'DisplayName', 'asymptote');
        plot(ax_r1, t, sqrt(msd_th_my), '--', 'DisplayName', 'theory');
        %plot(ax_r1, [0, max(t)], [0 0], 'Color', [0 0 0], 'DisplayName', 'zero', 'LineWidth', 1);
        plot(ax_r1, [tau tau], [0, sqrt(1.3*a)], '-', 'Color', [0 0 0], 'DisplayName', '\tau');
        x_draw = big_time * tau;
        plot(ax_r1, [x_draw x_draw], [0, sqrt(msd_th(x_draw, tau, T) + a/2)], '--', 'Color', [0 0 0], 'DisplayName', 'line time');
        x_draw = small_time * tau;
        plot(ax_r1, [x_draw x_draw], [0, sqrt(msd_th(x_draw, tau, T) + a/2)], '-.', 'Color', [0 0 0], 'DisplayName', 'ballistic time');
        grid(ax_r1, 'on');
        xlabel(ax_r1, 'time (L-J)');
        ylabel(ax_r1, 'MSD^{1/2} (L-J)');
        title(ax_r1, ['n = ' num2str(n) ';  T = ' num2str(T) ';  t_{lin} = ' num2str(big_time) '\tau;  c_{line} = ' num2str(eps_err(Lmd.a, Lmd.k)) '; c_{log} = ' num2str(eps_err(Lmd.log_t, Lmd.log_l))]);        
        legend(ax_r1, 'show', 'Location', 'best');        
        if(close_figs)
            savefig(fig_r1, fullfile(model_name, [fig_r1_name '.fig']));
            saveas(fig_r1, fullfile(model_name, [fig_r1_name '_whole.png']));            
        end
        x_draw = t_line(1) * 1.5;
        xlim(ax_r1, [0 x_draw]);
        ylim(ax_r1, [0 sqrt(msd_th(x_draw, tau, T))]);
        if(close_figs)
            legend(ax_r1, 'off');
            saveas(fig_r1, fullfile(model_name, [fig_r1_name '_main.png']));
            close(fig_r1);
        end        
        
        fig_r2err_name = 'r2_err';
        fig_r2err = figure;
        ax_r2err = axes;
        plot(ax_r2err, t, eps_err(msd, msd_th_my), 'DisplayName', 'error');
        hold(ax_r2err, 'on');
        grid(ax_r2err, 'on');
        xlabel(ax_r2err, 'time (L-J)');
        ylabel(ax_r2err, 'error');
        title(ax_r2err, ['|r^2 / r_{th}^2 - 1| | n = ' num2str(n) ';  T = ' num2str(T) ';  t_{lin} = ' num2str(big_time) '\tau;']);
        legend(ax_r2err, 'show', 'Location', 'best');        
        if(close_figs)
            savefig(fig_r2err, fullfile(model_name, [fig_r2err_name '.fig']));
            saveas(fig_r2err, fullfile(model_name, [fig_r2err_name '_whole.png']));            
        end
        
        set(0, 'DefaultLineLineWidth', line_w_was);        
    end
end

