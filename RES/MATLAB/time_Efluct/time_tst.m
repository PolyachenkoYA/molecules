close all;
clear;

% round(exp(linspace(log(128), log(10000), 10)))
% unique(M.','rows').' % not tested

[names, colors] = def_names;
names.dt_model.type = 1;
names.N_model.type = 2;
names.rc_model.type = 3;
names.dt_model.suff = 'dt';
names.N_model.suff = 'N';
names.rc_model.suff = 'rc';
names.group_names = {'LP_rc_big', 'LP_long', 'LP_N_Estd', 'LP_dt_t1e2', 'LP_dt_t1e3', 'LP_dt_t1e4'};
N_groups = length(names.group_names); % '' is the last, we don't need it
names.group_base_model_names = {names.rc_model.suff,...
                                names.dt_model.suff,...
                                names.N_model.suff,...
                                names.dt_model.suff,...
                                names.dt_model.suff,...
                                names.dt_model.suff};

if((length(colors) < max(N_groups, length(names.group_base_model_names))) ||...
   (length(names.group_base_model_names) ~= N_groups))
    disp('not enough colors/groups/base_names');
    return;
end

to_draw.Estd = 0;
to_draw.Estd = 1;
to_draw.Estd_N = 0;
to_draw.Estd_N = 1;
to_draw.dt_time = 0;
to_draw.dt_time = 1;
to_draw.N_time = 0;
to_draw.N_time = 1;
to_draw.rc_time = 0;
to_draw.rc_time = 1;

if(to_draw.dt_time)
    figure;
    time_ax = axes;
    hold(time_ax, 'on');
end
if(to_draw.N_time)
    figure;
    Ntime_ax = axes;
    hold(Ntime_ax, 'on');
end
if(to_draw.Estd)
    figure;
    std_ax = axes;
    hold(std_ax, 'on');
end
if(to_draw.Estd_N)
    figure;
    stdN_ax = axes;
    hold(stdN_ax, 'on');
end
if(to_draw.rc_time)
    figure;
    rc_ax = axes;
    hold(rc_ax, 'on');
    
    figure;
    rc_log_ax = axes;
    hold(rc_log_ax, 'on');    
end

comp_time_k = zeros(1, N_groups);
Estd_k = zeros(1, N_groups);
N_arr_k = zeros(1, N_groups);
for i = 1:N_groups
    [dt_arr{i}, comp_time{i}, time_fit{i}, Estd{i}, Estd_fit{i},...
     N_arr{i}, Ntime_fit{i}, Nstd_fit{i},...
     rc_arr{i}, rc_time_fit{i},...
     model{i}, params{i}] = ...
        group_proc_times(names, i, to_draw);
    if(isempty(model{i}))
        close all;
        return;
    end
    dt_real = 1 ./ dt_arr{i};
    
    if(to_draw.rc_time && (model{i}.type == names.rc_model.type))      
        R_a0 = params{i}.R/5;
        R_c0 = params{i}.R*sqrt(3)*(1 - eps*10);
        a0 = mean(comp_time{i}(rc_arr{i} <= R_a0));
        c0 = mean(comp_time{i}(rc_arr{i} >= R_c0));
                        
        alp = 1;
        l = 50000;
        min_std = 9999999;
        fit_points = (rc_arr{i} <= params{i}.R*(1 + 10*eps));
        r3 = rc_arr{i}(fit_points).^3;        
        r6 = r3.^2;
        y_exp = comp_time{i}(fit_points);
        b0 = (y_exp(end) - a0) / r6(end);
        b0 = rc_time_fit{i}(1);
        b_arr = linspace(b0*(1 - alp), b0*(1 + alp), l);          
        b1 = 0;        
        for j = 1:l            
            %std_bb = std(yexp_fit - (r6 * b_arr(j) + a0));
            std_bb = std(y_exp - (r3 * b_arr(j) + a0));
            if(min_std > std_bb)
                b1 = b_arr(j);
                min_std = std_bb;
            end
        end        
        %rc_time_fit{i} = [b1 0 0 0 0 0 a0]; 
        rc_time_fit{i} = [b1 0 0 a0]; 
        
        plot(rc_ax, rc_arr{i}, comp_time{i}, 'x', 'Color', colors(i, :),...
            'DisplayName', [names.group_names{i}...
                            '; (tf/tr)fit = ' num2str(round((params{i}.R^3 * 6/pi) * rc_time_fit{i}(1)/rc_time_fit{i}(4), 2)) ...
                            '; (tf/tr)c = ' num2str(c0 / a0 - 1)]);
                            %'; (tf/tr)fit = ' num2str(round((params{i}.R^3 * 6/pi)^2 * rc_time_fit{i}(1)/rc_time_fit{i}(7), 2)) ...                        
        y_exp = y_exp - a0;
        fit_points = (r3.^(1/3) >= params{i}.R *0.25).*(r3.^(1/3) < params{i}.R*(1 - 100*eps)) > 0;
        y_exp = y_exp(fit_points);
        x_exp = r3(fit_points).^(1/3);
        log_fit = polyfit(log(x_exp), log(y_exp), 1);
        plot(rc_log_ax, x_exp, y_exp, 'x', 'Color', colors(i, :), ...
            'DisplayName', [names.group_names{i} ...
                            '; p = ' num2str(log_fit(1))]);
        x = [x_exp(1) x_exp(end)];
        plot(rc_log_ax, x, exp(polyval(log_fit, log(x))), 'Color', colors(i, :), 'HandleVisibility', 'off');

        x = linspace(0, params{i}.R, 100);
        y = polyval(rc_time_fit{i}, x);
        plot(rc_ax, x, y, '-', 'Color', colors(i, :),'HandleVisibility','off');         
        plot(rc_ax, [0 params{i}.R/2], [a0 a0], 'DisplayName', 'a0');
        plot(rc_ax, [R_c0 max(rc_arr{i})], [c0 c0], 'DisplayName', 'c0');
        x = [params{i}.R params{i}.R];
        plot(rc_ax, x, [min(comp_time{i}) max(comp_time{i})], 'DisplayName', 'R cude', 'LineWidth', 2);
        x = [params{i}.R params{i}.R] * sqrt(3);
        plot(rc_ax, x, [min(comp_time{i}) max(comp_time{i})], 'DisplayName', 'Rmax', 'LineWidth', 2);
    end
    
    if(to_draw.dt_time && (model{i}.type == names.dt_model.type))
        plot(time_ax, dt_real, comp_time{i}, 'x', 'Color', colors(i, :),...
            'DisplayName', [names.group_names{i} ...
                           '; k = ' num2str(round(time_fit{i}(1), 2))...
                           '; std = ' num2str(round(std(comp_time{i} ./ exp(polyval(time_fit{i}, log(dt_real))) - 1),2))]);        
        x = log([min(dt_real) max(dt_real)]);
        plot(time_ax, exp(x), exp(polyval(time_fit{i}, x)), '-', 'Color', colors(i, :),'HandleVisibility','off'); 
    end
    
    if(to_draw.N_time && (model{i}.type == names.N_model.type))
        plot(Ntime_ax, N_arr{i}, comp_time{i}, 'x', 'Color', colors(i, :),...
            'DisplayName', [names.group_names{i} ...
                           '; k = ' num2str(round(Ntime_fit{i}(1), 2))...
                           '; std = ' num2str(round(std(comp_time{i} ./ exp(polyval(Ntime_fit{i}, log(N_arr{i}))) - 1),2))]);
        x = log([min(N_arr{i}) max(N_arr{i})]);
        plot(Ntime_ax, exp(x), exp(polyval(Ntime_fit{i}, x)), '-', 'Color', colors(i, :),'HandleVisibility','off');         
    end

    if(to_draw.Estd && (model{i}.type == names.dt_model.type))
        plot(std_ax, dt_real, Estd{i}, 'x', 'Color', colors(i, :),...
            'DisplayName', [names.group_names{i} ...
                           '; k = ' num2str(round(Estd_fit{i}(1), 2))...
                           '; std = ' num2str(round(std(Estd{i} ./ exp(polyval(Estd_fit{i}, log(dt_real))) - 1),2))]);
        x = log([min(dt_real) max(dt_real)]);    
        plot(std_ax, exp(x), exp(polyval(Estd_fit{i}, x)), '-', 'Color', colors(i, :),'HandleVisibility','off'); 
    end
    
    if(to_draw.Estd_N && (model{i}.type == names.N_model.type))
        plot(stdN_ax, N_arr{i}, Estd{i}, 'x', 'Color', colors(i, :),...
            'DisplayName', [names.group_names{i} ...
                           '; k = ' num2str(round(Nstd_fit{i}(1), 2))...
                           '; std = ' num2str(round(std(Estd{i} ./ exp(polyval(Nstd_fit{i}, log(N_arr{i}))) - 1),2))]);
        x = log([min(N_arr{i}) max(N_arr{i})]);    
        plot(stdN_ax, exp(x), exp(polyval(Nstd_fit{i}, x)), '-', 'Color', colors(i, :),'HandleVisibility','off'); 
    end    
end

if(to_draw.dt_time)
    time_leg = set_plot(time_ax, 'dt', 'real time', 'time(dt)', 'log', 'log');    
end

if(to_draw.N_time)
    Ntime_leg = set_plot(Ntime_ax, 'N', 'real time', 'time(N)', 'log', 'log');    
end

if(to_draw.Estd)
    std_leg = set_plot(std_ax, 'dt', 'E std', 'Estd(dt)', 'log', 'log');    
end

if(to_draw.Estd_N)
    Nstd_leg = set_plot(stdN_ax, 'N', 'E std', 'Estd(N)', 'log', 'log');
end

if(to_draw.rc_time)
    %rc_leg = set_plot(rc_ax, 'r_C', 'real time', 'time(r_C) | (t_r/t_f)_{th} = 0.68', 'linear', 'linear');
    rc_leg = set_plot(rc_ax, 'r_C', 'real time', 'time(r_C)', 'linear', 'linear');
    rc_log_leg = set_plot(rc_log_ax, 'r_C', 't - t0', 'time(r_C)', 'log', 'log');    
end
