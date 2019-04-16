close all;
%clear;

[names, colors] = def_names;

need_to_reload_data = 0;
names.group_name = 'time5_T0.44';
names.group_name = 'time5_T2_4M';
names.group_name = 'tdm_time10';
N_rep = 11;
models = [128 256 512 1024];
%models = [128,    256,    512,    640,   1024,   1280,   2048,   2560,   3200,   4096,   5120,   6400,   8192,  10240,  12800,  16000,  16384,  20480,  25600,  32000,  32768,  40960,  51200];
% 512, 2048, 6400, 8192
%models = [128,    256,    512,    640,   1024,   1280,   2048,   2560,   3200,   4096,   5120,   6400,   8192,  10240,  12800,  16000,  16384,  20480,  25600,  32000,  32768,  40960];
%models = [32000, 32768];

N_models = length(models);
model_suff = 'dt';
to_draw.r2 = 0;
to_draw.v2 = 1;
to_draw.v2_err = 1;
if((to_draw.r2 || to_draw.v2 || to_draw.v2_err) && (N_models > 4))
    disp(['too many (' num2str(N_models * (N_models - 1)/2) ') to draw']);
    return;
end

group_path = fullfile(names.data_path, names.group_name);
model_dirs = dir(group_path);
if(isempty(model_dirs))
    disp(['no "' fullfile(names.data_path, names.group_names{group_i}) '" found']);
    return;
end

% --------------------------------- load data ---------------------------
if(~exist('frames', 'var'))
    params = cell(N_models, N_rep);
    heads = cell(N_models, N_rep);
    model_names = cell(N_models, 1);
    Tmp = zeros(N_models, 1);
end
if(need_to_reload_data)
    for i1 = 1 : N_models
        model_names{i1} = [model_suff num2str(models(i1))];
        for i2 = 1:N_rep
            model_full_path = fullfile(group_path, [model_names{i1} '_' num2str(i2)]);
            [params{i1, i2}, heads{i1, i2}] = read_params(model_full_path, names, 'old_mode');

            if((i1 == 1) && (i2 == 1))
                time = dlmread(fullfile(model_full_path, names.time_filename));
                time(1) = [];
                % time has to be the same for all models because we analize each
                % pair
                
                framesN = heads{1,1}.framesN;
                Ntot = params{1,1}.Ntot;
                frames = zeros(N_models, N_rep, framesN, Ntot, 6);
                D = params{1,1}.R * 2;
            else
                if((heads{i1,i2}.framesN ~= framesN) || (params{i1,i2}.Ntot ~= Ntot) || (eps_err(params{i1,i2}.R, D/2) > 0.0001))
                    disp('wrong params of');
                    disp(model_full_path);
                    return;
                end
            end 

            for i3 = 1 : framesN 
                frames(i1, i2, i3, :, :) = dlmread(fullfile(model_full_path, names.frames_path, [num2str(i3 - 1) '.' names.frame_file_ext]), '', 2, 0);
                disp(['load data: ' num2str(100 * (i3 + ((i2 - 1) + (i1 - 1) * N_rep) * framesN) / (N_models * framesN * N_rep)) ' %']);
            end
        end
        Tmp(i1) = params{i1, 1}.Tmp;
    end
    Tmp_av = mean(Tmp);
    v2_th = 6 * Tmp_av;
end

% --------------------------------- average data -------------------------------
r2_av = zeros(N_rep, framesN, 1);
v2_av = zeros(N_rep, framesN, 1);
r2_final = zeros(framesN, 1);
v2_final = zeros(framesN, 1);
dq = zeros(framesN, Ntot, 6);
N_points = N_models * (N_models - 1) / 2;
stab_t = cell(N_models, 1);
td = cell(N_models, 1);
dt_ratio = cell(N_models, 1);
tdm_k = zeros(N_models - 2, 1);
tdm = zeros(N_models - 2, 1);
if((N_models - 1) > length(colors))
    disp('not enough colors');
    disp(['N_models = ' num2str(N_models) '; N_colors = ' num2str(length(colors))]);
    return;
end
for i1 = 1 : (N_models - 1)
    td{i1} = zeros(N_models - i1, 1);
    stab_t{i1} = zeros(N_models - i1, 1);
    dt_ratio{i1} = zeros(N_models - i1, 1);
    for i2 = (i1 + 1) : N_models
        for i3 = 1:N_rep      
            dq(:,:,:) = frames(i1, i3, :, :, :) - frames(i2, i3, :, :, :);
            unwrap_k = 0.8;
            unwrap_ind = dq(:,:,names.x_clmn : names.z_clmn) > unwrap_k*D;
            dq(unwrap_ind) = dq(unwrap_ind) - D;
            unwrap_ind = dq(:,:,names.x_clmn : names.z_clmn) < -unwrap_k*D;
            dq(unwrap_ind) = dq(unwrap_ind) + D;
            
            r2 = (dq(:,:,names.x_clmn).^2 + dq(:,:,names.y_clmn).^2 + dq(:,:,names.z_clmn).^2);
            v2 = (dq(:,:,names.vx_clmn).^2 + dq(:,:,names.vy_clmn).^2 + dq(:,:,names.vz_clmn).^2);
            for i_f = 1 : framesN
                r2_av(i3, i_f) = mean(r2(i_f, :));
                v2_av(i3, i_f) = mean(v2(i_f, :));
            end
            disp(['proc data: ' num2str(100 * ((i1 - 1) * (N_models - i1/2) + i2 - i1 - 1 + i3 / N_rep) / (N_models * (N_models - 1) / 2)) ' %']);    
        end
        for i_f = 1 : framesN                    
            r2_final(i_f) = mean(r2_av(:, i_f));
            v2_final(i_f) = mean(v2_av(:, i_f));
        end
        r2_final(1) = [];
        v2_final(1) = [];
        i2_0 = i2 - i1;
        plot_name = [num2str(models(i1)) ' / ' num2str(models(i2)) ' = ' num2str(models(i1)/models(i2))];
        
        v2_stab_err = eps_err(ones(framesN - 1, 1) * v2_th, v2_final);        
        stab_t{i1}(i2_0) = time(get_tdi(v2_stab_err, 0.5)); 
        dt_ratio{i1}(i2_0) = models(i1) / models(i2);
        if((i1 == 1) && (i2 == 2))
            [td_fig, td_ax, td_leg] = getFig('dt_1 / dt_2', 't_m');
        end
        fit_ind = ((time > stab_t{i1}(i2_0) * 0.1) .* (time < stab_t{i1}(i2_0) * 2/3)) > 0;
        fit_v2_t = time(fit_ind);
        log_fit = polyfit(fit_v2_t, log(v2_final(fit_ind)), 1);
        v2_tm_err = eps_err(exp(polyval(log_fit, time)), v2_final);
        td{i1}(i2_0) = (log(v2_th) - log_fit(2)) / log_fit(1);
        fit_ind = ((time > (td{i1}(i2_0) + 1)) .* (time < (td{i1}(i2_0) + 2))) > 0;
        fit_r2_t = time(fit_ind);
        r2_fit = polyfit(fit_r2_t, r2_final(fit_ind), 1);
        stab_t{i1}(i2_0) = - r2_fit(2)/r2_fit(1);
        
        if(to_draw.r2)
            getFig('time', '\Delta r^2', '\Delta r^2(time)', 'linear', 'log');
            plot(time, r2_final, 'DisplayName', plot_name);
            x_draw = [1 1] * stab_t{i1}(i2_0);
            y_draw = [min(r2_final) max(r2_final)*1.1];
            plot(x_draw, y_draw, 'DisplayName', ['t_{stab} = ' num2str(stab_t{i1}(i2_0))]);
            x_draw = [1 1] * td{i1}(i2_0);
            plot(x_draw, y_draw, 'DisplayName', ['t_m = ' num2str(td{i1}(i2_0))]);            
            
            getFig('time', '\Delta r^2', '\Delta r^2(time)', 'linear', 'linear');
            plot(time, r2_final, 'DisplayName', plot_name);
            x_draw = [1 1] * stab_t{i1}(i2_0);
            y_draw = [min(r2_final) max(r2_final)*1.1];
            plot(x_draw, y_draw, 'DisplayName', ['t_{stab} = ' num2str(stab_t{i1}(i2_0))]);
            x_draw = [1 1] * td{i1}(i2_0);
            plot(x_draw, y_draw, 'DisplayName', ['t_m = ' num2str(td{i1}(i2_0))]);
            x_draw = [stab_t{i1}(i2_0) fit_r2_t'];
            plot(x_draw, polyval(r2_fit, x_draw), 'DisplayName', 'linear fit');
        end

        if(to_draw.v2)
            getFig('time', '\Delta v^2', '\Delta v^2(time)', 'linear', 'log');
            plot(time, v2_final, 'DisplayName', plot_name);
            x_draw = [1 1] * stab_t{i1}(i2_0);
            y_draw = [min(v2_final) max(v2_final)*2];
            plot(x_draw, y_draw, 'DisplayName', ['t_{stab} = ' num2str(stab_t{i1}(i2_0))]);
            x_draw = [1 1] * td{i1}(i2_0);
            plot(x_draw, y_draw, 'DisplayName', ['t_m = ' num2str(td{i1}(i2_0))]);            
            x_draw = [td{i1}(i2_0) max(time)];
            y_draw = [1 1] * v2_th;
            plot(x_draw, y_draw, 'DisplayName', 'v^2_{th}');
            x_draw = [min(fit_v2_t) td{i1}(i2_0)];
            plot(x_draw, exp(polyval(log_fit, x_draw)), 'DisplayName', 'exp fit');
        end

        if(to_draw.v2_err)        
            getFig('time', '|v^2 / v_{th}^2 - 1|', 'V_{err}(time)', 'linear', 'log');
            plot(time, v2_stab_err, 'DisplayName', [plot_name '; v^2_{stab}']);           
            plot(time, v2_tm_err, 'DisplayName', [plot_name '; v^2_{exp}']);
            x_draw = [1 1] * stab_t{i1}(i2_0);
            y_draw = [min(v2_stab_err) max(v2_stab_err)*2];
            plot(x_draw, y_draw, 'DisplayName', ['t_{stab} = ' num2str(stab_t{i1}(i2_0))]);
            x_draw = [1 1] * td{i1}(i2_0);
            plot(x_draw, y_draw, 'DisplayName', ['t_m = ' num2str(td{i1}(i2_0))]);            
        end
    end     
    plot(td_ax, dt_ratio{i1}, td{i1}, 'x', 'DisplayName', num2str(1/models(i1), '%1.1e'), 'Color', colors(i1, :));
    if(length(dt_ratio{i1}) > 1)
        td_fit = polyfit(dt_ratio{i1}, td{i1}, 1);
        x_draw = [0, max(dt_ratio{i1})];
        plot(td_ax, x_draw, polyval(td_fit, x_draw), 'HandleVisibility', 'off', 'Color', colors(i1, :));
        tdm_k(i1) = td_fit(1);
        tdm(i1) = td_fit(2);
    end
end


x_tdm = 1./models(1:(end-2))';
tdm_fit = polyfit(log(x_tdm), tdm, 1);

getFig('dt', 't^d_m', ['t^d_m(dt) | tg = ' num2str(tdm_fit(1))], 'log', 'linear');
plot(x_tdm, tdm, 'x', 'DisplayName', 't^d_m(dt)');
x_draw = [min(x_tdm) max(x_tdm)];
plot(x_draw, polyval(tdm_fit, log(x_draw)), 'HandleVisibility', 'off');


