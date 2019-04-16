close all;

need_to_reload_data = 0;
if(need_to_reload_data)
    clear;
    need_to_reload_data = 1;
end
plot_all = 0;
n_arr = [0.8, 0.85, 0.88, 0.9, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1, 1.04, 1.06, 1.08, 1.1, 1.12, 1.14, 1.16, 1.18, 1.2, 1.25, 1.3, 1.35, 1.4];
T_arr = [1, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 3];

[names, colors] = def_names;
DATA_path = fullfile(pwd, names.data_path);
names.group_name = 'lindemann';
% --------------------------- load data ---------------------------------
N_T = length(T_arr);
N_n = length(n_arr);
if((N_T + N_n > 30) && plot_all)
    disp('Too many plots to draw');
    return;
end

% ----------------------- load data -------------------
group_path = fullfile(names.data_path, names.group_name);
if(~exist('x', 'var'))
    model_names = cell(N_T, N_n);    
    params = cell(N_T, N_n);
    heads = cell(N_T, N_n);
end
if(need_to_reload_data)
    for ni = 1:N_n
        for Ti = 1:N_T
            model_names{Ti, ni} = ['T' num2str(T_arr(Ti)) '_n' num2str(n_arr(ni))];
            model_full_path = fullfile(group_path, model_names{Ti, ni});
            [params{Ti, ni}, heads{Ti, ni}] = read_params(model_full_path, names);        

            if((Ti == 1) && (ni == 1))
                time = dlmread(fullfile(model_full_path, names.time_filename));

                framesN = heads{1,1}.framesN;
                Ntot = params{1,1}.Ntot;
                D = params{1,1}.R * 2;            
                data = zeros(framesN, Ntot, 6);
                x = zeros(N_n, N_T, framesN, Ntot, 3);
            else
                if(heads{Ti,ni}.framesN ~= framesN)
                    disp('wrong params of');
                    disp(model_full_path);
                    return;
                end            
            end
            
            params{Ti, ni}.framesN = framesN;
            for i = 1:framesN
                data(i, :, :) = dlmread(fullfile(model_full_path, names.frames_path, [num2str(i-1) '.' names.frame_file_ext]), '', 2, 0);
            end
            x(ni, Ti, :, :, :) = data(:, :, 1:3);
            disp(['load data: ' num2str(100 * (Ti + (ni-1)*N_T) / (N_T*N_n)) ' %']);
        end
    end
end

% --------------------------- proc data ---------------------------------
l = zeros(N_n, N_T);
a = zeros(N_n);
for ni = 1:N_n
    for Ti = 1:N_T
        l(ni, Ti) = sqrt(get_mean_disp(squeeze(x(ni, Ti, :, :, :)), params{Ti, ni}));
    end
    a(ni) = (4/params{Ti, ni}.n)^(1/3) / sqrt(2);
    l(ni, :) = l(ni, :) / a(ni);
end
if(plot_all)
    [T_fig, T_ax, T_leg] = getFig('n', 'l');
    [n_fig, n_ax, n_leg] = getFig('T', 'l');
    for Ti = 1:N_T
        plot(T_ax, n_arr, l(:, Ti), '-x', 'DisplayName', ['T = ' num2str(T_arr(Ti))]);
    end
    for ni = 1:N_n
        plot(n_ax, T_arr, l(ni, :), '-x', 'DisplayName', ['n = ' num2str(n_arr(ni))]);
    end
end

max_l = 10;
l(l > max_l) = max_l;
[fig, ax, leg] = getFig('n', 'T', 'l(T,n)');
surf(n_arr, T_arr, l', 'EdgeColor', 'interp', 'FaceColor', 'interp');
cb = colorbar(ax);
legend(ax, 'off');
xlim(ax, [min(n_arr) max(n_arr)]);
ylim(ax, [min(T_arr) max(T_arr)]);

