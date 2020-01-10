%clear;
close all;

need_to_load = 0;
reproc_data = 0;
r2_time_cut = 1;
cubo_time_cut = 20;
scale = 'linear';
[names, colors] = def_names;
names.model_name = 'nu0.1';
%names.model_name = 'tst';
%names.model_name = 'T0.5_n1.3';
DATA_path = fullfile(pwd, 'DATA', 'Anderson', '1');
model_path = fullfile(DATA_path, names.model_name);

% --------------------------- load data ---------------------------------
[params, head] = read_params(model_path, names);
time = dlmread(fullfile(model_path, names.time_filename));
lin_i = time > r2_time_cut;
cubo_time_i = time < cubo_time_cut;
N = head.framesN;

if(need_to_load)
    data = zeros(N, head.Ntot, 6);
    for i = 1:N
        data(i, :, :) = dlmread(fullfile(model_path, names.frames_path, [num2str(i-1) '.' names.frame_file_ext]), '', 2, 0);
        disp(['data load: ' num2str(i / head.framesN * 100) ' %']);
    end
    x = data(:, :, 1:3);
    v = data(:, :, 4:6);
    data = [];
end

% -------------------------- process data -------------------------------
if(reproc_data)
    [D_th, D_int, D_Einst, lin_fit, r2, v_cor] = get_diff_coef(params, time, lin_i, cubo_time_cut, x, v, N);
    % [D_th, D_int, D_Einst, lin_fit] = get_diff_coef(params, time, lin_i, x, v, N)
end

draw_diff(D_th, D_int, D_Einst, time, r2, lin_fit, v_cor, cubo_time_cut, scale);
