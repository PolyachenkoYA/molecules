%clear;
close all;

need_to_load = 1;
time_cut = 1;
scale = 'linear';
[names, colors] = def_names;
names.model_name = 'tst_diff_EvsC';
%names.model_name = 'tst';
%names.model_name = 'T0.5_n1.3';
DATA_path = fullfile(pwd, 'DATA');
model_path = fullfile(DATA_path, names.model_name);

% --------------------------- load data ---------------------------------
[params, head] = read_params(model_path, names);
time = dlmread(fullfile(model_path, names.time_filename));
lin_i = time > time_cut;
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
D_th = sqrt(params.Tmp / pi) * 2 / (3 * pi * params.n);
r2 = zeros(N,1);
v_cor = zeros(N,1);
dx = x(:, :, :) - x(1, :, :);
parfor i = 1:N
    r2(i) = mean(dx(i, :, 1) .* dx(i, :, 1) + dx(i, :, 2) .* dx(i, :, 2) + dx(i, :, 3) .* dx(i, :, 3));
    v_cor_c = 0;
    for j = 1:(N - i + 1)
        v_cor_c = v_cor_c + mean(v(j, :, 1) .* v(i+j-1, :, 1) + v(j, :, 2) .* v(i+j-1, :, 2) + v(j, :, 3) .* v(i+j-1, :, 3));
    end
    disp(i/N)
    v_cor(i) = v_cor_c / (N - i + 1);
end
D_int = sum(v_cor) * ((1/params.dumpDT) / 3);
lin_fit = polyfit(time(lin_i), r2(lin_i), 1);
D_Einst = lin_fit(1) / 6;

[fig, ax, leg] = getFig('time', '<r^2>',...
    ['D_{Einst} = ' num2str(D_Einst) ', D_{th} = ' num2str(D_th) ' | err = ' num2str(eps_err(D_Einst, D_th) * 100) ' %'],...
    scale, scale);
plot(time, r2, 'DisplayName', 'exp');
plot(time, polyval(lin_fit, time), '--', 'DisplayName', 'line fit');
%plot(time, r2 ./ (time * 6 * D_th), 'DisplayName', 'exp');

[fig2, ax2, leg2] = getFig('time', '<v(0)v(t)>',...
    ['D_{cubo} = ' num2str(D_int) ', D_{th} = ' num2str(D_th) ' | err = ' num2str(eps_err(D_int, D_th) * 100) ' %'],...
                           'linear', 'linear');
plot(time, v_cor, 'DisplayName', 'exp');
%plot(time, polyval(lin_fit, time), '--', 'DisplayName', 'line fit');

[fig3, ax3, leg3] = getFig('time', '|<v(0)v(t)>|',...
    ['D_{cubo} = ' num2str(D_int) ', D_{th} = ' num2str(D_th) ' | err = ' num2str(eps_err(D_int, D_th) * 100) ' %'],...
                           'linear', 'log');
plot(time, abs(v_cor), 'DisplayName', 'exp');
