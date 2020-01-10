%clear;
close all;

need_to_load = 1;
reproc_data = 0;
if(need_to_load)
    clear;
    need_to_load = 1;
    reproc_data = 1;
end
draw_all =  0;
to_draw = [0 0 0 1]; % [r2, cubo, |cubo|, Estd]

nu = [0.01, 0.1, 1];
nu = [0.0001, 0.0002, 0.0003, 0.0004, 0.0007, 0.0011, 0.0018, 0.0030, 0.0048, 0.0078, 0.0127, 0.0207, 0.0336, 0.0546, 0.0886, 0.1438, 0.2336, 0.3793, 0.6158];
N = length(nu);
r2_time_cut = 1;
cubo_time_cut = 20;
[names, colors] = def_names;
DATA_path = fullfile(pwd, 'DATA', 'Andersen');

% --------------------------- load data ---------------------------------
if(need_to_load)
    params = cell(1,N);
    head = cell(1,N);
    model_path = cell(1,N);
    for ni = 1:N
        model_path{ni} = fullfile(DATA_path, ['nu' num2str(nu(ni))]);
        [params{ni}, head{ni}] = read_params(model_path{ni}, names);
        if(ni == 1)
            time = dlmread(fullfile(model_path{ni}, names.time_filename));
            lin_i = time > r2_time_cut;
            Nfrm = head{1}.framesN;
            data = zeros(Nfrm, head{1}.Ntot, 6);
            x = zeros(N, Nfrm, head{1}.Ntot, 3);
            v = zeros(N, Nfrm, head{1}.Ntot, 3);
            E = zeros(N, Nfrm, 3);
        end

        for i = 1:Nfrm
            data(i, :, :) = dlmread(fullfile(model_path{ni}, names.frames_path, [num2str(i-1) '.' names.frame_file_ext]), '', 2, 0);
            disp(['data load: ' num2str(((ni - 1) * Nfrm + i) / Nfrm / N * 100) ' %']);
        end
        x(ni, :, :, :) = data(:, :, 1:3);
        v(ni, :, :, :) = data(:, :, 4:6);  
        
        E(ni, :, :) = dlmread(fullfile(model_path{ni}, names.E_filename));
    end
    data = [];
end

% -------------------------- process data -------------------------------
if(reproc_data)
    r2 = zeros(N, Nfrm);
    v_cor = zeros(N, Nfrm);
    D_th = zeros(1,N);
    D_int = zeros(1,N);
    D_Einst = zeros(1,N);
    lin_fit = cell(1,N);
    E_err = zeros(1,N);
    for ni = 1:N
        [D_th(ni), D_int(ni), D_Einst(ni), lin_fit{ni}, r2(ni, :), v_cor(ni, :)] = ...
            get_diff_coef(params{ni}, time, lin_i, cubo_time_cut, squeeze(x(ni,:,:,:)), squeeze(v(ni,:,:,:)), Nfrm);
        % [D_th, D_int, D_Einst, lin_fit, r2, v_cor] = get_diff_coef(params, time, lin_i, cubo_time_cut, x, v, Nfrm)
                
        E_err(ni) = std(E(ni, :, 3)) / mean(E(ni, :, 3));
    end
end

% -------------------------- draw res -------------------------------
if(draw_all)
    for ni = 1:N
        draw_diff(D_th(ni), D_int(ni), D_Einst(ni), ...
                  time, squeeze(r2(ni, :)), lin_fit{ni}, squeeze(v_cor(ni, :)), cubo_time_cut,...
                  'linear', to_draw);
        % draw_diff(D_th, D_int, D_Einst, time, r2, lin_fit, v_cor, cubo_time_cut, scale)
        
        if(to_draw(4))
            getFig('time', '\sigma_E/<E>', ['nu = ' num2str(nu(ni))]);
            plot(time, E(ni, :, 3) / mean(E(ni, :, 3)));
        end
    end
end

getFig('\nu', 'D', 'D(\nu)', 'log', 'log');
plot(nu, D_int, 'DisplayName', 'cubo');
plot(nu, D_Einst, 'DisplayName', 'Einst');

getFig('\nu', '\sigma_E/<E>', '(\sigma_E/<E>)(\nu)', 'log', 'log');
plot(nu, E_err, 'DisplayName', '\sigma_E/<E>');

