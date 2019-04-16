%clear;
close all;

need_to_reload_data = 1;
mode = 2;
plot_all = 0;

N_arr = [512 1024 2048];
[names, colors] = def_names;
names.group_name = 'time10';
names.group_name = 'time500';
names.group_name = 'T2_wrong';
names.group_name = 'T2';
names.group_name = 'T2_n0.6';

switch mode
    case 1
        N_T = length(T_arr);
        N_n = length(n_arr);
        if((N_T * N_n > 116) && plot_all)
            disp('Too many plots to draw');
            return;
        end

        % ----------------------- load data -------------------
        group_path = fullfile(names.data_path, names.group_name);
        if(~exist('pressure', 'var'))
            model_names = cell(N_T, N_n);    
            params = cell(N_T, N_n);
            heads = cell(N_T, N_n);
        end
        if(need_to_reload_data)
            for Ti = 1:N_T
                for ni = 1:N_n
                    model_names{Ti, ni} = ['T' num2str(T_arr(Ti)) '_n' num2str(n_arr(ni))];
                    model_full_path = fullfile(group_path, model_names{Ti, ni});
                    [params{Ti, ni}, heads{Ti, ni}] = read_params(model_full_path, names);        

                    if((Ti == 1) && (ni == 1))
                        time = dlmread(fullfile(model_full_path, names.time_filename));

                        framesN = heads{1,1}.framesN;
                        Ntot = params{1,1}.Ntot;
                        pressure = zeros(N_T, N_n, framesN);
                        D = params{1,1}.R * 2;            
                    else
                        if(heads{Ti,ni}.framesN ~= framesN)
                            disp('wrong params of');
                            disp(model_full_path);
                            return;
                        end            
                    end

                    pressure(Ti, ni, :) = dlmread(fullfile(model_full_path, names.pressure_filename));
                    disp(['load data: ' num2str(100 * (ni + (Ti-1)*N_n) / (N_T*N_n)) ' %']);
                end
            end
        end

        % -------------------------- proc data --------------------------
        P_av = zeros(N_T, N_n);
        P_std = zeros(N_T, N_n);
        P_err = zeros(N_T, N_n);
        rel_err = zeros(N_T, N_n);
        Pn_figs = cell(1,N_T);
        Pn_axes = cell(1,N_T);
        Pn_legs = cell(1,N_T);
        for Ti = 1:N_T
            for ni = 1:N_n        
                P_av(Ti, ni) = mean(pressure(Ti, ni, :));
                %P_av(Ti, ni) = 2 * n_arr(ni) * T_arr(Ti) - P_av(Ti, ni);
                P_std(Ti, ni) = std(pressure(Ti, ni, :));
                P_err(Ti, ni) = P_std(Ti, ni) / sqrt(framesN - 1);
                rel_err(Ti, ni) = P_err(Ti, ni) / P_av(Ti, ni);

                if(plot_all)
                    [fig, ax, leg] = getFig('time', 'P', ...
                                     ['T = ' num2str(T_arr(Ti)) '; n = ' num2str(n_arr(ni)) ' | err = ' num2str(rel_err(Ti, ni))]);

                    plot(time, squeeze(pressure(Ti, ni, :)), 'DisplayName', 'precise values');
                    x_draw = [min(time) max(time)];
                    plot(x_draw, [1 1] * P_av(Ti, ni), 'DisplayName', ['P_{av} = ' num2str(P_av(Ti, ni)) ' \pm ' num2str(P_std(Ti, ni))]);
                end
            end

            small_n_ind = 1:5;
            small_n_fit = polyfit(n_arr(small_n_ind), P_av(Ti, small_n_ind), 1);
            T_small_n = small_n_fit(1);

            [Pn_figs{Ti}, Pn_axes{Ti}, Pn_legs{Ti}] = getFig('\rho', 'P', ...
                                ['P(\rho) | T \approx ' num2str(T_small_n) '; err =  ' num2str(eps_err(T_arr(Ti), T_small_n))]);
            errorbar(Pn_axes{Ti}, n_arr, P_av(Ti, :), P_err(Ti, :), 'DisplayName', ['T = ' num2str(T_arr(Ti))]);
            %x_draw = [n_arr(small_n_ind(1)) n_arr(small_n_ind(end))];
            x_draw = [min(n_arr) max(n_arr)];
            plot(Pn_axes{Ti}, x_draw, polyval(small_n_fit, x_draw), 'DisplayName', 'linear fit');
            plot(Pn_axes{Ti}, x_draw, x_draw * T_arr(Ti), 'DisplayName', 'P = \rho T');
        end
    case 2
end  