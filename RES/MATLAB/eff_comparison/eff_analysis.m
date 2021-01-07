function mdl = eff_analysis(model_name, need_to_draw)
    mdl.model_name = model_name;
    draw_time = 0;
    save_pics = 0;
    draw_th = 0;
    draw_full = 1;
    need_to_comp = 0;
    p_normalize = 0;
    
    mdl.time_unit = 3600 * 1e9;
    filepath = 'DATA';
    filename = '';
    if(isempty(filename))
        filename = mdl.model_name;
    end
    mdl.data = load(fullfile(filepath, [filename '.dat']));
    
    if(strcmp(mdl.model_name, 'OMP'))
        mdl.p_arr = [1, 2, 4, 8, 16];
        mdl.N_arr = 2.^(4:17);
    else
        if(strcmp(mdl.model_name, 'MPI'))
            mdl.p_arr = 2.^(0:6);
            mdl.N_arr = 2.^(6:17);        
        else
            mdl.p_arr = 2.^(3:9);
            mdl.N_arr = 2.^(9:17);                
        end
    end
    [mdl.p_grid, mdl.N_grid] = meshgrid(mdl.p_arr, mdl.N_arr);    
    if(strcmp(mdl.model_name, 'OMP'))
        dim = 5;
        mdl.p_fit_ind = [];
        drawcut_ind = [];

        mdl.freq = 2.5;
        mdl.omp_fix = mdl.p_grid / 10;
        mdl.omp_fix = 0;
        %mdl.e_th = 1 ./ (24.5 ./ mdl.p_grid + 14 ./ mdl.N_grid + 1000 ./ mdl.N_grid ./ mdl.p_grid + 500 .* (mdl.p_grid ./ mdl.N_grid.^2) + mdl.omp_fix) * mdl.time_unit * mdl.freq;
        mdl.e_th = 1 ./ (25 ./ mdl.p_grid + 14 ./ mdl.N_grid + 1000 ./ mdl.N_grid ./ mdl.p_grid + ...
                         omp_lnch_t(mdl.p_grid) ./ mdl.N_grid.^2 / mdl.freq + mdl.omp_fix) * mdl.time_unit * mdl.freq;
    else
        if(strcmp(mdl.model_name, 'MPI'))
            dim = 3;
            mdl.p_fit_ind = [];        
            drawcut_ind = [];

            mdl.freq = 3.1;
            mdl.e_th = 1 ./ (35./mdl.p_grid + 1000 ./ mdl.N_grid ./ mdl.p_grid + ...
                             36./mdl.N_grid + 270 / mdl.freq * (mdl.p_grid./mdl.N_grid.^2) + ...
                             150 * (sqrt(mdl.p_grid) + 1)./mdl.N_grid.^(3/2)) * mdl.time_unit * mdl.freq;
            %mdl.e_th = 1 ./ (((30./mdl.p_grid + 13./mdl.N_grid + 100 * (mdl.p_grid./mdl.N_grid.^2)) / M) ... 
            %           + 8 ./ mdl.N_grid .* (0.1./bondw(mdl.N_grid ./ mdl.p_grid * 8) + 0.1./bondw(mdl.N_grid * 8)));            
        else  % CUDA            
            dim = 3;
            mdl.p_fit_ind = [];
            drawcut_ind = [];

            mdl.freq = 1.6 / 12;
            t_gm = 200;
            mdl.cuda_fix = 10 ./ mdl.N_grid ./ mdl.p_grid + mdl.p_grid ./ mdl.N_grid.^2 * 250 + mdl.p_grid  / 40000;
            %mdl.cuda_fix = 0;
            minN_grid = min(2048, mdl.N_grid);
            maxP_grid = max(32, mdl.p_grid);
            mdl.e_th = 1 ./ (72 ./ mdl.N_grid + t_gm ./ minN_grid ./ mdl.N_grid + t_gm ./ minN_grid ./ mdl.p_grid + ...
                             32 .* maxP_grid ./ minN_grid ./ mdl.p_grid + 1000 ./ minN_grid ./ mdl.N_grid + ...
                             mdl.cuda_fix) * mdl.time_unit * mdl.freq;
%            mdl.e_th = 1 ./ (72 ./ mdl.N_grid + t_gm ./ minN_grid ./ mdl.p_grid + ...
%                             32 .* maxP_grid ./ minN_grid ./ mdl.p_grid + (1000 + t_gm) ./ minN_grid ./ mdl.N_grid + mdl.cuda_fix) * mdl.time_unit * mdl.freq;                        
        end
    end    

    if(need_to_comp)
        max_iter = 400;
        N_fit_ind = 1:size(mdl.data, 1);
        if(isempty(mdl.p_fit_ind))
            mdl.p_fit_ind = 1:size(mdl.data, 2);
        end
        lmd = 0.01;    

        mdl.e_th = eff_fit(mdl.data, mdl.time_unit, dim, max_iter, N_fit_ind, mdl.p_fit_ind, lmd);
    end

    if(need_to_draw)
        eff_draw(mdl.data, mdl.e_th, filepath, mdl.model_name, mdl.p_arr, mdl.N_arr, drawcut_ind, draw_time, save_pics, draw_th, draw_full, p_normalize);
    end
end