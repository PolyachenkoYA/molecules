function res_file_create(base_path, n, T, filename, line_time, ballistic_time, sazerlend_time)
    res_file = fopen(fullfile(base_path, filename),'w');
    
    fprintf(res_file, '%f ', n);
    fprintf(res_file, '\n');
    fprintf(res_file, '%f ', T);
    fprintf(res_file, '\n');
    done_job = 0;
    whole_job = 0;
    for ni = n
        for Ti = T
            whole_job = whole_job + 1 / (ni * r_th(Ti, 3)^2);
        end
    end
    for ni = n
        for Ti = T
            [k_exp, Lmd, consist] = process_MSD(base_path, ni, Ti, ballistic_time, line_time, sazerlend_time, 0, 0, 0);       
            % [k_exp, Lmd, d_c] = draw_MSD(n, T, small_time, big_time, print_res, graph_mod, close_figs)        

            fprintf(res_file,'%f %f %f %f %f %f %f %f %f %f %f\n', Lmd.th_0, Lmd.th_S, Lmd.th_T, Lmd.a, Lmd.k, Lmd.log_t, Lmd.log_l, consist.th_line, consist.th_log, consist.exp_line, consist.exp_log);
            %disp([num2str(ni) ' ' num2str(Ti)  ' ' num2str(Lmd.th_0)  ' ' num2str(Lmd.th_S)  ' ' num2str(Lmd.th_T)  ' ' num2str(Lmd.a)  ' ' num2str(Lmd.k)  ' ' num2str(Lmd.log_t)  ' ' num2str(Lmd.log_l)]);
            
            done_job = done_job + 1 / (ni * r_th(Ti, 3)^2);
            disp([num2str(100 * done_job / whole_job) ' % done']);
        end
    end
    fclose(res_file);
end

