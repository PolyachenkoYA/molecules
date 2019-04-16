function [comp_time, Estd_abs, Estd_rel, N] = analize(full_path, names, params, head)
    comp_time = head.real_time;
    
    E = dlmread(fullfile(full_path, names.E_filename));
    Estd_abs = std(E(:,3))^2;
    Estd_rel = Estd_abs / mean(E(:,3))^2;
    N = params.Ntot;
end

