function [params, head] = read_params(full_path, names, mode)
    %disp(full_path);
    data = importdata(fullfile(full_path, names.params_filename));

    k = 1;
    if(exist('mode', 'var'))
        % ------- old one ----------    
        params.Ntot = data.data(k); k = k + 1;
        params.k12 = data.data(k); k = k + 1;
        params.k6 = data.data(k); k = k + 1;
        params.n = data.data(k); k = k + 1;
        params.r_max = data.data(k); k = k + 1;
        params.endT = data.data(k); k = k + 1;
        params.totalT = data.data(k); k = k + 1;
        params.dumpT = data.data(k); k = k + 1;
        params.dt = data.data(k); k = k + 1;
        params.Tmp = data.data(k); k = k + 1;
        params.mu = data.data(k); k = k + 1;
        params.kv = data.data(k); k = k + 1;
        params.dissipK = data.data(k); k = k + 1;
        params.TmpStabEps = data.data(k); k = k + 1;
        params.TmpStabGap = data.data(k); k = k + 1;
        params.compMode = data.data(k); k = k + 1;
        params.binOutF = data.data(k); k = k + 1;                
    else
        params.Ntot = data.data(k); k = k + 1;
        params.n = data.data(k); k = k + 1;
        params.n_cr = data.data(k); k = k + 1;
        params.r_cut = data.data(k); k = k + 1;
        params.endT = data.data(k); k = k + 1;
        params.dumpDT = data.data(k); k = k + 1;
        params.dt = data.data(k); k = k + 1;
        params.Tmp = data.data(k); k = k + 1;  
        params.dissipK = data.data(k); k = k + 1;
        params.TmpStabEps = data.data(k); k = k + 1;
        params.TmpStabGap = data.data(k); k = k + 1;
        params.compMode = data.data(k); k = k + 1;
        params.binOutF = data.data(k); k = k + 1;        
    end
            
    params.R = (params.Ntot / params.n)^(1/3) / 2;
    
    data = dlmread(fullfile(full_path, names.head_filename));
    k = 1;
    head.framesN = data(k); k = k + 1;
    head.R = data(k); k = k + 1;
    head.binOutF = data(k); k = k + 1;
    head.Ntot = data(k); k = k + 1;
    head.real_time = data(k); k = k + 1;
end

