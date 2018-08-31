function params = read_params(filename)
    params_data = importdata(filename); 

    i__ = 1;
    params.Ntot = params_data.data(i__); i__ = i__+1;
    params.k12 = params_data.data(i__); i__ = i__+1;
    params.k6 = params_data.data(i__); i__ = i__+1;
    params.R = params_data.data(i__); i__ = i__+1;
    params.r_max = params_data.data(i__); i__ = i__+1;
    params.endT = params_data.data(i__); i__ = i__+1;
    params.totalT = params_data.data(i__); i__ = i__+1;
    params.dumpDT = 1/params_data.data(i__); i__ = i__+1;
    params.dt = 1/params_data.data(i__); i__ = i__+1;
    params.Tmp = params_data.data(i__); i__ = i__+1;
    params.mu = params_data.data(i__); i__ = i__+1;
    params.kv = params_data.data(i__); i__ = i__+1;
    params.dissipK = params_data.data(i__); i__ = i__+1;
    params.TmpStabEps = params_data.data(i__); i__ = i__+1;
    params.TmpStabGap = params_data.data(i__); i__ = i__+1;
    params.compMode = params_data.data(i__); i__ = i__+1;
    params.binOutF = params_data.data(i__); i__ = i__+1;
end

