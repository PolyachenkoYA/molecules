function [dt_arr, comp_time, dt_time_fit, Estd_rel, Estd_fit, N_arr, Ntime_fit, Nstd_fit, rc_arr, rc_time_fit, model, params]...
            = group_proc_times(names, group_i, to_draw)
    dt_arr = [];
    comp_time = [];
    dt_time_fit = [];
    Estd_rel = [];
    Estd_fit = [];
    N_arr = [];
    Ntime_fit = [];
    Nstd_fit = [];
    rc_arr = [];
    rc_time_fit = [];
    model = [];
    params = [];

    %N_models = length(dt_arr);
    model_dirs = dir(fullfile(names.data_path, names.group_names{group_i}));
    if(isempty(model_dirs))
        disp(['no "' fullfile(names.data_path, names.group_names{group_i}) '" found']);
        return;
    end
    model = determine_model(names.group_base_model_names{group_i}, names);
    if(isempty(model))
        return;
    end
    if(~to_draw.Estd && ~to_draw.dt_time && (model.type == names.dt_model.type))
        return;
    end
    if(~to_draw.Estd_N && ~to_draw.N_time && (model.type == names.N_model.type))
        return;
    end
    if(~to_draw.rc_time && (model.type == names.rc_model.type))
        return;
    end
    
    [N_models, model_names, model_paths, model_prop_v] = getGroupProps(model_dirs, model.suff);
    
    comp_time = zeros(1,N_models);
    Estd_rel = zeros(1,N_models);
    Estd_abs = zeros(1,N_models);
    for i = 1:N_models        
        switch(model.type)
            case names.dt_model.type
                dt_arr(i) = model_prop_v(i);
            case names.N_model.type
                N_arr(i) = model_prop_v(i);
            case names.rc_model.type
                rc_arr(i) = model_prop_v(i);
        end        
        
        model_full_path = fullfile(model_paths{i}, model_names{i});        
        [params, head] = read_params(model_full_path, names);
        params.R = params.n;
        params.n = params.Ntot / (2 * params.R)^3;
        [comp_time(i), Estd_abs(i), Estd_rel(i)] = analize(model_full_path, names, params, head);
    end
    
    switch(model.type)
        case names.dt_model.type
            dt_real = 1 ./ dt_arr;          
            dt_time_fit = polyfit(log(dt_real), log(comp_time), 1);
            Estd_fit = polyfit(log(dt_real), log(Estd_rel), 1);  
        case names.N_model.type            
            N_points = 2;
            Ntime_fit = polyfit(log(N_arr((end - N_points + 1):end)), log(comp_time((end - N_points + 1):end)), 1);
            Nstd_fit =  polyfit(log(N_arr), log(Estd_rel), 1);   
        case names.rc_model.type
            fit_points = (rc_arr <= params.R*(1 + 10*eps));            
            %rc_time_fit = polyfit(rc_arr(fit_points), comp_time(fit_points), 6);            
            rc_time_fit = polyfit(rc_arr(fit_points), comp_time(fit_points), 3);            
    end
end

