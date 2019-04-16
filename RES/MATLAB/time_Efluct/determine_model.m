function model = determine_model(s, names)
    if(strcmp(s, names.dt_model.suff))
        model = names.dt_model;
    elseif(strcmp(s, names.N_model.suff))
        model = names.N_model;
    elseif(strcmp(s, names.rc_model.suff))
        model = names.rc_model;
    else
        disp('wrong model suff');
        model = [];
    end
end

