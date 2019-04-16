function [N_models, model_names, model_paths, arr_v] = getGroupProps(model_dirs, suff)
    model_dirs(1:2) = [];
    model_names = {model_dirs.name};
    model_paths = {model_dirs.folder};
    N_models = sum([model_dirs(~ismember({model_dirs.name},{'.','..'})).isdir]);   
    
    l = length(model_names);
    arr_v = zeros(1, l);
    for i = 1:l
        arr_v(i) = str2double(model_names{i}((length(suff)+1):end));
    end
    
    [arr_v, sort_i] = sort(arr_v);
    model_names_0 = model_names;
    model_path_0 = model_paths;
    for i = 1:l
        model_names{i} = model_names_0{sort_i(i)};
        model_paths{i} = model_path_0{sort_i(i)};
    end
end

