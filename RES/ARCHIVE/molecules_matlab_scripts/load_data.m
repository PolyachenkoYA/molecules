function [n, T, l_th, l_e] = load_data(filename)
    data = dlmread(filename);
    
    k = 1;
    n = data(:,k); k = k+1;
    T = data(:,k); k = k+1;
    l_th{1} = data(:,k); k = k+1;
    l_th{2} = data(:,k); k = k+1;
    l_th{3} = data(:,k); k = k+1;
    
    l_e{1} = data(:,k); k = k+1;
    l_e{2} = data(:,k); k = k+1;
    l_e{3} = data(:,k); k = k+1;
%    l_e{4} = data(:,k); k = k+1;
end

