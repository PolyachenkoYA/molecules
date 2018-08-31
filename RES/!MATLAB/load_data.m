function [n,T,th0,th1,th2,l1,l2] = load_data(filename)
    data = dlmread(filename);
    
    k = 1;
    n = data(:,k); k = k+1;
    T = data(:,k); k = k+1;
    th0 = data(:,k); k = k+1;
    th1 = data(:,k); k = k+1;
    th2 = data(:,k); k = k+1;
    
    l1 = data(:,k); k = k+1;
    l2 = data(:,k); k = k+1;
end

