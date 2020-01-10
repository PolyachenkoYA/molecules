clear;
close all;

draw_all =  0;

nu = [0.0001, 0.0002, 0.0003, 0.0004, 0.0007, 0.0011, 0.0018, 0.0030, 0.0048, 0.0078, 0.0127, 0.0207, 0.0336, 0.0546, 0.0886, 0.1438, 0.2336, 0.3793, 0.6158];
N = length(nu);
[names, colors] = def_names;
DATA_path = fullfile(pwd, 'DATA', 'Andersen');

params = cell(1,N);
head = cell(1,N);
model_path = cell(1,N);
for ni = 1:N
    model_path{ni} = fullfile(DATA_path, ['nu' num2str(nu(ni))]);
    disp(model_path{ni});
    [params{ni}, head{ni}] = read_params(model_path{ni}, names);
    if(ni == 1)
        Nfrm = head{1}.framesN;
        time = dlmread(fullfile(model_path{ni}, names.time_filename));
        E = zeros(N, Nfrm, 3);
    end

    E(ni, :, :) = dlmread(fullfile(model_path{ni}, names.E_filename));
end

E_err = zeros(1,N);
Cv = zeros(1,N);
x_draw = zeros(1,N);
for ni = 1:N
    E_mean = mean(E(ni, :, 3));
    E_err(ni) = std(E(ni, :, 3)) / E_mean;
    Cv(ni) = (mean(E(ni, :, 3).^2) - E_mean^2) / params{ni}.Ntot /params{ni}.Tmp^2;
    x_draw(ni) = nu(ni);
end

getFig('\nu', 'C_V', 'C_V(\nu)', 'log', 'log');
%getFig('\nu', 'C_V', 'C_V(\nu)');
plot(nu, Cv, 'DisplayName', 'C_V');




