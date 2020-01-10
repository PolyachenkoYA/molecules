clear;
close all;

model_name = 'tst';
[names, colors] = def_names;

E_data = dlmread(fullfile(names.data_path, model_name, names.E_filename));
time = dlmread(fullfile(names.data_path, model_name, names.time_filename));

H = E_data(:,4);
E = E_data(:,3);
Eerr = abs(E/mean(E) - 1);
Herr = abs(H/mean(H) - 1);

getFig('time', 'H', [model_name '; E_{std} = ' num2str(std(Eerr)) '; H_{std} = ' num2str(std(Herr))], 'linear', 'log');
plot(time, Eerr, 'DisplayName', 'E');
plot(time, Herr, 'DisplayName', 'H');

