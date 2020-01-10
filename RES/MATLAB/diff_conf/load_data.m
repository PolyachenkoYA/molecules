function  data = load_data(filename)
    raw_data = dlmread(filename);
    data.n = raw_data(1, raw_data(1,:) > 0);
    data.T = raw_data(2, raw_data(2,:) > 0);
    raw_data(1:2, :) = [];
    data.Lmd = raw_data;
end
