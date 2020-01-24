function w = bondw(N)
    sz = size(N);
    lin_ind = N < 3e3;
    w(lin_ind) = N(lin_ind).^0.8 * 1.32e7;
    w(~lin_ind) = 8e9;
    w = reshape(w, sz);
end

