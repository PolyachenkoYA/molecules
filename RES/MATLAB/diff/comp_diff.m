function [r2, v_cor] = comp_diff(N, x, v)
    r2 = zeros(N,1);
    v_cor = zeros(N,1);
    dx = x(:, :, :) - x(1, :, :);
    %parfor i = 1:N
    for i = 1:N
        r2(i) = mean(dx(i, :, 1) .* dx(i, :, 1) + dx(i, :, 2) .* dx(i, :, 2) + dx(i, :, 3) .* dx(i, :, 3));
        v_cor_c = 0;
        for j = 1:(N - i + 1)
            v_cor_c = v_cor_c + mean(v(j, :, 1) .* v(i+j-1, :, 1) + v(j, :, 2) .* v(i+j-1, :, 2) + v(j, :, 3) .* v(i+j-1, :, 3));
        end
        disp(['v cor proc: ' num2str((1 - (1-i/N)^2)*100) ' %']);
        v_cor(i) = v_cor_c / (N - i + 1);
    end
end

