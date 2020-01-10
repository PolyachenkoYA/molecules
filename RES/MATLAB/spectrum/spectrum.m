to_load = 1;  % switch here
if(to_load)
    clear;
    to_load = 1;
end
to_proc = 0; % switch here
if(to_load)
    to_proc = 1;
end
close all;
draw_all = 0;

[names, colors] = def_names;
names.model_name = 'tst100';
model_full_path = fullfile(names.data_path, 'spectrum', names.model_name);

% ------------------------- load data ------------------------------
if(to_load)
    [params, head] = read_params(model_full_path, names);
    N = head.framesN;
    Npart = params.Ntot;
    data = zeros(Npart, 6);
    time = dlmread(fullfile(model_full_path, names.time_filename));
    x = zeros(N, Npart, 3);
    for fi = 1:N
        data(:, :) = dlmread(fullfile(model_full_path, names.frames_path, [num2str(fi - 1) '.' names.frame_file_ext]), '', 2, 0);
        x(fi, :, :) = data(:, 4:6);
        disp(['load: ' num2str(100 * fi / N) ' %']);
    end
    data = [];
end

% --------------------------- proc data ---------------------------
if(to_proc)
    N1 = floor(N/2 + 1);
    A = zeros(Npart, N1, 3);
    for pi = 1:Npart
        for i = 1:3
            %x(:, pi, i) = x(:, pi, i) - mean(x(:, pi, i));
            [A(pi, :, i), w] = fft_m(x(:, pi, i), 1/params.dumpDT);
        end        
        disp(['proc: ' num2str(100 * pi / Npart) ' %']);        
    end
    A_av = zeros(N1, 3);
    x_av = zeros(N, 3);
    for fi = 1:N1        
        for i = 1:3
            A_av(fi, i) = mean(A(:, fi, i));
            x_av(fi, i) = mean(x(fi, :, i));
        end
    end
    A = [];
    for fi = 1:N1
        A_av(fi, 1) = mean(A_av(fi, :));
    end
    %A_av_tot = A_av(:, 1) / mean(A_av(:, 1));
    A_av_tot = A_av(:, 1) / sqrt(params.Tmp);
    %A_av(:, 1) = A_av(:, 1) / sqrt(params.Tmp);
    
    filt_N = 8;
    filt_coef = ones(1,filt_N)/filt_N;    
    A_av_filt = sym_filter(filt_coef, 1, A_av_tot);
    
    [r2, v_cor] = comp_diff(N, zeros(size(x)), x);
    [A_cor, w_cor] = fft_m(v_cor, 1/params.dumpDT);    
    %A_cor = A_cor / mean(A_cor);
    A_cor = A_cor / sqrt(8*params.Tmp/pi) / sqrt(params.Tmp);
    
    A_cor_filt = sym_filter(filt_coef, 1, A_cor);    
end

% ----------------------------------- draw res -------------------------

getFig('\omega', 'A_U', ['A_u(\omega); (\rho, T) = (' num2str(params.n) ';  ' num2str(params.Tmp) '); ' names.model_name]);
%getFig('\omega', 'A/<A>', ['A(\omega)/<A(\omega)>; (\rho, T) = (' num2str(params.n) ';  ' num2str(params.Tmp) '); ' names.model_name]);
plot(w, A_av_tot, 'DisplayName', '$F[v(t)](\omega) / \sqrt{T}$', 'Color', colors(1, :));
plot(w, A_av_filt, 'DisplayName', ['filtered $\pm ' num2str((filt_N-1)*(w(2)-w(1))) '$'], 'LineWidth', 2, 'Color', colors(1, :));
plot(w_cor, A_cor, 'DisplayName', '$F[V_{cor}(t)](\omega) / \sqrt{T} / \langle v \rangle$', 'Color', colors(2, :));
plot(w_cor, A_cor_filt, 'DisplayName', ['filtered $\pm ' num2str((filt_N-1)*(w_cor(2)-w_cor(1))) '$'], 'LineWidth', 2, 'Color', colors(2, :));
axis(gca, 'tight');

getFig('\omega', 'A_v/A_{cor}', ['A_v/A_{cor}; (\rho, T) = (' num2str(params.n) ';  ' num2str(params.Tmp) '); ' names.model_name]);
plot(w, A_av_filt./A_cor_filt, 'DisplayName', '$A_v(\omega)/A_{cor}(\omega)$');
%axis tight;

if(draw_all)
    for i = 1:3
        [fig1, ax1, leg1] = getFig('t', 'x', ['x_' num2str(i) '(t) | ' names.model_name]);
        plot(time,x(:, i, 1));

        [fig2, ax2, leg2] = getFig('w', 'A', ['A_' num2str(i) '(w) | ' names.model_name]);
        plot(w, A(i, :));
    end
end
