close all;
clear;

single_fig = 1;
do_th = [1 1 1];
N = 1000;
Tbnd = [1  5];
fig_names = ["th_0 Obchfiz" "th_1 (sazerland)" "th_2 T>1"];

xT = linspace(Tbnd(1), Tbnd(2), N);    
fig_abs = figure('Name', 'abolute theory comparison');
ax_abs = axes(fig_abs);
hold on;
grid on;
title('\lambda(T), n = 0.001');

fig_rel = figure('Name', 'relative theory comparison');
ax_rel = axes(fig_rel);
hold on;
grid on;
title('|\lambda(T) / \lambda_0(T) - 1|');
for i_do = 1:length(do_th)
    plot(ax_abs, xT, Lmd_th(r_th(xT, i_do), 0.001), 'DisplayName', fig_names(i_do));
    xlabel(ax_abs, 'T');
    ylabel(ax_abs, '\lambda (\sigma)');
    
    if(i_do > 1)
        %plot(ax_rel, xT, lmd_th(r_th(xT, i_do), 0.001) ./ lmd_th(r_th(xT, 1), 0.001) - 1, 'DisplayName', fig_names(i_do));
        plot(ax_rel, xT, eps_err(Lmd_th(r_th(xT, i_do), 0.001) , Lmd_th(r_th(xT, 1), 0.001)), 'DisplayName', fig_names(i_do));
        xlabel(ax_rel, 'T');
        ylabel(ax_rel, '|\lambda(T) / \lambda_0(T) - 1|');
    end
end

legend(ax_abs, 'show');
legend(ax_abs, 'Location','best');
legend(ax_rel, 'show');
legend(ax_rel, 'Location','best');
