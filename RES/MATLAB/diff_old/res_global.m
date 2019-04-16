clear;
close all;
[n, T, lmd_th, lmd_e] = load_data('./DATA/2048_new.txt');

do_th = [1 1 1];
%do_exp = [0 1 1 1];
do_exp = [0 1 0];
cmp_th = 0;
cmp_exp = 0;

brd = [-0.5, 0.5];
ticks = [0.3 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.2, 1.4, 1.7, 2, 2.5, 3];
%tit = 'log_{10}(err)';
tit = '|\lambda / \lambda_0 - 1|';
log_scale = 0;

if(cmp_th)
    cmp_Lmd(T, n, lmd_th{1}, lmd_th{2}, brd, [tit ' | THEORY_0  vs  THEORY_1'], ticks, log_scale);
    cmp_Lmd(T, n, lmd_th{1}, lmd_th{3}, brd, [tit ' | THEORY_0  vs  THEORY_2'], ticks, log_scale);
    cmp_Lmd(T, n, lmd_th{2}, lmd_th{3}, brd, [tit ' | THEORY_1  vs  THEORY_2'], ticks, log_scale);
end

if(cmp_exp)
    cmp_Lmd(T, n, lmd_e{1}, lmd_e{2}, brd, [tit ' | EXP_1  vs  EXP_2'], ticks, log_scale);
    cmp_Lmd(T, n, lmd_e{1}, lmd_e{3}, brd, [tit ' | EXP_1  vs  EXP_3'], ticks, log_scale);
    cmp_Lmd(T, n, lmd_e{2}, lmd_e{3}, brd, [tit ' | EXP_2  vs  EXP_3'], ticks, log_scale);
end

for i_th = 1:length(do_th)
    for i_exp = 1:length(do_exp)
        if(do_th(i_th) && do_exp(i_exp))
            cmp_Lmd(T, n, lmd_e{i_exp}, lmd_th{i_th}, brd, [tit ' | THEORY_' num2str(i_th - 1) '  vs  EXP_' num2str(i_exp)], ticks, log_scale);
        end
    end
end

clear;
