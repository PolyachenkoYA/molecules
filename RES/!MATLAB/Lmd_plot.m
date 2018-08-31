[n, T, lmd_th0, lmd_th1, lmd_th2, lmd_1, lmd_2] = load_data('./DATA/Lmd_data2048.txt');

do_th = [1 1 0];
cmp_th = 0;
cmp_exp = 1;

brd = [-0.5, 0.5];
log_scale = 0;

if(cmp_th)
    cmp_Lmd(T, n, lmd_th0, lmd_th1, brd, 'log10(err) | THEORY_0  vs  THEORY_1', log_scale);
    cmp_Lmd(T, n, lmd_th0, lmd_th2, brd, 'log10(err) | THEORY_0  vs  THEORY_2', log_scale);
end

if(cmp_exp)
    cmp_Lmd(T, n, lmd_1, lmd_2, brd, 'log10(err) | EXP_1  vs  EXP_2', log_scale);
end

if(do_th(1))
    cmp_Lmd(T, n, lmd_1, lmd_th0, brd, 'log10(err) | THEORY_0  vs  EXP_1', log_scale);
    cmp_Lmd(T, n, lmd_2, lmd_th0, brd, 'log10(err) | THEORY_0  vs  EXP_2', log_scale);
end

if(do_th(2))
    cmp_Lmd(T, n, lmd_1, lmd_th1, brd, 'log10(err) | THEORY_1  vs  EXP_1', log_scale);
    cmp_Lmd(T, n, lmd_2, lmd_th1, brd, 'log10(err) | THEORY_1  vs  EXP_2', log_scale);
end

if(do_th(3))
    cmp_Lmd(T, n, lmd_1, lmd_th2, brd, 'log10(err) | THEORY_2  vs  EXP_1', log_scale);
    cmp_Lmd(T, n, lmd_2, lmd_th2, brd, 'log10(err) | THEORY_2  vs  EXP_2', log_scale);
end

clear;
