clear;
close all;

compr_eff(32, []);
compr_eff([], 32000);
compr_eff(64, 65000);

getFig('p', 'launch time (s)', 't(p)', 'log', 'log');
p = 2.^(0:4);
plot(p, omp_lnch_t(p) * 1e-9);
