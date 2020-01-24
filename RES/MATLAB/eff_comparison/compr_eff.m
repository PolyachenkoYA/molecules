function compr_eff(p, N)
    N = 2^round(log(N) / log(2));
    [mpi_model, omp_model, cuda_model] = accum_eff_data;

    if(~isempty(p))
        [~, ax_pfix, ~] = ...
            getFig('N', '$\eta$', ['$p = ' num2str(p) '; \hspace{10pt} \eta = T/dt \cdot N^2 / t_{wallclock}$'], 'log', 'log');
    else
        ax_pfix = [];
    end
    if(~isempty(N))
        [~, ax_Nfix, ~] = ...
            getFig('p', '$\eta$', ['$N = ' num2str(N) '; \hspace{10pt} \eta = T/dt \cdot N^2 / t_{wallclock}$'], 'log', 'log');
    else
        ax_Nfix = [];    
    end
    add_model_to_compar(p, N, mpi_model, ax_pfix, ax_Nfix, getMyColor(2));
    add_model_to_compar(p, N, omp_model, ax_pfix, ax_Nfix, getMyColor(3));
    add_model_to_compar(p, N, cuda_model, ax_pfix, ax_Nfix, getMyColor(4));
end