function [mpi_model, omp_model, cuda_model] = accum_eff_data(file_path)
    if(~exist('file_path', 'var'))
        file_path = [];
    end
    mpi_model = eff_analysis('MPI', 0);
    omp_model = eff_analysis('OMP', 0);
    cuda_model = eff_analysis('CUDA', 0);
    if(~isempty(file_path))
        save('DATA/eff_data', 'cuda_model', 'mpi_model', 'omp_model');
    end
end