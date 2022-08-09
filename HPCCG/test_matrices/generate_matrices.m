dim = 2000;
n_samples = 1024;

cond_start = 10^2;
cond_end = 10^7;
density = 0.1;
generation_method = 2;
cond_space = linspace(cond_start, cond_end, n_samples);

for k = 1:numel(cond_space)
    A = sprandsym(dim, density, 1/cond_space(k), generation_method);    
    condition_number = condest(A);
    dm2hb(sprintf('A_%06d_%d.dat', k, condition_number), A);
end