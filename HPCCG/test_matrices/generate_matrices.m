dim = 1000;
n_samples = 1024;

cond_start = 10;
cond_end = 10000;
cond_space = linspace(cond_start, cond_end, n_samples);

for k = 1:numel(cond_space)
    A = sprandsym(dim, 0.2, 1/cond_space(k), 2);
    condition_number = condest(A);
    dm2hb(sprintf('A_%06d_%d.dat', k, condition_number), A);
end