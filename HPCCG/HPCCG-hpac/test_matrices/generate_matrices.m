dim = 1000;
n_samples = 1024;

cond_start = 10;
cond_end = 10000;
cond_space = linspace(cond_start, cond_end, n_samples);

A = sprandsym(dim, 0.2, 1/cond_space(1), 2);
dm2hb('A_000001.dat', 'A')
for k = 2:numel(cond_space)
    A = sprandsym(dim, 0.2, 1/cond_space(k), 2);
    dm2hb(sprintf('A_%06d.dat', k), 'A')
end