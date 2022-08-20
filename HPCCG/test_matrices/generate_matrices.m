dim = 200;
density = 0.1;
generation_method = 2;
cond_space = readmatrix("./cond_space.txt");

for k = 1:numel(cond_space)
    A = sprandsym(dim, density, 1/cond_space(k), generation_method);    
    condition_number = condest(A);
    dm2hb(sprintf('./A_%06d_%d.dat', k, condition_number), A);
end