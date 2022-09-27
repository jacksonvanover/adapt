dim = 200;
density = 0.1;
generation_method = 1;
cond_space = readmatrix("./cond_space.txt");

if license('checkout','Distrib_Computing_Toolbox') == 0
    for k = 1:numel(cond_space)
        A = sprandsym(dim, density, 1/cond_space(k), generation_method);    
        dm2hb(sprintf('./A_%06d_%d.dat', k, cond_space(k)), A);
    end
else
    p = parpool;

    parfor k = 1:numel(cond_space)
        A = sprandsym(dim, density, 1/cond_space(k), generation_method);    
        dm2hb(sprintf('./A_%06d_%d.dat', k, cond_space(k)), A);
    end
end