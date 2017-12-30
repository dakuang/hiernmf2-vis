function WHinit = init_nmf_hier5(X, k)

[m, n] = size(X);

WHinit.Winit = rand(m, 4*k);
WHinit.Hinit = rand(4*k, n);
