function run_hiernmf2(data_file, num_leafs)

load(data_file);

alg_matrix = cell(1, 1);
alg_matrix{1} = @matrix_orig;

alg_init = cell(1, 1);
alg_init{1} = @init_nmf_hier5;

alg = cell(1, 1);

alg{1}.matrix = 1;
alg{1}.init = 1;
alg{1}.func = @(X, k, WHinit, random_run)alg_nmfsh_comb_hier8(X, k, WHinit, random_run, 3, 2, true, 'projgrad', 'ndcg_part', true, @anls_entry_rank2_precompute);
alg{1}.hier = true;

i = 1;
fix_k = [num_leafs];
label = {};
label_names = {};
clustering_hier8(i, X, fix_k(i), label, label_names, voc, 'hier8', [data_file, '_', num2str(i)], '', 1, [1], alg_matrix, alg_init, alg{:});
