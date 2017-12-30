function run_hiernmf2(data_file, num_leafs)

load(data_file);

alg_matrix = cell(1, 1);
alg_matrix{1} = @matrix_orig;

alg_init = cell(1, 1);
alg_init{1} = @init_nmf_hier5;

alg = cell(1, 1);

alg{1}.matrix = 1;
alg{1}.init = 1;
alg{1}.func = @(X, k, WHinit, computeobj, dataset, run)alg_nmfsh_comb_flat8(X, k, WHinit, 2, true, 'projgrad', computeobj, @anls_entry_rank2_precompute, @anls_entry_blockpivot_precompute);
alg{1}.hier = false;
alg{1}.init_for_hier = @init_flatnmf8_for_hier;
alg{1}.computeobj = true;

i = 1;
fix_k = [num_leafs];
label = {};
label_names = {};
clustering_hier8_large_WH(i, X, fix_k(i), label, label_names, voc, 'flat2', [data_file, '_', num2str(i)], '', 1, [1], alg_matrix, alg_init, alg{:});

load(['flat2_data_', data_file, '_', num2str(i), '.mat']);
W = final_WH_summary{:, :, end}.W;
H = final_WH_summary{:, :, end}.H;
top_keywords = final_voc_summary{:, :, end};
save(['NMF_result_', data_file, '_', num2str(num_leafs), '.mat'], 'W', 'H', 'top_keywords', '-v7.3');
