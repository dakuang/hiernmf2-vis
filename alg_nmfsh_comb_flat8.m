function [cluster, iter, obj, timing, W, H] = alg_nmfsh_comb_flat8(X, k, WHinit, vec_norm, normW, conv, computeobj, anls_alg_rank2, anls_alg, tol, maxiter)

if ~exist('computeobj', 'var')
	computeobj = true;
end

if ~exist('tol', 'var')
	tol = 1e-4;
end

if ~exist('maxiter', 'var')
	maxiter = 10000;
end

t0 = tic;

trial_allowance = 3;
[clusters, objs, timings, Ws, priorities] = alg_nmfsh_comb_hier8_large(X, k, WHinit, trial_allowance, vec_norm, normW, conv, 'ndcg_part', anls_alg_rank2, tol, maxiter);
W = Ws{end};
clear clusters objs timings Ws priorities;
H = nnlsm_blockpivot(W, X, 0);

timing = toc(t0);

obj.obj = norm(X, 'fro')^2 - 2 * trace(W' * (X * H')) + trace((W' * W) * (H * H'));
obj.grad = 0;
iter = 0;

[max_val, cluster] = max(H);
