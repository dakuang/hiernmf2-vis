function [clusters, objs, timings, Ws, priorities] = alg_nmfsh_comb_hier8_large(X, k, WHinit, trial_allowance, vec_norm, normW, conv, scoring, anls_alg, tol, maxiter)

t0 = tic;

unbalanced = 0.1;

if ~exist('tol', 'var')
	tol = 1e-4;
end

if ~exist('maxiter', 'var')
	maxiter = 10000;
end

if ~exist('anls_alg', 'var')
	anls_alg = @anls_entry_rank2_precompute;
end

[m, n] = size(X);
clusters = cell(1, k-1);
objs = zeros(1, k-1);
timings = zeros(1, k-1);
Ws = cell(1, k-1);
priorities = cell(1, k-1);
init_used = 0;
term_subset = find(sum(X, 2) ~= 0);
X_subset = X(term_subset, :);
obj = 1e308;
init_used = init_used + 1;

[W, H] = nmfsh_comb_rank2(X_subset, 2, WHinit.Winit(term_subset, 2*(init_used-1)+1:2*init_used), WHinit.Hinit(2*(init_used-1)+1:2*init_used, :), vec_norm, normW, conv, tol, maxiter, anls_alg);
[max_val, cluster_subset] = max(H);
cluster = zeros(1, n);
cluster(cluster_subset == 1) = 0;
cluster(cluster_subset == 2) = 0.5;
levels = 0;
priority = 1e308;
W_buffer = cell(1, 1);
W_buffer{1} = zeros(m, 2);
W_buffer{1}(term_subset, :) = W;

for i = 1 : k-1
	clusters{i} = floor(cluster);
	unique_cluster = unique(clusters{i}(clusters{i} >= 0));
	[max_priority, unique_cluster_pos] = max(priority);
	if max_priority < 0 
		for j = i : k-1 
			objs(j) = objs(i-1);
			timings(j) = timings(i-1);
			Ws{j} = Ws{i-1};
			priorities{j} = priorities{i-1};
		end 
		break;
	end
	label = unique_cluster(unique_cluster_pos);
	split = 2^levels(unique_cluster_pos);
	clusters{i}(cluster == label + 0.5) = label + split;
	cluster(cluster == label + 0.5) = label + split;
	unique_cluster = unique(clusters{i}(clusters{i} >= 0));
	num_cols_move = length(unique(clusters{i}(clusters{i} > label+split)));
	label_pos = find(unique_cluster == label);
	label_split_pos = find(unique_cluster == label+split);
	levels(end+1) = 0;
	levels(end-num_cols_move+1:end) = levels(end-num_cols_move:end-1);
	levels(label_pos) = levels(label_pos) + 1;
	levels(label_split_pos) = levels(label_pos);
	Ws{i} = zeros(m, i+1);
	if i > 1
		Ws{i}(:, 1:i) = Ws{i-1};
		Ws{i}(:, end-num_cols_move+1:end) = Ws{i}(:, end-num_cols_move:end-1);
	end
	Ws{i}(:, label_pos) = W_buffer{unique_cluster_pos}(:, 1);
	Ws{i}(:, label_split_pos) = W_buffer{unique_cluster_pos}(:, 2);
	W_buffer{end+1} = [];
	W_buffer(end-num_cols_move+1:end) = W_buffer(end-num_cols_move:end-1);
	priority(end+1) = 0;
	priority(end-num_cols_move+1:end) = priority(end-num_cols_move:end-1);
	timings(i) = toc(t0);

	subset = find(cluster == label);
	init_used = init_used + 1;
	[cluster, W_buffer_one, priority_one] = trial_split(cluster, trial_allowance, unbalanced, priority, X, subset, label, Ws{i}(:, label_pos), scoring, WHinit, init_used, vec_norm, normW, conv, tol, maxiter,  anls_alg);
	W_buffer{label_pos} = W_buffer_one;
	priority(label_pos) = priority_one;

	subset = find(cluster == label+split);
	init_used = init_used + 1;
	[cluster, W_buffer_one, priority_one] = trial_split(cluster, trial_allowance, unbalanced, priority, X, subset, label+split, Ws{i}(:, label_split_pos), scoring, WHinit, init_used, vec_norm, normW, conv, tol, maxiter, anls_alg);
	W_buffer{label_split_pos} = W_buffer_one;
	priority(label_split_pos) = priority_one;
	
	priorities{i} = priority;
end

%--------------------------------------

function [cluster, W_buffer_one, priority_one] = trial_split(cluster, trial_allowance, unbalanced, priority, X, subset, new_label, W_parent, scoring, WHinit, init_used, vec_norm, normW, conv, tol, maxiter, anls_alg)

[m, n] = size(X);

trial = 0;
subset_backup = subset;
while trial < trial_allowance
	[cluster_subset, W_buffer_one, priority_one] = potential_split(X, subset, new_label, W_parent, scoring, WHinit, init_used, vec_norm, normW, conv, tol, maxiter, anls_alg);
	if priority_one < 0
		break;
	end
	unique_cluster_subset = unique(cluster_subset);
	if length(unique_cluster_subset) ~= 2
		error('Invalid number of unique sub-clusters!');
	end
	length_cluster1 = length(find(cluster_subset == unique_cluster_subset(1)));
	length_cluster2 = length(find(cluster_subset == unique_cluster_subset(2)));
	if min(length_cluster1, length_cluster2) < unbalanced * length(cluster_subset)
		[min_val, idx_small] = min([length_cluster1, length_cluster2]);
		subset_small = find(cluster_subset == unique_cluster_subset(idx_small));
		subset_small = subset(subset_small);
		[cluster_subset_small, W_buffer_one_small, priority_one_small] = potential_split(X, subset_small, new_label, W_buffer_one(:, idx_small), scoring, WHinit, init_used, vec_norm, normW, conv, tol, maxiter, anls_alg);
		if priority_one_small < min(priority(priority > 0))
			trial = trial + 1;
			if trial < trial_allowance
				disp(['Drop ', num2str(length(subset_small)), ' documents ...']);
				cluster(subset_small) = -1;
				subset = setdiff(subset, subset_small);
			end
		else
			break;
		end
	else
		break;
	end
end

if trial == trial_allowance
	W_buffer_one = zeros(m, 2);
	priority_one = -2;
	disp(['Recycle ', num2str(length(find(cluster(subset_backup) == -1))), ' documents ...']);
	cluster(subset_backup) = new_label;
else
	cluster(subset) = cluster_subset;
end

%--------------------------------------

function [cluster_subset, W_buffer_one, priority_one] = potential_split(X, subset, new_label, W_parent, scoring, WHinit, init_used, vec_norm, normW, conv, tol, maxiter, anls_alg)

[m, n] = size(X);
if length(subset) <= 3
	cluster_subset = new_label * ones(1, length(subset));
	W_buffer_one = zeros(m, 2);
	priority_one = -1;
else
	term_subset = find(sum(X(:, subset), 2) ~= 0);
	X_subset = X(term_subset, subset);
	obj = 1e308;
	[W, H] = nmfsh_comb_rank2(X_subset, 2, WHinit.Winit(term_subset, 2*(init_used-1)+1:2*init_used), WHinit.Hinit(2*(init_used-1)+1:2*init_used, subset), vec_norm, normW, conv, tol, maxiter, anls_alg);
	[max_val, cluster_subset] = max(H);
	cluster_subset(cluster_subset == 2) = new_label + 0.5;
	cluster_subset(cluster_subset == 1) = new_label;
	W_buffer_one = zeros(m, 2);
	W_buffer_one(term_subset, :) = W;
	if length(unique(cluster_subset)) > 1
		priority_one = compute_priority(W_parent, W_buffer_one, scoring);
	else
		priority_one = -2;
	end
end

%--------------------------------------

function priority = compute_priority(W_parent, W_child, scoring)

n = length(W_parent);
[sorted_parent, idx_parent] = sort(W_parent, 'descend');
[sorted_child1, idx_child1] = sort(W_child(:, 1), 'descend');
[sorted_child2, idx_child2] = sort(W_child(:, 2), 'descend');

switch scoring

case 'ndcg_topk'
k = 10;
idx_parent = idx_parent(1 : min(n, k));
idx_child1 = idx_child1(1 : min(n, k));
idx_child2 = idx_child2(1 : min(n, k));
intersection = intersect(idx_child1, idx_child2);
if length(intersection) ~= 0
	for i = 1 : length(idx_child1)
		if length(find(intersection == idx_child1(i))) > 0
			idx_child1(i) = 0;
		end
	end
	for i = 1 : length(idx_child2)
		if length(find(intersection == idx_child2(i))) > 0
			idx_child2(i) = 0;
		end
	end
end
priority = NDCG_topk(idx_parent, idx_child1) * NDCG_topk(idx_parent, idx_child2);

case 'ratio_topk'
k = 10;
idx_parent = idx_parent(1 : min(n, k));
idx_child1 = idx_child1(1 : min(n, k));
idx_child2 = idx_child2(1 : min(n, k));
intersection = intersect(idx_child1, idx_child2);
if length(intersection) ~= 0
	for i = 1 : length(idx_child1)
		if length(find(intersection == idx_child1(i))) > 0
			idx_child1(i) = 0;
		end
	end
	for i = 1 : length(idx_child2)
		if length(find(intersection == idx_child2(i))) > 0
			idx_child2(i) = 0;
		end
	end
end
priority = length(intersect(idx_parent, idx_child1)) * length(intersect(idx_parent, idx_child2)) / n / n;

case 'ndcg_full'
weight = log(n:-1:1)';
first_zero = find(sorted_parent == 0, 1);
if length(first_zero) > 0
	weight(first_zero:end) = 1;
end
[sorted, idx1] = sort(idx_child1);
[sorted, idx2] = sort(idx_child2);
max_pos = max(idx1, idx2);
discount = log(n-max_pos(idx_parent)+1);
discount(discount == 0) = log(2);
weight = weight ./ discount;
priority = NDCG_full(idx_parent, idx_child1, weight) * NDCG_full(idx_parent, idx_child2, weight);

case 'ndcg_part'
n_part = length(find(W_parent ~= 0));
if n_part <= 1
	priority = -3;
else
	weight = log(n:-1:1)';
	first_zero = find(sorted_parent == 0, 1); 
	if length(first_zero) > 0 
		weight(first_zero:end) = 1;
	end
	weight_part = zeros(n, 1);
	weight_part(1:n_part) = log(n_part:-1:1)';
	[sorted, idx1] = sort(idx_child1);
	[sorted, idx2] = sort(idx_child2);
	max_pos = max(idx1, idx2);
	discount = log(n-max_pos(idx_parent)+1);
	discount(discount == 0) = log(2);
	weight = weight ./ discount;
	weight_part = weight_part ./ discount;
	priority = NDCG_part(idx_parent, idx_child1, weight, weight_part) * NDCG_part(idx_parent, idx_child2, weight, weight_part);
end

case 'ndcg_fullpart'
weight = log(n:-1:1)';
first_zero = find(sorted_parent == 0, 1); 
if length(first_zero) > 0 
	weight(first_zero:end) = 1;
end
n_part = length(find(W_parent ~= 0));
weight_part = zeros(n, 1);
weight_part(1:n_part) = log(n:-1:n-n_part+1)';
[sorted, idx1] = sort(idx_child1);
[sorted, idx2] = sort(idx_child2);
max_pos = max(idx1, idx2);
discount = log(n-max_pos(idx_parent)+1);
discount(discount == 0) = log(2);
weight = weight ./ discount;
weight_part = weight_part ./ discount;
priority = NDCG_part(idx_parent, idx_child1, weight, weight_part) * NDCG_part(idx_parent, idx_child2, weight, weight_part);

end

%--------------------------------------

function score = NDCG_full(ground, test, weight)

[sorted, seq_idx] = sort(ground);
weight = weight(seq_idx);

n = length(test);
uncum_score = weight(test);
uncum_score(2:n) = uncum_score(2:n) ./ log2(2:n)';
cum_score = cumsum(uncum_score);

ideal_score = sort(weight, 'descend');
ideal_score(2:n) = ideal_score(2:n) ./ log2(2:n)';
cum_ideal_score = cumsum(ideal_score);

score = cum_score ./ cum_ideal_score;
score = score(end);

%--------------------------------------

function score = NDCG_part(ground, test, weight, weight_part)

[sorted, seq_idx] = sort(ground);
weight_part = weight_part(seq_idx);

n = length(test);
uncum_score = weight_part(test);
uncum_score(2:n) = uncum_score(2:n) ./ log2(2:n)';
cum_score = cumsum(uncum_score);

ideal_score = sort(weight, 'descend');
ideal_score(2:n) = ideal_score(2:n) ./ log2(2:n)';
cum_ideal_score = cumsum(ideal_score);

score = cum_score ./ cum_ideal_score;
score = score(end);

%--------------------------------------

function score = NDCG_topk(ground, test)

n = length(test);
uncum_score = zeros(1, n);
for i = 1 : n
	uncum_score(i) = (length(find(ground == test(i))) > 0);
end
uncum_score(2:n) = uncum_score(2:n) ./ log2(2:n);
cum_score = cumsum(uncum_score);

ideal_score = ones(1, n);
ideal_score(2:n) = ideal_score(2:n) ./ log2(2:n);
cum_ideal_score = cumsum(ideal_score);

score = cum_score ./ cum_ideal_score;
score = score(end);
