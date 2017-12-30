function [X, term_subset, doc_subset] = tfidf(X, min_term_count, min_doc_count)

[m, n] = size(X);
term_subset = 1 : m;
doc_subset = 1 : n;
while true
	subset = (sum(X, 2) >= min_term_count & sum(X~=0, 2) < size(X, 2));
	X = X(subset, :);
	term_subset = term_subset(subset);
	pattern = (X ~= 0);
	subset = find(sum(pattern) >= min_doc_count);
	if length(subset) == size(X, 2)
		break;
	else
		X = X(:, subset);
		doc_subset = doc_subset(subset);
	end
end
df = full(sum(pattern, 2));
n = size(X, 2);
idf = log(n./df);

[idx, jdx, vals] = find(X);
X = sparse(idx, jdx, log(vals) + 1);

X = bsxfun(@times, X, idf);

D = full(1./sqrt(sum(X.^2)));
X = bsxfun(@times, X, D);
