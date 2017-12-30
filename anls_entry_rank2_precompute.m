function H = anls_entry_rank2_precompute(left, right, H)

% left: 2*2
% right: n*2
% Returning H of size n*2 also

n = size(right, 1);

solve_either = zeros(n, 2);
solve_either(:, 1) = right(:, 1) ./ left(1,1);
solve_either(:, 2) = right(:, 2) ./ left(2,2);
cosine_either = bsxfun(@times, solve_either, [sqrt(left(1,1)), sqrt(left(2,2))]);
choose_first = (cosine_either(:, 1) >= cosine_either(:, 2));
solve_either(choose_first, 2) = 0;
solve_either(~choose_first, 1) = 0;

H = (left \ right')';
use_either = ~all(H>0, 2);
H(use_either, :) = solve_either(use_either, :);
