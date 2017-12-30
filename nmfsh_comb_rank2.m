function [W, H, iter, grad] = nmfsh_comb_rank2(A, k, Winit, Hinit, vec_norm, normW, conv, tol, maxiter, anls_alg)

[m, n] = size(A);
W = Winit;
H = Hinit;

switch conv
case 'projgrad'
left = H * H';
right = A * H';
for iter = 1 : maxiter
	if left(1,1) == 0 | left(2,2) == 0
		grad = -1;
		disp('H is rank-deficient!');
		return;
	end
        W = anls_alg(left, right, W);
	left = W' * W;
	right = A' * W;
	if left(1,1) == 0 | left(2,2) == 0
		grad = -1;
		disp('W is rank-deficient!');
		return;
	end
        H = anls_alg(left, right, H')';
        gradH = left * H - right';
	left = H * H';
	right = A * H';
	gradW = W * left - right;
	if iter == 1
	        initgrad = sqrt(norm(gradW(gradW<=0|W>0))^2 + norm(gradH(gradH<=0|H>0))^2);
		continue;
	else
        	projnorm = sqrt(norm(gradW(gradW<=0|W>0))^2 + norm(gradH(gradH<=0|H>0))^2);
	end
	if projnorm < tol * initgrad
		break;
	end
end
grad = projnorm / initgrad;
end

if vec_norm ~= 0
	if normW
        	norms = sum(W.^vec_norm) .^ (1/vec_norm);
	        W = bsxfun(@rdivide, W, norms);
        	H = bsxfun(@times, H, norms');
	else    
        	norms = sum(H.^vec_norm, 2) .^ (1/vec_norm);
	        W = bsxfun(@times, W, norms');
        	H = bsxfun(@rdivide, H, norms);
	end
end
