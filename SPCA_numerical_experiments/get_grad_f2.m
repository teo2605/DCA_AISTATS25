function df2 = get_grad_f2(x, Sigma, lambda)
% evaluate the gradient of
%   f2(x) = 0.5 x^T Sigma x - 0.5 lambda ||x||^2
    if nargin < 3; lambda = 0; end
    df2 = Sigma * x - lambda * x;
end