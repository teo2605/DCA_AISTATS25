function f2_val = get_f2(x,Sigma,lambda)
% evaluate the function value of 
%   f2(x) = 0.5 x^T Sigma x - 0.5 lambda ||x||^2
    if nargin < 3; lambda = 0; end
    f2_val = 1/2 * x.'*Sigma*x - lambda*norm(x,2)^2/2;
end