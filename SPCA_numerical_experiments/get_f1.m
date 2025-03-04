function f1_val = get_f1(x,kappa,eta,lambda)
% evaluate the function value of 
%   f1(x) = kappa ||x||_1 + 0.5 (eta-lambda) ||x||_2^2 + id(B(0,1)) 
    if nargin < 4; lambda = 0; end
    if norm(x,2) > 1
        f1_val = Inf;
    else
        f1_val = kappa*norm(x,1) + (eta-lambda)*norm(x,2)^2/2;
    end
end