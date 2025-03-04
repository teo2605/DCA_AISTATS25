function df1 = get_grad_f1(x,kappa,eta,lambda)
% evaluate a subgradient of
%   f1(x) = kappa ||x||_1 + 0.5 (eta-lambda) ||x||_2^2 + id(B(0,1))
    if nargin < 4; lambda = 0; end
    if norm(x,2) > 1 + eps
        % compute the normal cone to unit ball, i.e., a subgradient of
        % indicator of unit ball
        error("Empty set for ||x||>1!");
    end
    df1 = kappa*sign(x) + (eta-lambda)*x; % a subgradient of f1
end