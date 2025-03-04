function [lambda_opt, p_opt, id_p_opt] = get_opt_lambda(L1,mu1,L2,mu2)
% Compute optimal curvature shift, given all four curvature parameters
    if nargin == 0   
        mu1 = 1; L1 = 2; mu2 = .5; L2 = 1.5;
    end    
    tol = 1e-12; % accuracy on span of lambda values
    lambda_max = mu1 + min(0,(mu2-mu1)/2) - tol;
    lambda_vec = sort([0, linspace(-1 , lambda_max, 1+1e5)]);
    if mu1 + mu2 <= 0 && (mu1 ~= 0 && mu2 ~= 0)
        lambda_max = (mu1+mu2)/2 - tol;
        lambda_vec = sort([0, linspace(-1 , lambda_max, 1+1e5)]);
    end
    p_lambda = NaN(size(lambda_vec));
    id_p = NaN(size(lambda_vec));
    for ll = 1 : length(lambda_vec)
        lambda = lambda_vec(ll);     
        [p_lambda(ll), id_p(ll)] = get_p_lambda(L1,mu1,L2,mu2,lambda);
    end    
    [p_opt, ll_best] = max(p_lambda);
    lambda_opt = lambda_vec(ll_best);
    id_p_opt = id_p(ll_best);
end