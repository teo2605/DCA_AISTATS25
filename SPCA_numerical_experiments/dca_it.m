function xp = dca_it(x,kappa,eta,Sigma,lambda)
    if nargin <= 4; lambda = 0; end
    xp = get_grad_f1_conj( get_grad_f2(x,Sigma,lambda) ,kappa,eta,lambda);
end