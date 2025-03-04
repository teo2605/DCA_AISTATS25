function [p_lambda, id_p] = get_p_lambda(L1,mu1,L2,mu2,lambda)
% Compute the corresponding denominator pi when the functions are shifted
% by -lambda x^2/2

    % shifted curvatures
     L1_shif =  L1 - lambda;
    mu1_shif = mu1 - lambda;
     L2_shif =  L2 - lambda;
    mu2_shif = mu2 - lambda;

    [p_lambda, id_p] = compute_DCA_rate_6_regimes(L1_shif, mu1_shif, L2_shif, mu2_shif);
end