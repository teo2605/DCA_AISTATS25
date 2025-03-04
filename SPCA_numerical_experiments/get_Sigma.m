function Sigma = get_Sigma(n,sparse_density)
% Script to generate the sampple covariance matrix Sigma = A^T*A,
% where A is sparse matrix with randomly normal distributed entries
    A = sprandn(20*n,n,sparse_density);    
    Sigma_ = A.'*A;
    Sigma = Sigma_ / norm(full(Sigma_));    
end