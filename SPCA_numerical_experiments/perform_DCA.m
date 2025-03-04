function [G_nrm,F_vals,xs,num_it] = perform_DCA(x0,N,kappa,eta,Sigma,lambda)
% Input:
    % x0     -- starting point
    % N      -- number of iterations
    % kappa  -- weight of l1 norm
    % eta    -- weight of l2-regularization
    % Sigma  -- Covariance matrix; Sigma = A'*A, where A has dimension (20*n) x n
    % lambda -- curvature shift in f1 and f2
    % 
    % f1(x) = kappa*norm(x,1) + id_{unit ball} + eta/2*norm(x,2)^2 -
    %           lambda/2*norm(x,2)^2
    % f2(x) = 1/2*x'*Sigma*x - lambda/2*norm(x,2)^2
    % F(x) = f1(x) - f2(x)
%
% Output: 
%   G_nrm = \|\partial f1(x_k) - \nabla f2(x_k)\|^2, k = 0,...,N
%   F_vals = f1(x_k) - f2(x_k), k = 0,...,N
%   xs    = argmin F(x_k), k=0,...,N
%   num_it = argmin (G_nrm < tol)
%
tol = eps; % Default tol 
n = length(Sigma);
x_its = NaN(n,N+1);
x_its(:,1) = x0;
g1_vals = NaN(n,N+1); g1_vals(:,1) = get_grad_f1(x0,kappa,eta,lambda);
g2_vals = NaN(n,N+1); g2_vals(:,1) = get_grad_f2(x0,Sigma,lambda);
f1_vals = NaN(1,N+1); f1_vals(:,1) = get_f1(x0,kappa,eta,lambda);
f2_vals = NaN(1,N+1); f2_vals(:,1) = get_f2(x0,Sigma,lambda);
num_it = N+1;
for i = 1 : N
    x_its(:,i+1)   = dca_it     (x_its(:,i),kappa,eta,Sigma,lambda);
    g1_vals(:,i+1) = g2_vals(:,i);
    g2_vals(:,i+1) = get_grad_f2(x_its(:,i+1),Sigma,lambda);
    f1_vals(:,i+1) = get_f1     (x_its(:,i+1),kappa,eta,lambda);
    f2_vals(:,i+1) = get_f2     (x_its(:,i+1),Sigma,lambda);
    if norm(g1_vals(:,i+1) - g2_vals(:,i+1))^2 <= tol
        num_it = i+1;
        break;
    end
end
F_vals = f1_vals - f2_vals;
G_vals = g1_vals - g2_vals;
G_nrm = vecnorm(G_vals,2,1).^2;
xs = x_its(:,num_it);
end