clear; clc; close all;
% We compare (in terms of worst-case rates):
% - rate of DCA for a given splitting
% - rate of DCA for the best splitting wrt shifting curvature lambda
% - rate of Gradient Descent (GD) with optimal stepsize
% for 2 cases:
%   - smooth nonconvex: mu_F = -L_F
%   - smooth convex:    mu_F = 0
%
%% smooth nonconvex case (mu_F = -L_F)
mu1 = 1.5; L1 = 2; mu2 = 1; L2 = 2.5; % curvatures of f1 and f2
mu_F = mu1 - L2; L_F = L1 - mu2; % curvatures of F
fprintf("Smooth nonconvex case \n\n")
%% DCA rates
%% initial rate (1/denominator*N)
lambda = 0; % 
[p_dca_init, id_p_init] = get_p_lambda(L1,mu1,L2,mu2,lambda); % initial denominator
%% best rate (1/denominator*N)
lambda_opt = get_opt_lambda(L1,mu1,L2,mu2); % 1.0991; best curvature shift
lambda = lambda_opt;
[p_dca_opt, id_p_opt] = get_p_lambda(L1,mu1,L2,mu2,lambda); % best denominator 
 L1_opt =  L1 - lambda_opt;
mu1_opt = mu1 - lambda_opt;
 L2_opt =  L2 - lambda_opt;
mu2_opt = mu2 - lambda_opt;
%% GD rates
gamma_opt_ncvx = 2/sqrt(3) * 1/L_F; % optimal stepsize for smooth nonconvex functions
gamma = gamma_opt_ncvx/L_F; N = 1;
p_GD_opt = gamma * ( 2 + gamma*L_F .* gamma*mu_F ./ ( 2-gamma*mu_F - gamma*L_F ) ) * N; % optimal denominator for smooth nonconvex functions
%%
fprintf("DCA rates: \n")
fprintf("\t\t Initial splitting: \t\t\t\t Regime p%d = %.4f: \t mu1 = %.4f, L1 = %.4f, mu2 = %.4f, L2 = %.4f \n", id_p_init, p_dca_init, mu1, L1, mu2, L2 )
fprintf("\t\t Best splitting (lambda = %.4f): \t Regime p%d = %.4f: \t mu1 = %.4f, L1 = %.4f, mu2 = %.4f, L2 = %.4f \n\n", lambda_opt, id_p_opt, p_dca_opt, mu1_opt, L1_opt, mu2_opt, L2_opt )
fprintf("GD rate: \t\t p_GD = %.4f: \t muF = %.4f, LF = %.4f \n\n\n", p_GD_opt, mu_F, L_F )


%% smooth convex case (mu_F = -L_F)
% mu1 = 1.5; L1 = 2; mu2 = 1; L2 = 1.5; % lambda_opt = 7/6;
mu1 = 1; L1 = 2; mu2 = .1; L2 = 1; % lambda_opt = .4;
mu_F = mu1 - L2; L_F = L1 - mu2; 
fprintf("\nSmooth convex case \n\n")
%% DCA rates
%% initial rate (1/denominator*N)
lambda = 0; % 
[p_dca_init, id_p_init] = get_p_lambda(L1,mu1,L2,mu2,lambda); % initial denominator
%% best rate (1/denominator*N)
lambda_opt = get_opt_lambda(L1,mu1,L2,mu2); % 1.0991; best curvature shift
lambda = lambda_opt;
[p_dca_opt, id_p_opt] = get_p_lambda(L1,mu1,L2,mu2,lambda); % best denominator 
 L1_opt =  L1 - lambda_opt;
mu1_opt = mu1 - lambda_opt;
 L2_opt =  L2 - lambda_opt;
mu2_opt = mu2 - lambda_opt;
%% GD rates
gamma_opt_cvx = 3/2 * 1/L_F; % optimal stepsize for smooth convex functions after one iteration
gamma = gamma_opt_cvx/L_F; N = 1;
p_GD_opt = 2*gamma*N; % optimal denominator for smooth nonconvex functions
%%
fprintf("DCA rates: \n")
fprintf("\t\t Initial splitting: \t\t\t\t Regime p%d = %.4f: \t mu1 = %.4f, L1 = %.4f, mu2 = %.4f, L2 = %.4f \n", id_p_init, p_dca_init, mu1, L1, mu2, L2 )
fprintf("\t\t Best splitting (lambda = %.4f): \t Regime p%d = %.4f: \t mu1 = %.4f, L1 = %.4f, mu2 = %.4f, L2 = %.4f \n\n", lambda_opt, id_p_opt, p_dca_opt, mu1_opt, L1_opt, mu2_opt, L2_opt )
fprintf("GD rate: \t\t p_GD = %.4f: \t muF = %.4f, LF = %.4f \n\n\n", p_GD_opt, mu_F, L_F )