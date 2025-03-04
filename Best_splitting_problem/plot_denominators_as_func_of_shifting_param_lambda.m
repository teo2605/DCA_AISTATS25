clear; clc; close all;
% This script plots the values of possible denominators as function of
% curvature shifting parameter lambda. 
% In all these cases, the four curvature parameters are assumed available.

mu1 = 1.5; L1 = 3; mu2 = -1; L2 = 4; % best is p3!
mu1 = .2; L1 = 3; mu2 = .1; L2 = 1000; % best is p3!!!
mu1 = .2; L1 = 1000; mu2 = .1; L2 = 3; % best is p1 with mu2 < 0 !!!
mu1 = .2; L1 = 100; mu2 = .1; L2 = 100; % best is p3 when mu1 > mu2 > 0!!!
mu1 = .02; L1 = 200; mu2 = .01; L2 = 100; % best is p3!!!
mu1 = 0; L1 = 1; mu2 = 0; L2 = 1; % both convex => already best splitting
mu1 = 0.1; L1 = 3; L2 = 4; mu2 = -0.1; % not converging unless lambda<0 shift
mu1 = 0.1; L1 = 3; L2 = 4; mu2 = 0.1; % start from both strongly convex with same mu => best is to shift both to convex
mu1 = 0.2; L1 = 3; L2 = 4; mu2 = 0.1; % start from both strongly convex with mu1>mu2 => best is to make mu2 weakly convex
mu1 = 0.1; L1 = 3; L2 = 4; mu2 = 0.2; % start from both strongly convex with mu1<mu2 => best is to make mu1 convex
mu1 = 0.1; L1 = 400; L2 = 300; mu2 = 0.2; % start from both strongly convex with mu1<mu2 => best is to make mu1 convex
mu1 = 1e-3; L1 = 4; L2 = 3; mu2 = 2e-3; % start from both strongly convex with mu1<mu2 => best is to make mu1 convex
mu1 = 0.1; L1 = 4; L2 = 3; mu2 = 0e-3; % start from both strongly convex with mu1<mu2 => best is to make mu1 convex
mu1 = 2; mu2 = -1.75; L1 = 4; L2 = 3; % start in p4, best in p3
mu1 = 2.99; mu2 = -2.9; L1 = 4; L2 = 3; % start from p4, best in p3
mu1 = 0.1; L1 = 3; L2 = 4; mu2 = -0.2; % not converging unless shifted with some lambda<0

tol = 1e-6;
lambda_max = mu1+min(0,(mu2-mu1)/2) - tol;

lambda_vec = sort([0, linspace(-1 , lambda_max, 1+1e5)]);

if mu1 + mu2 <= 0 && (mu1 ~= 0 && mu2 ~= 0)
    lambda_max = (mu1+mu2)/2 - tol;
    lambda_vec = sort([0, linspace(-1 , lambda_max, 1+1e5)]);
end

p_lambda = 0 * lambda_vec;
id_p = 0 * lambda_vec;
for ll = 1 : length(lambda_vec)
    lambda = lambda_vec(ll);
    [p_lambda(ll), id_p(ll)] = get_p_lambda(L1,mu1,L2,mu2,lambda);
end

[p_max, ll_best] = max(p_lambda);
lambda_opt = lambda_vec(ll_best);

%% figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure(1); clf;
subplot(2,1,1)
plot(lambda_vec, p_lambda, 'LineWidth', 3); grid on; hold on;
plot(lambda_opt, p_lambda(ll_best),'r*', 'MarkerSize', 10)
xlim([min(lambda_vec), max(max(lambda_vec), 0)])
xlabel('$\lambda$')
ylabel('Value $p(\lambda)$')

ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 12;
ax.GridLineWidth = 2;


subplot(2,1,2)
plot(lambda_vec, id_p, '*', 'LineWidth', 1, 'MarkerSize', 5); grid on; hold on;
plot(lambda_opt, id_p(ll_best),'r*', 'MarkerSize', 10)
xlim([min(lambda_vec), max(max(lambda_vec), 0)])
ylim([0 6])
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;
ax.GridLineWidth = 2;

xlabel('$\lambda$')
ylabel('Regime id $p_i(\lambda)$')

title_str = sprintf('Best curvature shift $-\\lambda \\|\\cdot\\|^2/2$ for \n setup: $\\mu_1 = %.3f$, $L_1 = %.2f$, $\\mu_2 = %.3f$, $L_2 = %.2f$', ...
                    mu1, L1, mu2, L2);
S = sgtitle(title_str,'Interpreter','latex');
set(gcf,'Position', [680   248   560   750])
%% display messages
mu_F = mu1 - L2; L_F = L1 - mu2;
fprintf("Objective F: mu_F = %.2f, L_F = %.2f \n\n", mu_F, L_F)
p_init = 0;
if (mu1 == 0 && mu2 == 0) || (mu1+mu2 > 0)
    p_init = p_lambda(id_p(lambda_vec == 0));
end

% shifted cvx
 L1_opt =  L1 - lambda_opt;
mu1_opt = mu1 - lambda_opt;
 L2_opt =  L2 - lambda_opt;
mu2_opt = mu2 - lambda_opt;

fprintf("Initial splitting: \t Regime p%d: \t mu1 = %.4f, L1 = %.4f, mu2 = %.4f, L2 = %.4f \n\n", id_p(lambda_vec == 0), mu1, L1, mu2, L2 )
fprintf("Initial: Regime p%d = %.4f for lambda = %.4f \n\n", id_p(lambda_vec == 0), p_init, 0) ; % Regime p0 means no convergence as p0=0!

fprintf("Best splitting: \t Regime p%d: \t mu1 = %.4f, L1 = %.4f, mu2 = %.4f, L2 = %.4f \n\n",  id_p(ll_best), mu1_opt, L1_opt, mu2_opt, L2_opt )
fprintf("Best: Regime p%d = %.4f for lambda = %.4f \n\n", id_p(ll_best), p_max, lambda_vec(ll_best))

fprintf("Improvement: \t %.4f  \n\n", abs(p_init - p_max)/p_init*100 )

fig_name = sprintf("Best_splitting_L1=%.4f_mu1=%.4f_L2=%.4f_mu2=%.4f_regime_p%d.pdf", L1, mu1, L2, mu2, id_p(ll_best));
sppi = get(groot,"ScreenPixelsPerInch");
exportgraphics(gcf,fig_name,"Resolution",sppi)