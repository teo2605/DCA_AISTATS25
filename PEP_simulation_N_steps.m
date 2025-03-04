function [] = PEP_simulation_N_steps(L1,mu1,L2,mu2,N,Delta)
%% To exemplify the rates from Corollary 3.2 and their exactness 
%%
if nargin == 0
    %% setup
    clear; clc; close all;
    N   =  20;                 % number of iterations
    Delta = 1;                 % initial gap F0-FN
    % some examples with parameters of specific regimes
    L1 = 4; mu1 = 1.5; L2 = 3; mu2 =  1;    % regime p1 (f2 convex)
    L1 = 2; mu1 = 1; L2 = 1.5; mu2 =  -0.25;    % regime p1 (f2 weakly convex)
    % L2 = 4; mu1 = 1.5; L1 = 3; mu2 =  1;    % regime p2
    % L2 = 4; mu1 = 1.5; L1 = 3; mu2 = -1;    % regime p3
    % L1 = 2; mu1 = 1; L2 = 0.5; mu2 = -0.5;  % regime p4
    L1 = 2; mu1 = 1; L2 = -0.5; mu2 = -0.75;  % regime p4
    % L1 = 2; mu1 = 1; L2 = 0.9; mu2 =  -0.25; % regime p5
    % L1 = 2; mu1 = 1; L2 = 2.5; mu2 =  2.25; % regime p6
end
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fntsz = 16;
%% compute worst-case values with PEP
its_vec = 1 : N;
wc_vals = 0 * its_vec;
parfor k = 1 : length(its_vec)
    wc_vals(k) = PEP_DCA_no_opt(its_vec(k),L1,mu1,L2,mu2,Delta);
end
%% get theoretical bound after one iteration (Theorem 3.1)
[p_th, id_p] = compute_DCA_rate_6_regimes(L1,mu1,L2,mu2);
%% compare numerical findings with theory
figure(Color='white'); legend_names = {};
plot(its_vec, wc_vals, '--*', 'LineWidth', 1,'MarkerSize', 7); hold on; grid on;
legend_names{end+1} = "PEP result";

plot(its_vec, Delta./(p_th*its_vec), '--', 'LineWidth', 1, 'Marker','diamond','MarkerSize', 7);

legend_names{end+1} = sprintf("Theory: $\\frac{1}{p_{%d} \\cdot N}$", id_p);
leg = legend(legend_names, 'Interpreter','latex', 'FontSize', fntsz, 'Location','northeast');
leg.AutoUpdate = "off";

ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 16;
ax.GridLineWidth = 2;
xlabel('Iteration $k$','Interpreter', 'latex', 'FontSize', fntsz)
ylabel('Worst-case values wc$(k)$','Interpreter', 'latex', 'FontSize', fntsz)
title_str = sprintf('Setup: $L_1=%.2f$, $\\mu_1=%.2f$, $L_2=%.2f$, $\\mu_2=%.2f$', L1, mu1, L2, mu2);
title(title_str, 'Interpreter', 'latex', 'FontSize', fntsz)

plot_name = sprintf('PEP_Sim_Regime_p%d.pdf', id_p);
exportgraphics(gcf, plot_name);
end