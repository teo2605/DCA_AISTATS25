clear; clc; close all;
%% parameters for the plots
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
load('rng_1.mat'); rng(rngstate);
fntsz = 18;
color_vec = [230 25 75 % red
             60 180 75 % green
             240 50 230 % magenta
             0 130 200 % blue
             255 225 25 % yellow             
             128 0 0 % maroon
             245 130 48 % orange             
             ] ...
             / 256; % maroon
%% Setup of parameters
N = 3000; % # iterations
n = 200; % size of A
num_d = 1000; % #initializations inside the unit ball
eps_vec = 10.^(-1:-1:-12); % error vector for N_epsilon
sparse_density = 0.1; % 10%
%% hyperparameters of function f1
kappa = 0.02;
eta = 0.2; % 10*kappa
%% covariance matrix for function f2
load("Sigma_vals.mat"); % this loads the Sigma matrix used in the experiments; use get_Sigma.m to generate another Sigma
Sigma = Sigma_200;
%% generate num_d initial points within the unit ball
x0_unit_ball = generate_x0_unit_ball(n,num_d);
%% Theoretical rate and Optimal curvature shift 
mu2 = min(eig(Sigma)); L2 = max(eig(Sigma)); 
mu1 = eta; L1 = Inf;
[p_th, id] = compute_DCA_rate_6_regimes(L1,mu1,L2,mu2);
[lambda_opt, p_lambda_opt, id_p_opt] = get_opt_lambda(L1,mu1,L2,mu2); % theoretical optimal splitting
[p_lambda_0, id_p_0] = get_p_lambda(L1,mu1,L2,mu2,0); % initial problem
fprintf("\t Setup: mu1=%.2f, L1=%.2f,  mu2=%.4f, L2=%.4f; p%d = %f\n\n", mu1, L1, mu2, L2, id, p_th);
%% Generate vector of splittings lambda
lambda_max = min([mu1, (mu1+mu2)/2]);
lambda_vec = [0, lambda_opt, -lambda_opt/2, -lambda_opt, lambda_opt/2];
%% Get F_vals at some iterate x
F_vals = @(x) ( get_f1(x,kappa,eta,0) - get_f2(x,Sigma,0) );
%% Generate cells for num_d initializations
sim_all = cell(num_d,1);
%% Simulation for various initializations x0
parfor dd = 1 : num_d
    x0 = x0_unit_ball(:,dd);
    xs_vec = NaN(length(lambda_vec),n);
    G_nrm_all = NaN(length(lambda_vec),N+1);
    F_vals_all = NaN(length(lambda_vec),N+1);
    good_sol = zeros(length(lambda_vec),1);
    num_it_all = zeros(length(lambda_vec),1);
    tols_its_all = zeros(length(lambda_vec), length(eps_vec));
    for ll = 1 : length(lambda_vec)
        [G_nrm_all(ll,:), F_vals_all(ll,:), xs_vec(ll,:), num_it_all(ll)] = perform_DCA(x0,N,kappa,eta,Sigma,lambda_vec(ll));
        if nnz(xs_vec(ll,:)) <= 2*sparse_density*n && nnz(xs_vec(ll,:)) >= 0.25*sparse_density*n % select solutions of desired sparisty and nontrivial
            good_sol(ll) = 1;
        end
    end
    nnz_xs = zeros(1,length(lambda_vec));
    for ll = 1 : length(lambda_vec)
        nnz_xs(ll) = nnz(xs_vec(ll,:));
    end
    sim_all{dd}.Fs_all = min(F_vals_all,[],2);
    if ( sum(good_sol) < length(lambda_vec) ) || norm(nnz_xs - mean(nnz_xs), "inf") > 1e-8 || norm( sim_all{dd}.Fs_all - mean(sim_all{dd}.Fs_all) ) > 1e-8
        sim_all{dd}.flag_bad_xs = true; % remove runs if not the same solutions obtain for all lambdas
    else
        sim_all{dd}.flag_bad_xs = false;
        %% get number of iterates to reach certain accuracy of ||G||^2
        for ll = 1 : length(lambda_vec)
            for tt = 1 : length(eps_vec)               
               N_tol = find( get_min_along_iterates(G_nrm_all(ll,:)) <= eps_vec(tt),1); 
               if isempty(N_tol); N_tol = N; end               
               tols_its_all(ll,tt) = N_tol;               
            end
        end
    end
    % save data
    sim_all{dd}.num_it_all = num_it_all;
    sim_all{dd}.F_vals_all = F_vals_all;
    sim_all{dd}.Fs_all = min(F_vals_all,[],2);
    sim_all{dd}.G_nrm_all = G_nrm_all;
    sim_all{dd}.xs_vec = xs_vec;
    sim_all{dd}.good_sol = good_sol;
    sim_all{dd}.x0 = x0;
    sim_all{dd}.tols_its_all = tols_its_all;
    sim_all{dd}.nnz_xs = nnz_xs;
end
%% process data
%% remove simulations not converging to the same optimum point
sim_all_good = {};
for dd = 1 : num_d
    if ~sim_all{dd}.flag_bad_xs
        sim_all_good{end+1} = sim_all{dd}; 
    end
end
num_d_good = length(sim_all_good);
%% get average number of iterations for each lambda
N_eps_d = 0 * sim_all_good{1}.tols_its_all;
parfor dd = 1 : num_d_good
    N_eps_d = N_eps_d + sim_all_good{dd}.tols_its_all;
end
N_eps = N_eps_d / num_d_good;
x0_all_vals = NaN(n, num_d_good);
parfor dd = 1 : num_d_good
    x0_all_vals(:,dd) = sim_all_good{dd}.x0; 
end
%% set the legends
leg_default = cell(1,length(lambda_vec));
for ll = 1 : length(lambda_vec)
    switch lambda_vec(ll)
        case lambda_opt
            leg_default{ll} = sprintf('$\\lambda^*=%.1f$',lambda_opt);
        case 0
            leg_default{ll} = sprintf('$\\lambda=0$');
        otherwise
            leg_default{ll} = sprintf('$\\lambda=%.1f$',lambda_vec(ll));
    end    
end
leg_th = [sprintf('wc rate: $\\lambda=0$'), sprintf('wc rate: $\\lambda^*=%.1f$',lambda_opt), leg_default]; % including th. plots
%%

%% plots for some initial condition x0
id = randi(num_d_good,1); 
figure(1); clf;
semilogy(1:N, ( F_vals(sim_all_good{id}.x0) - sim_all_good{id}.Fs_all(1) )./(p_lambda_0*(1:N)), 'LineStyle','--','linewidth',3,'Color',color_vec(1,:)); hold on;
semilogy(1:N, ( F_vals(sim_all_good{id}.x0) - sim_all_good{id}.Fs_all(2) )./(p_lambda_opt*(1:N)), 'LineStyle','--','linewidth',3,'Color',color_vec(2,:)); hold on;
for ll = 1 : length(lambda_vec)
    lambda = lambda_vec(ll);
    semilogy(0:N, get_min_along_iterates(sim_all_good{id}.G_nrm_all(ll,:)/2), '-', 'LineWidth', 3,'Color',color_vec(ll,:)); grid on; hold on;    
end
legend(leg_th, 'Interpreter', 'latex', 'FontSize', fntsz, 'Location','best')
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 16;
ax.GridLineWidth = 1;
xlabel('Iteration $N$','Interpreter', 'latex', 'FontSize', fntsz)
ylabel('$\frac{1}{2}\min_{0 \leq k \leq N} \{\|g_1^k-\nabla f_2(x^k)\|^2\}$','Interpreter', 'latex', 'FontSize', fntsz)
title_str = sprintf('Setup: $\\eta=%.1f$, $\\kappa=%.2f$, $n=%d$; $\\mu_2=%.4f$, $L_2=%.0f$', eta, kappa, n, mu2, L2);
title(title_str, 'Interpreter', 'latex', 'FontSize', fntsz)
ylim([eps 1])
xlim([0 1000])
%% figure of ||G_N|| (not the min) for a given sample = id
figure(2); clf;
for ll = 1 : length(lambda_vec)
    semilogy(0:N,sim_all_good{id}.G_nrm_all(ll,:),'LineWidth',3,'color',color_vec(ll,:)); grid on; hold on;    
end
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 16;
ax.GridLineWidth = 1;
xlabel('Iteration $N$','Interpreter', 'latex', 'FontSize', fntsz)
ylabel('$\frac{1}{2}\|g_1^N-\nabla f_2(x^N)\|^2$','Interpreter', 'latex', 'FontSize', fntsz)
title(title_str, 'Interpreter', 'latex', 'FontSize', fntsz);
legend(leg_default, 'Interpreter', 'latex', 'FontSize', fntsz, 'Location','best')
ylim([eps 1])
xlim([0 1000])
%% figure of N_tols vs. tols
figure(3); clf; leg_names_N_tols = {};
semilogy(1:N, ( F_vals(sim_all_good{id}.x0) - sim_all_good{id}.Fs_all(1) )./(p_lambda_0*(1:N)), 'LineStyle','--','linewidth',3,'Color',color_vec(1,:)); hold on;
semilogy(1:N, ( F_vals(sim_all_good{id}.x0) - sim_all_good{id}.Fs_all(2) )./(p_lambda_opt*(1:N)), 'LineStyle','--','linewidth',3,'Color',color_vec(2,:)); hold on;
for ll = 1 : length(lambda_vec)
    semilogy(N_eps(ll,:),eps_vec,'-*','LineWidth',3,'color',color_vec(ll,:)); grid on; hold on;
end
xlabel('Iteration $N$','Interpreter', 'latex', 'FontSize', fntsz)
ylabel('$\min_{0 \leq k \leq N} \{\|g_1^k - \nabla f_2(x^k)\|^2\}$','Interpreter', 'latex', 'FontSize', fntsz)
title(title_str, 'Interpreter', 'latex', 'FontSize', fntsz);
legend(leg_th, 'Interpreter', 'latex', 'FontSize', fntsz, 'Location','east')
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 16;
ax.GridLineWidth = .75;
ylim([1e-12 1])
xlim([0 1000])
% exportgraphics(ax, 'Sim_eta=0.2_kappa=0.02_several_splittings_avg_num_its.pdf');
%% Generate table of N_epsilon
fprintf('\t %d / %d successful runs \n\n', num_d_good, num_d);
T = array2table(N_eps);
row_names = cell(length(lambda_vec),1);
% Rows: {'lambda=0','lambda_*','-lambda_*','lambda_*/2'}
for ll = 1 : length(lambda_vec);  row_names{ll} = sprintf('lambda = %.1f',lambda_vec(ll)); end
T.Properties.RowNames = row_names;
for ee = 1 : length(eps_vec); T.Properties.VariableNames{ee} = sprintf('eps=1e%d',log10(eps_vec(ee))); end
T.Properties.Description = 'Case 2';
disp(T)