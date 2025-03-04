function [wc, unbounded, result] = PEP_DCA_with_opt(N,L1,mu1,L2,mu2,Delta)
%% Default run
if nargin == 0
    clear; clc;
    N   =  1;                  % number of iterations
    Delta = 1;                 % initial gap F0-FN
    L2 = 4; L1 = 3; mu1 = 1.5; mu2 = -1;    
end
%% Initialize an empty PEP
P = pep();
%% Set up the objective function
param1.L  = L1;     % Smoothness parameter
param1.mu = mu1;    % Curvature parameter
param2.L  = L2;     % Smoothness parameter
param2.mu = mu2;    % Curvature parameter
f1 = P.DeclareFunction('Hypoconvex', param1);
f2 = P.DeclareFunction('Hypoconvex', param2);
F = f1 - f2;
%% Declare data structures
% N steps of DCA -- x_2 to x_{N+1}, starting from x_1
x_        = cell(N+1,1);
f1_       = cell(N+1,1);
f2_       = cell(N+1,1);
g1_       = cell(N+1,1);
%% tags
tag_xf1     = @(i)(sprintf('f1_x_%d', i-1));
tag_xf2     = @(i)(sprintf('f2_x_%d', i-1));
%% Define a stationary point x_s_
[~, F_s_] = F.OptimalPoint();
q = max(0, L1-mu2); % additional term in denominator due to using F* for initial condition
%% Define the starting point
x_{1}  = P.StartingPoint();		 % x0=x_{1} is some starting point
[g1_{1}, f1_{1}, x_{1}] = f1.oracle(x_{1}, tag_xf1(1));
%% Algorithm -- Iterations of DCA
for i = 1 : N
    % g1_{i+1} \in \partial f2(x_{i})
    [g1_{i+1}, f2_{i}, x_{i}]     = f2.oracle(x_{i}, tag_xf2(i));
    % Generate x_{i+1} such that it is satisfied the necessary condition:
    % \partial f1(x_{i+1}) = g1_{i+1} = \partial f2(x_{i}) (" = g2_{i}")
    x_{i+1}  = Point('Point');
    f1_{i+1} = Point('Function value');
    f1.AddComponent(x_{i+1}, g1_{i+1}, f1_{i+1}, tag_xf1(i+1)); % 
end
[g2_Np1_, f2_{N+1}, ~] = f2.oracle(x_{N+1}, tag_xf2(N+1));
%% Set up the initial condition
P.InitialCondition( (f1_{1} - f2_{1}) - F_s_ <= Delta );
%% Add Performance measure
subgrad_nrm_ = cell(N+1, 1);
for i = 1 : N
    subgrad_nrm_{i} = ( g1_{i} - g1_{i+1} ) ^2;
    %% Extra condition for F* responsible for the additional term in the denominator
    P.AddConstraint( f1_{i} - f2_{i} - F_s_ - 1/(2*q) * subgrad_nrm_{i} >= 0 ); % Optimality constraint
end
subgrad_nrm_{N+1} = ( g1_{N+1} - g2_Np1_ ) ^2;
P.AddConstraint( f1_{N+1} - f2_{N+1} - F_s_ - 1/(2*q) * subgrad_nrm_{N+1} >= 0 ); % Optimality constraint
%%
for i = 1 : N+1    
    P.PerformanceMetric( subgrad_nrm_{i} / 2 );
end
%% Solve the PEP
P.TraceHeuristic(0);
result = P.solve(0);
wc     = result.WCperformance;
%% Evaluate the output
unbounded = false;
if strcmp(result.solverDetails.info, 'Unbounded objective function (MOSEK)')
    fprintf("\n\n \t\t UNBOUNDED -- Setup: L1 = %.4f, L2 = %.4f, mu1 = %.4f, mu2 = %.4f \n\n",...
            L1, L2, mu1, mu2);
    unbounded = true;
    return;
end
%% Compare with the theoretical findings (our six regimes)
if nargin == 0    
    p_pep = (1/wc - 1/q) / N; % if sublinear rate
    [p_th, id_p] = compute_DCA_rate_6_regimes(L1,mu1,L2,mu2);
    if id_p == 0
        disp("Regime not found");
    else    
        fprintf("Regime p%d \n \t L1 = %.2f \t mu1 = %.2f \t L2 = %.2f \t mu2 = %.2f\n", ...
                    id_p, L1, mu1, L2, mu2);
        fprintf(" Coefficient of N in denominator: \n " + ...
            "\t with PEP: %.4f \n " + ...
            "\t theory: %.4f \n", p_pep, p_th);        
    end
end