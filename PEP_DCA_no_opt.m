function [wc, unbounded, result] = PEP_DCA_no_opt(N,L1,mu1,L2,mu2,Delta)
%% Used to verify Theorem 3.1
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
% F=f1-f2 is the objective function (not explicitly defined)
f1 = P.DeclareFunction('Hypoconvex', param1);
f2 = P.DeclareFunction('Hypoconvex', param2);
%% Declare data structures
% N steps of DCA -- x_2 to x_{N+1}, starting from x_1
x_        = cell(N+1,1);
f1_       = cell(N+1,1);
f2_       = cell(N+1,1);
g1_       = cell(N+1,1);
%% tags
tag_xf1     = @(i)(sprintf('f1_x_%d', i-1));
tag_xf2     = @(i)(sprintf('f2_x_%d', i-1));
%% Define the starting point
x_{1}  = P.StartingPoint();		 % x0=x_{1} is some starting point
[g1_{1}, f1_{1}, x_{1}] = f1.oracle(x_{1}, tag_xf1(1));
%% (5) Algorithm -- Iterations of DCA
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
P.InitialCondition( (f1_{1} - f2_{1}) - (f1_{N+1}-f2_{N+1}) <= Delta );
%% Add Performance measure
subgrad_nrm_ = cell(N+1, 1);
for i = 1 : N
    subgrad_nrm_{i} = ( g1_{i} - g1_{i+1} ) ^2;
end
subgrad_nrm_{N+1} = ( g1_{N+1} - g2_Np1_ ) ^2;
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
    p_pep = 1 / (wc * N); % if sublinear rate
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
    %% Display multipliers of interpolation inequalities
    fprintf("Multipliers of interpolation inequalities:")
    [result.dualnames', num2cell(result.dualvalues)']
    %% Coefficients from Table 1 (pi = sigma + sigma^+;) + multipliers alpha in the proofs     
    %% p1
    if nargin == 0 && id_p == 1
        sigma = 1/(L2) * (L2-mu1)/(L1-mu1);
        sigma_plus = 1/(L2) * (1 + (1/L2 - 1/L1) / (1/mu1 - 1/L1));
        alpha_th      = (mu1 / L2) * (L1 - L2) / (L1 - mu1);
        alpha_proof   = alpha_th / (alpha_th+1);
        [alpha_th; alpha_proof; sigma; sigma_plus]'
    end
    %% p2
    if nargin == 0 && id_p == 2
        sigma_plus = 1/(L1) * (L1-mu2)/(L2-mu2);
        sigma = 1/(L1) * (1 + (1/L1 - 1/L2) / (1/mu2 - 1/L2));
        alpha_th      = (mu2 / L1) * (L2 - L1) / (L2 - mu2);
        [alpha_th; sigma; sigma_plus]'
    end  
    %% p3
    if nargin == 0 && id_p == 3
        sigma = 1/L1 * (1 + 1/L1 / (1/mu1 + 1/mu2 + 1/L2 - 1/L1));
        sigma_plus = 1/(L2 + mu2);        
        alpha_th      = -mu2 / (L2 + mu2);
        [alpha_th; sigma; sigma_plus]'
    end
    %% p4
    if nargin == 0 && id_p == 4
        sigma = 0;
        sigma_plus = (mu1+mu2) / (-mu2)^2;
        alpha_th      = (mu1+mu2) / (-mu2);
        [alpha_th; sigma; sigma_plus]'
    end
    %% p5
    if nargin == 0 && id_p == 5
        sigma = (mu1+L2) / (L2)^2;
        sigma_plus = 0;        
        alpha_th      = (mu1) / (L2);
        [alpha_th; sigma; sigma_plus]'
    end
    %% p6
    if nargin == 0 && id_p == 6
        sigma = (mu2+L1) / (L1)^2;
        sigma_plus = 0;        
        alpha_th      = (mu2) / (L1);
        [alpha_th; sigma; sigma_plus]'
    end
end