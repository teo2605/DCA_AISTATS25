function L_NoLips_takeI
% In this example, we use a Bregman gradient method for
% solving the constrained smooth strongly convex minimization problem
%   min_x { F(x) = f_1(x) + f_2(x) }
%   for notational convenience we denote xs=argmin_x F(x);
% where f_1(x) is h-smooth and convex and where f_2(x) is a closed convex
% indicator function.
%
% We show how to compute the worst-case value of F(xN)-F(xs) when xN is
% obtained by doing N steps of the method starting with an initial
% iterate satisfying Dh(x*,x0)<=1. (Dh is the Bregman distance generated by
% h, between x* and x0)
%
% [1] Heinz H. Bauschke, Jérôme Bolte, and Marc Teboulle. "A Descent Lemma
%     Beyond Lipschitz Gradient Continuity: First-Order Methods Revisited
%     and Applications." (2017)
%
% [2] Radu-Alexandru Dragomir, Adrien B. Taylor, Alexandre d’Aspremont, and
%     Jérôme Bolte. "Optimal Complexity and Certification of Bregman
%     First-Order Methods". (2019)
%
% DISCLAIMER: This example requires some experience with PESTO and PEPs
% (see Section 4 in [2]).

% (0) Initialize an empty PEP
P = pep();

L = 1; % d = Lh - f1 is convex
d  = P.DeclareFunction('Convex');
f1 = P.DeclareFunction('Convex');
h  = (d + f1)/L;
f2 = P.DeclareFunction('ConvexIndicator');

F  = f1 + f2;
% (2) Set up the starting point and initial condition
x0        = P.StartingPoint();        % x0 is some starting point
[ xs, Fs] = F.OptimalPoint('opt');    % xs is an optimal point, and fs=F(xs)
[gfs, fs] = f1.oracle('opt');

[gh0, h0] = h.oracle(x0,'x0');
[gf0, f0] = f1.oracle('x0');
[ghs, hs] = h.oracle(xs,'opt');

P.InitialCondition( hs - h0 - gh0 * (xs - x0) <= 1);    % Initial condition Dh(x*,x0)<=1

% (3) Algorithm
gamma = 1/L/2;        % stepsize
N     = 3;          % number of iterations

x = x0;
gfx = cell(N+1,1); gfx{1} = gf0;
ffx = cell(N+1,1); ffx{1} = f0;
ghx = cell(N+1,1); ghx{1} = gh0;

for i = 1:N
    name = sprintf('x%d',i);
    [x, ghx{i+1}, hx] = mirror(gfx{i}, ghx{i}, h+f2, gamma, name);
    [gfx{i+1}, ffx{i+1}] = F.oracle(name);
end

% (4) Set up the performance measure
P.PerformanceMetric(ffx{N+1}-fs);

% (5) Solve the PEP
P.solve()

% (6) Evaluate the output
double(ffx{N+1}-fs)   % worst-case objective function accuracy

% Result should match 1/gamma/N (see [2])
end