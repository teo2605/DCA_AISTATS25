clc; clear; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fntsz = 16;
% Proposition D.1 -- plot the function example reaching the worst-case in
% regime p1
%
xmax = 4;

L1 = 2; mu1 = .1; L2 = .5; mu2 = -.01;    % regime p1 (f2 may be weakly convex)
p1 = compute_DCA_rate_6_regimes(L1,mu1,L2,mu2);

N = 3; DeltaF = 1;
U = -sqrt(2*DeltaF/(p1*N));

f1_0 = DeltaF;
f2_0 = 0;
g2_0 = U;
x0   = 0;
f2_x = @(x) ( L2 * (x - x0).^2/2 + g2_0 * (x - x0) + f2_0 );
g2_x = @(x) ( L2 * (x - x0) + g2_0 );

k_vec = 0 : N;
xk_vec = x0 - k_vec * U/L2;
g1k_vec = g2_0 - (k_vec-1)*U;
g2k_vec = g1k_vec - U;
f2k_vec = f2_x(xk_vec);
f1k_vec = f2k_vec + (N-k_vec)/N * DeltaF;
xbark_vec = xk_vec(1:end-1) - (L2-mu1)/(L1-mu1) * U/L2;

Fk_vec = f1k_vec - f2k_vec;
x_vec = sort(unique([-1, xk_vec, xbark_vec, linspace(-.5, xmax, 10001)]));

f1_xbars = [];
f2_xvec = f2_x(x_vec);
f1_xvec = 0 * x_vec;

g1_xbars = [];
g2_xvec = g2_x(x_vec);
g1_xvec = 0 * x_vec;

parfor id_x = 1 : length(x_vec)
    x = x_vec(id_x);
    if x <= x0
        f1_xvec(id_x) = L1 * (x-x0)^2/2 + g1k_vec(1)*(x-x0) + f1_0;
        g1_xvec(id_x) = L1 * (x-x0) + g1k_vec(1);
    else
        if x >= xk_vec(end)
            f1_xvec(id_x) = L1 * (x-xk_vec(end))^2/2 + g1k_vec(end)*(x-xk_vec(end)) + f1k_vec(end);
            g1_xvec(id_x) = L1 * (x-xk_vec(end)) + g1k_vec(end);
        else
            id_k = find(x > xk_vec, 1, "last");
            if x < xbark_vec(id_k)
                f1_xvec(id_x) = L1 * (x-xk_vec(id_k))^2/2 + g1k_vec(id_k)*(x-xk_vec(id_k)) + f1k_vec(id_k);
                g1_xvec(id_x) = L1 * (x-xk_vec(id_k)) + g1k_vec(id_k);
            else
                f1_xvec(id_x) = mu1 * (x-xk_vec(id_k+1))^2/2 + g1k_vec(id_k+1)*(x-xk_vec(id_k+1)) + f1k_vec(id_k+1);
                g1_xvec(id_x) = mu1 * (x-xk_vec(id_k+1)) + g1k_vec(id_k+1);
            end
            if x == xbark_vec(id_k)
                f1_xbars = [f1_xbars, f1_xvec(id_x)];
                g1_xbars = [g1_xbars, g1_xvec(id_x)];
            end
        end
        
    end
end

%% plot for function values
F_xvec = f1_xvec - f2_xvec;
figure(Color='white'); legend_names = {};
plot(x_vec, f1_xvec, 'LineWidth', 2); hold on; legend_names{end+1} = '$f_1(x)$';
plot(x_vec, f2_xvec, 'LineWidth', 2); hold on; legend_names{end+1} = '$f_2(x)$';
plot(x_vec, f1_xvec - f2_xvec, 'LineWidth', 2); hold on;  legend_names{end+1} = '$F(x)$';
xlabel('$x$','Interpreter', 'latex', 'FontSize', fntsz)
ylabel('Function values','Interpreter', 'latex', 'FontSize', fntsz)
leg = legend(legend_names, 'Interpreter','latex', 'FontSize', fntsz, 'AutoUpdate','off');

plot(xk_vec, f1k_vec, 'k*', 'MarkerSize', 8);
plot(xbark_vec, f1_xbars, 'mo', 'Marker','pentagram', 'MarkerSize', 8);

ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;
ax.GridLineWidth = 1;
grid on;
xlim([-0.5 xmax])

title_str = sprintf('Setup: $L_1=%.2f$, $\\mu_1=%.2f$, $L_2=%.2f$, $\\mu_2=%.2f$', L1, mu1, L2, mu2);
leg.Location = 'northwest';
plot_name = sprintf('Worst_case_function_Regime_p1_function_values.pdf');
exportgraphics(gcf, plot_name);
%% plot for gradients
G_xvec = g1_xvec - g2_xvec;
figure(Color='white'); clf; legend_names = {};
plot(x_vec, g1_xvec, 'LineWidth', 2); hold on; legend_names{end+1} = '$\nabla f_1(x)$';
plot(x_vec, g2_xvec, 'LineWidth', 2); hold on; legend_names{end+1} = '$\nabla f_2(x)$';
plot(x_vec, G_xvec, 'LineWidth', 2); hold on;  legend_names{end+1} = '$\nabla F(x)$';
xlabel('$x$','Interpreter', 'latex', 'FontSize', fntsz)
ylabel('Gradient values','Interpreter', 'latex', 'FontSize', fntsz)
legend(legend_names, 'Interpreter','latex', 'FontSize', fntsz, 'Location','northwest','AutoUpdate','off');

plot(xk_vec, g1k_vec, 'k*', 'MarkerSize', 8);
plot(xbark_vec, g1_xbars, 'mo', 'Marker','pentagram', 'MarkerSize', 8);



ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;
ax.GridLineWidth = 1;
grid on;

ylim([-1.5, 1.5])
xlim([-0.5 xmax])

%% Illustrate the DCA iterations
for k = 1 : N
    plot([xk_vec(k) xk_vec(k)], [g1k_vec(k), g2k_vec(k)], 'k--', 'LineWidth', 1);
    plot([xk_vec(k) xk_vec(k+1)], [g2k_vec(k), g2k_vec(k)], '--', 'Color', "#A2142F", 'LineWidth', 1);
end
%%
plot_name = sprintf('Worst_case_function_Regime_p1_gradients_DCA_steps.pdf');
exportgraphics(gcf, plot_name);
disp(title_str)