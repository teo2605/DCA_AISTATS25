clc; clear; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fntsz = 16;
% Proposition D.2 -- plot the function example reaching the worst-case in
% regime p2
%
xmax = 5;
L2 = 2; mu2 = .1; L1 = .5; mu1 = .1;    % regime p2 (f2 convex; L2>L1)
% L2 = 2; mu2 = 1; L1 = 1.5; mu1 = 0.25;
p2 = compute_DCA_rate_6_regimes(L1,mu1,L2,mu2);

N = 3; DeltaF = 1;
U = -sqrt(2*DeltaF/(p2*N));

f10 = DeltaF;
f20 = 0;
g20 = 0;
g10 = U;
x0   = 0;
f1_x = @(x) ( L1 * (x - x0).^2/2 + g10 * (x - x0) + f10 );
g1_x = @(x) ( L1 * (x - x0) + g10 );

k_vec = 0 : N;
xk_vec = x0 - k_vec * U/L1;
g1k_vec = g1_x(xk_vec);
g2k_vec = g1k_vec - U;
f1k_vec = f1_x(xk_vec);
f2k_vec = f1k_vec - (N-k_vec)/N * DeltaF;
xbark_vec = xk_vec(1:end-1) - (L2-L1)/(L2-mu2) * U/L1;
%%
Fk_vec = f1k_vec - f2k_vec;
x_vec = sort(unique([-1, xk_vec, xbark_vec, linspace(-.5, xmax, 10001)]));

f2_xbars = [];
f1_xvec = f1_x(x_vec);
f2_xvec = 0 * x_vec;

g2_xbars = [];
g1_xvec = g1_x(x_vec);
g2_xvec = 0 * x_vec;

parfor id_x = 1 : length(x_vec)
    x = x_vec(id_x);
    if x <= x0
        f2_xvec(id_x) = mu2 * (x-x0)^2/2 + g2k_vec(1)*(x-x0) + f20;
        g2_xvec(id_x) = mu2 * (x-x0)     + g2k_vec(1);
    else
        if x >= xk_vec(end)
            f2_xvec(id_x) = mu2 * (x-xk_vec(end))^2/2 + g2k_vec(end)*(x-xk_vec(end)) + f2k_vec(end);
            g2_xvec(id_x) = mu2 * (x-xk_vec(end))     + g2k_vec(end);
        else
            id_k = find(x > xk_vec, 1, "last");
            if x < xbark_vec(id_k)
                f2_xvec(id_x) = mu2 * (x-xk_vec(id_k))^2/2 + g2k_vec(id_k)*(x-xk_vec(id_k)) + f2k_vec(id_k);
                g2_xvec(id_x) = mu2 * (x-xk_vec(id_k))     + g2k_vec(id_k);
            else
                f2_xvec(id_x) = L2 * (x-xk_vec(id_k+1))^2/2 + g2k_vec(id_k+1)*(x-xk_vec(id_k+1)) + f2k_vec(id_k+1);
                g2_xvec(id_x) = L2 * (x-xk_vec(id_k+1))     + g2k_vec(id_k+1);
            end
            if x == xbark_vec(id_k)
                f2_xbars = [f2_xbars, f2_xvec(id_x)];
                g2_xbars = [g2_xbars, g2_xvec(id_x)];
            end
        end
    end
end

%% plot for function values
F_xvec = f1_xvec - f2_xvec;
figure(Color='white'); legend_names = {};
plot(x_vec, f1_xvec, 'LineWidth', 2); hold on; legend_names{end+1} = '$f_1(x)$';
plot(x_vec, f2_xvec, 'LineWidth', 2); hold on; legend_names{end+1} = '$f_2(x)$';
plot(x_vec, F_xvec, 'LineWidth', 2); hold on;  legend_names{end+1} = '$F(x)$';
xlabel('$x$','Interpreter', 'latex', 'FontSize', fntsz)
ylabel('Function values','Interpreter', 'latex', 'FontSize', fntsz)
leg = legend(legend_names, 'Interpreter','latex', 'FontSize', fntsz, 'AutoUpdate','off');

plot(xk_vec, f2k_vec, 'k*', 'MarkerSize', 8);
plot(xbark_vec, f2_xbars, 'mo', 'Marker','pentagram', 'MarkerSize', 8);

ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;
ax.GridLineWidth = 1;
grid on;
xlim([-1 xmax])

title_str = sprintf('Setup: $L_1=%.2f$, $\\mu_1=%.2f$, $L_2=%.2f$, $\\mu_2=%.2f$', L1, mu1, L2, mu2);
% title(title_str, 'Interpreter', 'latex', 'FontSize', fntsz)
leg.Location = 'northwest';
plot_name = sprintf('Worst_case_function_Regime_p2_function_values.pdf');
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

plot(xk_vec, g2k_vec, 'k*', 'MarkerSize', 8);
plot(xbark_vec, g2_xbars, 'mo', 'Marker','pentagram', 'MarkerSize', 8);

ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;
ax.GridLineWidth = 1;
grid on;

ylim([-.75, 2])
xlim([-1, xmax])

%% Illustrate the DCA iterations
for k = 1 : N
    plot([xk_vec(k) xk_vec(k)], [g1k_vec(k), g2k_vec(k)], 'k--', 'LineWidth', 1);
    plot([xk_vec(k) xk_vec(k+1)], [g2k_vec(k), g2k_vec(k)], '--', 'Color', "#A2142F", 'LineWidth', 1);
end
%%
plot_name = sprintf('Worst_case_function_Regime_p2_gradients_DCA_steps.pdf');
exportgraphics(gcf, plot_name);
disp(title_str)