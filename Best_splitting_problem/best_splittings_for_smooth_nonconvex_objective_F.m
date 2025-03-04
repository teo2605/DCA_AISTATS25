clear; clc; close all;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
fntsz = 18;
%%
LF = 1.5; muF = -0.5; % parameters of objective F

tol = 0;
num_mu2_vec = 3001;
mu2_vec = sort(unique([linspace(-LF, 2*LF, num_mu2_vec)]));
num_mu2_vec = length(mu2_vec);

num_L2_vec = 3001;
L2_vec = sort(unique([linspace(-muF, 2*LF, num_L2_vec)]));
num_L2_vec = length(L2_vec);

id_p = zeros(num_mu2_vec, num_L2_vec);
p_vals = zeros(num_mu2_vec, num_L2_vec);
sel_vec = [];
parfor LL2 = 1 : num_L2_vec
    L2 = L2_vec(LL2);
    mu1 = L2 + muF;
    id_vec = 0 * mu2_vec';
    p_vec = 0 * mu2_vec' - inf;
    for mm2 = 1 : num_mu2_vec
        mu2 = mu2_vec(mm2);
        L1 = mu2 + LF;
        if mu1 + mu2 > 0 && mu2 < L2
            [p_vec(mm2), id_vec(mm2)] = compute_DCA_rate_6_regimes(L1,mu1,L2,mu2);
        end
    end
    id_p(:, LL2) = id_vec;
    p_vals(:, LL2) = p_vec;
end
L2_vec_sel = L2_vec(L2_vec>=0);
mu1_vec = muF + L2_vec;
thr_S = -mu1_vec .* L2_vec_sel ./ (mu1_vec + L2_vec_sel);
mu2_lim = -mu1_vec .* (1+0*L2_vec);
%% add fake point of p2 for scaling the colorbar
id_L2 = find(L2_vec == LF);
id_mu2 = find(mu2_vec == muF);
id_p(id_mu2,id_L2) = 2;

%% figure
fig1 = figure(Color='white'); clf;
fig1.WindowState = 'maximized';
pause(1)
imagesc(L2_vec, mu2_vec, id_p); hold on; grid on;
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 16;
ax.GridLineWidth = 2;
xlabel('${L}_2$','Interpreter', 'latex', 'FontSize', fntsz)
ylabel('${\mu}_2$','Interpreter', 'latex', 'FontSize', fntsz)
C = [  .75 0 .75 % 4 cyclam
        0 .75 .75  % 3 cyan
        1, .647, 0 % 2 orange
        .5, .5, .5 % 1 gray
        ];
C(end+1,:)= 1;
colormap(flipud(C))
colorbar('Ticks',linspace(0,6,7),...
    'TickLabels',["--" "p1" "p2" "p3" "p4" % "p5" "p6" ...
    ],'FontSize', fntsz);
set(gca,'YDir','normal') 
%%
xlim([-muF Inf])

legend_names = {};

plot( L2_vec_sel, thr_S, 'r', 'Linewidth', 5 )
legend_names{end+1} = '$B=0$';

plot( [min(L2_vec) max(L2_vec)], [min(L2_vec) max(L2_vec)], 'k', 'Linewidth', 3 )
legend_names{end+1} = '$L_2 = \mu_2$';


%%
leg = legend(legend_names, 'Interpreter','latex', 'FontSize', fntsz, 'Location','northwest');
leg.AutoUpdate = "off";

x_ticks = gca().XTick;
x_ticks_new = sort(unique([x_ticks(1), -muF, x_ticks(2:end)]));
xticks(x_ticks_new)

%% example: best splitting when muF = -0.5; LF = 1.5; mu2 = 0.5; L2 = 1.5
mu2_ex = 0.5; L2_ex = 1.5;
mu1_ex = muF + L2_ex; L1_ex = LF + mu2_ex;
lambda_star = get_opt_lambda(L1_ex,mu1_ex,L2_ex,mu2_ex); % 0.6733
mu2_ex_1_opt = mu2_ex - lambda_star;
L2_ex_1_opt = L2_ex - lambda_star;
plot([L2_ex L2_ex_1_opt], [mu2_ex mu2_ex_1_opt], 'w--', 'Linewidth', 1)
plot(L2_ex, mu2_ex, 'r*', 'Linewidth', 3, 'MarkerSize',8, 'Marker','x')
plot(L2_ex_1_opt, mu2_ex_1_opt, 'g*', 'Linewidth', 3, 'MarkerSize',8)

%% best splitting when muF = -0.5; LF = 1.5; mu2 = 0.75; L2 = 1
mu2_ex = .75; L2_ex = 1;
mu1_ex = muF + L2_ex; L1_ex = LF + mu2_ex;
lambda_star = get_opt_lambda(L1_ex,mu1_ex,L2_ex,mu2_ex); % 0.5
mu2_ex_1_opt = mu2_ex - lambda_star;
L2_ex_1_opt = L2_ex - lambda_star;
% figure;
plot([L2_ex L2_ex_1_opt], [mu2_ex mu2_ex_1_opt], 'w--', 'Linewidth', 1)
plot(L2_ex, mu2_ex, 'r*', 'Linewidth', 3, 'MarkerSize',8, 'Marker','x')
plot(L2_ex_1_opt, mu2_ex_1_opt, 'g*', 'Linewidth', 3, 'MarkerSize',8)


%% best splitting when muF = -0.5; LF = 1.5; mu2 = 0.5; L2 = 2.5
mu2_ex = .5; L2_ex = 2.5;
mu1_ex = muF + L2_ex; L1_ex = LF + mu2_ex;
lambda_star = get_opt_lambda(L1_ex,mu1_ex,L2_ex,mu2_ex); % 0.8703
mu2_ex_1_opt = mu2_ex - lambda_star;
L2_ex_1_opt = L2_ex - lambda_star;
% figure;
plot([L2_ex L2_ex_1_opt], [mu2_ex mu2_ex_1_opt], 'w--', 'Linewidth', 1)
plot(L2_ex, mu2_ex, 'r*', 'Linewidth', 3, 'MarkerSize',8, 'Marker','x')
plot(L2_ex_1_opt, mu2_ex_1_opt, 'g*', 'Linewidth', 3, 'MarkerSize',8)

%% best splitting when muF = -0.5; LF = 1.5; mu2 = 0.5; L2 = 2.5
mu2_ex = -1; L2_ex = 2;
mu1_ex = muF + L2_ex; L1_ex = LF + mu2_ex;
lambda_star = get_opt_lambda(L1_ex,mu1_ex,L2_ex,mu2_ex); % -0.6749
mu2_ex_1_opt = mu2_ex - lambda_star;
L2_ex_1_opt = L2_ex - lambda_star;

plot([L2_ex L2_ex_1_opt], [mu2_ex mu2_ex_1_opt], 'w--', 'Linewidth', 1)
plot(L2_ex, mu2_ex, 'r*', 'Linewidth', 3, 'MarkerSize',8, 'Marker','x')
plot(L2_ex_1_opt, mu2_ex_1_opt, 'g*', 'Linewidth', 3, 'MarkerSize',8)


%%
p_vals_new = p_vals;
axis xy; 
pause(1)
%% add contour plot
contour_vals = unique(sort([0, 0.25, 0.5, .75, 1, 1.25, 1.5, 1.75, 2, 2.5]));
h_ax = gca;
h_ax_c = axes('position', get(h_ax, 'position'), 'Color', 'none');
[h_contour, Cntr] = contour(h_ax_c, p_vals_new, ...
               contour_vals , 'ShowText', 'on', 'LineWidth', 2, 'TextStep', 2, 'LineColor', 'w');   
hold on;
h_ax_c.Color = 'none';
h_ax_c.XTick = [];
h_ax_c.YTick = [];

x0 = min(gca().XLim); y0 = min(gca().YLim); 
ymax = max(gca().YLim); xmax = max(gca().XLim); m = 1;
plot([x0, num_mu2_vec], [num_mu2_vec*((LF-muF)/(LF+max(mu2_vec))) num_mu2_vec], 'k', 'LineWidth', 5)
ax = gca;
ax.FontWeight = 'bold';
%%
sppi = get(groot,"ScreenPixelsPerInch");
exportgraphics(gcf,"contour_plot_all_pi_fixed_muF_LF.pdf","Resolution",sppi)