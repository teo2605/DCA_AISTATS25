clear; clc; close all;
fntsz = 18;

mu1 = 1; L1 = 2;
tol = 0;
num_mu2_vec = 2001;
mu2_vec = sort(unique([0, linspace(-mu1, 1.25*L1, num_mu2_vec)]));
num_mu2_vec = length(mu2_vec);

num_L2_vec = 2001;
L2_vec = sort(unique([0, linspace(-mu1, 1.5*L1, num_L2_vec)]));
num_L2_vec = length(L2_vec);

id_p = zeros(num_mu2_vec, num_L2_vec);
p_vals = zeros(num_mu2_vec, num_L2_vec);
parfor LL2 = 1 : num_L2_vec
    L2 = L2_vec(LL2);
    id_vec = 0 * mu2_vec';
    p_vec = 0 * mu2_vec' - inf;
    for mm2 = 1 : num_mu2_vec
        if mu1 + mu2_vec(mm2) > 0 && mu2_vec(mm2) < L2
            [p_vec(mm2), id_vec(mm2)] = compute_DCA_rate_6_regimes(L1,mu1,L2,mu2_vec(mm2));            
        end
    end
    id_p(:, LL2) = id_vec;
    p_vals(:, LL2) = p_vec;
end
L2_vec_sel = L2_vec(L2_vec>=0);
thr_B = -mu1 * L2_vec_sel ./ (mu1 + L2_vec_sel); % from condition B=0
mu2_lim = -mu1*(1+0*L2_vec);


%% figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

figure(Color='white');
imagesc(L2_vec, mu2_vec, id_p); hold on; grid on;
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;
ax.GridLineWidth = 2;
xlabel('${L}_2$','Interpreter', 'latex', 'FontSize', fntsz)
ylabel('${\mu}_2$','Interpreter', 'latex', 'FontSize', fntsz)
    
C = [
    .0 .5 .0 % 6 green
    1, 1, 0 % 5 yellowish
    .75 0 .75 % 4 cyclam
    0 .75 .75  % 3 cyan
    1, .647, 0 % 2 orange
    .5, .5, .5 % 1 gray
    ];
C(end+1,:)= 1;
colormap(flipud(C))
colorbar('Ticks',linspace(0,6,7),...
    'TickLabels',["--" "p1" "p2" "p3" "p4" "p5" "p6"],'FontSize', fntsz);

set(gca,'YDir','normal') 

xlim([-mu1 max(L2_vec)])
ylim([-mu1 max(mu2_vec)])

legend_names = {};

plot( L2_vec_sel, thr_B, 'r', 'Linewidth', 6 )
legend_names{end+1} = '$B=0$';

plot( [min(L2_vec) max(L2_vec)], [min(L2_vec) max(L2_vec)], 'k', 'Linewidth', 6 )
legend_names{end+1} = '$L_2 = \mu_2$';

%%
leg = legend(legend_names, 'Interpreter','latex', 'FontSize', fntsz, 'Location','northwest');
leg.AutoUpdate = "off";

%% overlay the contour plot
pause(.5)
fig1 = figure(1);
fig1.WindowState = 'maximized';
axis xy; 
pause(2)
h_ax = gca;
h_ax_c = axes('position', get(h_ax, 'position'), 'Color', 'none');
contour_vals = unique(sort([0,.2, .4, .6, .75, .85, .9, .95, 1., 1.1, 1.4, 2, 5, 10, 50]));
[h_contour, Cntr] = contour(h_ax_c, p_vals, ...
               contour_vals , 'ShowText', 'on', 'LineWidth', 2.5, 'TextStep', 1);
hold on;
h_ax_c.Color = 'none';
h_ax_c.XTick = [];
h_ax_c.YTick = [];

x0 = min(gca().XLim); y0 = min(gca().YLim); 
ymax = max(gca().YLim); xmax = max(gca().XLim); m = 1;
plot([x0, 1753], [y0 ymax], 'k', 'LineWidth', 6.5); % the separation line L2=mu2 is ticker to hide multiple overlayed contour values

ax = gca;
ax.FontWeight = 'bold';

%%
sppi = get(groot,"ScreenPixelsPerInch");
exportgraphics(gcf,"contour_plot_all_pi_mu1_1_L1_2.pdf","Resolution",sppi)