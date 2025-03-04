clear; clc; close all;

L1 = 20; L2 = 2;
num_mu1_vec = 1001;
mu1_vec = sort(unique([0, L2, linspace(0, L1, num_mu1_vec)]));
num_mu1_vec = length(mu1_vec);
num_mu2_vec = 1001;
mu2_vec = sort(unique([0, linspace(-L1, L2, num_mu2_vec)]));
num_mu2_vec = length(mu2_vec);

id_p = zeros(num_mu2_vec, num_mu1_vec);
p_all = 0 * id_p;
for mm1 = 1 : num_mu1_vec
    mu1 = mu1_vec(mm1);
    id_vec = 0 * mu2_vec';
    p_vec = 0 * id_vec';
    for mm2 = 1 : num_mu2_vec
        if mu1 + mu2_vec(mm2) > 0
            [p_vec(mm2), id_vec(mm2)] = compute_DCA_rate_6_regimes(L1,mu1,L2,mu2_vec(mm2));
        end
    end
    id_p(:, mm1) = id_vec;
    p_all(:, mm1) = p_vec;
end
thr_B = -L2 * mu1_vec ./ (L2 + mu1_vec);
mu2_lim = -mu1_vec;


% set p6 (fake) to adjust the colorbar
id_mu1_p6 = find(mu1_vec == 0);
id_mu2_p6 = find(mu2_vec == 0);
id_p(id_mu2_p6, id_mu1_p6) = 6;


%% figure
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');


clf; gcf = imagesc(mu1_vec, mu2_vec, id_p); hold on; grid on;
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;
ax.GridLineWidth = 2;
xlabel('${\mu}_1$','Interpreter', 'latex', 'FontSize', 14)
ylabel('${\mu}_2$','Interpreter', 'latex', 'FontSize', 14)
title(sprintf('${L}_1=\\infty$, ${L}_2=%d$', L2),'Interpreter', 'latex', 'FontSize', 14)
    
C = [
    % .75 .75 0  % olive                
    % .0 .0 1 % blue
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
        'TickLabels',["--" "p1" "p2" "p3" "p4" "p5" "p6"]);
set(gca,'YDir','normal') 

xlim([0 L1])
ylim([-L1 L2])


legend_names = {};

plot( mu1_vec, thr_B, 'r', 'Linewidth', 3 )
legend_names{end+1} = '$B=0$';

plot( mu1_vec, mu2_lim, 'k', 'Linewidth', 3 )
legend_names{end+1} = '$\mu_1 + \mu_2 = 0$';

legend(legend_names, 'Interpreter','latex', 'FontSize', 14, 'Location','southwest','AutoUpdate','off');
%% change labels
xtick_labels = ax.XTickLabel; xtick_vals = ax.XTick;
xtick_labels(2) = {'$L_2$'}; xtick_vals(2)=L2;
xtick_labels(3) = {'$\dots$'};
xtick_labels(4) = {'$\infty$'};
xtick_labels = xtick_labels(1:4);
xtick_vals = xtick_vals([1:3,5]);
xticks(xtick_vals); xticklabels(xtick_labels);

ytick_labels = ax.YTickLabel; ytick_vals = ax.YTick;
ytick_labels(end+1) = {'$L_2$'}; ytick_vals(end+1)=L2;
ytick_labels(end-2) = {'$\vdots$'};
ytick_labels = ytick_labels(end-3 : end);
ytick_labels(1) = {'$-\infty$'};
ytick_vals = ytick_vals([1, end-2 : end]);
yticks(ytick_vals); yticklabels(ytick_labels);

%%
ax = gca;
ax.FontWeight = 'bold';
% return
exportgraphics(ax, 'All_regimes_1_step_L1=inf_L2=2_mu1_geq_0.pdf')