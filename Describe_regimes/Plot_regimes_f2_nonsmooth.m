clear; clc; close all;

L2 = 20; L1 = 2;
num_mu1_vec = 1001;
mu1_vec = sort(unique([0, linspace(0, L1, num_mu1_vec)]));
num_mu1_vec = length(mu1_vec);
num_mu2_vec = 1001;
mu2_vec = sort(unique([0, L1, linspace(-L1, L2, num_mu2_vec)]));
num_mu2_vec = length(mu2_vec);

id_p = zeros(num_mu2_vec, num_mu1_vec);
for mm1 = 1 : num_mu1_vec
    mu1 = mu1_vec(mm1);
    id_vec = 0 * mu2_vec';
    for mm2 = 1 : num_mu2_vec
        mu2 = mu2_vec(mm2);
        if mu1 + mu2 > 0
            if (mu2 >= 0)
                if (mu2 >= L1)
                    id_vec(mm2) = 6;
                else
                    id_vec(mm2) = 2;
                end
            else
                id_vec(mm2) = 3;
            end
        end
    end
    id_p(:, mm1) = id_vec;
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
figure(1); clf;

gcf = imagesc(mu1_vec, mu2_vec, id_p); hold on; grid on;
ax = gca;
ax.FontWeight = 'bold';
ax.FontSize = 13;
ax.GridLineWidth = 2;
xlabel('${\mu}_1$','Interpreter', 'latex', 'FontSize', 14)
ylabel('${\mu}_2$','Interpreter', 'latex', 'FontSize', 14)
title(sprintf('${L}_1=%d$, ${L}_2=\\infty$', L1),'Interpreter', 'latex', 'FontSize', 14)
    

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

% plot( mu1_vec, mu2_lim, 'r', 'Linewidth', 3 )
% legend_names{end+1} = 'threshold $B$';

plot( mu1_vec, mu2_lim, 'k', 'Linewidth', 3 )
legend_names{end+1} = '$\mu_1 + \mu_2 = 0$';


% plot([L2 L2], [-L2/2 L2],'r', 'Linewidth', 3 )

leg = legend(legend_names, 'Interpreter','latex', 'FontSize', 14, 'Location','southwest','AutoUpdate','off');

% set(leg,...
%     'Position',[0.116904761904762 0.120158730158731 0.260006815987136 0.117380954924084],...
%     'FontSize',14,...
%     'AutoUpdate','off');

%% change labels
xtick_labels = [{'0'}, {'$\dots$'}, {'$L_1$'}];
xtick_vals = [0, L1/2, L1];
xticks(xtick_vals); xticklabels(xtick_labels);


ytick_labels = [{'$-L_1$'}, {'0'}, {'$L_1$'}, {'$\vdots$'}, {'$\infty$'}];
ytick_vals = [-L1, 0, mu2_vec(find(mu2_vec == L1) + 2), 2*L1, 3*L1];
yticks(ytick_vals); yticklabels(ytick_labels);

ylim([-L1, 3*L1])
%%
ax = gca;
ax.FontWeight = 'bold';
exportgraphics(ax, 'All_regimes_1_step_L2=inf_L1=2_mu1_geq_0.pdf')