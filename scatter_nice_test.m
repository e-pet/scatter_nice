N = 10000;
mu = [2;3];
dat = [1, -.3; -.3, 1.5] * randn(2, N) + mu;
color = exp(sqrt(sum(dat.^2))) + exp(10-(1:N)/20);

figure;
tiledlayout(2, 2);
nexttile;
scatter(dat(1, :), dat(2, :), [], color);
title('scatter');
colorbar;
nexttile;
scatter_nice(dat(1, :), dat(2, :), color);
title('scatter\_nice');

% plotting two groups
N2 = 2000;
mu2 = [3; 2];
dat2 = [1, .4; .4, 0.7] * randn(2, N2) + mu2;

nexttile;
scatter(dat2(1, :), dat2(2, :));
hold on;
scatter(dat(1, :), dat(2, :));
legend();
title('scatter');

nexttile;
x_data = {dat(1, :), dat2(1, :)};
y_data = {dat(2, :), dat2(2, :)};
scatter_nice(x_data, y_data);
title('scatter\_nice');