clear;clc;close all;

az0 = 125.61; % Direction of I-10
ex0 = sin((az0-90)/180*pi);
ey0 = cos((az0-90)/180*pi);

% location of station 1001
stlo0 = -116.322817210;
stla0 = 33.816908394;
d0 = 4.357; %km
f = 5; %Hz

%% event 15
event_dir = '../event15/';
stnm = load([event_dir 'stidx_6.txt']);
n_window = 13;
t_b = 390; dt = 30;
loc = [7.0031, 6.0443, 4.9083, 3.0139, 1.9768, 1.4158, 0.74124,...
    -0.0238, -0.65992, -1.7066, -2.5711, -3.2843, -4.3724];
for j = 1:length(stnm)
    stidx = stnm(j);
    fnm = [event_dir 'data/R' num2str(stidx) '.SAC'];
    [hdr,data] = load_sac(fnm);
    y = (hdr.stla - stla0)/180*pi*6371;
    x = (hdr.stlo - stlo0)*cos(stla0/180*pi)/180*pi*6371;
    dist_0(j) = sqrt(x^2+y^2);
    for i = 1:n_window
        n_1 = round((t_b + dt*(i-1))/hdr.delta);
        n_2 = round((t_b + dt*(i+1))/hdr.delta);
        A_mdn(i, j) = median(abs(data(n_1:n_2)));
    end
end

%% event22
event_dir = 'event22/';
n_window = 12;
t_b_0 = 150; dt = 30;
loc_22 = [-0.88488, -0.61877, -0.52079, -0.081876, 0.18516, 0.5258,...
    1.0083, 1.3178, 1.4179, 1.8582, 2.7111, 2.9692];

for j = 1:length(stnm)
    stidx = stnm(j);
    fnm = [event_dir 'data/R' num2str(stidx) '.SAC'];
    [hdr,data] = load_sac(fnm);
    for i = 1:n_window
        n_1 = round((t_b + dt*(i-1))/hdr.delta);
        n_2 = round((t_b + dt*(i+1))/hdr.delta);
        A_mdn_22(i, j) = median(abs(data(n_1:n_2)));
    end
end
%%
figure(1)
dist_all = (dist_0+d0)';
for i = 1:13
    A_mdn(i, :) = log(A_mdn(i,:)/A_mdn(i,1));
    A_mdn(i, :) = A_mdn(i, :)./(dist_all./sqrt(dist_all.^2+loc(i)^2))';
    plot(dist_all, A_mdn(i, :), 'LineWidth', 1, 'Color', [200 200 200]/255);
    hold on
end
for i = 1:12
    A_mdn_22(i, :) = log(A_mdn_22(i,:)/A_mdn_22(i,1));
    A_mdn_22(i, :) = A_mdn_22(i, :)./(dist_all./sqrt(dist_all.^2+loc_22(i)^2))';
    plot(dist_all, A_mdn_22(i, :), 'LineWidth', 1, 'Color', [200 200 200]/255);
    hold on
end
%%
stacked = smooth(mean([A_mdn;A_mdn_22], 1), 10);
A_std = std([A_mdn;A_mdn_22], 0, 1)';
plot(dist_all, stacked, 'LineWidth', 2, 'color', 'red'); hold on
plot(dist_all, stacked-A_std, 'LineWidth', 1, 'color', 'k'); hold on
plot(dist_all, stacked+A_std, 'LineWidth', 1, 'color', 'k'); hold on
patch([dist_all;flipud(dist_all)],...
    [stacked-A_std;flipud(stacked+A_std)],'r','FaceA',.2,'EdgeA',0)

set(gca,'FontSize',15);
xlabel('Distance to R1001')
ylabel('lnA')
xlim([min(dist_all) max(dist_all)])

%%
load('event15/velo.mat')
velo_model = interp1(dist_ctr, V0_ctr, dist_0);
velo_std = interp1(dist_ctr, V0_std, dist_0);

for i = 1:length(dist_all)-1
    Q(i) = 1000*pi*f*(dist_0(i+1)-dist_0(i))/((stacked(i)- ...
        stacked(i+1))*velo_model(i));
    Q_std(i) = Q(i)*sqrt((A_std(i)/stacked(i))^2+(velo_std(i)/velo_model(i))^2);
end

figure(2)
Q([6,11,16,17]) = [];
Q_std([6,11,16,17]) = [];
dist_all([6,11,16,17]) = [];
Q_mean = mean(Q);

plot([4,8], [Q_mean Q_mean], 'r', 'LineWidth', 2); hold on
errorbar(dist_all(1:end-1), Q, Q_std, 'o', 'linewidth', ...
    2, 'color', [100 100 100]/255); hold on
scatter(dist_all(1:end-1), Q, 100, 'filled', 'k');
legend('Mean Q-value')
set(gca,'FontSize',20);
xlabel('Distance to I-10(km)')
ylabel('Q-value')
xlim([min(dist_all) max(dist_all)])