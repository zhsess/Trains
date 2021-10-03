clear;clc;close all;

az0 = 125.61;
ex0 = sin((az0-90)/180*pi);
ey0 = cos((az0-90)/180*pi);
stlo0 = -116.322817210; % location of station 1001
stla0 = 33.816908394;

%% event 15
event_dir = '../event15/';
stnm = load([event_dir 'stidx_lin.txt']);
n_window = 13;
t = 310:10:910;
dt = 30;
for j = 1:length(stnm)
    stidx = stnm(j);
    fnm = [event_dir 'data/R' num2str(stidx) '.SAC'];
    [hdr,data] = load_sac(fnm);
    y = (hdr.stla - stla0)/180*pi*6371;
    x = (hdr.stlo - stlo0)*cos(stla0/180*pi)/180*pi*6371;
    dist_0(j) = ex0*x + ey0*y;
    for i = 1:length(t)
        n_1 = round((t(i)-dt)/hdr.delta);
        n_2 = round((t(i)+dt)/hdr.delta);
        A_mdn(i, j) = median(abs(data(n_1:n_2)));
    end
end

figure(1)
c_rgb = colormap(jet);
idx_dist = round((dist_0-dist_0(1))/(dist_0(end)-dist_0(1))*length(c_rgb));
idx_dist(1) = 1;
for j = 1:length(stnm)
    plot(t-t(1), A_mdn(:,j),'linewidth',1.5,'color',c_rgb(idx_dist(j),:));
    hold on
end
xlim([0,600])
xlabel('Time(s)')
ylabel('Amplitude')
set(gca,'FontSize',20);
colorbar
colorbar('Ticks',[0, 53/256, 210/256, 1],...
        'TickLabels',{'0','BF','MCF','3km'})