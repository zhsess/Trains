clear; clc; close all;

event_dir = '../event15/';
t_b = 300; t_e = 900;

%% BF subarray
stnm = load([event_dir 'stidx_BF.txt'])';
stnm = stnm(:); I = ~isnan(stnm);
stnm = sort(stnm(I));

[hdr0, ~] = load_sac([event_dir 'data/R' num2str(stnm(1)) '.SAC']);
n1 = round(t_b/hdr0.delta);
n2 = round(t_e/hdr0.delta);

All_x = [];
All_y = [];
All_t = [];
All_CC = [];
for l = 1:length(stnm)
    disp(l)
    stidx1 = stnm(l);
    fnm1 = [event_dir 'data/R' num2str(stidx1) '.SAC'];
    [hdr1, data1] = load_sac(fnm1);
    for j = 1:length(stnm)
        stidx2 = stnm(j);
        fnm2 = [event_dir 'data/R' num2str(stidx2) '.SAC'];
        [hdr2,data2] = load_sac(fnm2);
        
        x = (hdr2.stlo - hdr1.stlo)*cos(hdr1.stla/180*pi)/180*pi*6371;
        y = (hdr2.stla - hdr1.stla)/180*pi*6371;
        All_x = [All_x; x];
        All_y = [All_y; y];
        
        [r, lags]= xcorr(data1(n1:n2),data2(n1:n2),500,'coeff');
        [~,I] = max(r);
        All_t = [All_t; lags(I)*hdr1.delta];
        All_CC = [All_CC; r(I)];
    end
end
save([event_dir 'S1_BF.mat'], 'All_x', 'All_y', 'All_CC', 'All_t')
clearvars -except event_dir n1 n2

%% MCF subarray
stnm = load([event_dir 'stidx_MCF.txt'])';
stnm = stnm(:); I = ~isnan(stnm);
stnm = sort(stnm(I));

All_x = [];
All_y = [];
All_t = [];
All_CC = [];
for l = 1:length(stnm)
    disp(l)
    stidx1 = stnm(l);
    fnm1 = [event_dir 'data/R' num2str(stidx1) '.SAC'];
    [hdr1, data1] = load_sac(fnm1);
    for j = 1:length(stnm)
        stidx2 = stnm(j);
        fnm2 = [event_dir 'data/R' num2str(stidx2) '.SAC'];
        [hdr2,data2] = load_sac(fnm2);
        
        x = (hdr2.stlo - hdr1.stlo)*cos(hdr1.stla/180*pi)/180*pi*6371;
        y = (hdr2.stla - hdr1.stla)/180*pi*6371;
        All_x = [All_x; x];
        All_y = [All_y; y];
        
        [r, lags]= xcorr(data1(n1:n2),data2(n1:n2),500,'coeff');
        [~,I] = max(r);
        All_t = [All_t; lags(I)*hdr1.delta];
        All_CC = [All_CC; r(I)];
    end
end
save([event_dir 'S1_MCF.mat'], 'All_x', 'All_y', 'All_CC', 'All_t')

%% Plot
clearvars -except event_dir; clc; close all; figure(1)

load([event_dir 'S1_BF.mat'])
subplot(2,2,3)
I = All_CC>0.7;
scatter(All_x(I),All_y(I),100,-All_t(I),'filled'); axis equal; hold on;
colormap(turbo); colorbar; box on;
b = -All_t(I);
G = [All_x(I) All_y(I)];
m = G\b;
ex_1 = m(1)/(norm(m));
ey_1 = m(2)/(norm(m));
plot([-ex_1 ex_1],[-ey_1 ey_1],'linewidth',2);
xlim([-max(All_x(I)*1.5) max(All_x(I)*1.5)])
ylim([-max(All_x(I)*1.5) max(All_x(I)*1.5)])
xlabel('East'); ylabel('North')
title('BF')
set(gca,'fontsize',20)

load([event_dir 'S1_MCF.mat'])
subplot(2,2,1)
I = All_CC>0.7;
scatter(All_x(I),All_y(I),100,-All_t(I),'filled'); axis equal; hold on;
colormap(turbo); colorbar; box on;
b = -All_t(I);
G = [All_x(I) All_y(I)];
m = G\b;
ex_2 = m(1)/(norm(m));
ey_2 = m(2)/(norm(m));
plot([-ex_2 ex_2],[-ey_2 ey_2],'linewidth',2);
xlim([-max(All_x(I)*1.5) max(All_x(I)*1.5)])
ylim([-max(All_x(I)*1.5) max(All_x(I)*1.5)])
title('MCF')
ylabel('North')
set(gca,'fontsize',20)

subplot(1,2,2)
az0 = 125.61; % I-10 heading in degrees;
ex = sin((az0-90)/180*pi);
ey = cos((az0-90)/180*pi);
plot([0 ex],[0 ey],'-k','linewidth',2); axis equal; hold on;
plot([0 ex_1],[0 ey_1],'-r','linewidth',2); hold on
plot([0 ex_2],[0 ey_2],'-r','linewidth',2);
legend('I-10', 'BF', 'MCF')
set(gca,'fontsize',20)
