clear; clc; close all;

az0 = 125.61; % I-10 heading in degrees;
ex0 = sin(az0/180*pi);
ey0 = cos(az0/180*pi);
stlo0 = -116.322817210; % location of station 1001
stla0 = 33.816908394;

event_dir = '../event15/';
load([event_dir 'mw/dists.mat']);
stnm = load([event_dir 'stidx_lin.txt']);
for j = 1:length(stnm)-2
    t_b = 300; t_e = 900;
    stidx1 = stnm(j);
    stidx2 = stnm(j+1);
    fnm1 = [event_dir 'data/R' num2str(stidx1) '.SAC'];
    fnm2 = [event_dir 'data/R' num2str(stidx2) '.SAC'];
    
    [hdr1,data1_0] = load_sac(fnm1);
    [hdr2,data2_0] = load_sac(fnm2);

    % Median amplitude
    n1 = round(t_b/hdr1.delta);
    n2 = round(t_e/hdr1.delta);
    data1 = data1_0(n1: n2);
    data2 = data2_0(n1: n2);
    A_mdn(j) = sqrt(mean(data1.^2));

    y0 = (hdr2.stla+hdr1.stla - 2*stla0)/360*pi*6371;
    x0 = (hdr2.stlo+hdr1.stlo - 2*stlo0)*cos(stla0/180*pi)/360*pi*6371;
    dist_0(j) = ex0*y0 - ey0*x0;
    
    % Rayleigh phase velocity
    dy = (hdr2.stla - hdr1.stla)/180*pi*6371;
    dx = (hdr2.stlo - hdr1.stlo)*cos(hdr1.stla/180*pi)/180*pi*6371;

    for i  = 6:10
        t_b = t(i)-30; t_e = t(i)+30;
        n1 = round(t_b/hdr1.delta);
        n2 = round(t_e/hdr1.delta);
        data1 = data1_0(n1: n2);
        data2 = data2_0(n1: n2);

        x_tr = -dists(i)*ex0;
        y_tr = -dists(i)*ey0;
        x_st = x0 - 4.5*ey0;
        y_st = y0 + 4.5*ex0;
        ex = (x_st-x_tr)/sqrt((x_st-x_tr)^2+(y_st-y_tr)^2);
        ey = (y_st-y_tr)/sqrt((x_st-x_tr)^2+(y_st-y_tr)^2);
        dist = ex*dx + ey*dy;
        
        [r, lags]= xcorr(data1,data2,500,'coeff');
        [~,I] = max(r);
        dt0 = lags(I)*hdr1.delta;
        V0(j,i) = -dist*1e3/dt0;
        CC(j,i) = r(I);
    end
    
end
save([event_dir 'S2.mat'],'CC','dist_0','V0','A_mdn','stnm')

%%
clearvars -except event_dir; clc; close all;
load([event_dir 'S2.mat'])
d0 = 4.5; %km

figure(1)
A_mdn = smooth(A_mdn, 4, 'moving');
plot(dist_0,A_mdn/(max(A_mdn))*3,'r','linewidth',4); hold on;
for i = 1:1:length(dist_0)
    if mod(i,3) == 1
        text(dist_0(i), 0.05, ['R' num2str(stnm(i))],'rotation', 45,'fontsize',12)
    end
    if stnm(i) == 1028
        plot([dist_0(i) dist_0(i)],[0 3],'-k','linewidth',2)
    elseif  stnm(i) == 1098
        plot([dist_0(i) dist_0(i)],[0 3],'--k','linewidth',2)
    end
end

All_dist = repmat(dist_0,5,1)';
All_V0 = V0(:, 6:10);
All_CC = CC(:, 6:10);

I = All_CC > 0.6;
scatter(All_dist(I), All_V0(I)/1e3,300, All_CC(I), 'filled'); hold on;

dist_ctr = 0:0.1:3;
for i = 1:length(dist_ctr)
    dist1 = dist_ctr(i) - 0.1;
    dist2 = dist_ctr(i) + 0.1;
    I1 = All_dist>dist1;
    I2 = All_dist<=dist2;
    V0_ctr(i) = mean(All_V0(I1&I2&I));
    V0_std(i) = std(All_V0(I1&I2&I));
end
I = V0_ctr~=Inf;
dist_ctr = dist_ctr(I);
V0_ctr = V0_ctr(I);
V0_std = V0_std(I);
V0_ctr = smooth(V0_ctr, 4, 'moving');
V0_std = smooth(V0_std, 4, 'moving');

plot(dist_ctr,V0_ctr/1e3,'-k','linewidth',4);
plot(dist_ctr,(V0_ctr-V0_std)/1e3,'-k','linewidth',2);
plot(dist_ctr,(V0_ctr+V0_std)/1e3,'-k','linewidth',2);

colorbar; box on; grid on
text(6.4,2.8,'NE','fontsize',40)
text(0.2+d0,2.8,'SW','fontsize',40)

legend('Median Amplitude','Banning strand','Mission Creek strand','CC')
title('Rayleigh wave velocity derived from Train Signal')
xlabel('Distance to I-10 (km)')
ylabel('Rayleigh Velocity (km/s)')
xlim([0 max(dist_0)])
ylim([0 3])
colormap(turbo)
caxis([0.6 1])
set(gca,'fontsize',20)
set(gcf,'renderer','painters')

save([event_dir 'velo.mat'], 'dist_ctr', 'V0_ctr', 'V0_std')