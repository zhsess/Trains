clear; clc; close all;

az0 = 125.61; % I-10 heading in degrees;
ex0 = sin((az0-90)/180*pi);
ey0 = cos((az0-90)/180*pi);
stlo0 = -116.322817210; % location of station 1001
stla0 = 33.816908394;

src_la = 33.790690;
src_lo = -116.363377;

event_dir = '../event15/';
stnm = load([event_dir 'stidx_lin.txt']);

t_b = 300; t_e = 900;

for j = 1:length(stnm)-1
    disp(j)
    stidx1 = stnm(j);
    stidx2 = stnm(j+1);
    fnm1 = [event_dir 'data/R' num2str(stidx1) '.SAC'];
    fnm2 = [event_dir 'data/R' num2str(stidx2) '.SAC'];
    
    [hdr1,data1] = load_sac(fnm1);
    [hdr2,data2] = load_sac(fnm2);
    n1 = round(t_b/hdr1.delta);
    n2 = round(t_e/hdr1.delta);
    data1 = data1(n1: n2);
    data2 = data2(n1: n2);
    
    % Median amplitude
    A_mdn(j) = median(abs(data1));
    
    % Rayleigh phase velocity
    y = (hdr2.stla - hdr1.stla)/180*pi*6371;
    x = (hdr2.stlo - hdr1.stlo)*cos(hdr1.stla/180*pi)/180*pi*6371;
    az = azimuth(src_la,src_lo,hdr1.stla,hdr1.stlo);
    ex = sin(az/180*pi);
    ey = cos(az/180*pi);
    dist = ex*x + ey*y;
    
    [r, lags]= xcorr(data1,data2,500,'coeff');
    [~,I] = max(r);
    dt0 = lags(I)*hdr1.delta;
    V0(j) = dist*1e3/abs(lags(I)*hdr1.delta);
    CC(j) = r(I);
    
%     % plot pdf
%     f = figure('visible','off');
%     clf;
%     set(gcf,'PaperPositionMode','auto');
%     set(gcf,'Units','inches');
%     afFigurePosition = [0 0 12 13.5];
%     set(gcf,'Position',afFigurePosition);
%     set(gcf,'PaperSize',[afFigurePosition(1)*2+afFigurePosition(3) afFigurePosition(2)*2+afFigurePosition(4)]);
% 
%     subplot(2,1,1)
%     t = (n1:n2)*hdr1.delta - t_b;
%     plot(t, data1, 'k'); hold on;
%     plot(t, data2, 'r');
%     set(gca,'fontsize',15)
%     title(['Station R' num2str(stidx1) ' & R' num2str(stidx2)]);
%     xlim([0 t_e-t_b])
%     xlabel('Time/s')
%     ylabel('Amplitude')
%     
%     subplot(2,1,2)
%     plot(lags*hdr1.delta,r,'-k','linewidth',2); hold on;
%     plot([dt0 dt0],[-1 1],'r--','linewidth',2);
%     plot([0 0],[-1 1],'--','color',[0.5 0.5 0.5],'linewidth',2)
%     ylabel('Xcorr')
%     xlabel('Time/s')
%     xlabel('time (s)')
%     title(['Apparent Velocity:' num2str(V0(j)) ' m/s'])
%     set(gca,'fontsize',15)
%     ylim([-1 1])
%     xlim([-0.5 0.5])
%     
%     saveas(f,[event_dir '/R' num2str(stidx1)  '_Xcc.pdf'])
    
    y = (hdr2.stla+hdr1.stla - 2*stla0)/360*pi*6371;
    x = (hdr2.stlo+hdr1.stlo - 2*stlo0)*cos(stla0/180*pi)/360*pi*6371;
    dist_all(j) = ex0*x + ey0*y;
end
save([event_dir 'S2_CC.mat'],'CC','dist_all','V0','A_mdn','stnm')

%%
clearvars -except event_dir; clc; close all;
load([event_dir 'S2_CC.mat'])
d0 = 4.357; %km

figure(1)
plot(dist_all+d0,A_mdn/(max(A_mdn))*3,'r','linewidth',4); hold on;
for i = 1:1:length(dist_all)
    if mod(i,3) == 1
        text(dist_all(i)+d0,0.05,['R' num2str(stnm(i))],'rotation', 45,'fontsize',12)
    end
    if stnm(i) == 1028
        plot([dist_all(i)+d0 dist_all(i)+d0],[0 3],'-k','linewidth',2)
    elseif  stnm(i) == 1098
        plot([dist_all(i)+d0 dist_all(i)+d0],[0 3],'--k','linewidth',2)
    end
end

All_dist = dist_all;
All_V0 = V0;
All_CC = CC;

% event_dir = 'event22/';
% load([event_dir 'V0_CC.mat'])
% All_dist = horzcat(All_dist, dist_all);
% All_V0 = horzcat(All_V0, V0);
% All_CC = horzcat(All_CC, CC);

I = All_CC > 0.6;
scatter(All_dist(I)+d0, All_V0(I)/1e3,300, All_CC(I), 'filled'); hold on;

dist_ctr = 0:0.1:3;
for i = 1:length(dist_ctr)
    dist1 = dist_ctr(i) - 0.3;
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

plot(dist_ctr+d0,V0_ctr/1e3,'-k','linewidth',4);
plot(dist_ctr+d0,(V0_ctr-V0_std)/1e3,'-k','linewidth',2);
plot(dist_ctr+d0,(V0_ctr+V0_std)/1e3,'-k','linewidth',2);

colorbar; box on; grid on
text(6.4+d0,2.8,'NE','fontsize',40)
text(0.2+d0,2.8,'SW','fontsize',40)

legend('Median Amplitude','Banning strand','Mission Creek strand','CC')
title('Rayleigh wave velocity derived from Train Signal')
xlabel('Distance to I-10 (km)')
ylabel('Rayleigh Velocity (km/s)')
xlim([d0 max(dist_all)+d0])
ylim([0 3])
colormap(turbo)
caxis([0.6 1])
set(gca,'fontsize',20)
set(gcf,'renderer','painters')

save([event_dir 'velo.mat'], 'dist_ctr', 'V0_ctr', 'V0_std')