clear; clc; close all;

 % Location of I-10
az_r = 125.61;
ex_r = sin(az_r/180*pi);
ey_r = cos(az_r/180*pi);
lo0 = -116.3549;
la0 = 33.7862;

% location of station 1001
az0 = 39;
ex0 = sin(az0/180*pi);
ey0 = cos(az0/180*pi);
stlo0 = -116.322817210; 
stla0 = 33.816908394;

[disto_0, azo_0] = distance(la0, lo0, stla0, stlo0);
disto_0 = disto_0*pi*6371/180;
xo_0 = disto_0*sin(azo_0/180*pi);
yo_0 = disto_0*cos(azo_0/180*pi);


% % %% Linear Array
% % for evid = 1:29
% %     event_id = ['event' num2str(evid, '%02d')];
% %     event_dir = ['../event/' event_id '/'];
% %     load([event_dir 'mw/dists.mat']);
% %     stnm = load([event_dir 'stidx_lin.txt']);
% %     for j = 1:length(stnm)-2
% %         t_b = 300; t_e = 900;
% %         stidx1 = stnm(j);
% %         stidx2 = stnm(j+1);
% %         fnm1 = [event_dir 'data/R' num2str(stidx1) '.SAC'];
% %         fnm2 = [event_dir 'data/R' num2str(stidx2) '.SAC'];
% %         
% %         [hdr1,data1_0] = load_sac(fnm1);
% %         [hdr2,data2_0] = load_sac(fnm2);
% %     
% %         % Median amplitude
% %         n1 = round(t_b/hdr1.delta);
% %         n2 = round(t_e/hdr1.delta);
% %         data1 = data1_0(n1: n2);
% %         data2 = data2_0(n1: n2);
% %         A_mdn(j) = sqrt(mean(data1.^2));
% %     
% %         y0 = (hdr2.stla+hdr1.stla - 2*stla0)/360*pi*6371;
% %         x0 = (hdr2.stlo+hdr1.stlo - 2*stlo0)*cos(stla0/180*pi)/360*pi*6371;
% %         dist_0(j) = ex0*x0 + ey0*y0;
% %         
% %         % Rayleigh phase velocity
% %         dy = (hdr2.stla - hdr1.stla)/180*pi*6371;
% %         dx = (hdr2.stlo - hdr1.stlo)*cos(hdr1.stla/180*pi)/180*pi*6371;
% %     
% %         for i  = 5:10
% %             t_b = t(i)-30; t_e = t(i)+30;
% %             n1 = round(t_b/hdr1.delta);
% %             n2 = round(t_e/hdr1.delta);
% %             data1 = data1_0(n1: n2);
% %             data2 = data2_0(n1: n2);
% %     
% %             x_tr = dist_tr_BFs(i)*ex_r;
% %             y_tr = dist_tr_BFs(i)*ey_r;
% %             x_st = x0 + xo_0;
% %             y_st = y0 + yo_0;
% %             ex = (x_st-x_tr)/sqrt((x_st-x_tr)^2+(y_st-y_tr)^2);
% %             ey = (y_st-y_tr)/sqrt((x_st-x_tr)^2+(y_st-y_tr)^2);
% %             dist = ex*dx + ey*dy;
% %             
% %             [r, lags]= xcorr(data1,data2,500,'coeff');
% %             [~,I] = max(r);
% %             dt0 = lags(I)*hdr1.delta;
% %             V0(j,i) = -dist*1e3/dt0;
% %             CC(j,i) = r(I);
% %         end
% %         
% %     end
% %     save([event_dir 'S2.mat'],'CC','dist_0','V0','A_mdn','stnm','disto_0')
% % end


% % %% MCF 2D array
% % stnm_MCF = load([event_dir 'stidx_MCF.txt']);
% % for i = 2:5
% %     for j = 1:18
% %         t_b = 300; t_e = 900;
% %         stidx1 = stnm_MCF(i, j);
% %         stidx2 = stnm_MCF(i, j+1);
% %         fnm1 = [event_dir 'data/R' num2str(stidx1) '.SAC'];
% %         fnm2 = [event_dir 'data/R' num2str(stidx2) '.SAC'];
% %     
% %         [hdr1,data1_0] = load_sac(fnm1);
% %         [hdr2,data2_0] = load_sac(fnm2);
% % 
% %         y0 = (hdr2.stla+hdr1.stla - 2*stla0)/360*pi*6371;
% %         x0 = (hdr2.stlo+hdr1.stlo - 2*stlo0)*cos(stla0/180*pi)/360*pi*6371;
% %         dist_0_MCF(i,j) = ex0*x0 + ey0*y0;
% %     
% %         % Rayleigh phase velocity
% %         dy = (hdr2.stla - hdr1.stla)/180*pi*6371;
% %         dx = (hdr2.stlo - hdr1.stlo)*cos(hdr1.stla/180*pi)/180*pi*6371;
% % 
% %         for k  = 6:10
% %             t_b = t(k)-30; t_e = t(k)+30;
% %             n1 = round(t_b/hdr1.delta);
% %             n2 = round(t_e/hdr1.delta);
% %             data1 = data1_0(n1: n2);
% %             data2 = data2_0(n1: n2);
% %     
% %             x_tr = dist_tr(k)*ex_r;
% %             y_tr = dist_tr(k)*ey_r;
% %             x_st = x0 + xo_0;
% %             y_st = y0 + yo_0;
% %             ex = (x_st-x_tr)/sqrt((x_st-x_tr)^2+(y_st-y_tr)^2);
% %             ey = (y_st-y_tr)/sqrt((x_st-x_tr)^2+(y_st-y_tr)^2);
% %             dist = ex*dx + ey*dy;
% %             
% %             [r, lags]= xcorr(data1,data2,500,'coeff');
% %             [~,I] = max(r);
% %             dt0 = lags(I)*hdr1.delta;
% %             V0_MCF(i,j,k) = -dist*1e3/dt0;
% %             CC_MCF(i,j,k) = r(I);
% %         end
% %     end
% % end
% % save([event_dir 'S2_MCF.mat'],'CC_MCF','dist_0_MCF','V0_MCF', 'stnm_MCF')


%%
for evid = 1:29
    event_id = ['event' num2str(evid, '%02d')];
    event_dir = ['../event/' event_id '/'];
    load([event_dir 'S2.mat'])
    if evid==1
        All_dist = repmat(dist_0,6,1)';
        All_V0 = V0(:, 5:10);
        All_CC = CC(:, 5:10);
        A_mm = A_mdn/max(A_mdn);
    else
        All_dist = [All_dist; repmat(dist_0,6,1)'];
        All_V0 = [All_V0; V0(:, 5:10)];
        All_CC = [All_CC; CC(:, 5:10)];
        A_mm = A_mm + A_mdn/max(A_mdn);
    end
end

I = (All_CC>0.6) & (All_V0>0) & (All_V0<4000);

All_dist_CC = All_dist(I);
All_V0_CC = All_V0(I);
All_CC_CC = All_CC(I);

% All_dist_MCF = repmat(dist_0_MCF(2:5,:),1,1,5);
% All_V0_MCF = V0_MCF(2:5, :, 6:10);
% All_CC_MCF = CC_MCF(2:5, :, 6:10);
% I_MCF = (All_CC_MCF>0.7) & (All_V0_MCF>0) & (All_V0_MCF<3000);


%% Velocity along linear array
% stnm = stnm(1:112);
% v_mean = zeros(112, 1);
% v_std = zeros(112, 1);
% for i = 1:112
%     v_ls = [];
%     for j = 1:5
%         if I(i,j) == 1
%             v_ls = [v_ls, All_V0(i,j)];
%         end
%     end
%     v_mean(i) = mean(v_ls);
%     v_std(i) = std(v_ls);
% end
% v_mean = smooth(v_mean, 5, 'moving');
% v_std = smooth(v_std, 5, 'moving');
% save([event_dir 'velo_sta.mat'], 'stnm', 'v_mean', 'v_std')

%%

fig = figure(1);
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);

d0 = 4.5; %km
All_V0_CC = sqrt(800*All_V0_CC);

yyaxis left;
plot(All_dist_CC, All_V0_CC/1e3,'.', MarkerSize=20, Color=[0.6,0.6,0.6]); hold on;

dist_ctr = 0:0.05:3;
for i = 1:length(dist_ctr)
    dist1 = dist_ctr(i) - 0.05;
    dist2 = dist_ctr(i) + 0.05;
    I1 = All_dist_CC>dist1;
    I2 = All_dist_CC<=dist2;
    V0_ctr(i) = mean(All_V0_CC(I1&I2));
    V0_std(i) = std(All_V0_CC(I1&I2));
end
I = V0_ctr~=Inf;
dist_ctr = dist_ctr(I);
V0_ctr = V0_ctr(I);
V0_std = V0_std(I);
V0_ctr = smooth(V0_ctr, 3, 'moving');
V0_std = smooth(V0_std, 3, 'moving');

plot(dist_ctr,V0_ctr/1e3,'-k','linewidth',4);
plot(dist_ctr,(V0_ctr-V0_std)/1e3,'-k','linewidth',2);
plot(dist_ctr,(V0_ctr+V0_std)/1e3,'-k','linewidth',2);

ylabel('Rayleigh wave velocity (km/s)')


% Amplitude
yyaxis right;

A_corrected = A_mm.*sqrt(dist_0+disto_0);
A_plot = smooth(A_corrected, 3, 'moving');
plot(dist_0,A_plot/(max(A_plot)),'r','linewidth',4); hold on;
for i = 1:1:length(dist_0)
%     if mod(i,3) == 1
%         text(dist_0(i), 0.05, ['R' num2str(stnm(i))],'rotation', 45,'fontsize',20)
%     end
    if stnm(i) == 1028
        plot([dist_0(i) dist_0(i)],[0 2],'-k','linewidth',2)
    elseif  stnm(i) == 1098
        plot([dist_0(i) dist_0(i)],[0 2],'--k','linewidth',2)
    end
end
ylabel('Relative Amplitude')


legend('','Velocity','1-\sigma of Velocity','','Amplitude (corrected)')
xlabel('Distance to R1001 (km)')
xlim([0 max(dist_0)])
ylim([0 1])
set(gca,'fontsize',20)
set(gcf,'renderer','painters')


save('../info/velo_amp.mat', 'dist_ctr', 'V0_ctr', 'V0_std', 'A_corrected')