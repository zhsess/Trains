clc; close all; clear
% 1 for BF, 2 for MCF

event_dir = '../event15/';

% event15
dist_bs = zeros(100,13);
th_mcfs = zeros(100,13);
th_bfs = zeros(100,13);

az0 = 125.61; % I-10 heading in degrees;
ex = sin(az0/180*pi); ey = cos(az0/180*pi);

lo0 = -116.3549; la0 = 33.7862; % Center of I-10
BF_lo = -116.320897; BF_la = 33.820251; % BF
MCF_lo = -116.305676; MCF_la = 33.834889; % MCF

[dist_BF, az_BF] = distance(la0, lo0, BF_la, BF_lo);
[dist_MCF, az_MCF] = distance(la0, lo0, MCF_la, MCF_lo);
dist_BF = dist_BF*pi*6371/180;
dist_MCF = dist_MCF*pi*6371/180;
x_BF = dist_BF*sin(az_BF/180*pi);
y_BF = dist_BF*cos(az_BF/180*pi);
x_MCF = dist_MCF*sin(az_MCF/180*pi);
y_MCF = dist_MCF*cos(az_MCF/180*pi);


for i = 1:100 % times of boot-strap
    randI_BF = rand(5151,1)<0.6;
    randI_MCF = rand(4851,1)<0.6;
    disp(i)
    for fig_idx = 1:13
        load([event_dir 'mw/' num2str(fig_idx) '_BF.mat'])
        I =(All_CC>0.7)&(abs(All_t)<0.1)&randI_BF;
        b = -All_t(I);
        G = [All_x(I) All_y(I)];
        m = G\b;
        ex_BF = m(1)/(norm(m));
        ey_BF = m(2)/(norm(m));
    
        load([event_dir 'mw/' num2str(fig_idx) '_MCF.mat'])
        I = (All_CC>0.7)&(abs(All_t)<0.1)&randI_MCF;
        b = -All_t(I);
        G = [All_x(I) All_y(I)];
        m = G\b;
        ex_MCF = m(1)/(norm(m));
        ey_MCF = m(2)/(norm(m));
         
        dist_tr_BF = (ey_BF*x_BF-ex_BF*y_BF)/(ex*ey_BF-ey*ex_BF);
        dist_tr_MCF = (ey_MCF*x_MCF-ex_MCF*y_MCF)/(ex*ey_MCF-ey*ex_MCF);
        dist = (dist_tr_BF+dist_tr_MCF)/2;
        dist_bs(i, fig_idx) = dist;

        th_BF = 90 - atan(ey_BF/ex_BF)*180/pi;
        th_MCF = 90 - atan(ey_MCF/ex_MCF)*180/pi;
        th_BFs(i, fig_idx) = th_BF;
        th_MCFs(i, fig_idx) = th_MCF;

    end
end

save([event_dir 'mw/dists_bs.mat'], 'dist_bs', 'th_BFs', 'th_MCFs')



