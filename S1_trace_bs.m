clc; close all; clear
% 1 for BF, 2 for MCF

event_dir = '../event15/';
if exist([event_dir 'mw'],'dir')==0.
    mkdir([event_dir 'mw']);
end
% event15
t_0 = 390; dt = 30;

dist_bs = zeros(100,13);
th_mcfs = zeros(100,13);
th_bfs = zeros(100,13);
e2 = [];
t = 390:30:750;

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
        ex1 = m(1)/(norm(m));
        ey1 = m(2)/(norm(m));
    
        load([event_dir 'mw/' num2str(fig_idx) '_MCF.mat'])
        I = (All_CC>0.7)&(abs(All_t)<0.1)&randI_MCF;
        b = -All_t(I);
        G = [All_x(I) All_y(I)];
        m = G\b;
        ex2 = m(1)/(norm(m));
        ey2 = m(2)/(norm(m));
    
        az0 = 125.61; % I-10 heading in degrees;
        ex = sin(az0/180*pi); ey = cos(az0/180*pi);
    
        lo0 = -116.355645; la0 = 33.786583; % Center of I-10
        stlo1 = -116.320897; stla1 = 33.820251; % BF
        stlo2 = -116.305676; stla2 = 33.834889; % MCF
    
        [dist1, az1] = distance(la0, lo0, stla1, stlo1);
        [dist2, az2] = distance(la0, lo0, stla2, stlo2);
        dist1 = dist1*pi*6371/180;
        dist2 = dist2*pi*6371/180;
        x1 = dist1*sin(az1/180*pi); y1 = dist1*cos(az1/180*pi);
        x2 = dist2*sin(az2/180*pi); y2 = dist2*cos(az2/180*pi);
        
        dist_1 = -(ey1*x1-ex1*y1)/(ex*ey1-ey*ex1);
        dist_2 = -(ey2*x2-ex2*y2)/(ex*ey2-ey*ex2);
        dist = (dist_1+dist_2)/2;
        dist_bs(i, fig_idx) = dist;

        th_bf = 90 - atan(ey1/ex1)*180/pi;
        th_mcf = 90 - atan(ey2/ex2)*180/pi;
        th_bfs(i, fig_idx) = th_bf;
        th_mcfs(i, fig_idx) = th_mcf;

    end
end

save([event_dir 'mw/dists_bs.mat'], 'dist_bs', 'th_mcfs', 'th_bfs')



