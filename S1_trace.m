clc; close all; clear
% 1 for BF, 2 for MCF

event_dir = '../event15/';
if exist([event_dir 'mw'],'dir')==0.
    mkdir([event_dir 'mw']);
end
% event15
t_0 = 420; dt = 30;

t = [];
dist_tr = [];
e_tr1 = [];
e_tr2 = [];
for fig_idx = 1:13
    % BF
    clc; close all;
    clearvars -except t_0 dt fig_idx event_dir dist_tr e_tr1 e_tr2 t; 
    
    disp(['BF: fig_idx=' num2str(fig_idx)])
    t = [t t_0+dt*(fig_idx-1)];
    t_b = t_0+dt*(fig_idx-2);
    t_e = t_0+dt*fig_idx;
    stnm = load([event_dir 'stidx_BF.txt'])';
    stnm = stnm(:); I = ~isnan(stnm);
    stnm = stnm(I);
    stnm = sort(stnm);

    All_x = [];
    All_y = [];
    All_t = [];
    All_CC = [];
    for l = 1:length(stnm)-1
        stidx0 = stnm(l);
        disp(l)
        for j = l+1:length(stnm)
            stidx1 = stnm(j);
            fnm1 = [event_dir 'data/R' num2str(stidx0) '.SAC'];
            fnm2 = [event_dir 'data/R' num2str(stidx1) '.SAC'];

            [hdr1,data1] = load_sac(fnm1);
            [hdr2,data2] = load_sac(fnm2);

            y = (hdr2.stla - hdr1.stla)/180*pi*6371;
            x = (hdr2.stlo - hdr1.stlo)*cos(hdr1.stla/180*pi)/180*pi*6371;
            All_x = [All_x; x];
            All_y = [All_y; y];

            n1 = round(t_b/hdr1.delta);
            n2 = round(t_e/hdr1.delta);

            [r, lags]= xcorr(data1(n1:n2),data2(n1:n2),500,'coeff');
            [~,I] = max(r);
            dt0 = lags(I)*hdr1.delta;

            All_CC = [All_CC; r(I)];
            All_t = [All_t; lags(I)*hdr1.delta];
        end
    end
    save([event_dir 'mw/' num2str(fig_idx) '_BF.mat'], 'All_CC', 'All_t', 'All_x', 'All_y')
    
    % MCF
    clearvars -except t_0 dt fig_idx event_dir dist_tr e_tr1 e_tr2 t; clc;
    
    disp(['MCF: fig_idx=' num2str(fig_idx)])
    t_b = t_0+dt*(fig_idx-1);
    t_e = t_0+dt*(fig_idx+1);
    stnm = load([event_dir 'stidx_MCF.txt'])';
    stnm = stnm(:); I = ~isnan(stnm);
    stnm = stnm(I);
    stnm = sort(stnm);

    All_x = [];
    All_y = [];
    All_t = [];
    All_CC = [];
    for l = 1:length(stnm)-1
        stidx0 = stnm(l);
        disp(l)
        for j = l+1:length(stnm)
            stidx1 = stnm(j);
            fnm1 = [event_dir 'data/R' num2str(stidx0) '.SAC'];
            fnm2 = [event_dir 'data/R' num2str(stidx1) '.SAC'];

            [hdr1,data1] = load_sac(fnm1);
            [hdr2,data2] = load_sac(fnm2);

            y = (hdr2.stla - hdr1.stla)/180*pi*6371;
            x = (hdr2.stlo - hdr1.stlo)*cos(hdr1.stla/180*pi)/180*pi*6371;
            All_x = [All_x; x];
            All_y = [All_y; y];

            n1 = round(t_b/hdr1.delta);
            n2 = round(t_e/hdr1.delta);

            [r, lags]= xcorr(data1(n1:n2),data2(n1:n2),500,'coeff');
            [~,I] = max(r);
            dt0 = lags(I)*hdr1.delta;

            All_CC = [All_CC; r(I)];
            All_t = [All_t; lags(I)*hdr1.delta];
        end
    end
    save([event_dir 'mw/' num2str(fig_idx) '_MCF.mat'], 'All_CC', 'All_t', 'All_x', 'All_y')
    %%
    clearvars -except t_0 dt fig_idx event_dir dist_tr e_tr1 e_tr2 t; clc;

    figure(1)
    load([event_dir 'mw/' num2str(fig_idx) '_BF.mat'])
    
    subplot(2,3,4)
    I =(All_CC>0.7)&(abs(All_t)<0.1);
    scatter(All_x(I),All_y(I),20,-All_t(I),'filled'); axis equal; hold on;
    colormap(turbo); colorbar; box on;
    b = -All_t(I);
    G = [All_x(I) All_y(I)];
    m = G\b;
    ex1 = m(1)/(norm(m));
    ey1 = m(2)/(norm(m));
    plot([-ex1 ex1],[-ey1 ey1],'linewidth',2);
    xlim([-0.12 0.12])
    ylim([-0.12 0.12])
    xlabel('East')
    ylabel('North')
    title('BF')
    set(gca,'fontsize',15)

    load([event_dir 'mw/' num2str(fig_idx) '_MCF.mat'])
    subplot(2,3,1)
    I = All_CC>0.7;
    scatter(All_x(I),All_y(I),20,-All_t(I),'filled'); axis equal; hold on;
    colormap(turbo); colorbar; box on;
    b = -All_t(I);
    G = [All_x(I) All_y(I)];
    m = G\b;
    ex2 = m(1)/(norm(m));
    ey2 = m(2)/(norm(m));
    plot([-ex2 ex2],[-ey2 ey2],'linewidth',2);
    xlim([-0.12 0.12])
    ylim([-0.12 0.12])
    ylabel('North')
    title('MCF')
    set(gca,'fontsize',15)

    az0 = 125.61; % I-10 heading in degrees;
    ex = sin(az0/180*pi); ey = cos(az0/180*pi);

    lo0 = -116.3549; la0 = 33.7862; % Center of I-10
    stlo1 = -116.320897; stla1 = 33.820251; % BF
    stlo2 = -116.305676; stla2 = 33.834889; % MCF

    [dist1, az1] = distance(la0, lo0, stla1, stlo1);
    [dist2, az2] = distance(la0, lo0, stla2, stlo2);
    dist1 = dist1*pi*6371/180;
    dist2 = dist2*pi*6371/180;
    x1 = dist1*sin(az1/180*pi); y1 = dist1*cos(az1/180*pi);
    x2 = dist2*sin(az2/180*pi); y2 = dist2*cos(az2/180*pi);
    
    dist_1 = (ey1*x1-ex1*y1)/(ex*ey1-ey*ex1);
    dist_2 = (ey2*x2-ex2*y2)/(ex*ey2-ey*ex2);
    dist = (dist_1+dist_2)/2;
    dist_tr = [dist_tr dist];
    e_tr1 = [e_tr1 [ex1;ey1]];
    e_tr2 = [e_tr2 [ex2;ey2]];

    subplot(2,3,[2,3,5,6])
    plot([-5*ex 5*ex],[-5*ey 5*ey],'-.k','linewidth',5); hold on;axis equal;
    plot(0,0,'r.','MarkerSize', 50)
    plot([x1 x1-10*ex1], [y1 y1-10*ey1], 'r-', 'LineWidth',3); hold on
    % plot([x2 x2-10*ex2], [y2 y2-10*ey2], 'r-', 'LineWidth',3); hold on
    plot(x1, y1, 'k^', 'MarkerSize', 20, 'MarkerFaceColor', 'r'); hold on
    plot(x2, y2, 'k^', 'MarkerSize', 20, 'MarkerFaceColor', 'b'); 
    legend('','','','Array1','Array 2','Location','northwest')
    box on;
    xlim([-3 6])
    ylim([-3 6])
    xlabel('East')
    ylabel('North')
    title([num2str(t(end)) 's Distance:' num2str(dist) 'km'])
    set(gca,'fontsize',15)
    
    F = getframe(gcf);
    I=frame2im(F);
    [I,map] = rgb2ind(I,256);
    
    if fig_idx == 1
        imwrite(I,map,[event_dir 'mv.gif'],'gif','Loopcount',inf,'DelayTime',0.5);
    else
        imwrite(I,map,[event_dir 'mv.gif'],'gif','WriteMode','append','DelayTime',0.5);
    end
end

e_MCF = e_tr1;
e_BF = e_tr2;
save([event_dir 'mw/dists.mat'], 'dist_tr','e_MCF','e_BF', 't')



