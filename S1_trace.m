clc; close all; clear
% 1 for BF, 2 for MCF

n_events = 29;
d = struct( ...
    'event01', 420, ...
    'event02', 420, ...
    'event03', 420, ...
    'event04', 420, ...
    'event05', 420, ...
    'event06', 420, ...
    'event07', 420, ...
    'event08', 420, ...
    'event09', 420, ...
    'event10', 420, ...
    'event11', 420, ...
    'event12', 420, ...
    'event13', 390, ...
    'event14', 420, ...
    'event15', 420, ...
    'event16', 420, ...
    'event17', 420, ...
    'event18', 420, ...
    'event19', 420, ...
    'event20', 420, ...
    'event21', 420, ...
    'event22', 420, ...
    'event23', 420, ...
    'event24', 420, ...
    'event25', 420, ...
    'event26', 420, ...
    'event27', 420, ...
    'event28', 420, ...
    'event29', 420 ...
    );

for i = 15:29

    clearvars -except n_events d i;

    event_id = ['event' num2str(i, '%02d')];
    disp(event_id)
    event_dir = ['../event/' event_id '/'];
    if exist([event_dir 'mw'],'dir')==0.
        mkdir([event_dir 'mw']);
    end
    t_0 = d.(event_id); 
    dt = 30;

    dist_tr_BFs = [];
    dist_tr_MCFs = [];
    e_BF = [];
    e_MCF = [];
    
    t = [];
    for fig_idx = 1:13
        
        disp(['fig_idx: ' num2str(fig_idx)])
        close all;

        t = [t t_0+dt*(fig_idx-1)];
        t_b = t_0+dt*(fig_idx-2);
        t_e = t_0+dt*fig_idx;
        n1 = round(t_b*500);
        n2 = round(t_e*500);

        % BF
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
    
                [r, lags]= xcorr(data1(n1:n2),data2(n1:n2),500,'coeff');
                [~,I] = max(r);
                dt0 = lags(I)*hdr1.delta;
    
                All_CC = [All_CC; r(I)];
                All_t = [All_t; lags(I)*hdr1.delta];
            end
        end
        save([event_dir 'mw/' num2str(fig_idx) '_BF.mat'], 'All_CC', 'All_t', 'All_x', 'All_y')
        
        % MCF
        t_b = t_0+dt*(fig_idx-2);
        t_e = t_0+dt*fig_idx;
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

        %% Plot

        figure(1)
        % Plot BF array
        load([event_dir 'mw/' num2str(fig_idx) '_BF.mat'])
        subplot(2,3,4)
        I =(All_CC>0.7)&(abs(All_t)<0.1);
        scatter(All_x(I),All_y(I),20,-All_t(I),'filled'); axis equal; hold on;
        colormap(turbo); colorbar; box on;

        b = -All_t(I);
        G = [All_x(I) All_y(I)];
        m = G\b;
        ex_BF = m(1)/(norm(m));
        ey_BF = m(2)/(norm(m));
        e_BF = [e_BF [ex_BF;ey_BF]];

        plot([-ex_BF ex_BF],[-ey_BF ey_BF],'linewidth',2);
        xlim([-0.12 0.12])
        ylim([-0.12 0.12])
        xlabel('East')
        ylabel('North')
        title('BF')
        set(gca,'fontsize',15)

        % Plot MCF array
        load([event_dir 'mw/' num2str(fig_idx) '_MCF.mat'])
        subplot(2,3,1)
        I = All_CC>0.7;
        scatter(All_x(I),All_y(I),20,-All_t(I),'filled'); axis equal; hold on;
        colormap(turbo); colorbar; box on;

        b = -All_t(I);
        G = [All_x(I) All_y(I)];
        m = G\b;
        ex_MCF = m(1)/(norm(m));
        ey_MCF = m(2)/(norm(m));
        e_MCF = [e_MCF [ex_MCF;ey_MCF]];

        plot([-ex_MCF ex_MCF],[-ey_MCF ey_MCF],'linewidth',2);
        xlim([-0.12 0.12])
        ylim([-0.12 0.12])
        ylabel('North')
        title('MCF')
        set(gca,'fontsize',15)

        % Plot geometry
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
        
        dist_tr_BF = (ey_BF*x_BF-ex_BF*y_BF)/(ex*ey_BF-ey*ex_BF);
        dist_tr_MCF = (ey_MCF*x_MCF-ex_MCF*y_MCF)/(ex*ey_MCF-ey*ex_MCF);
        dist_tr_BFs = [dist_tr_BFs dist_tr_BF];
        dist_tr_MCFs = [dist_tr_MCFs dist_tr_MCF];

        subplot(2,3,[2,3,5,6])
        plot([-5*ex 5*ex],[-5*ey 5*ey],'-.k','linewidth',5); hold on;axis equal;
        plot(0,0,'r.','MarkerSize', 50)
        plot([x_BF x_BF-10*ex_BF], [y_BF y_BF-10*ey_BF], 'r-', 'LineWidth',3); hold on
        plot(x_BF, y_BF, 'k^', 'MarkerSize', 20, 'MarkerFaceColor', 'r'); hold on
        plot(x_MCF, y_MCF, 'k^', 'MarkerSize', 20, 'MarkerFaceColor', 'b'); 
        legend('','','','BF','MCF','Location','northwest')
        box on;
        xlim([-3 6])
        ylim([-3 6])
        xlabel('East')
        ylabel('North')
        title([num2str(t(end)-600) 's Distance:' num2str(dist_tr_BF) 'km'])
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

    th_BF = 90 - atan(e_BF(2,:)./e_BF(1,:))*180/pi;
    th_MCF = 90 - atan(e_MCF(2,:)./e_MCF(1,:))*180/pi;
    save([event_dir 'mw/dists.mat'], 'dist_tr_BFs', 'dist_tr_MCFs', 'e_MCF','e_BF', 'th_BF','th_MCF', 't')

end





