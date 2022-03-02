clear;clc;close all;

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

f = 4; %Hz

%% event 15
event_dir = '../event15/';
load([event_dir 'mw/dists.mat']);
stnm = load([event_dir 'stidx_q.txt']);

n_window = 13;
t_b = 390; dt = 30;

x_tr = dists*cos(az_r-az0-90);
y_tr = dists*cos(az_r-az0-90);


for j = 1:length(stnm)
    stidx = stnm(j);
    fnm = [event_dir 'data/R' num2str(stidx) '.SAC'];
    [hdr,data] = load_sac(fnm);

    [dist, az] = distance(la0, lo0, hdr.stla, hdr.stlo);
    dist = dist*pi*6371/180;
    x = dist*sin((az-az0)/180*pi);
    y = dist*cos((az-az0)/180*pi);
    disp(x)
end
%     dist = dist*pi*6371/180;
%     for i = 6:10
%         n_1 = round((t_b + dt*(i-1))/hdr.delta);
%         n_2 = round((t_b + dt*(i+1))/hdr.delta);
%         A(i, j) = sqrt(mean(data(n_1:n_2).^2));
%     end
% end
% 
% for i = 6:10
%     plot(smooth(A(i,:),10,'moving'))
%     hold on
% end