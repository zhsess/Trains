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
load([event_dir 'velo.mat'])
stnm = load([event_dir 'stidx_lin.txt']);

n_window = 13;
t_b = 420; dt = 30;

x_tr = dist_tr*cos(az0+90-az_r);
y_tr = dist_tr*sin(az0+90-az_r);

dist_st = [];
y_st = [];
for j = 1:length(stnm)
    stidx = stnm(j);
    fnm = [event_dir 'data/R' num2str(stidx) '.SAC'];
    [hdr,data] = load_sac(fnm);

    [dist, az] = distance(la0, lo0, hdr.stla, hdr.stlo);
    dist = dist*pi*6371/180;
    
    x = dist*sin((az-az0)/180*pi);
    y = dist*cos((az-az0)/180*pi);

    dist_st(j) = dist;
    y_st(j) = y;
    
    for i = 1:n_window
        n_1 = round((t_b + dt*(i-2))/hdr.delta);
        n_2 = round((t_b + dt*i)/hdr.delta);
        A(i, j) = sqrt(mean(abs(data(n_1:n_2)).^2));
        cost(i, j) = (y-y_tr(i))/sqrt((x-x_tr(i))^2+(y-y_tr(i))^2);
    end
end

for i = 1:n_window
    A(i,:) = smooth(A(i,:), 20, 'moving');
end

n_b = 25; n_e = 80;
d_b = 8; d_e = 12;
len = length(stnm(n_b:n_e));
lnA = zeros((d_e-d_b+1)*(len-1), 1);
G = zeros((d_e-d_b+1)*(len-1), len);
for i = 1:d_e-d_b+1
    for j = 1:len-1
        lnA(j+(i-1)*len) = -log(A(i+d_b-1,j+n_b+1)/A(i+d_b-1,j+n_b))/(pi*f);
        for k = 1:j
            G(j+(i-1)*len, k) = 1/cost(i+d_b-1,j+31) - 1/cost(i+d_b-1,j+n_b);
        end
        G(j+(i-1)*len, j+1) = 1/cost(i+d_b-1,j+31);
    end
end
t = (G\lnA)';

y_st = y_st - y_st(1);
v_section = interp1(dist_ctr, V0_ctr, y_st, 'cubic');
v_section_std = interp1(dist_ctr, V0_std, y_st, 'cubic');
d_y = y_st(n_b+1:n_e) - y_st(n_b:n_e-1);

Q = smooth(1000*d_y./(t(2:end).*v_section(n_b:n_e-1)), 5, 'moving');
Q_ub = smooth(1000*d_y./(t(2:end).*(v_section(n_b:n_e-1)-v_section_std(n_b:n_e-1))), 5, 'moving');
Q_lb = smooth(1000*d_y./(t(2:end).*(v_section(n_b:n_e-1)+v_section_std(n_b:n_e-1))), 5, 'moving');

plot(y_st(n_b:n_e-1), Q, 'black', 'linewidth', 2)
hold on
area = fill([y_st(n_b:n_e-1) fliplr(y_st(n_b:n_e-1))], [Q_ub' fliplr(Q_lb')], 'r', 'edgealpha', '0', 'facealpha', '.2');
xlabel('Distance to R1001 (km)')
ylabel('Q-value')
set(gca,'FontSize',20)





% len = length(stnm);
% lnA = zeros(4*(len-1), 1);
% G = zeros(4*(len-1), len);
% for i = 1:4
%     for j = 1:len-1
%         lnA(j+(i-1)*len) = log((A(8+i,j+1)*A(8,j))/(A(8+i,j)*A(8,j+1)))/(pi*f);
%         for k = 1:j
%             G(j+(i-1)*len, k) = 1/cost(8+i,j+1) - 1/cost(8+i,j);
%         end
%         G(j+(i-1)*len, j+1) = 1/cost(8+i,j+1) - 1;
%     end
% end
% t = G\lnA;
% plot(t)

% Q = l/(ct)