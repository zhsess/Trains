%%
% The dataset contains vertical amplitude of xxx train events detected using the continuous waveforms recorded by the SoSAFZ nodal array. 
% The movement of each train event is determined by resolving the Rayleigh wavefront at two 2D subarrays. For each location of the train, 
% we estimate the median absolution value of vertical component of ground motions recorded by each node in a xx s window. Trains generally 
% generate stronger ground motions when approaching a seismograph and weaker ground motions when leaving due to geometrical spreading and attenuation.
% 
% The data are documented in Matlab structures for each train event. For example, train event T001 saved in "T001.mat" has a data structure:
% 
% T.evid = 'T001'; % event ID
% T.b = datetime(Y,M,D,H,MI,S,MS); % start time of this event
% T.e = datetime(Y,M,D,H,MI,S,MS); % end time of this event
% 
% T.train_loc = [-3.00 -2.00 -1.00 0.00 1.00 2.00 3.00]; % 1*n train location, unit: km
% T.train_loc_lat = [33.1323 33.1323 33.1323 33.1323 33.1323 33.1323 33.1323]; % 1*n latitude of train locations
% T.train_loc_lon = [-117.1323 -117.1323 -117.1323 -117.1323 -117.1323 -117.1323 -117.1323]; % 1*n longitude of train locations
% 
% % i ranges from 1 to m corresponding to all stations in the array.
% T.st{i}.lat = 33.1323; % station latitude
% T.st{i}.lon = -117.1323; % station longitude
% T.st{i}.flag = 1; % = 1: station i has good amplitude measurements; = 0: no measurements
% 
% % if T.st{i}.flag = 0; no need to assign the following variables
% T.st{i}.A = [1.00 2.00 3.00 4.00 3.00 2.00 1.00]; % 1*n amplitude, unit: count
%%

clear;clc;close all;
event_dir = '../event15/';
load([event_dir 'mw/dists.mat'])
T.evid = 'T015';
t_b = 390; dt=30; n_window = 13;
time_center = datetime(2020,03,07,17,34,40);

T.b = time_center + seconds(t_b - 600);
T.e = time_center + seconds(t_b + dt*n_window +dt - 600);

T.train_loc = dist_tr;

lon0 = -116.355645; lat0 = 33.786583; % Center of I-10
az0 = 125.61; % I-10 heading in degrees;
ex = sin(az0/180*pi); ey = cos(az0/180*pi);
T.train_loc_lat = lat0 - T.train_loc*ey/111.19;
T.train_loc_lon = lon0 - T.train_loc*ex/(111.19*cos(lat0*pi/180));

st_all = load('stations_YA.txt');
st_used = load([event_dir 'stidx_lin.txt']);
for i = 1:132
    T.st{i}.id = st_all(i, 1);
    T.st{i}.lat = st_all(i, 3);
    T.st{i}.lon = st_all(i, 2);
    if ismember(st_all(i,1), st_used)
        T.st{i}.flag = 1;
        T.st{i}.A = [];
        fnm = [event_dir 'data/R' num2str(T.st{i}.id) '.SAC'];
        [hdr,data] = load_sac(fnm);
        for n = 1:n_window
            n_1 = round((t_b + dt*(n-1))/hdr.delta);
            n_2 = round((t_b + dt*(n+1))/hdr.delta);
            A_m = median(abs(data(n_1:n_2)));
            T.st{i}.A = [T.st{i}.A A_m];
        end
    else
        T.st{i}.flag = 0;
    end
end
save([event_dir 'T15.mat'], 'T');
