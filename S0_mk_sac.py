from obspy.core import *
import matplotlib.pyplot as plt
import numpy as np
import os
import glob

# load time table
with open('../pick/train_gmt.txt', 'r') as time_table:
    records = time_table.readlines()
records = [UTCDateTime(i) for i in records]

# load stations
with open('../info/stations_YA.txt', 'r') as station_list:
    stations = station_list.readlines()
st_dict = {}
for i in stations:
    idx, lon, lat, ele = i.split()
    idx, lon, lat, ele = int(idx),float(lon),float(lat),float(ele)
    st_dict[idx] = [lon, lat, ele]

for i in range(len(records)):
    print('Processing event: ' + str(i+1).zfill(2))
    time = records[i]
    event_dir = '../event/event' + str(i+1).zfill(2)
    if not os.path.exists(event_dir + '/data'):
        os.makedirs(event_dir + '/data')


    mon, day = time.month, time.day
    for i in range(1001,9999):
        ls = glob.glob('/Volumes/Hao/SoSAF_data/*/YA/R'+str(i)+'/*Z__2020'+str(mon).zfill(2)+str(day).zfill(2)+'T*__*')
        t_start, t_end = time-600, time+600
        if len(ls)>0:
            st = read(ls[0])
            tr = st[0].slice(t_start, t_end)
            if len(tr.data)>=500000:
                tr = tr.detrend(type='linear').filter('bandpass',freqmin=1,freqmax=5)
                tr.write(event_dir+'/data/R'+str(i)+'.SAC', format='SAC')
                
                tr = read(event_dir+'/data/R'+str(i)+'.SAC')[0]
                tr.stats['sac']['stlo'] = st_dict[i][0]
                tr.stats['sac']['stla'] = st_dict[i][1]
                tr.stats['sac']['stel'] = st_dict[i][2]
                tr.write(event_dir+'/data/R'+str(i)+'.SAC', format='SAC')

    sac_ls = glob.glob(event_dir+'/data/*.SAC')
    idx_ls = [int(i[-8:-4]) for i in sac_ls]
    idx_ls.sort()

    # make stidx_linear
    stidx_lin = open(event_dir+'/stidx_lin.txt', 'w')
    for i in idx_ls:
        if 1000<i<1125:
            stidx_lin.write(str(i)+'\n')
    stidx_lin.close()

    # make stidx_BF
    stidx_BF = open(event_dir+'/stidx_BF.txt', 'w')
    k=0
    for idx in idx_ls:
        if 1008<idx<1030:
            k += 1
            stidx_BF.write(str(idx)+' ')
    for j in range(21-k):
        stidx_BF.write('nan ')
    stidx_BF.write('\n')
    for i in range(2,6):
        k = 0
        for idx in idx_ls:
            if i*1000<idx<i*1000+30:
                k += 1
                stidx_BF.write(str(idx)+' ')
        for j in range(21-k):
            stidx_BF.write('nan ')
        stidx_BF.write('\n')
    stidx_BF.close()

    # make stidx_MCF
    stidx_MCF = open(event_dir+'/stidx_MCF.txt', 'w')
    k=0
    for idx in idx_ls:
        if 1089<idx<1110:
            k += 1
            stidx_MCF.write(str(idx)+' ')
    for j in range(21-k):
        stidx_MCF.write('nan ')
    stidx_MCF.write('\n')
    for i in range(6,10):
        k = 0
        for idx in idx_ls:
            if i*1000<idx<i*1000+30:
                k += 1
                stidx_MCF.write(str(idx)+' ')
        for j in range(21-k):
            stidx_MCF.write('nan ')
        stidx_MCF.write('\n')
    stidx_MCF.close()
