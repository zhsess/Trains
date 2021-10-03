from obspy.core import *
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

time = UTCDateTime('2020-03-08T12:10:00Z')
event_dir = '../event26'

st_lst = open('./stations_YA.txt')
hdl = st_lst.readline()
stas = st_lst.readlines()
st_lst.close()

st_dict = {}
for i in stas:
    idx, lon, lat, ele = i.split()
    idx, lon, lat, ele = int(idx),float(lon),float(lat),float(ele)
    st_dict[idx] = [lon, lat, ele]

mon, day = time.month, time.day
for i in range(1001,9999):
    ls = glob.glob('/Volumes/Data/SoSAF/2020/*/*/R'+str(i)+'/*Z__2020' \
                        +str(mon).zfill(2)+str(day).zfill(2)+'T*__*')
    print(ls)
    if len(ls)>0:
        st = read(ls[0])
        tr = st[0].slice(t_start, t_end)
        if len(tr.data)>=500000:
            print(i)
            tr = tr.detrend(type='linear')
            tr.write(event_dir+'/data/R'+str(i)+'.SAC', format='SAC')
            
            tr = read('./data/R'+str(i)+'.SAC')[0]
            tr.stats['sac']['stlo'] = st_dict[i][0]
            tr.stats['sac']['stla'] = st_dict[i][1]
            tr.stats['sac']['stel'] = st_dict[i][2]
            tr.write(event_dir+'/data/R'+str(i)+'.SAC', format='SAC')
