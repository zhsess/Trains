import glob

evdir = '../event22/'
sac_ls = glob.glob(evdir+'data/*.SAC')
idx_ls = [int(i[-8:-4]) for i in sac_ls]
idx_ls.sort()

# make stidx_linear
stidx_lin = open(evdir+'stidx_lin.txt', 'w')
for i in idx_ls:
    if 1000<i<1125:
        stidx_lin.write(str(i)+'\n')
stidx_lin.close()

# make stidx_BF
stidx_BF = open(evdir+'stidx_BF.txt', 'w')
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
stidx_MCF = open(evdir+'stidx_MCF.txt', 'w')
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
