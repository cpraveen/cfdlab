from matplotlib import rcParams
rcParams['font.size'] = 14
rcParams['font.family'] = 'serif'
rcParams['figure.autolayout'] = True
rcParams['lines.linewidth'] = 2
rcParams['lines.markersize'] = 6
rcParams['axes.titlesize'] = 14
rcParams['axes.labelsize'] = 14
rcParams['text.usetex'] = False    # This will use Latex fonts (slow)

import numpy as np
import matplotlib.pyplot as plt
import os

dirs = ['tgv_kep','tgv_kg','tgv_kepec']
for dir in dirs:
    c="cd "+dir+" && grep it,t,ke,ent= log.txt |awk '{print $2, $3, $4, $5}' > global.txt"
    os.system(c)
    d = np.loadtxt(dir+'/global.txt')
    plt.figure(1)
    plt.plot(d[:,1],d[:,2]/d[0,2])
    plt.figure(2)
    plt.plot(d[:,1],d[:,3]/np.abs(d[0,3]))
    #plt.ticklabel_format(axis="y", style="sci", scilimits=(0,0))

plt.figure(1)
plt.xlabel('Time')
plt.ylabel('Relative Kinetic energy')
plt.legend(dirs)
plt.savefig('tgv_ke.pdf')

plt.figure(2)
plt.xlabel('Time')
plt.ylabel('Relative Entropy')
plt.legend(dirs)
plt.savefig('tgv_ent.pdf')
