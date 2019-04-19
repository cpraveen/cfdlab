import numpy as np
import argparse

from matplotlib import rcParams
rcParams['figure.autolayout'] = True

import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('-x', type=int, help='x axis', default=1)
parser.add_argument('-y', type=int, help='y axis', default=2)
parser.add_argument('-f', type=str, help='files', nargs='+', required=True)
parser.add_argument('-m', type=int, help='mark every', default=0)
args = parser.parse_args()

ix = args.x - 1
iy = args.y - 1

markers = ['s','*','v','o','d']

plt.figure()
m = 0
for f in args.f:
    d = np.loadtxt(f)
    if args.m > 0:
        plt.plot(d[:,ix], d[:,iy], marker=markers[m], markevery=args.m)
    else:
        plt.plot(d[:,ix], d[:,iy])
    m += 1

plt.title('Showing column '+str(ix+1)+' vs '+str(iy+1))
plt.legend(args.f)
plt.show()
