#!/usr/bin/env python
import os
import subprocess
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('integers', metavar='INTS', type=int, nargs='+',
                    help='Values for J and N steps.')

args = parser.parse_args()

if len(args.integers)==2:
    J = args.integers[0]
    N = args.integers[1]
else:
    raise ValueError("J and N step values not specified.")

slowdown = 5

data = np.loadtxt( "./DATA/Heat_equation_J_"+ str(J) + "_N_" + str(N) + ".dat" )

x = data[:,0]
t = data[:,1]
u = data[:,2]

x_axis = x[0:J+1]
time = t[0:(J+1)*(N+1):J+1]
dt = time[1] - time[0]
Tmax = time[-1]

u.shape = (N+1,J+1)

files = []
print("Saving frames ...")

for i in range(N+1):
    plt.cla()
    axes = plt.gca()
    axes.set_xlim([0,1])
    axes.set_ylim([0,110])
    plt.xlabel('x')
    plt.ylabel('u(x)')
    title = 'Heat equation'
    plt.title(title)
    time_str = 't=%1.2f' %time[i]
    plt.text(0.07, 0.94, time_str,
         horizontalalignment='center',
         transform = axes.transAxes)
    plt.plot(x_axis, u[0,:], "k--", label="u(x,0)")
    plt.plot(x_axis, u[i,:], "k-", label="u(x,t)")
    plt.legend(loc='upper right')
    fname = '_tmp%03d.png' % i
    plt.savefig(fname)
    files.append(fname)

print("Saving complete.")
print("Making movie Heat_equation.avi ...")
fps = N / (Tmax * slowdown)
FNULL = open(os.devnull, 'w')  # For supressing subprocess output
subprocess.call("mencoder 'mf://_tmp*.png' -mf type=png:fps=" + str(fps) +
                " -ovc lavc -lavcopts vcodec=msmpeg4v2 -oac mp3lame -o " +
                "./DATA/Heat_equation.avi", shell=True, stdout=FNULL,
                stderr=subprocess.STDOUT)
#vcodec = wmv2 for .mpg and -oac copy -o
print("Movie complete.")

# cleanup
for fname in files:
    os.remove(fname)
