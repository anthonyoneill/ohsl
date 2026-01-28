#!/usr/bin/env python
import os
import subprocess
import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('integers', metavar='INTS', type=int, nargs='+',
                    help='Value for N = number of nodes.')

args = parser.parse_args()

if len(args.integers)==1:
    N = args.integers[0]
else:
    raise ValueError("Number of nodes not specified.")

data = np.loadtxt( "./DATA/Blasius_equation_N_"+ str(N) + ".dat" )

eta = data[:,0]
f = data[:,1]
fd = data[:,2]
fdd = data[:,3]

print("Plotting solution ...")

fig, ax = plt.subplots()
ax.plot( eta, f, "k", linewidth=2, label=r"$f$" )
ax.plot( eta, fd, "k--", linewidth=2, label=r"$f'$" )
ax.plot( eta, fdd, "k:", linewidth=2, label=r"$f''$" )
ax.set_xlim( [0, 10.0] )
ax.set_ylim( [0, 2.0] )
plt.xlabel(r"$\eta$")
title = 'Blasius equation solution with N=' + str(N) + ' nodes'
plt.title(title)
plt.legend(loc='upper right')
plt.grid()
plt.savefig( "DATA/Blasius_equation_N_" + str(N) + ".png" )
print("Plot saved as DATA/Blasius_equation_N_" + str(N) + ".png")
#plt.show()    
