#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import sys
from MDAnalysis import Universe
from MDAnalysis.analysis.rms import *

# reference structure
reference_file = '../trjconv/ref.gro'
ref = Universe(reference_file)

# use not original trr, but modified (PBC) xtc
structure_file = reference_file
trajectory = '../trjcat/all.xtc'
uni = Universe(structure_file, trajectory)

R = RMSD(uni, ref, select='name CA')
print(R.run())

# write
# with open('rmsd.csv', 'w+') as f:
#    for i in range(l):
#        f.write("%8d, %3.4f\n" % (i, ary[i]))
np.savetxt('rmsd.csv', R.rmsd, delimiter=', ')

# plotting
import matplotlib as mpl
mpl.use('Agg')
import seaborn as sbn

df = R.rmsd.T
sbn.tsplot(df[2], time=df[1])
sbn.plt.title('RMSD of all CA')
sbn.plt.xlabel('Time (ps)')
sbn.plt.ylabel('RMSD (A)')
sbn.plt.savefig('rmsd.png')
