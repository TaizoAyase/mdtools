#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sbn
import sys
from tqdm import tqdm

import mdtraj as md


### parameters ###

reference_file = '../test_files/chainB.gro'
structure_file = reference_file
trajectory = '../test_files/chainB.xtc'

target_selection = 'name CA'

iterload = True
chunk_size = 100

##################

# reference structure
ref = md.load(reference_file)

target_idx = ref.topology.select(target_selection)

# use not original trr, but modified (PBC) xtc
if not iterload:
    print('iterload OFF.')
    print('Load traj file.')
    traj = md.load(trajectory, top=ref)
    print('Calc RMSD ...')
    rmsd = md.rmsd(traj, ref, atom_indices=target_idx)
    t = traj.time / 1000 # ns
else:
    rmsd = []
    t = []
    print('iterload ON.')
    print('Calc RMSD ...')
    for chunk in tqdm(md.iterload(trajectory, top=ref, chunk=chunk_size)):
        rmsd.append(md.rmsd(chunk, ref, atom_indices=target_idx))
        t.append(chunk.time)
    rmsd = np.concatenate(rmsd)
    t = np.concatenate(t) / 1000 # ns

# write
# with open('rmsd.csv', 'w+') as f:
#    for i in range(l):
#        f.write("%8d, %3.4f\n" % (i, ary[i]))
np.savetxt('rmsd.dat', rmsd)

# plotting
print('Plotting.')
plt.plot(t, rmsd)
plt.title('RMSD of "%s"' % target_selection)
plt.xlabel('Time (ps)')
plt.ylabel('RMSD (nm)')
plt.savefig('rmsd.png')
