#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import sys
import mdtools.constants as const
import mdtraj as md
import itertools
from tqdm import tqdm

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns


### setup ###

# example...
target_pairs = [
    [278, 348],
    [165, 280],
    [128, 195],
    [128, 192],
    [271, 350],
    [271, 348],
]

include_mainchain = True
plotting = True

structure_file = '/path/to/structure/file'
trajectory_file = '/path/to/trajectory/file'

chunk_size = 100

#############

if include_mainchain:
    # atom list (including O and N of backbone)
    resid_dict = const.RESID_DICT
else:
    # atom list (sidechain only)
    resid_dict = const.RESID_DICT_SIDECHAIN

# load the structure file data
sys.stderr.write('Loading data ...\n')
ref = md.load(structure_file)

# prepare selection
sys.stderr.write('Preparing the selection array...\n')
calc_selections = []
for lst in target_pairs:
    sel = ref.topology.select("name CA and (residue %d or residue %d)" % tuple(lst))
    ca_atom1 = ref.topology.atom(sel[0])
    ca_atom2 = ref.topology.atom(sel[1])

    # raise error if selection was invalid
    if len(sel) != 2:
        raise RuntimeError

    # get atom names list from the dictionary
    try:
        atomset1 = resid_dict[ca_atom1.residue.name]
        atomset2 = resid_dict[ca_atom2.residue.name]
    except KeyError as e:
        sys.stderr.write('!!! Hydrogen bonding atom is not defined in selected residues. !!!')
        sys.stderr.write('If you want to calc distance between main-chain, set "include_mainchain = True"')
        raise e

    atomset_iter = itertools.product(atomset1, atomset2)

    selections_tmp = []
    for atoms in atomset_iter:
        inputs = (ca_atom1.residue.resSeq, atoms[0], ca_atom2.residue.resSeq, atoms[1])
        s = ref.topology.select('(residue %d and name %s) or (residue %d and name %s)' % inputs)
        selections_tmp.append(s)

    calc_selections.append(selections_tmp)

# calc
sys.stderr.write('Start distance calculation...\n')
results = [[] for _ in calc_selections]
itertraj = md.iterload(trajectory_file, top=ref, chunk=chunk_size)
for chunk in tqdm(itertraj):
    for i, sel in enumerate(calc_selections):
        coor = chunk.xyz[:, sel, :] # 4-dimensional ary: (time, #of selection pair, 2atoms, xyz)
        dist_vectors = coor[:, :, 0, :] - coor[:, :, 1, :] # 3-dimensional ary: (time, #of pair, xyz)
        dists = np.linalg.norm(dist_vectors, axis=2) # calc norm over xyz
        min_dist = np.min(dists, axis=1) # select the minimum value in the pairs -> 1-dimensional ary
        results[i].append(min_dist)

# convert into 2-dim numpy array
t_ary = np.array([np.concatenate(e) for e in results]).T

sys.stderr.write('\n')

# write-out to csv file
sys.stderr.write('Saving raw-distance file in rawdata.csv. \n')
np.savetxt('rawdata.csv', t_ary, delimiter=', ')

if not plotting:
    sys.stderr.write('Exit without plotting.\n')
    sys.exit(0)

# plotting
n_frames = t_ary.shape[0]
t = np.arange(n_frames)
for i in range(len(target_pairs)):
    plt.clf()
    title = 'mindist_%d-%d' % tuple(target_pairs[i])
    plt.plot(t, t_ary[:, i])
    plt.title(title)
    plt.xlabel('Time step')
    plt.ylabel('Minimum dist (nm)')
    plt.savefig('%s.png' % title)
    sys.stderr.write('%s.png saved !\n' % title)

sys.stderr.write('Exit normally.\n')
