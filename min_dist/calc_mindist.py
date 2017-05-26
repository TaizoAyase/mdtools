#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import sys
from MDAnalysis import *
from itertools import product
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

#structure_file = '/path/to/structure/file'
#trajectory_file = '/path/to/trajectory/file'
structure_file = '../test_files/chainB.gro'
trajectory_file = '../test_files/chainB.xtc'

#############
if include_mainchain:
    # atom list (O/N only)
    resid_dict = {
        'ALA': ['N', 'O'],
        'CYS': ['N', 'O', 'SG'],
        'ASP': ['N', 'O', 'OD1', 'OD2'],
        'ASPP': ['N', 'O', 'OD1', 'OD2'],
        'GLU': ['N', 'O', 'OE1', 'OE2'],
        'GLUP': ['N', 'O', 'OE1', 'OE2'],
        'PHE': ['N', 'O'],
        'GLY': ['N', 'O'],
        'HIS': ['N', 'O', 'ND1', 'NE2'],  # normal His
        'HSD': ['N', 'O', 'ND1', 'NE2'],  # neutral His proton on ND1
        'HDE': ['N', 'O', 'ND1', 'NE2'],  # neutral His proton on NE2
        'HSP': ['N', 'O', 'ND1', 'NE2'],  # protonated His
        'ILE': ['N', 'O'],
        'LYS': ['N', 'O', 'NZ'],
        'LEU': ['N', 'O'],
        'MET': ['N', 'O', 'SD'],
        'ASN': ['N', 'O', 'OD1', 'ND2'],
        'PRO': ['N', 'O'],
        'GLN': ['N', 'O', 'OE1', 'NE2'],
        'ARG': ['N', 'O', 'NE', 'NH1', 'NH2'],
        'SER': ['N', 'O', 'OG'],
        'THR': ['N', 'O', 'OG1'],
        'VAL': ['N', 'O'],
        'TRP': ['N', 'O', 'NE1'],
        'TYR': ['N', 'O', 'OH'],
    }
else:
    # atom list (sidechain only)
    resid_dict = {
        'CYS': ['SG'],
        'ASP': ['OD1', 'OD2'],
        'ASPP': ['OD1', 'OD2'],
        'GLU': ['OE1', 'OE2'],
        'GLUP': ['OE1', 'OE2'],
        'HIS': ['ND1', 'NE2'],  # normal His
        'HSD': ['ND1', 'NE2'],  # neutral His proton on ND1
        'HDE': ['ND1', 'NE2'],  # neutral His proton on NE2
        'HSP': ['ND1', 'NE2'],  # protonated His
        'LYS': ['NZ'],
        'MET': ['SD'],
        'ASN': ['OD1', 'ND2'],
        'GLN': ['OE1', 'NE2'],
        'ARG': ['NE', 'NH1', 'NH2'],
        'SER': ['OG'],
        'THR': ['OG1'],
        'TRP': ['NE1'],
        'TYR': ['OH'],
    }

# load data
sys.stderr.write('Loading data ...\n')
uni = Universe(structure_file, trajectory_file)

# prepare selection
sys.stderr.write('Preparing the selection array...\n')
calc_selections = []
for lst in target_pairs:
    sel = uni.select_atoms(
        *['protein and resid %d and name CA' % i for i in lst])

    # raise error if selection was invalid
    if sel.n_atoms != 2:
        raise RuntimeError

    try:
        atomset1 = resid_dict[sel.resnames[0]]
        atomset2 = resid_dict[sel.resnames[1]]
    except KeyError as e:
        sys.stderr.write(
            '!!! Hydrogen bonding atom is not defined in selected residues. !!!')
        sys.stderr.write(
            'If you want to calc distance between main-chain, set "include_mainchain = True"')
        raise e

    atomset_iter = product(atomset1, atomset2)

    selections_tmp = []
    for atoms in atomset_iter:
        s = uni.select_atoms(
            *['protein and resid %d and name %s' % (i, a) for (i, a) in zip(lst, atoms)])
        selections_tmp.append(s)

    calc_selections.append(selections_tmp)

# calc
t_tot = uni.trajectory.n_frames
t_ary = np.zeros((t_tot, len(target_pairs)))

for ts in uni.trajectory:
    sys.stderr.write('Calc for frame %d/%d ...\r' % (ts.frame, t_tot))
    # iterate over selection_pairs
    for i, sel in enumerate(calc_selections):
        d_ary = [s.bond.length() for s in sel]
        min_d = min(d_ary)
        t_ary[ts.frame, i] = min_d

sys.stderr.write('\n')

# write-out to csv file
sys.stderr.write('Saving raw-distance file in rawdata.csv. \n')
np.savetxt('rawdata.csv', t_ary, delimiter=', ')

if not plotting:
    sys.stderr.write('Exit without plotting.\n')
    sys.exit(0)

# plotting
t = np.arange(t_tot)
for i in range(len(target_pairs)):
    plt.clf()
    title = 'mindist_%d-%d' % tuple(target_pairs[i])
    plt.plot(t, t_ary[:, i])
    plt.title(title)
    plt.xlabel('Time step')
    plt.ylabel('Minimum dist (A)')
    plt.savefig('%s.png' % title)
    sys.stderr.write('%s.png saved !\n' % title)

sys.stderr.write('Exit normally.\n')
