#!/usr/bin/env python

from __future__ import print_function
from mdtools.rms_tools import centroid, kabsch, rmsd
import numpy as np
import sys
from MDAnalysis import *

### setup ###

reference_file = '/path/to/reference/pdb_or_gro'
structure_file = reference_file
trajectory_file = '/path/to/trajectory/file'

fit = 'protein and name CA'
calc = {
    'TM1': '101:133',
    'TM2': '137:157',  # etc...
}
selection_tmpl = 'name CA and resid %s'

rmsd_out = 'rmsd.csv'

plotting = True

##################

### main ###
# 1. Load Universe
sys.stderr.write('Loading files ...\n')
ref_uni = Universe(reference_file)
uni = Universe(structure_file, trajectory_file)


# 2. prepare selection
# 2.1. reference
sys.stderr.write('Preparing selections ...\n')

sel_ref = ref_uni.select_atoms(fit)  # for overall fitting

ref_selections = {}  # for each calc elements reference
for k, v in calc.items():
    sel = ref_uni.select_atoms(selection_tmpl % v)
    ref_selections[k] = sel

# 2.2. mobile
sel_mob = uni.select_atoms(fit)  # for overall fitting

selections = {}  # for each calc elements
for k, v in calc.items():
    sel = uni.select_atoms(selection_tmpl % v)
    selections[k] = sel


# 3. prepare for result
t_tot = uni.trajectory.n_frames
z = np.zeros(t_tot)
results = {}
for k in calc.keys():
    results[k] = z.copy()


# 4. reference coordinate sets
ref_coord = sel_ref.coordinates()
com_ref = centroid(ref_coord)
ref_coord -= com_ref
ref_coord_sets = {}
for k, v in ref_selections.items():
    coord = v.coordinates() - com_ref
    ref_coord_sets[k] = coord


# 5. start calc
for ts in uni.trajectory:
    sys.stderr.write("Calc for %d/%d ... \r" % (ts.frame, t_tot))

    mob_coord = sel_mob.coordinates()
    com_mob = centroid(mob_coord)
    mob_coord -= com_mob

    rot_matrix = kabsch(mob_coord, ref_coord)

    for k, sel in selections.items():
        cal_coord = sel.coordinates()
        cal_coord -= com_mob
        moved_coord = np.dot(cal_coord, rot_matrix)
        val = rmsd(moved_coord, ref_coord_sets[k])
        results[k][ts.frame] = val

sys.stderr.write('\n')

# 6. saving in csv file
with open(rmsd_out, 'w+') as f:
    tmpl = '%s, ' * (len(results) - 1) + '%s\n'
    string = tmpl % tuple(results.keys())
    f.write(string)
    for i in xrange(t_tot):
        tmpl = '%2.4f, ' * (len(results) - 1) + '%2.4f\n'
        d = [results[k][i] for k in results.keys()]
        string = tmpl % tuple(d)
        f.write(string)

if not plotting:
    sys.stderr.write('End program without plotting.')
    sys.exit(0)

# 7. plotting
sys.stderr.write('Start plotting ...\n')
import matplotlib as mpl
mpl.use('Agg')
import seaborn as sbn

t = np.arange(t_tot)
sbn.plt.hold(False)
for k, v in results.items():
    sbn.tsplot(v, time=t)
    sbn.plt.title('RMSD of %s; fitting %s' % (k, fit))
    sbn.plt.xlabel('Time step')
    sbn.plt.ylabel('RMSD (A)')
    sbn.plt.savefig('RMSD_%s.png' % k)

sys.stderr.write('Exitting normally.\n')
