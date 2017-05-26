#!/usr/bin/env python

from __future__ import print_function, absolute_import
from mdtools.rms_tools import fit, centroid, kabsch, rmsd
from tqdm import tqdm
import mdtraj as md
import numpy as np
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


### input files ###

reference_file = '/path/to/your/reference/file'
topology_file = '/path/to/your/topology/file'
trajectory_file = '/path/to/your/trajectroy/file'

chunk_size = 100

###################

def output(ary):
    l = len(ary)
    fmt = "%.3f, " * (l - 1) + "%.3f\n"
    return fmt % tuple(ary)


def average(gen, fit_idx, write_idx, ref_coord):
    # gen is the iterload generator object

    coord_ave = np.zeros((len(write_idx), 3))
    frame_tot = 0

    com_ref = centroid(ref_coord)
    ref_coord -= com_ref

    for chunk in tqdm(gen):
        for j in range(chunk.n_frames):
            mobile_coord = chunk.xyz[j, fit_idx, :]
            target_coord = chunk.xyz[j, write_idx, :]

            com_mobile = centroid(mobile_coord)
            trans_vect = com_ref - com_mobile
            mobile_coord -= com_mobile

            rotation_matrix = kabsch(mobile_coord, ref_coord)

            moved_coord = np.dot(target_coord - com_mobile, rotation_matrix) + com_ref
            coord_ave += moved_coord

        frame_tot += chunk.n_frames

    coord_ave /= frame_tot
    return coord_ave


#####
### main ###

sys.stderr.write('Load files ...\n')
ref = md.load(reference_file)

sel_ca = ref.topology.select('name CA')
sel_all = ref.topology.select('protein')

ref_idx = np.in1d(sel_all, sel_ca)

init_frame = md.load_frame(trajectory_file, 0, top=topology_file)
init_coord = init_frame.xyz[0, sel_all]

ref_ca_coord = ref.xyz[0, ref_idx, :]

itertraj = md.iterload(trajectory_file, top=topology_file, chunk=chunk_size)
coord_ave = average(itertraj, sel_ca, sel_all, ref_ca_coord)
rmsd_1 = rmsd(coord_ave, init_coord) * 10.0 # nm to A
sys.stderr.write("1st RMSD: %2.8f\n" % rmsd_1)

rmsd_ary = [rmsd_1]

# calculate average structure
# iterate averaging until converge to 0.001 A
i = 1
while rmsd_ary[-1] > 0.001:
    i += 1
    itertraj = md.iterload(trajectory_file, top=topology_file, chunk=chunk_size)
    new_coord_ave = average(itertraj, sel_ca, sel_all, coord_ave[ref_idx, :])
    rmsd_val = rmsd(new_coord_ave, coord_ave) * 10.0 # nm to A
    rmsd_ary.append(rmsd_val)
    sys.stderr.write("%dth RMSD: %2.4f\n" % (i, rmsd_val))
    coord_ave = new_coord_ave

sys.stderr.write(str(rmsd_ary))
sys.stderr.write("\n")
sys.stderr.write("Finish averaging.\n")

### calc RMSF ###
ref_average = coord_ave[ref_idx, :]
ref_average -= centroid(ref_average)

frame_tot = 0
sumsquare = 0

itertraj = md.iterload(trajectory_file, top=topology_file, chunk=chunk_size)
for chunk in tqdm(itertraj):
    for j in range(chunk.n_frames):
        mobile_coord = chunk.xyz[j, sel_ca, :]
        com_mobile = centroid(mobile_coord)

        mobile_coord -= com_mobile
        rotation_matrix = kabsch(mobile_coord, ref_average)
        moved_coord = np.dot(mobile_coord, rotation_matrix)

        sumsquare += np.sum((moved_coord - ref_average) ** 2, axis=1)

    frame_tot += chunk.n_frames

sumsquare /= frame_tot
rmsf = np.sqrt(sumsquare)
np.savetxt('rmsf.txt', rmsf)

# plotting
plt.plot(rmsf)
plt.xlabel("Resid number")
plt.ylabel("RMSF (nm)")
plt.title('RMSF for all CA atoms')
plt.savefig('rmsf.png')

# save to pdb file
# set RMSF to b-factor col.
idx = ref.topology.select('all')
b_fac = []
j = 0
for i in idx:
    if i in sel_ca:
        b_fac.append(rmsf[j] * 10.0)
        j += 1
    else:
        b_fac.append(0)
ref.save_pdb('average.pdb', bfactors=b_fac)
