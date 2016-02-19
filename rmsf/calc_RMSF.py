#!/usr/bin/env python

from __future__ import print_function
from MDAnalysis import *
import numpy as np
import sys
import matplotlib as mpl
mpl.use('Agg')
import seaborn as sbn


### input files ###

reference_file = '/path/to/your/reference/file'
topology_file = '/path/to/your/topology/file'
trajectory_file = '/path/to/your/trajectroy/file'

###################

def centroid(X):
    C = np.sum(X, axis=0) / len(X)
    return C


# implement of Kabsch algorithm for calc rotation matrix
def kabsch(P, Q):
    A = np.dot(P.T, Q)
    V, S, W = np.linalg.svd(A)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        # S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # U:rotation Matrix
    U = np.dot(V, W)
    return U


def rmsd(calc, ref):
    # rmsd for all residues => return a single rmsd value, mean for all resid
    msd = np.sum((calc - ref) ** 2)
    result = np.sqrt(msd / len(calc))
    return(result)


def output(ary):
    l = len(ary)
    fmt = "%.3f, " * (l - 1) + "%.3f\n"
    return fmt % tuple(ary)


def average(uni, ref_coord):
    # reset timestep to 0
    uni.trajectory[0]

    sel_cal = uni.select_atoms("name CA")
    sel_write = uni.select_atoms("protein")

    mat_shape = sel_write.coordinates().shape
    coord_ave = np.zeros(mat_shape)
    frame_tot = uni.trajectory.n_frames

    com_ref = centroid(ref_coord)
    ref_coord -= com_ref

    for ts in uni.trajectory:
        sys.stderr.write("Calc for %d/%d ... \r" % (ts.frame, frame_tot))

        mobile_coord = sel_cal.coordinates()
        target_coord = sel_write.coordinates()

        com_mobile = centroid(mobile_coord)
        trans_vect = com_ref - com_mobile
        mobile_coord -= com_mobile

        rotation_matrix = kabsch(mobile_coord, ref_coord)

        moved_coord = np.dot(target_coord - com_mobile, rotation_matrix) + com_mobile + trans_vect
        coord_ave += moved_coord

    coord_ave /= frame_tot
    sys.stderr.write("\n")
    return coord_ave


#####
### main ###

sys.stderr.write('Load files ...\n')
uni = Universe(topology_file, trajectory_file)
ref = Universe(reference_file)

sel_ca = uni.select_atoms('name CA')
ref_ca = ref.select_atoms('name CA')
sel_target = uni.select_atoms('protein')

ref_idx = ref_ca.indices

init_coord = sel_target.coordinates()

coord_ave = average(uni, ref_ca.coordinates())
rmsd_1 = rmsd(coord_ave, init_coord)
sys.stderr.write("1st RMSD: %2.8f\n" % rmsd_1)

rmsd_ary = [rmsd_1]

# calculate average structure
# iterate averaging until converge to 0.001 A
i = 1
while rmsd_ary[-1] > 0.001:
    i += 1
    new_coord_ave = average(uni, coord_ave[ref_idx])
    rmsd_val = rmsd(new_coord_ave, coord_ave)
    rmsd_ary.append(rmsd_val)
    sys.stderr.write("%dth RMSD: %2.4f\n" % (i, rmsd_val))
    coord_ave = new_coord_ave

sys.stderr.write(str(rmsd_ary))
sys.stderr.write("\n")

### calc RMSF ###
uni.trajectory[0]
ref_average = coord_ave[ref_idx]
frame_tot = uni.trajectory.n_frames

sumsquare = 0
com_ref = centroid(ref_average)
ref_average -= com_ref
for ts in uni.trajectory:
    sys.stderr.write("Calc for %d/%d ... \r" % (ts.frame, frame_tot))

    mobile_coord = sel_ca.coordinates()
    com_mobile = centroid(mobile_coord)
    trans_vect = com_ref - com_mobile

    mobile_coord -= com_mobile
    rotation_matrix = kabsch(mobile_coord, ref_average)
    # moved_coord = np.dot(mobile_coord - com_mobile, rotation_matrix) + com_mobile + trans_vect
    moved_coord = np.dot(mobile_coord, rotation_matrix)

    sumsquare += np.sum((moved_coord - ref_average) ** 2, axis=1)

sys.stderr.write("\n")

sumsquare /= frame_tot
rmsf = np.sqrt(sumsquare)
np.savetxt('rmsf.txt', rmsf)

# plotting
sbn.tsplot(rmsf, time=ref_idx)
sbn.plt.xlabel("Resid number")
sbn.plt.ylabel("RMSF (A)")
sbn.plt.title('RMSF for all CA atoms')
sbn.plt.savefig('rmsf.png')

# save to pdb file
f = open("average.pdb", "w+")

i = 0
for (at, coord) in zip(sel_target, coord_ave):
    i += 1
    resn = int(at.resnum) - 100 - 1
    if "CA" in at.name:
        b_fac = rmsf[resn]
    else:
        b_fac = 0
    lst = (i, at.name, at.resname, at.resnum, coord[0], coord[1], coord[2], b_fac)
    f.write("ATOM  %5d %-4s %3s   %3d    %8.3f%8.3f%8.3f  1.00  %1.2f\n" % lst)
f.close()
