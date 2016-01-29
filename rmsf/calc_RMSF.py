#!/usr/bin/env python

from __future__ import print_function
from MDAnalysis import *
import numpy as np
import sys
import matplotlib.pyplot as plt

#####
def fit(P, Q):
    # P:mobile, Q:reference
    P -= centroid(P)
    Q -= centroid(Q)
    rot = kabsch(P, Q)
    moved = np.dot(P, rot)
    return(moved)

def centroid(X):
    C = np.sum(X, axis = 0) / len(X)
    return C

# implement of Kabsch algorithm for calc rotation matrix
def kabsch(P, Q):
    A = np.dot(P.T, Q)
    V, S, W = np.linalg.svd(A)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        #S[-1] = -S[-1]
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

    sel_cal = uni.selectAtoms("name CA")
    sel_write = uni.selectAtoms("protein")

    mat_shape = sel_write.coordinates().shape
    coord_ave = np.zeros(mat_shape)
    frame_tot = uni.trajectory.numframes
    
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

#uni = Universe("../../npt.gro", "../../p1.trr")
#ref = Universe("../../npt.gro")
uni = Universe("ref.gro", "./trajout.xtc")
ref = Universe("ref.gro")

sel_ca = uni.selectAtoms('name CA')
ref_ca = ref.selectAtoms('name CA')
sel_target = uni.selectAtoms('protein')

ref_idx = ref_ca.indices()

init_coord = sel_target.coordinates()

coord_ave = average(uni, ref_ca.coordinates())
rmsd_1 = rmsd(coord_ave, init_coord)
print("1st RMSD: %2.8f" % rmsd_1)

rmsd_ary = [rmsd_1]

# iterate averaging until converge
i = 1
while rmsd_ary[-1] > 0.001:
    i += 1
    new_coord_ave = average(uni, coord_ave[ref_idx])
    rmsd_val = rmsd(new_coord_ave, coord_ave)
    rmsd_ary.append(rmsd_val)
    print("%dth RMSD: %2.4f" % (i, rmsd_val))
    coord_ave = new_coord_ave

print(rmsd_ary)

### calc RMSF ###
uni.trajectory[0]
ref_average = coord_ave[ref_idx]
frame_tot = uni.trajectory.numframes

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
    #moved_coord = np.dot(mobile_coord - com_mobile, rotation_matrix) + com_mobile + trans_vect
    moved_coord = np.dot(mobile_coord, rotation_matrix)

    sumsquare += np.sum((moved_coord - ref_average) ** 2, axis=1)

sys.stderr.write("\n")

sumsquare /= frame_tot
rmsf = np.sqrt(sumsquare)
np.savetxt('rmsf.txt', rmsf)
print(rmsf)

#plt.plot(ref_idx, rmsf)
#plt.xlabel("Resid number")
#plt.ylabe("RMSF (A)")
#plt.title('RMSF for all CA atoms')
#plt.savefig('rmsf.png')

f = open("average.pdb", "w+")

i = 0
for (at, coord) in zip(sel_target, coord_ave):
    i += 1
    resn = int(at.resnum) - 100 - 1
    if "CA" in at.name:
        b_fac = rmsf[resn]
    else:
        b_fac = 0
    list = (i, at.name, at.resname, at.resnum, coord[0], coord[1], coord[2], b_fac)
    f.write("ATOM  %5d  %-3s %3s   %3d      %2.3f  %2.3f  %2.3f  1.00  %1.2f\n" % list)
f.close()
