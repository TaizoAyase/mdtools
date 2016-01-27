#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import math
import sys
from MDAnalysis import *

'''
program for calc Ramachandran distribution for each MD frame

               (r4)
         v2    / v3
    (r2) --> (r3)
v1 /
(r1)
'''


def calc_dihedral(v1, v2, v3):
    ax1 = np.cross(v3, v2)
    ax2 = np.cross(v2, v1)

    # calc angle
    n1 = np.linalg.norm(ax1, axis=1)
    n2 = np.linalg.norm(ax2, axis=1)
    cos_th = np.sum(ax1 * ax2, axis=1) / n1 / n2
    # theta = np.arccos(cos_th) / np.pi * 180.0

    # replace outlier cos value
    # to avoid numpy.arccos value error
    # e.g., error when cos_th = 1.00000012, theta = nan
    cos_th[cos_th > 1.0] = 1.0
    cos_th[cos_th < -1.0] = -1.0
    theta = np.arccos(cos_th) / np.pi * 180.0
    ### debug ###
    # if np.sum(theta != theta): # check nan
    #  print("cosine value:", output(cos_th, "%.8f"))
    #  print("theta  value:", output(theta, "%.1f"))

    # add sign
    v = np.cross(ax1, ax2)
    tf_mat = np.sum(v * v2, axis=1) > 0
    dot_sign = (tf_mat - 1) + tf_mat
    theta *= dot_sign
    return(theta)


def save(theta, out):
    # format = "%.1f, " * (len(theta) - 1) + "%.1f\n"
    # np.savetxt(out, theta, fmt = "%.1f ", newline = "")
    # np.savetxt(out, theta, delimiter = ", ", newline = ", ")
    # out.write("\n")
    str = output(theta)
    out.write(str)
    return


def output(ary, format="%.2f"):
    l = len(ary)
    fmt = (format + ", ") * (l - 1) + format + "\n"
    return(fmt % tuple(ary))


def calc(uni):
    res_start = 1
    res_end = 44

    sel_C = uni.selectAtoms(*["resid %d and name C" %
                              x for x in xrange(res_start, res_end + 1)])
    sel_N = uni.selectAtoms(*["resid %d and name N" %
                              x for x in xrange(res_start, res_end + 1)])
    sel_CA = uni.selectAtoms(
        *["resid %d and name CA" % x for x in xrange(res_start, res_end + 1)])

    phiout = open("phi_time.csv", "w")
    psiout = open("psi_time.csv", "w")
    i = 0

    for ts in uni.trajectory:
        sys.stderr.write("Processing frame %d/%d...\r" %
                         (ts.frame, len(uni.trajectory)))

        mat_C = sel_C.coordinates()
        mat_N = sel_N.coordinates()
        mat_CA = sel_CA.coordinates()

        nres = mat_C.shape[0]

        # calc dihedral:phi
        v1 = mat_C[1:] - mat_CA[1:]
        v2 = mat_CA[1:] - mat_N[1:]
        v3 = mat_N[1:] - mat_C[:-1]

        theta = calc_dihedral(v1, v2, v3)
        save(theta, phiout)

        ###
        # calc dihedral:psi
        v1 = mat_N[:-1] - mat_CA[:-1]
        v2 = mat_CA[:-1] - mat_C[:-1]
        v3 = mat_C[:-1] - mat_N[1:]

        theta = calc_dihedral(v1, v2, v3)
        save(theta, psiout)

    phiout.close()
    psiout.close()
    sys.stderr.write('\nDone.\n')
    return

# main
if __name__ == '__main__':
    if len(sys.argv) < 3:
        usage = "usage: python %s reference trajectory" % sys.argv[0]
        print(usage)
        sys.exit(1)

    ref = sys.argv[1]
    xtc = sys.argv[2]
    trj = Universe(ref, xtc)
    calc(trj)
