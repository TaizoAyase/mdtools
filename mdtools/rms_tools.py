#!/usr/bin/env python

from __future__ import print_function
import numpy as np
from MDAnalysis import *


def fit(P, Q):
    # P:mobile, Q:reference
    P -= centroid(P)
    Q -= centroid(Q)
    rot = kabsch(P, Q)
    moved = np.dot(P, rot)
    return(moved)


def centroid(X):
    C = np.sum(X, axis=0) / len(X)
    return(C)


# implementation of Kabsch algorithm for calc rotation matrix

def kabsch(P, Q):
    A = np.dot(P.T, Q)
    V, S, W = np.linalg.svd(A)
    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

    if d:
        # S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    # U:rotation Matrix
    U = np.dot(V, W)
    return(U)


def rmsd(calc, ref):
    # rmsd for all residues => return a single rmsd value, mean for all
    msd = np.sum((calc - ref) ** 2)
    result = np.sqrt(msd / len(calc))
    return(result)


def superpose(coord, ref):
    ref -= centroid(ref)
    coord -= centroid(coord)
    rot_matrix = kabsch(coord, ref)
    moved_coord = np.dot(coord, rot_matrix)
    return moved_coord
