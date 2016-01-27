#!/usr/bin/env python

from __future__ import print_function
from MDAnalysis import *
import numpy as np
import sys

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
  return(C)

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
  return(U)

def average(uni, rot):
  pass

def output(ary):
  l = len(ary)
  fmt = "%.3f, " * (l - 1) + "%.3f\n"
  return(fmt % tuple(ary))

#####

uni = Universe("../trjconv/ref.gro", "../trjcat/all.xtc")
ref = Universe("../../4_eq2/npt.gro")

sel_ca = uni.selectAtoms('name CA')
ref_ca = ref.selectAtoms('name CA')
sel_target = uni.selectAtoms('protein')

mat_shape = sel_target.coordinates().shape
coord_ave = np.zeros(mat_shape)
frame_tot = uni.trajectory.numframes

for ts in uni.trajectory:
  sys.stderr.write("Calc for %d/%d ... \r" % (ts.frame, frame_tot))

  mobile_coord = sel_ca.coordinates()
  ref_coord    = ref_ca.coordinates()
  target_coord = sel_target.coordinates()
  
  com_mobile = centroid(mobile_coord)
  com_ref    = centroid(ref_coord)

  trans_vect = com_ref - com_mobile

  mobile_coord -= com_mobile
  ref_coord    -= com_ref

  rotation_matrix = kabsch(mobile_coord, ref_coord)

  moved_coord = np.dot(target_coord - com_mobile, rotation_matrix) + com_mobile + trans_vect

  coord_ave += moved_coord

coord_ave /= frame_tot
sys.stderr.write("\n")

#print(coord_ave)

i = 0
for (at, coord) in zip(sel_target, coord_ave):
  i += 1
  list = (i, at.name, at.resname, at.resnum, coord[0], coord[1], coord[2])
  print("ATOM  %5d  %-3s %3s   %3d      %2.3f  %2.3f  %2.3f  1.00  0.00" % list)

