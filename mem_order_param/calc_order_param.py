#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import sys
from MDAnalysis import *

'''
Script for calculate (first) membrane order
This script is originally written for Charmm36 POPC/POPE lipids
Edit the atom selections for other types of lipid

ref: https://www.mpibpc.mpg.de/276911/Vermeer_EBJ_2007.pdf 
'''

### set up ###

structure_file = '/path/to/structure/file'
trajectory_files = '/path/to/trajectory/files/ or list_of_trajectory_files'

output1 = "ole_par.csv"
output2 = "par_par.csv"

##############


def calc(uni, file1, file2):

    print("Select hydrogen and carbon atoms of lipid tail ...")

    # selection for oleoyl tail
    # number of carbons
    r = np.arange(2, 19)
    i_range = r[(r!= 9) * (r!= 10)] # remove i = 9, 10; these two form carbon double bond
    ole_atom_num = len(i_range)
    ole_carbon_sel = uni.select_atoms(*["resname POPE and name C2%d" % i for i in i_range])
    ole_hR_sel     = uni.select_atoms(*["resname POPE and name H%dR" % i for i in i_range])
    ole_hS_sel     = uni.select_atoms(*["resname POPE and name H%dS" % i for i in i_range])

    # selection for parmitoil tail
    # number of carbons
    r = np.arange(2, 17)
    par_atom_num = len(r)
    par_carbon_sel = uni.select_atoms(*["resname POPE and name C3%d" % i for i in r])
    par_hX_sel     = uni.select_atoms(*["resname POPE and name H%dX" % i for i in r])
    par_hY_sel     = uni.select_atoms(*["resname POPE and name H%dY" % i for i in r])

    out1 = open(file1, 'w+')
    out2 = open(file2, 'w+')

    sys.stderr.write("Counting residue number ...\n")
    resnum = uni.select_atoms("resname POPE").numberOfResidues()

    frame_total = uni.trajectory.n_frames
    for ts in uni.trajectory:
        sys.stderr.write("Calc for frame %d/%d ...\r" % (ts.frame, frame_total))

        # calc for oleoyl
        ole_carbon_coor = ole_carbon_sel.coordinates()
        ole_hR_coor     = ole_hR_sel.coordinates()
        ole_hS_coor     = ole_hS_sel.coordinates()
        ole_param = calc_order_param(ole_carbon_coor, ole_hR_coor, ole_hS_coor)
        #print(ole_param.shape)
        averaged = average_all_resid(ole_param, ole_atom_num, resnum)
        out1.write(output(averaged))

        # calc for parmitoil
        par_carbon_coor = par_carbon_sel.coordinates()
        par_hX_coor     = par_hX_sel.coordinates()
        par_hY_coor     = par_hY_sel.coordinates()
        par_param = calc_order_param(par_carbon_coor, par_hX_coor, par_hY_coor)
        #print(par_param.shape)
        averaged = average_all_resid(par_param, par_atom_num, resnum)
        out2.write(output(averaged))

    out1.close()
    out2.close()
    sys.stderr.write('\nDone.\n')
    return()

def calc_order_param(carbon, hydrogen1, hydrogen2):
    v1 = hydrogen1 - carbon
    v2 = hydrogen2 - carbon

    cd1 = v1[:, 2] # dot product to z-axis
    cd2 = v2[:, 2] 
    #cd_r1 = np.sqrt(np.sum(cd1 ** 2, axis = -1))
    #cd_r2 = np.sqrt(np.sum(cd2 ** 2, axis = -1))
    v1_norm = np.linalg.norm(v1, axis = 1)
    v2_norm = np.linalg.norm(v2, axis = 1)

    #cos1 = cd1 / cd_r1
    #cos2 = cd2 / cd_r2
    cos1 = cd1 / v1_norm
    cos2 = cd2 / v2_norm

    S_cd1 = -0.5 * (3. * np.square(cos1) - 1)
    S_cd2 = -0.5 * (3. * np.square(cos2) - 1)
    S_cd = (S_cd1 + S_cd2) / 2
    return(S_cd)

# average order params along all residues
def average_all_resid(ary, num_atoms, resnum):
    new_ary = ary.reshape((num_atoms, resnum))
    averaged_ary = np.average(new_ary, axis = 1)
    return(averaged_ary)

def output(ary):
    l = len(ary)
    fmt = "%.3f, " * (l - 1) + "%.3f\n"
    return(fmt % tuple(ary))


#############################################
  
if __name__ == '__main__':
    print("Loading universe ...")
    uni = Universe(structure_file, trajectory_files)

    param = calc(uni, output1, output2)
