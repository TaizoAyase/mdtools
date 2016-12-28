#!/usr/bin/env python

from __future__ import print_function
from mdtools.rms_tools import fit, centroid, kabsch, rmsd, superpose
from tqdm import tqdm
from MDAnalysis import *
import numpy as np
import scipy.linalg
import sys
import os

#
# script for time-structure based independent component analysis(tICA)
# ref: doi:10.1063/1.3554380
#
# finally, cuemol (http://www.cuemol.org/) scene file
# with protein PDB and eigen vectors will be written
#

##### setup #####
#
reference_file = '/path/to/reference/pdb_or_gro'
structure_file = reference_file
trajectory_file = '/path/to/trajectory/file'

fit = 'name CA'

cov_out = 'covar.csv'
vec_out = 'vec_all.csv'
proj_out = 'proj_out.csv'
pdb_out = 'prot.pdb'
qsc_out = 'showvec.qsc'

icomp = 5  # max num of components to write
scale_factor = 5  # for showing vectors in qsc

lag_step = 100  # steps = 1ns
lag_time = 1  # ns, used for calc of decay time-constant

##################

# loading trajectory
sys.stderr.write('Loading Universe ...\n')
uni = Universe(structure_file, trajectory_file)
uni_ref = Universe(reference_file)

sel_ca = uni.select_atoms(fit)
sel_ca_ref = uni_ref.select_atoms(fit)

if sel_ca.n_atoms != sel_ca_ref.n_atoms:
    raise RuntimeError('Atom number is not matched.')

ref_coord = sel_ca_ref.coordinates()
n_frames = uni.trajectory.n_frames
n_atoms = sel_ca_ref.n_atoms
dof = n_atoms * 3

# make average coordinate
coord_ave = np.zeros(dof)
sys.stderr.write('Calc average coordinate...\n')
for ts in tqdm(uni.trajectory):
    coord = superpose(sel_ca.coordinates(), ref_coord).flatten()
    coord_ave += coord
coord_ave /= n_frames


# build covariance matrix
# sum(deviation x deviation) / frames
sys.stderr.write('Calc deviation from averaged structure ...\n')
deviation_all = np.zeros((n_frames, dof))
uni.trajectory[0]
for ts in tqdm(uni.trajectory):
    coord = superpose(sel_ca.coordinates(), ref_coord).flatten()
    deviation = coord - coord_ave
    deviation_all[ts.frame] = deviation

sys.stderr.write('Calc time-lagged correlation matrix ...\n')
offset_correl_tmp = np.zeros((dof, dof))
for i in tqdm(xrange(0, n_frames - lag_step)):
    sys.stderr.write('Building %d/%d ...\r' % (i, n_frames - lag_step))
    offset_correl_tmp += np.outer(deviation_all[i], deviation_all[i+lag_step])
sys.stderr.write('\n')

offset_correl_tmp /= (n_frames - lag_step)
# symmetrization
offset_correl = (offset_correl_tmp + offset_correl_tmp.T)/2

# np.cov result is defferent by n_frames order
sys.stderr.write('Building covariance matrix ...\n')
cov = np.cov(deviation_all.T)

np.savetxt(cov_out, cov, delimiter=', ')

# calc generalized eigen values
sys.stderr.write('Calc eigen vectors and eigen values ...\n')
values, vectors = scipy.linalg.eigh(offset_correl, b=cov, turbo=True)

sys.stderr.write('Write out the eigen vector/value file ...\n')
np.savetxt('eigen_val.csv', values, delimiter=', ')
np.savetxt('eigen_vec.csv', vectors, delimiter=', ')


# print the decay time-constant of autocorrelation function
# estimated from each independent component
sys.stderr.write('Estimated decay time-constant ...\n')

# loop with reverse of values array
time_constant = -1 * lag_time / np.log(values[::-1])
for i, val in enumerate(time_constant):
    mes = 'IC%-2d: %+5.2f ns' % (i, val)
    print(mes)
    if abs(val) < lag_time:
        break

# projection eigen vectors to real space
# by g = C * eig_vec
proj_vec = np.dot(cov, vectors)

# print independent components in 3D-vectors format
# loop with reverse vector matrix
sys.stderr.write('Write out the vector file ...\n')
f = open(vec_out, 'w+')

for i, vec in enumerate(proj_vec[::-1]):
    f.write('@%d comp\n' % i)
    for j in xrange(n_atoms):
        vec_x = vec[3*j + 0]
        vec_y = vec[3*j + 1]
        vec_z = vec[3*j + 2]
        vec_len = np.linalg.norm(vec)
        f.write('%6d, %9.5f, %9.5f, %9.5f, %9.5f\n' % (j, vec_x, vec_y, vec_z, vec_len))
    if i > icomp:
        break
f.close()


# projection
# projection must be projected to g vector
sys.stderr.write('Write out the projection file ...\n')
projected_deviation = np.dot(deviation_all, vectors.T)
np.savetxt(proj_out, projected_deviation[:, :icomp], delimiter=', ')

# write out qsc-file to visualize
sys.stderr.write('Write out the QSC file ...\n')
prot_sel = uni_ref.select_atoms('protein')
prot_sel.write(pdb_out, format='PDB')

f = open(qsc_out, 'w+')

header = '''<?xml version="1.0" encoding="utf-8"?>
<scene>
    <qsc_opts base64="false" compress="none" version="QDF0"/>
    <object type="MolCoord" alt_src="%s" name="prot" sel="*" srctype="pdb" src="%s">
        <coloring type="PaintColoring">
            <paint sel="sheet" color="SteelBlue"/>
            <paint sel="helix" color="khaki"/>
            <paint sel="nucleic" color="yellow"/>
            <paint sel="*" color="FloralWhite"/>
        </coloring>
        <ropts/>
        <renderer type="tube" group="" name="tube1" style="DefaultHSCPaint"/>
        <renderer type="*selection" visible="false"/>
        <renderer type="*namelabel" style="DefaultAtomLabel"/>
''' % (os.path.abspath(pdb_out), pdb_out)

f.write(header)

rend_line = '\t\t<renderer type="atomintr" color="#BF0000" mode="fancy" showlabel="false" stipple0="1000.0" name="dom%d" visible="false" width="0.3">\n'
vect_line = '\t\t<line pos1="(%3.5f, %3.5f, %3.5f)" pos2="(%3.5f, %3.5f, %3.5f)"/>\n'
tail_line = '\t\t</renderer>\n'
ref_coord = sel_ca_ref.coordinates()

for i in xrange(icomp):
    f.write(rend_line % i)
    vec = proj_vec[:, -(i+1)].reshape(n_atoms, 3)
    for j, v in enumerate(vec):
        plus_vec  = ref_coord[j] + v * scale_factor
        minus_vec = ref_coord[j] - v * scale_factor
        f.write(vect_line % tuple(np.array((plus_vec, minus_vec)).flatten()))
    f.write(tail_line)

tail = '''    </object>
<camera center="(66.727969,61.483686,86.935145)" centerMark="crosshair" distance="200" name="__current" perspec="false" rotation="(0.054345,-0.623993,-0.479576,0.614562)" slab="146.069874" stereoDist="1" stereoMode="none" zoom="120.876559"/>
</scene>'''
f.write(tail)
f.close()
