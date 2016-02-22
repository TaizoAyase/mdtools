#!/usr/bin/env python

from __future__ import print_function
from mdtools.rms_tools import fit, centroid, kabsch, rmsd, superpose
from MDAnalysis import *
import numpy as np
import scipy
import sys
import os

#
# script for performing Principal component analysis (PCA) for MD trajectory
# performed by diagonalizing the covariance matrix
# covariance matrix is built based on deviation
# from avaraged structure of protein C-alhpa atoms
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

eigval_threshold = 90  # % of total was displayed
prcomp = 5  # max num of components to write
scale_factor = 5  # for showing vectors in qsc


##################

sys.stderr.write('')

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
for ts in uni.trajectory:
    sys.stderr.write('Loading trajectory %d/%d ...\r' % (ts.frame, n_frames))
    coord = superpose(sel_ca.coordinates(), ref_coord).flatten()
    coord_ave += coord
sys.stderr.write('\n')
coord_ave /= n_frames


# build covariance matrix
# sum(deviation x deviation) / frames
sys.stderr.write('Calc covariance matrix ...\n')
deviation_all = np.zeros((n_frames, dof))
uni.trajectory[0]
for ts in uni.trajectory:
    sys.stderr.write('Loading trajectory %d/%d ...\r' % (ts.frame, n_frames))
    coord = superpose(sel_ca.coordinates(), ref_coord).flatten()
    deviation = coord - coord_ave
    deviation_all[ts.frame] = deviation

sys.stderr.write('\n')

# np.cov result is defferent by n_frames ** 2 order
sys.stderr.write('Building covariance matrix ...\n')
deviation_all /= n_frames
cov = np.cov(deviation_all.T) * n_frames * n_frames

np.savetxt(cov_out, cov, delimiter=', ')

# calc eigen values
sys.stderr.write('Calc eigen vectors and eigen values ...\n')
values, vectors = scipy.linalg.eigh(cov, turbo=True)

sys.stderr.write('Write out the eigen vector/value file ...\n')
np.savetxt('eigen_val.csv', values, delimiter=', ')
np.savetxt('eigen_vec.csv', vectors, delimiter=', ')

# print principal component contributions
eigval_tot = np.sum(values)
sys.stderr.write('Calc contributions ...\n')

sum_val = 0
# loop with reverse of values array
for i, val in enumerate(values[::-1]):
    mes = 'Comp%d: %5.2f %%' % (i, val/eigval_tot * 100)
    print(mes)
    sum_val += val/eigval_tot * 100
    if sum_val > eigval_threshold:
        print('----- %4.2f %% total -----' % sum_val)
        break

# print components in 3D-vectors format
# loop with reverse vector matrix
sys.stderr.write('Write out the vector file ...\n')
f = open(vec_out, 'w+')

for i, vec in enumerate(vectors[::-1]):
    f.write('@%d comp\n' % i)
    for j in xrange(n_atoms):
        vec_x = vec[3*j + 0]
        vec_y = vec[3*j + 1]
        vec_z = vec[3*j + 2]
        vec_len = np.linalg.norm(vec)
        f.write('%6d, %9.5f, %9.5f, %9.5f, %9.5f\n' % (j, vec_x, vec_y, vec_z, vec_len))
    if i > prcomp:
        break
f.close()


# projection
sys.stderr.write('Write out the projection file ...\n')
proj_vec = np.dot(deviation_all, vectors)
np.savetxt('proj_all.csv', proj_vec, delimiter=', ')

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

for i in xrange(prcomp):
    f.write(rend_line % i)
    vec = vectors[:, -(i+1)].reshape(n_atoms, 3)
    for j, v in enumerate(vec):
        plus_vec  = ref_coord[j] + v * np.sqrt(values[-(i+1)]) * scale_factor
        minus_vec = ref_coord[j] - v * np.sqrt(values[-(i+1)]) * scale_factor
        f.write(vect_line % tuple(np.array((plus_vec, minus_vec)).flatten()))
    f.write(tail_line)

tail = '''    </object>
<camera center="(66.727969,61.483686,86.935145)" centerMark="crosshair" distance="200" name="__current" perspec="false" rotation="(0.054345,-0.623993,-0.479576,0.614562)" slab="146.069874" stereoDist="1" stereoMode="none" zoom="120.876559"/>
</scene>'''
f.write(tail)
f.close()
