#!/usr/bin/env python

from __future__ import print_function
from mdtools.rms_tools import fit, centroid, kabsch, rmsd, superpose
from tqdm import tqdm
import mdtraj as md
import numpy as np
import scipy.linalg
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

reference_file = '../test_files/chainB.gro'
structure_file = reference_file
trajectory_file = '../test_files/chainB.xtc'

fit = 'name CA'

cov_out = 'covar.csv'
vec_out = 'vec_all.csv'
proj_out = 'proj_out.csv'
pdb_out = 'prot.pdb'
qsc_out = 'showvec.qsc'

eigval_threshold = 90  # % of total was displayed
prcomp = 5  # max num of components to write
scale_factor = 5  # for showing vectors in qsc

chunk_size = 1

##################

# loading trajectory
sys.stderr.write('Loading Universe ...\n')
ref = md.load(reference_file)

sel_ca = ref.topology.select(fit)
sel_ca_ref = ref.topology.select(fit)

if len(sel_ca) != len(sel_ca_ref):
    raise RuntimeError('Atom number is not matched.')

ref_coord = ref.xyz[0, sel_ca_ref, :] * 10. # nm to A
n_frames = 0
n_atoms = len(sel_ca)
dof = n_atoms * 3

# make average coordinate
coord_ave = np.zeros(dof)
sys.stderr.write('Calc average coordinate...\n')
for chunk in tqdm(md.iterload(trajectory_file, top=structure_file, chunk=chunk_size)):
    coord = superpose(chunk.xyz[0, sel_ca_ref, :] * 10, ref_coord).flatten()
    coord_ave += coord
    n_frames += chunk.n_frames
coord_ave /= n_frames


# build covariance matrix
# sum(deviation x deviation) / frames
sys.stderr.write('Calc covariance matrix ...\n')
deviation_all = np.zeros((n_frames, dof))
i = 0
for chunk in tqdm(md.iterload(trajectory_file, top=structure_file, chunk=chunk_size), total=n_frames):
    coord = superpose(chunk.xyz[0, sel_ca_ref, :] * 10, ref_coord).flatten()
    deviation = coord - coord_ave
    deviation_all[i] = deviation
    i += 1


sys.stderr.write('Building covariance matrix ...\n')
cov = np.cov(deviation_all.T)

np.savetxt(cov_out, cov, delimiter=', ')

# calc eigen values
sys.stderr.write('Calc eigen vectors and eigen values ...\n')
values, vectors = scipy.linalg.eigh(cov)

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
    for j in range(n_atoms):
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
np.savetxt(proj_out, proj_vec[:, :prcomp:-1], delimiter=', ')

# write out qsc-file to visualize
sys.stderr.write('Write out the QSC file ...\n')
ref.save_pdb(pdb_out)


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
#ref_coord = sel_ca_ref.positions
ref_coord = ref.xyz[0, sel_ca_ref, :] * 10 # nm to A

for i in range(prcomp):
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
