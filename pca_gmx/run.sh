#!/bin/tcsh

# preparation
gmx make_ndx -f ../../5_prod1/npt.gro
gmx convert-tpr -s ../../5_prod1/prod1_00000.tpr -n index.ndx

# run
gmx covar -f ../trjcat2/all.xtc -s ../trjconv/ref.gro
gmx anaeig -v eigenvec.trr -eig eigenval.xvg -f ../trjcat2/all.xtc -s tpxout.tpr -n index.ndx -rmsf -comp
