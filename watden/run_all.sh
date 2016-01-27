#!/bin/tcsh

##### 
# !!!!! 
# first of all, make index.ndx with "gmx make_ndx"
# which includes both 'C-alpha and Water' group
# !!!!!
####

gmx convert-tpr -s ../../5_prod1/prod1_00000.tpr -n index.ndx

ruby ./trjconv_pbc.rb

ruby ./all_trjcat.rb
rm ./prod1*.xtc

gmx trjconv -f all.xtc -o all_fit.xtc -s ./tpxout.tpr -fit rot+trans <<EOS
3
0
EOS
rm all.xtc

python ./calc_density.py

ruby ./convert.rb density.dx > reformat_density.dx

### make averaged strucutre of protein
python ./make_average.py > ./average.pdb
