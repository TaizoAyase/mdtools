#!/bin/bash

trajectory_file=$1
tprfile=$2

gmx convert-tpr -s $tprfile <<EOS
2
EOS

gmx rmsf -f $trajectory_file -s tprout.tpr 
