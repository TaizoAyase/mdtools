#!/usr/bin/env python

from MDAnalysis import Universe
from MDAnalysis.analysis.density import density_from_Universe

uni = Universe("ref.gro", "./all_fit.xtc")
D = density_from_Universe(uni, delta = 1.0, atomselection = "name OW")
D.export("density.dx")
