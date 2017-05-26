#!/usr/bin/env python

from __future__ import print_function
import numpy as np


BOLTZMANN = 0.0019872041 # kcal / mol

RESID_DICT = {
    'ALA': ['N', 'O'],
    'CYS': ['N', 'O', 'SG'],
    'ASP': ['N', 'O', 'OD1', 'OD2'],
    'ASPP': ['N', 'O', 'OD1', 'OD2'],
    'GLU': ['N', 'O', 'OE1', 'OE2'],
    'GLUP': ['N', 'O', 'OE1', 'OE2'],
    'PHE': ['N', 'O'],
    'GLY': ['N', 'O'],
    'HIS': ['N', 'O', 'ND1', 'NE2'],  # normal His
    'HSD': ['N', 'O', 'ND1', 'NE2'],  # neutral His proton on ND1
    'HDE': ['N', 'O', 'ND1', 'NE2'],  # neutral His proton on NE2
    'HSP': ['N', 'O', 'ND1', 'NE2'],  # protonated His
    'ILE': ['N', 'O'],
    'LYS': ['N', 'O', 'NZ'],
    'LEU': ['N', 'O'],
    'MET': ['N', 'O', 'SD'],
    'ASN': ['N', 'O', 'OD1', 'ND2'],
    'PRO': ['N', 'O'],
    'GLN': ['N', 'O', 'OE1', 'NE2'],
    'ARG': ['N', 'O', 'NE', 'NH1', 'NH2'],
    'SER': ['N', 'O', 'OG'],
    'THR': ['N', 'O', 'OG1'],
    'VAL': ['N', 'O'],
    'TRP': ['N', 'O', 'NE1'],
    'TYR': ['N', 'O', 'OH'],
}


RESID_DICT_SIDECHAIN = {
    'CYS': ['SG'],
    'ASP': ['OD1', 'OD2'],
    'ASPP': ['OD1', 'OD2'],
    'GLU': ['OE1', 'OE2'],
    'GLUP': ['OE1', 'OE2'],
    'HIS': ['ND1', 'NE2'],  # normal His
    'HSD': ['ND1', 'NE2'],  # neutral His proton on ND1
    'HDE': ['ND1', 'NE2'],  # neutral His proton on NE2
    'HSP': ['ND1', 'NE2'],  # protonated His
    'LYS': ['NZ'],
    'MET': ['SD'],
    'ASN': ['OD1', 'ND2'],
    'GLN': ['OE1', 'NE2'],
    'ARG': ['NE', 'NH1', 'NH2'],
    'SER': ['OG'],
    'THR': ['OG1'],
    'TRP': ['NE1'],
    'TYR': ['OH'],
}