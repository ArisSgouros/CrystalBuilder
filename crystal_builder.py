###############################################################################
# MIT License
#
# Copyright (c) 2023 ArisSgouros
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
###############################################################################

import sys
import os
import math as m
import numpy as np
import argparse
import copy as cp
from module.aux import UniquePairType
from module.export import ExportLammpsDataFile, ExportLammpsDumpFile
from module.read import ReadBasis
from module.net_types import Atom, Bond, AtomType, BondType
from module.calculate_bond import MinImag, CalculateBonds

parser = argparse.ArgumentParser(description='Generate crystal structure')
parser.add_argument('path_basis', type=str, help='Path of basis file (Lammps datafile)')
parser.add_argument('cells', type=str, help='delimited list of replications')
parser.add_argument('-periodic', '--periodic', type=str, default="1,1,1", help='delimited list of periodicity directions')
parser.add_argument('-rc', '--rc', type=str, default="", help='list of cutoff radii of bonded interactions')
parser.add_argument('-drc', '--drc', type=float, default=0.1, help='tolerance of bonded interactions')
parser.add_argument('-file_pos', '--file_pos', type=str, default='pos.dat', help='Name of the Lammps data file')
parser.add_argument('-file_dump', '--file_dump', type=str, default='dump.lammpstrj', help='Name of the Lammps dump file')

# constants
kDim = 3
kTol = 1.e-8

if __name__ == "__main__":
   args = parser.parse_args()
   path_basis = args.path_basis
   cells = [int(item) for item in args.cells.split(',')]
   periodicity = [int(item) for item in args.periodic.split(',')]
   rc_list = [float(item) for item in args.rc.split(',')]
   drc = args.drc
   calc_bond = any(abs(i) > kTol for i in rc_list)

   file_pos = args.file_pos
   file_dump = args.file_dump

   print('path        : ', path_basis)
   print('cells       : ', cells)
   print('periodicity : ', periodicity)
   print('calc_bond   : ', calc_bond)
   print('rc(s)       : ', rc_list)
   print('drc         : ', drc)

   print("Reading the datafile:", path_basis)
   basis_vectors, basis_orig, basis_atoms, atom_types = ReadBasis(path_basis)
   print()

   print("basis origin:")
   print('   ', basis_orig)
   print("basis vectors:")
   for basis_vector in basis_vectors:
      print('   ', basis_vector)
   print("atom types:")
   for atom_type in atom_types.values():
      print('   ',atom_type.type, atom_type.mass, atom_type.name)
   print("basis atoms:")
   print("   %-8s %-8s %-16s %-16s %-16s %-16s" % ("id", "type", "charge", "x", "y", "z"))
   for atom in basis_atoms:
      print("   %-8d %-8d %-16.9f %-16.9f %-16.9f %-16.9f" % (atom.aid, atom.type, atom.qq, atom.rr[0], atom.rr[1], atom.rr[2]))
   print("box vectors:")
   box_vectors = []
   box_orig = []
   for ii in range(kDim):
      box_vectors.append(basis_vectors[ii]*cells[ii])
      box_orig.append(basis_orig[ii]*cells[ii])
      print('   ', box_vectors[ii])
   print("box vectors (lammps constraints):")
   box_lmp = {'xlo':box_orig[0], 'xhi':box_vectors[0][0], \
              'ylo':box_orig[1], 'yhi':box_vectors[1][1], \
              'zlo':box_orig[2], 'zhi':box_vectors[2][2], \
              'xy':box_vectors[1][0], 'xz':box_vectors[2][0], 'yz':box_vectors[2][1]}
   for key in box_lmp:
      print("   %-3s : %-16.9f" % (key, box_lmp[key]))

   print()
   print("Generating the atomic coordinates..")
   atoms = []
   aid = 1
   imol = 1
   for iz in range(cells[2]):
      for iy in range(cells[1]):
         for ix in range(cells[0]):
            for atom in basis_atoms:
               rr = atom.rr + basis_vectors[0]*ix + basis_vectors[1]*iy +basis_vectors[2]*iz
               atoms.append(Atom(aid, imol, atom.type, atom.qq, rr))
               aid+=1
   bonds = []
   bond_types = []
   if calc_bond:
      print("Generating bonds..")
      bonds, bond_types = CalculateBonds(atoms, box_lmp, periodicity, rc_list, drc)

   print()
   print("Generating a lammps datafile with name " + file_pos + " ..")
   ExportLammpsDataFile(file_pos, box_lmp, atom_types, bond_types, atoms, bonds)
   print()
   print("Generating a lammps dump file with name " + file_dump + " ..")
   ExportLammpsDumpFile(file_dump, box_lmp, atoms)
