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

def ExportLammpsDataFile(filename, box_lmp, atom_types, bond_types, atoms, bonds):
   xlo = box_lmp['xlo']
   xhi = box_lmp['xhi']
   ylo = box_lmp['ylo']
   yhi = box_lmp['yhi']
   zlo = box_lmp['zlo']
   zhi = box_lmp['zhi']
   xy  = box_lmp['xy']
   xz  = box_lmp['xz']
   yz  = box_lmp['yz']
   is_triclinic = (xy != 0.0 or xz != 0.0 or yz != 0.0)

   f = open(filename, 'w')

   f.write('# LAMMPS data file\n')
   f.write('\n')
   f.write('%d atoms\n' % len(atoms))
   f.write('%d bonds\n' % len(bonds))
   #f.write('%d angles\n' % len(angles))
   #f.write('%d dihedrals\n' % len(dihedrals))
   f.write('\n')
   f.write('%d atom types\n' % len(atom_types))
   f.write('%d bond types\n' % len(bond_types))
   #f.write('%d angle types\n' % len(angle_types))
   #f.write('%d dihedrals types\n' % len(dihedral_types))
   f.write('\n')
   f.write('%-16.9f %-16.9f xlo xhi\n' % (xlo, xhi))
   f.write('%-16.9f %-16.9f ylo yhi\n' % (ylo, yhi))
   f.write('%-16.9f %-16.9f zlo zhi\n' % (zlo, zhi))
   if is_triclinic:
      f.write('%-16.9f %-16.9f %-16.9f xy xz yz\n' % (xy, xz, yz))
   f.write('\n')
   f.write('Masses\n')
   f.write('\n')
   for atom_type in atom_types.values():
      f.write("%d %-16.9f # %s\n" % (atom_type.type, atom_type.mass, atom_type.name))
   f.write('\n')
   f.write('Atoms\n')
   f.write('\n')
   for iat in atoms:
      f.write("%d %d %d %-16.9f %-16.9f %-16.9f %-16.9f # %s\n" % ( iat.aid, iat.molid, iat.type, iat.qq, iat.rr[0], iat.rr[1], iat.rr[2], atom_types[iat.type].name ))
   if bonds:
      f.write('\n')
      f.write('Bonds\n')
      f.write('\n')
      for bond in bonds:
         iatom = bond.iatom
         jatom = bond.jatom
         btype_str = UniquePairType(iatom.type, jatom.type)
         f.write("%d %d %d %d # %s\n" % (bond.bid, bond.type, iatom.aid, jatom.aid, btype_str))
   #if n_angles:
   #   f.write('\n')
   #   f.write('Angles\n')
   #   f.write('\n')
   #   for id in range(n_angles):
   #      f.write("%d %d %d %d %d\n" % (id+1, 1, angles[id][0]+1, angles[id][1]+1, angles[id][2]+1 ))
   #if n_dihedrals:
   #   f.write('\n')
   #   f.write('Dihedrals\n')
   #   f.write('\n')
   #   for id in range(n_dihedrals):
   #      f.write("%d %d %d %d %d %d\n" % (id+1, 1, dihedrals[id][0]+1, dihedrals[id][1]+1, dihedrals[id][2]+1, dihedrals[id][3]+1 ))

   f.close()

def ExportLammpsDumpFile(filename, box_lmp, atoms):
   xlo = box_lmp['xlo']
   xhi = box_lmp['xhi']
   ylo = box_lmp['ylo']
   yhi = box_lmp['yhi']
   zlo = box_lmp['zlo']
   zhi = box_lmp['zhi']
   xy  = box_lmp['xy']
   xz  = box_lmp['xz']
   yz  = box_lmp['yz']
   is_triclinic = (xy != 0.0 or xz != 0.0 or yz != 0.0)

   fout = open(filename, 'w')
   fout.write("ITEM: TIMESTEP\n")
   fout.write("0\n")
   fout.write("ITEM: NUMBER OF ATOMS\n")
   fout.write("%-8d\n" % len(atoms))

   if is_triclinic:
      xlo_bound = xlo + min(0.0,xy,xz,xy+xz)
      xhi_bound = xhi + max(0.0,xy,xz,xy+xz)
      ylo_bound = ylo + min(0.0,yz)
      yhi_bound = yhi + max(0.0,yz)
      zlo_bound = zlo
      zhi_bound = zhi
      fout.write("ITEM: BOX BOUNDS xy xz yz\n")
      fout.write("%-16.9f %-16.9f %-16.9f\n" % (xlo_bound, xhi_bound, xy) )
      fout.write("%-16.9f %-16.9f %-16.9f\n" % (ylo_bound, yhi_bound, xz) )
      fout.write("%-16.9f %-16.9f %-16.9f\n" % (zlo_bound, zhi_bound, yz) )
   else:
      fout.write("ITEM: BOX BOUNDS pp pp pp\n")
      fout.write("%-16.9f %-16.9f\n" % (xlo, xhi) )
      fout.write("%-16.9f %-16.9f\n" % (ylo, yhi) )
      fout.write("%-16.9f %-16.9f\n" % (zlo, zhi) )
   fout.write("ITEM: ATOMS id mol type xu yu zu\n")

   for iat in atoms:
      fout.write("%-d %-d %-d %-f %-f %-f\n" % (iat.aid, iat.molid, iat.type, iat.rr[0], iat.rr[1], iat.rr[2]))

   fout.close()
