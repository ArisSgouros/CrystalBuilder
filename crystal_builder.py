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

class AtomType:
   def __init__(self, type, mass, name):
      self.type = type
      self.mass = mass
      self.name = name

class Atom:
   def __init__(self, aid, mol, type, qq, rr):
      self.aid = aid
      self.molid = mol
      self.type = type
      self.qq = qq
      self.rr = rr

class BondType:
   def __init__(self, type):
      self.type = type

class Bond:
   def __init__(self, bid, iatom, jatom):
      self.bid = bid
      self.iatom = iatom
      self.jatom = jatom
      self.type = -1

# Parse basis files subject to the constraint coordinate system used in Lammps [see https://docs.lammps.org/Howto_triclinic.html]
def ReadBasis(path_basis):
   n_atom = 0
   n_atom_types = 0
   xlo = 0.0
   xhi = 0.0
   ylo = 0.0
   yhi = 0.0
   zlo = 0.0
   zhi = 0.0
   xy = 0.0
   xz = 0.0
   yz = 0.0

   g = open(path_basis, 'r')
   lines = []
   line_num = 0
   while True:
      line = g.readline()
      # Check if this is the end of file; if yes break the loop.
      if line == '':
         break

      line_split = line.split()

      if "atoms" in line:
         n_atom = int(line_split[0])
      if "atom types" in line:
         n_atom_types = int(line_split[0])
      if "xlo xhi" in line:
         xlo = float(line_split[0])
         xhi = float(line_split[1])
      if "ylo yhi" in line:
         ylo = float(line_split[0])
         yhi = float(line_split[1])
      if "zlo zhi" in line:
         zlo = float(line_split[0])
         zhi = float(line_split[1])
      if "xy xz yz" in line:
         xy = float(line_split[0])
         xz = float(line_split[1])
         yz = float(line_split[2])
      if "Masses" in line:
         Masses_start = line_num + 2
      if "Atoms" in line:
         Atoms_start = line_num + 2
      lines.append(line)
      line_num += 1
   g.close()
   #
   # Read the Masses section
   atom_types    = {}
   try:
      for ii in range(n_atom_types):
         line  = lines[Masses_start+ii].split()
         itype = int(line[0])
         mass  = float(line[1])
         name  = ""
         if len(line) == 4:
            name = line[3]
         atom_types[itype] = AtomType(itype, mass, name)
   except:
      print("Problem with reading masses section/not existing")
      sys.exit()
   #
   # Read the Atoms section
   imol = 1
   basis_atoms   = []
   try:
      for ii in range(n_atom):
         line  = lines[Atoms_start+ii].split()
         aid   = int(line[0])
         itype = int(line[2])
         qq    = float(line[3])
         xx    = float(line[4])
         yy    = float(line[5])
         zz    = float(line[6])
         rr    = np.array([xx, yy, zz])
         basis_atoms.append(Atom(aid, imol, itype, qq, rr))
   except:
      print("Problem with reading atom section")
      sys.exit()

   # Calculate the origin
   basis_orig = np.array([xlo, ylo, zlo])

   # Calculate the basis vector [see https://docs.lammps.org/Howto_triclinic.html]
   basis_vectors = []
   basis_vectors.append(np.array([xhi-xlo, 0.0    , 0.0]))
   basis_vectors.append(np.array([xy     , yhi-ylo, 0.0]))
   basis_vectors.append(np.array([xz     , yz     , zhi-zlo]))
   return [basis_vectors, basis_orig, basis_atoms, atom_types]

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

def MinImag(rr, lx, ly, lz, xy, xz, yz, periodicity):
   rr_imag = np.array([ rr[0] - xy*round(rr[1]/ly) - xz*round(rr[2]/lz), \
                        rr[1] - yz*round(rr[2]/lz), \
                        rr[2] ])
   rr_imag[0] -= int(periodicity[0])*lx*round(rr_imag[0]/lx)
   rr_imag[1] -= int(periodicity[1])*ly*round(rr_imag[1]/ly)
   rr_imag[2] -= int(periodicity[2])*lz*round(rr_imag[2]/lz)
   return rr_imag

def UniquePairType(itype, jtype):
   sort_types = sorted(set({itype, jtype}))
   return str(sort_types[0])+" "+str(sort_types[1])

def CalculateBonds(atoms, box_lmp, periodicity, rc_list, drc):
   lx = box_lmp['xhi'] - box_lmp['xlo']
   ly = box_lmp['yhi'] - box_lmp['ylo']
   lz = box_lmp['zhi'] - box_lmp['zlo']
   xy  = box_lmp['xy']
   xz  = box_lmp['xz']
   yz  = box_lmp['yz']

   bonds = []
   bond_types = {}

   natom = len(atoms)
   bid = 1

   for rc in rc_list:
      rmin2 = pow(rc - drc, 2)
      rmax2 = pow(rc + drc, 2)
      for ii in range(natom):
         ri = atoms[ii].rr
         iid = atoms[ii].aid
         for jj in range(ii+1, natom):
            rj = atoms[jj].rr
            jid = atoms[ii].aid
            rij = rj - ri
            rij = MinImag(rij, lx, ly, lz, xy, xz, yz, periodicity)
            rij2 = np.dot(rij, rij)
            if rmin2 < rij2 < rmax2:
               bonds.append(Bond(bid, atoms[ii], atoms[jj]))
               bid += 1

   # determine the individual bond types
   itype = 1
   for bond in bonds:
      btype_str = UniquePairType(bond.iatom.type, bond.jatom.type)
      if not bond_types.get(btype_str):
         bond_types[btype_str] = BondType(itype)
         itype += 1
      bond.type = bond_types[btype_str].type

   return bonds, bond_types

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
