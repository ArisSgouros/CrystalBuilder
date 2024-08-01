import numpy as np
import copy as cp
import sys

class Atom:
   def __init__(self, aid, mol, type, qq, rr):
      self.aid = aid
      self.molid = mol
      self.type = type
      self.qq = qq
      self.rr = rr

def ReadDump(filename):
   foo = open(filename, "r")

   foo.readline() # ITEM: TIMESTEP
   foo.readline() #
   foo.readline() # NUMBER OF ATOMS
   nat = int(foo.readline())
   foo.readline() # ITEM: BOX BOUNDS pp pp pp
   xlo, xhi = [float(val) for val in foo.readline().split()]
   ylo, yhi = [float(val) for val in foo.readline().split()]
   zlo, zhi = [float(val) for val in foo.readline().split()]
   foo.readline() # ITEM: ATOMS id mol type xu yu zu

   box = [xlo, xhi, ylo, yhi, zlo, zhi]

   lx = xhi - xlo
   ly = yhi - ylo
   lz = zhi - zlo

   atoms = []

   for iat in range(nat):
      data = foo.readline().split()
      id_ = int(data[0])
      mol = int(data[1])
      type_ = int(data[2])
      xx = float(data[3])
      yy = float(data[4])
      zz = float(data[5])

      qq = 0.0
      rr = [xx, yy, zz]
      atoms.append(Atom(id_, mol, type_, qq, rr))

   return [box, atoms]

def ExportDump(box, atoms, filename):
   foo = open(filename, "w")

   foo.write("ITEM: TIMESTEP\n")
   foo.write("0\n")
   foo.write("ITEM: NUMBER OF ATOMS\n")
   foo.write("%d\n" % (len(atoms)))
   foo.write("ITEM: BOX BOUNDS pp pp pp\n")
   foo.write("%f %f\n" % (box[0], box[1]))
   foo.write("%f %f\n" % (box[2], box[3]))
   foo.write("%f %f\n" % (box[4], box[5]))
   foo.write("ITEM: ATOMS id mol type xu yu zu\n")

   id_ = 1
   for atom in atoms:
      foo.write("%d %d %d %f %f %f\n" % (id_, atom.molid, atom.type, atom.rr[0], atom.rr[1], atom.rr[2]))
      id_ += 1
   foo.close()
 


#
# Interface
#

#file_in = "o.ti3c2_30x6x1.lammpstrj"
file_in = "o.ti3c2_62x12x1.lammpstrj"
box, atoms = ReadDump(file_in)


atoms_out = cp.deepcopy(atoms)

file_out = "o.mod.lammpstrj"



lx = box[1] - box[0]
ly = box[3] - box[2]
pore_radius = 5.2
pore_xy_rel = [[ix, iy, pore_radius] for ix in [0.083333333,0.25,0.416666667,0.583333333,0.75,0.916666667] for iy in [0.25, 0.75]]

def IsInsidePore(rr_at, rr_cyl, radius):
   dx = rr_cyl[0] - rr_at[0]
   dy = rr_cyl[1] - rr_at[1]
   dr = np.sqrt(dx**2 + dy**2)
   return (dr < radius)


for [ix, iy, irad] in pore_xy_rel:

   pore_xy = [ix*lx, iy*ly]
   pore_radius = irad

   atoms_out = [atom for atom in atoms_out if not IsInsidePore(atom.rr, pore_xy, pore_radius)]


ExportDump(box, atoms_out, file_out)

sys.exit()

