###############################################################################
# MIT License
#
# Copyright (c) 2024 ArisSgouros
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

parser = argparse.ArgumentParser(description='Convert lattice constants to lammps dimensions')
parser.add_argument('-lattice', '--lattice', type=str, default="", help='delimited list of lattice constants (alpha, beta, gamma, ang_alpha, ang_beta, ang_gamma')

if __name__ == "__main__":
   args = parser.parse_args()

   arg_latt = [float(item) for item in args.lattice.split(',')]

   is_lattice = (len(arg_latt) == 6)

   [aa, bb, cc, anga, angb, angc] = arg_latt

   cosa = np.cos(np.radians(anga))
   cosb = np.cos(np.radians(angb))
   cosc = np.cos(np.radians(angc))

   lx = aa
   xy = bb*cosc
   xz = cc*cosb
   ly = np.sqrt(bb*bb - xy*xy)
   yz = (bb*cc*cosa - xy*xz) / ly
   lz = np.sqrt(cc*cc - xz*xz - yz*yz)

   print(aa, bb, cc, cosa, cosb, cosc)
   print(lx, ly, lz, xy, xz, yz)
 
