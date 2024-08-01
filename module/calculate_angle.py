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
from module.net_types import Atom, AtomType, Angle

def is_list_unique(list_in):
   return len(set(list_in)) == len(list_in)

def CalculateAngles(atoms):
   angles = []

   # populate the angle list
   id_ = 1
   for atom in atoms:
      atom.exist = False
   for iatom in atoms:
      for jatom in iatom.neigh:
         # angle ijk
         if jatom.exist: continue
         for katom in jatom.neigh:
            angle = [iatom, jatom, katom]
            if katom.exist or not is_list_unique(angle): continue
            angles.append(Angle(id_, iatom, jatom, katom))
            id_ += 1
         # angle jik
         for katom in iatom.neigh:
            angle = [jatom, iatom, katom]
            if katom.exist or jatom.aid >= katom.aid or not is_list_unique(angle): continue
            angles.append(Angle(id_, jatom, iatom, katom))
            id_ += 1
      iatom.exist = True

   return angles
