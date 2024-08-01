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
from module.net_types import Atom, AtomType, Dihed

def is_list_unique(list_in):
   return len(set(list_in)) == len(list_in)

def CalculateDiheds(atoms):
   diheds = []

   # populate the dihed list
   id_ = 1
   for atom in atoms:
      atom.exist = False
   for iatom in atoms:
      for jatom in iatom.neigh:
         # dihed ijkl
         if jatom.exist: continue
         for katom in jatom.neigh:
            if katom.exist or katom == iatom: continue
            for latom in katom.neigh:
               dihed = [iatom, jatom, katom, latom]
               if latom.exist or not is_list_unique(dihed): continue
               diheds.append(Dihed(id_, iatom, jatom, katom, latom))
               id_ += 1
         # dihed jikl
         for katom in iatom.neigh:
            if katom.exist or katom == jatom: continue
            for latom in katom.neigh:
               dihed = [jatom, iatom, katom, latom]
               if latom.exist or not is_list_unique(dihed): continue
               diheds.append(Dihed(id_, jatom, iatom, katom, latom))
               id_ += 1
      iatom.exist = True

   return diheds
