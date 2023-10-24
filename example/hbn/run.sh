#!/bin/bash
python ../../crystal_builder.py in.basis_hbn.dat "4,4,2" --file_pos="o.hbn.dat" --file_dump="o.hbn.lammpstrj" --rc="1.44568507405082" --drc=0.1 > o.log
