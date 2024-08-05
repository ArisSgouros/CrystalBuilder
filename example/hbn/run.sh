#!/bin/bash
python ../../crystal_builder.py \
  in.basis_hbn.dat "4,4,2" \
  --file_pos="o.hbn.dat" \
  --file_dump="o.hbn.lammpstrj" \
  --file_xyz="o.hbn.xyz" \
  --rc="1.44568507405082" \
  --drc=0.1 \
  --bond=1 \
  --angle=1 \
  --dihed=1 \
   > o.log
