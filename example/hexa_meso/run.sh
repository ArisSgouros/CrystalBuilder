#!/bin/bash
python ../../crystal_builder.py \
  in.basis_hexa.dat "10,2,1" \
  --file_pos="o.hexa.dat" \
  --file_dump="o.hexa.lammpstrj" \
  --file_xyz="o.hexa.xyz" \
  --rc="1.0" \
  --drc=0.1 \
  --bond=1 \
  --angle=1 \
  --dihed=1 \
  > o.log
