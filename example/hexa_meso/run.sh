#!/bin/bash
python ../../crystal_builder.py in.basis_hexa.dat "10,2,1" --file_pos="o.hexa.dat" --file_dump="o.hexa.lammpstrj" --file_xyz="o.hexa.xyz" --rc="0.0" --drc=0.1 > o.log
