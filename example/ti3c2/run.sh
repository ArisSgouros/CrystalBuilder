#!/bin/bash
nx=5
ny=2
nz=1
python ../../crystal_builder.py in.basis_ti3c2.dat "$nx,$ny,$nz" --file_pos="o.ti3c2_"$nx"x"$ny"x"$nz".dat" --file_dump="o.ti3c2_"$nx"x"$ny"x"$nz".lammpstrj" --file_xyz="o.ti3c2_"$nx"x"$ny"x"$nz".xyz" --rc="0.0" --bond=0 > o.log
#nx=31
#ny=6


nx=62
ny=12
nz=1
python ../../crystal_builder.py in.basis_ti3c2.dat "$nx,$ny,$nz" --file_pos="o.ti3c2_"$nx"x"$ny"x"$nz".dat" --file_dump="o.ti3c2_"$nx"x"$ny"x"$nz".lammpstrj" --file_xyz="o.ti3c2_"$nx"x"$ny"x"$nz".xyz" --rc="0.0" --bond=0 > o.log
python pypore.py
#nx=31
#ny=6
nz=2
python ../../crystal_builder.py in.basis_ti3c2.dat "$nx,$ny,$nz" --file_pos="o.ti3c2_"$nx"x"$ny"x"$nz".dat" --file_dump="o.ti3c2_"$nx"x"$ny"x"$nz".lammpstrj" --file_xyz="o.ti3c2_"$nx"x"$ny"x"$nz".xyz" --rc="0.0" --bond=0 > o.log
