#!/bin/bash
python ../../crystal_builder.py in.basis_graphene.dat "20,10,1" --file_pos="o.20_10.dat" --file_dump="o.20_10.lammpstrj" --file_xyz="o.20_10.xyz" --rc="0.0" --drc=0.1 > o.log
python ../../crystal_builder.py in.basis_graphene.dat "30,4,1" --file_pos="o.30_4.dat" --file_dump="o.30_4.lammpstrj" --file_xyz="o.30_4.xyz" --rc="0.0" --drc=0.1 > o.log
python ../../crystal_builder.py in.basis_graphene.dat "200,4,1" --file_pos="o.200_4.dat" --file_dump="o.200_4.lammpstrj" --file_xyz="o.200_4.xyz" --rc="0.0" --drc=0.1 > o.log
