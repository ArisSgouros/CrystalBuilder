#!/bin/bash
python ../../crystal_builder.py in.basis_nacl.dat "5,5,5" --file_pos="o.nacl.dat" --file_dump="o.nacl.lammpstrj" --file_xyz="o.nacl.xyz" --rc="2.8201" > o.log
