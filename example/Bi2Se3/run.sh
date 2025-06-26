#!/bin/bash

# Monolayer w/ 3x3 replications, bond type does not depend on bond length
xrep=3
yrep=3
name="o.type_not_depend_blen"$xrep"x"$yrep
python ../../crystal_builder.py \
  i.basis_Bi2Se3_1layer.dat "3,3,1" \
  --file_pos=$name".dat" \
  --file_dump=$name".lammpstrj" \
  --file_xyz=$name".xyz" \
  --rc="2.88,3.11" \
  --drc=0.1 \
  --bond=1 \
  --angle=0 \
  --dihed=0 \
  --cis_trans=0 \
  --diff_bond_len=0 \
   > $name".log"


# Monolayer w/ 3x3 replications, bond type depends on bond length
xrep=3
yrep=3
name="o.type_depends_blen"$xrep"x"$yrep
python ../../crystal_builder.py \
  i.basis_Bi2Se3_1layer.dat "3,3,1" \
  --file_pos=$name".dat" \
  --file_dump=$name".lammpstrj" \
  --file_xyz=$name".xyz" \
  --rc="2.88,3.11" \
  --drc=0.1 \
  --bond=1 \
  --angle=0 \
  --dihed=0 \
  --cis_trans=0 \
  --diff_bond_len=1 \
  --diff_bond_fmt="%.2f" \
   > $name".log"
