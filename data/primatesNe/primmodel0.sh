#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1gb
#SBATCH -o primmodel0.out
#SBATCH -e primmodel0.err

../coevol -f -d prim.phy -t prim.rootedtree -cal prim.calib 80 80 -c prim.lht -dsom primmodel0
