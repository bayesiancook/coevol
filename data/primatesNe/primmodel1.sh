#!/bin/bash
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=1gb
#SBATCH -o primmodel1.out
#SBATCH -e primmodel1.err

../nearlyneutral -f -d prim.phy -t prim.rootedtree -cal prim.calib 80 80 -c prim.lht primmodel1
