#!/bin/bash
#SBATCH --partition=long
#SBATCH --time=5-12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem-per-cpu=1gb
#SBATCH -o primbd.out
#SBATCH -e primbd.err
module purge
module use /softs/modulefiles
module load gcc/4.7.3 mpi/openmpi/1.8.4/openmpi-gcc47

../coevol -f -d prim.phy -t prim.rootedtree -cal prim.calib 90 20 -bd -c hyperprim.lht -dsom primcovbd1 &
../coevol -f -d prim.phy -t prim.rootedtree -cal prim.calib 90 20 -bd -c hyperprim.lht -dsom primcovbd2 &
../coevol -f -d prim.phy -t prim.rootedtree -cal prim.calib 90 20 -bd -c hyperprim.lht -dsom -diag primdiagbd1 &
../coevol -f -d prim.phy -t prim.rootedtree -cal prim.calib 90 20 -bd -c hyperprim.lht -dsom -diag primdiagbd2 &
../nearlyneutral -f -d prim.phy -t prim.rootedtree -cal prim.calib 90 20 -bd -c hyperprim.lht primmechbd1 &
../nearlyneutral -f -d prim.phy -t prim.rootedtree -cal prim.calib 90 20 -bd -c hyperprim.lht primmechbd2 &
../nearlyneutral -f -d prim.phy -t prim.rootedtree -cal prim.calib 90 20 -bd -c hyperprim.lht -ancpol primmechbdancpol1 &
../nearlyneutral -f -d prim.phy -t prim.rootedtree -cal prim.calib 90 20 -bd -c hyperprim.lht -ancpol primmechbdancpol2 &
wait
