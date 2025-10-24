#!/bin/bash
#PBS -N Sn3O4
#PBS -e job.err
#PBS -o job.o
#PBS -q workq
#PBS -l nodes=1:ppn=16

echo "Loading Siesta..."
module load siesta-5.2.2

echo "Checking Siesta path..."
which siesta
ls -l $(which siesta)

echo "Setting stack size..."
ulimit -s unlimited

echo "Setting MPI fabrics..."
export I_MPI_FABRICS=shm
export OMP_NUM_THREADS=1
#export I_MPI_HYDRA_TOPOLIB=ipl

echo "Changing to work directory..."
cd $PBS_O_WORKDIR
pwd

mpirun -hostfile $PBS_NODEFILE -np $NCPUS siesta calc.fdf --o calc.out





