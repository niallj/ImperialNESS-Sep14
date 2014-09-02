#PBS -q tng
#PBS -l select=1:ncpus=8
#PBS -l walltime=06:00:00

module load intel-suite mpi lammps/4Feb14

cp data.soret $TMPDIR
cp soret.in $TMPDIR

mpiexec lammps -in soret.in

mv $TMPDIR $PBS_O_WORKDIR/run
