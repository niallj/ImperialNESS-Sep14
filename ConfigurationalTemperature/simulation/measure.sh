#PBS -q pqchem
#PBS -N Measure
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=8

module load intel-suite lammps/4Feb14 mpi

cp $PBS_O_WORKDIR/measure.in $TMPDIR
cp $PBS_O_WORKDIR/settings.in $TMPDIR
cp $PBS_O_WORKDIR/equilibration_run/equilibrated.restart $TMPDIR

mpiexec lammps -in measure.in

mv $TMPDIR $PBS_O_WORKDIR/measure
