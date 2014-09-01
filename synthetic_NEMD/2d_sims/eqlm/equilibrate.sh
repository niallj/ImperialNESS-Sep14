#PBS -N equilibrate
#PBS -l select=1:ncpus=1
#PBS -q tng

EXEC=$HOME/ImperialNESS-Sep14/synthetic_NEMD/bin/sllod_2d.exe

cp $PBS_O_WORKDIR/input.dat $TMPDIR

$EXEC

cp -r $TMPDIR $PBS_O_WORKDIR
