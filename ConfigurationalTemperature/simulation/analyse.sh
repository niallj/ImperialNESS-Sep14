#PBS -q tng
#PBS -l select=1:ncpus=1
#PBS -N analyse
#PBS -l walltime=01:00:00

#uncomment as appropriate
PROGRAM=$HOME/ImperialNESS-Sep14/ConfigurationalTemperature/analysis/C++/analyse 
#PROGRAM=$HOME/ImperialNESS-Sep14/ConfigurationalTemperature/analysis/FORTRAN/analyse

cp $PBS_O_WORKDIR/out.lammpstrj.bin $TMPDIR

$PROGRAM out.lammpstrj.bin

mv $TMPDIR/temperatures.dat $PBS_O_WORKDIR/temperatures.dat
