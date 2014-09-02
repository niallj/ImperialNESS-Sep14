#PBS -q tng
#PBS -l select=1:ncpus=1
#PBS -N analyse
#PBS -l walltime=01:00:00

#uncomment as appropriate
PROGRAM=$HOME/ImperialNESS-Sep14/ConfigurationalTemperature/analysis/C++/analyse 
#PROGRAM=$HOME/ImperialNESS-Sep14/ConfigurationalTemperature/analysis/FORTRAN/analyse

#copy the trajectory file to the temporary directory
#(make sure the trajectory file is located in the same folder as this script!)
cp $PBS_O_WORKDIR/out.lammpstrj.bin $TMPDIR

#execute the C++/FORTRAN analysis code
$PROGRAM out.lammpstrj.bin

#copy the temperature output back to the initial directory
mv $TMPDIR/temperatures.dat $PBS_O_WORKDIR/temperatures.dat
