#PBS -N equilibrate
#PBS -l select=1:ncpus=1
#PBS -q tng

EXEC=$HOME/ImperialNESS-Sep14/synthetic_NEMD/bin/sllod_2d.exe

#copy the input file from the directory in which the script was executed to the temporary storage
cp $PBS_O_WORKDIR/input.dat $TMPDIR

#run the sllod executable
$EXEC

#copy the files generated from the temportary directory back to wherever we started from
cp -r $TMPDIR $PBS_O_WORKDIR/equilibrated
