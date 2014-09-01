#makefile for surface_tension_2d
#dont forget to use tabs!!
#spaces are significant!

# define default settings
SHELL = /bin/sh
BINROOT =$(PWD)
EX = sllod_3d.exe

#---- intel F95 compiler
FC = gfortran -O3  
#FC = gfortran   
  
SWITCH = -g  

OBJECTS = header_file.o sllod_3d.o read_inputdata.o read_config.o rhs.o rk4.o\
	write_config.o cellvect.o cartfrac.o fraccart.o initial.o cellmap.o cellpairs.o\
	averager.o pair_energy.o  pair_force.o setup_linkcell.o snap.o pbcs.o\
	analyse_properties.o invert.o maxwell.o 

# Makefile
$(EX): $(OBJECTS)
	$(FC) -o $(EX) $(SWITCH) $(OBJECTS)
header_file.mod: header_file.o header_file.f90

	$(FC) -c $(SWITCH) header_file.f90
header_file.o: header_file.f90

	$(FC) -c $(SWITCH) header_file.f90
sllod_3d.o: header_file.mod sllod_3d.f90

	$(FC) -c $(SWITCH) sllod_3d.f90
read_inputdata.o: header_file.mod read_inputdata.f90 

	$(FC) -c $(SWITCH) read_inputdata.f90

read_config.o: header_file.mod read_config.f90 
	$(FC) -c $(SWITCH) read_config.f90

write_config.o: header_file.mod write_config.f90 
	$(FC) -c $(SWITCH) write_config.f90

rhs.o: header_file.mod rhs.f90 
	$(FC) -c $(SWITCH) rhs.f90

cellmap.o: header_file.mod cellmap.f90 
	$(FC) -c $(SWITCH) cellmap.f90

cellpairs.o: header_file.mod cellpairs.f90 
	$(FC) -c $(SWITCH) cellpairs.f90

setup_linkcell.o: header_file.mod setup_linkcell.f90 
	$(FC) -c $(SWITCH) setup_linkcell.f90

initial.o: header_file.mod initial.f90 
	$(FC) -c $(SWITCH) initial.f90

averager.o: header_file.mod averager.f90 
	$(FC) -c $(SWITCH) averager.f90

pbcs.o: header_file.mod pbcs.f90 
	$(FC) -c $(SWITCH) pbcs.f90

pair_energy.o: header_file.mod pair_energy.f90 
	$(FC) -c $(SWITCH) pair_energy.f90

pair_force.o: header_file.mod pair_force.f90 
	$(FC) -c $(SWITCH) pair_force.f90

fraccart.o: header_file.mod fraccart.f90 
	$(FC) -c $(SWITCH) fraccart.f90

cartfrac.o: header_file.mod cartfrac.f90 
	$(FC) -c $(SWITCH) cartfrac.f90

analyse_properties.o: header_file.mod analyse_properties.f90 
	$(FC) -c $(SWITCH) analyse_properties.f90

cellvect.o: header_file.mod cellvect.f90 
	$(FC) -c $(SWITCH) cellvect.f90

snap.o: header_file.mod snap.f90 
	$(FC) -c $(SWITCH) snap.f90

rk4.o: rk4.f90 
	$(FC) -c $(SWITCH) rk4.f90

invert.o: invert.f90 
	$(FC) -c $(SWITCH) invert.f90

maxwell.o: maxwell.f90 
	$(FC) -c $(SWITCH) maxwell.f90
%: %.o
	$(FC) -c $(SWITCH) $< 

clean:
	rm -f *.o *.mod *.MOD
