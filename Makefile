# ----------------------------------------------------------------
# Project: Complete One loop Predictions for B->mu+mu- in the MSSM
# Please cite: Dedes, Rosiek, Tanedo (2007)
# ----------------------------------------------------------------
# Main file for B->mu+mu- calculation
# filename: Makefile
# Current version: 1.00 (SM result)
# By Flip Tanedo (philip.tanedo@gmail.com)
# Date: 29 November 2007
# Latest version available at: n/a
# ----------------------------------------------------------------
# This is the README FILE
# ----------------------------------------------------------------
# RECORD OF REVISIONS
#	Date	 Programmer	Description
#	====	 ==========	==================================
#	18/11/07 Flip Tanedo	First version (running)
#	29/11/07 Flip Tanedo	SM result
# ----------------------------------------------------------------
# LIST OF FILES (on/off)
# Go?   file
# ===   ==========================================================
# on	driver.f95			Runs everything
# on	program_parameters.f95		Non-physics parameters
# on	physics_data.f95		Physics parameters
# on	tuning_data.f95			"Unphysical" parameters
# on	integrals.f95			Passarino-Veltman 
# on	sorter.f95			For 'Integrals'
# on 	errorlog.f95			For error tracking
# on	diagrams_box.f95		Box diagrams
# on	diagrams_dselfenergy.f95	Self energies
# on	diagrams_zpenguin.f95		Z-penguins
# on	amplitude.f95			Wilson Coefficients
# on	branchingratio.f95		Branching Ratio
# ----------------------------------------------------------------
# This is my makefile
# Things to do: Map dependencies, clean up
# http://www.tlug.org.za/old/csslug/writing_makefiles.html
#
# OLD MAKEFILE:
#
#all: 
#	gfortran driver.f95 program_parameters.f95 physics_data.f95 integrals.f95  sorter.f95 tuning_data.f95 diagrams_box.f95 diagrams_zpenguin.f95 diagrams_dselfenergy.f95 amplitude.f95 branchingratio.f95 errorlog.f95 temporary.f95 -o Bmumu
#	./Bmumu
#
#
#
all:  Bmumu

FOPT = -g
F95 = gfortran
PP = program_parameters
PD = physics_data
TD = tuning_data
DB = diagrams_box
DD = diagrams_dselfenergy
DZ = diagrams_zpenguin
BR = branchingratio
COMP = $(F95) $(FOPT)

program_parameters.o: $(PP).f95
	$(COMP) -c $(PP).f95

physics_data.o: $(PD).f95 $(PP).o
	$(COMP) -c $(PD).f95

tuning_data.o: $(TD).f95 $(PP).o
	$(COMP) -c $(TD).f95

integrals.o: integrals.f95 $(PP).o $(TD).o sorter.o errorlog.o
	$(COMP) -c integrals.f95 

sorter.o: sorter.f95 $(PP).o $(TD).o
	$(COMP) -c sorter.f95

errorlog.o: errorlog.f95 $(PP).o
	$(COMP) -c errorlog.f95

diagrams_box.o: $(DB).f95 $(PP).o $(PD).o integrals.o $(TD).o 
	$(COMP) -c $(DB).f95

diagrams_dselfenergy.o: $(DD).f95 $(PP).o $(PD).o integrals.o $(TD).o 
	$(COMP) -c $(DD).f95

diagrams_zpenguin.o: $(DZ).f95 $(PP).o $(PD).o integrals.o $(TD).o 
	$(COMP) -c $(DZ).f95

amplitude.o: amplitude.f95 $(PP).o $(PD).o $(DD).o $(DB).o $(DZ).o
	$(COMP) -c amplitude.f95

branchingratio.o: branchingratio.f95 $(PP).o $(PD).o amplitude.o
	$(COMP) -c branchingratio.f95

Bmumu: driver.f95 $(PP).o $(PD).o $(TD).o integrals.o errorlog.o $(DD).o $(DB).o $(DZ).o amplitude.o $(BR).o
	$(COMP) driver.f95 $(PP).f95 $(PD).f95 $(TD).f95 integrals.f95 sorter.f95 errorlog.f95 $(DD).f95 $(DB).f95 $(DZ).f95 amplitude.f95 $(BR).f95 -o Bmumu
	./Bmumu

clean:
	rm -rf *.o *.log *.mod
# type in "make clean" to remove objects, modules, log files (leave behind source only)
