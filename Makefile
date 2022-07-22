############## VARIABLES #############
# Compiler settings
FC = gfortran 
FFLAGS = -Wall -std=f2018 -Wextra -Werror -pedantic-errors -g -fopenmp

#Name of the default output excutable
OUT = simp.out

# Path to main folders
SRC_PATH  = $(shell pwd)/fortran/src/
OBJ_PATH  = $(shell pwd)/fortran/obj/
BIN_PATH  = $(shell pwd)/fortran/bin/
DATA_PATH = $(shell pwd)/fortran/data/

# Variable with source file names
SRC_ext_force      = $(wildcard $(SRC_PATH)ext_force/*.f90)
SRC_in_out         = $(wildcard $(SRC_PATH)in_out/*.f90)
SRC_modules        = $(wildcard $(SRC_PATH)modules/*.f90)
SRC_int_force      = $(wildcard $(SRC_PATH)int_force/*.f90)
SRC_integration    = $(wildcard $(SRC_PATH)integration/*.f90)
SRC_neighbor       = $(wildcard $(SRC_PATH)neighbor/*.f90)
SRC_phy_quantities = $(wildcard $(SRC_PATH)phy_quantities/*.f90)
SRC_sph            = $(wildcard $(SRC_PATH)sph/*.f90)
SRC_main           = $(SRC_PATH)simp.f90

# Variable with obj file names
OBJ_modules        = $(patsubst $(SRC_PATH)modules/%.f90,$(OBJ_PATH)%.o,$(SRC_modules))
OBJ_ext_force      = $(patsubst $(SRC_PATH)ext_force/%.f90,$(OBJ_PATH)%.o,$(SRC_ext_force))
OBJ_in_out         = $(patsubst $(SRC_PATH)in_out/%.f90,$(OBJ_PATH)%.o,$(SRC_in_out))
OBJ_int_force      = $(patsubst $(SRC_PATH)int_force/%.f90,$(OBJ_PATH)%.o,$(SRC_int_force))
OBJ_integration    = $(patsubst $(SRC_PATH)integration/%.f90,$(OBJ_PATH)%.o,$(SRC_integration))
OBJ_neighbor       = $(patsubst $(SRC_PATH)neighbor/%.f90,$(OBJ_PATH)%.o,$(SRC_neighbor))
OBJ_phy_quantities = $(patsubst $(SRC_PATH)phy_quantities/%.f90,$(OBJ_PATH)%.o,$(SRC_phy_quantities))
OBJ_sph            = $(patsubst $(SRC_PATH)sph/%.f90,$(OBJ_PATH)%.o,$(SRC_sph))
OBJ_main           = $(OBJ_PATH)simp.o

# Object variable
ALL_OBJ  = $(OBJ_modules)
ALL_OBJ += $(OBJ_ext_force)
ALL_OBJ += $(OBJ_in_out)
ALL_OBJ += $(OBJ_int_force)
ALL_OBJ += $(OBJ_integration)
ALL_OBJ += $(OBJ_neighbor)
ALL_OBJ += $(OBJ_phy_quantities)
ALL_OBJ += $(OBJ_sph)

# Executable variable
OUT = simp.out
EXEC = $(BIN_PATH)$(OUT)
NPROC = 4
export OMP_NUM_THREADS=$(NPROC)



############## General commands #############
all: create_exe

	mkdir $(shell pwd)/fortran/data/$(OUT)
	rsync -r $(SRC_PATH) $(shell pwd)/fortran/data/$(OUT)/model/

	./fortran/bin/$(OUT)

all_background: create_exe 

	mkdir $(shell pwd)/fortran/data/$(OUT)
	rsync -r $(SRC_PATH) $(shell pwd)/fortran/data/$(OUT)/model/
	nohup ./fortran/bin/$(OUT) > $(OUT).txt &

create_exe : $(EXEC)

#Compilation call for individual uses
compile_modules: $(OBJ_modules)
compile_ext_force: $(OBJ_ext_force)
compile_in_out: $(OBJ_in_out)
compile_int_force: $(OBJ_int_force)
compile_integration: $(OBJ_integration)
compile_neighbor: $(OBJ_neighbor)
compile_phy_quantities: $(OBJ_phy_quantities)
compile_sph: $(OBJ_sph)
compile_main: $(OBJ_main)

#Compile all files
compile_all: compile_modules compile_ext_force compile_in_out compile_int_force compile_integration compile_neighbor compile_phy_quantities compile_sph compile_main

#Create plots
graphics:
	python3 $(shell pwd)/python/src/data_visualization.py

#Cleanup
clean_obj:
	rm -f $(OBJ_PATH)/*

clean_exe:
	rm -f $(BIN_PATH)/*

clean_data:
	rm -f -r $(DATA_PATH)/*

clean_plots:
	rm -r -f $(shell pwd)/python/plots/*

clean_program: clean_obj clean_exe clean_data

clean_all: clean_program clean_plots




############## DEPENDENCIES #############
$(EXEC) : $(OBJ_main) $(ALL_OBJ)
	$(FC) $(FFLAGS) -o $@ $+

#Compile module in this specific order
$(OBJ_modules): 
	$(FC) $(FFLAGS) -J$(OBJ_PATH) -c $(SRC_PATH)modules/globals.f90 -o $(OBJ_PATH)globals.o
	$(FC) $(FFLAGS) -J$(OBJ_PATH) -c $(SRC_PATH)modules/data_variables.f90 -o $(OBJ_PATH)data_variables.o
	$(FC) $(FFLAGS) -J$(OBJ_PATH) -c $(SRC_PATH)modules/domain.f90 -o $(OBJ_PATH)domain.o
	$(FC) $(FFLAGS) -J$(OBJ_PATH) -c $(SRC_PATH)modules/boundary_domain.f90 -o $(OBJ_PATH)boundary_domain.o
	$(FC) $(FFLAGS) -J$(OBJ_PATH) -c $(SRC_PATH)modules/interfaces.f90 -o $(OBJ_PATH)interfaces.o
	$(FC) $(FFLAGS) -J$(OBJ_PATH) -c $(SRC_PATH)modules/tree_module.f90 -o $(OBJ_PATH)tree_module.o
	$(FC) $(FFLAGS) -J$(OBJ_PATH) -c $(SRC_PATH)modules/mpi_variables.f90 -o $(OBJ_PATH)mpi_variables.o

#Compile external force folder
$(OBJ_PATH)%.o : $(SRC_PATH)ext_force/%.f90
	$(FC) $(FFLAGS) -I$(OBJ_PATH) -o $@ -c $<

#Compile input/output folder
$(OBJ_PATH)%.o : $(SRC_PATH)in_out/%.f90
	$(FC) $(FFLAGS) -I$(OBJ_PATH) -c $< -o $@

#Compile internal force folder
$(OBJ_PATH)%.o : $(SRC_PATH)int_force/%.f90
	$(FC) $(FFLAGS) -I$(OBJ_PATH) -c $< -o $@

#Compile integration folder
$(OBJ_PATH)%.o : $(SRC_PATH)integration/%.f90
	$(FC) $(FFLAGS) -I$(OBJ_PATH) -c $< -o $@

#Compile neighbor folder
$(OBJ_PATH)%.o : $(SRC_PATH)neighbor/%.f90
	$(FC) $(FFLAGS) -I$(OBJ_PATH) -c $< -o $@

#Compile physical quantities folder
$(OBJ_PATH)%.o : $(SRC_PATH)phy_quantities/%.f90
	$(FC) $(FFLAGS) -I$(OBJ_PATH) -c $< -o $@

#Compile sph quantities folder
$(OBJ_PATH)%.o : $(SRC_PATH)sph/%.f90
	$(FC) $(FFLAGS) -I$(OBJ_PATH) -c $< -o $@

#Compile main program
$(OBJ_main) : $(SRC_main) $(ALL_OBJ)
	$(FC) $(FFLAGS) -I$(OBJ_PATH) -c $< -o $@
