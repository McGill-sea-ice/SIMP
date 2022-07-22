# Sea Ice Modelling Particles

The Sea Ice Modelling Particles (SIMP) is a FORTRAN sea-ice model based on
Smoothed Particle Hydrodynamics developped at McGill University under 
supervision of Bruno Tremblay by Oreste Marquis for its master degree. The 
current version of the model can simulate the Viscous-Plastic sea-ice 
dynamics in varying domain. A description of the program subroutine tested 
are showed in (Marquis et al. 2022). However, some of the subroutine in the
present model are still under development or have been developped for
reasearch purposes (and remain in the model for flexibility in the 
computation), but are no more used. 


## Getting Started

These instructions will give you a copy of the project up and running on
your local machine for development and testing purposes.

### Prerequisites

- Program is written in FORTRAN2018 and has been developed and compiled using the GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.1) 9.4.0.
- The simulation visualization plots have been developped using Python 3.8.8.

Make sure to have those installed or compatible versions on your machine to
use SIMP. 

### Installing

Installing gfortran :
https://fortran-lang.org/learn/os_setup/install_gfortran

Installing Python:
https://docs.python.org/3/using/index.html

Installing SIMP:
"git clone https://github.com/OresteMarquis/SIMP.git"

## Running SIMP
A Makefile should be present in the parent directory of SIMP. Make sure to be
is this directory when running the program. There 3 main commands and 
multiple more specific command use by the main ones (see the Makefile for
more details on the specific).

The first command is to run the simulation and is the "all" command (it can 
also be used to start the with the no hang up function using 
"all_background"):

`make all_background OUT=< name of the run > NPROC=< # Cores >`.

There is the visualization command:

`make graphics`

The last command is to clean up the various folder:

`make clean_all`.


### Test disclaimer
No standardize test have been made to ensure that the programme is not
broken when new commit where push to the main branch. Therefore, 
some function may have been created and use correctly in the past but are now broken without knowing. This issue is important for certain subroutine that have been abandonned in favor of other used for the run in Marquis et Al. 2022
paper. The subroutines that may have issues are highlighted using *** in the 
subroutine classification sections.

### Subroutine classification
Here is a brief description of the folder and files. Note that to change run
paramters only the following files need to be changed : globals.f90, data_variables.f90, domain.f90, boundary_domain.f90, 
initialize.f90, wind_forcing.f90 and water_drag.f90.

#### Programme architecture

	- SIMP/
	    - Makefile
	    - fortran/
		- template.txt
		- bin/
		- data/  
		- obj/  
		- src/  
		    - simp.f90
		    - ext_force/  
			- boundary_force.f90 
			- water_drag.f90 
			- wind_forcing.f90
		    - in_out/  
			- boundaries.f90 
			- initialize.f90 
			- output.f90
		    - integration/ 
			- single_step.f90 
			- time_integration.f90 
			- time_step.f90
		    - int_force/ 
			- bulk_viscosity.f90 
			- eos.f90 
			- internal_stress.f90 
			- shear_viscosity.f90 
			- stress_tensor.f90
		    - modules/ 
			- boundary_domain.f90 
			- data_variables.f90 
			- domain.f90 
			- globals.f90 
			- interfaces.f90 
			- mpi_variables.f90  
			- tree_module.f90
		    - neighbor/ 
			- NNPS.f90
		    - phy_quantities/ 
			- density.f90 
			- mass_evolution.f90 
			- sound_speed.f90
		    - sph/
			- boundary_kernel.f90 
			- kernel.f90 
			- regularization.f90 
			- smoothing_length.f90 
			- vector_divergence.f90 
			- vector_gradient.f90
	    - python/
		- plots/
		- src/
		    - data_visualization.py 
		    - plotting_functions.py

#### Folders

	- fortran/        : Folder containing all the fortran related folder and files.
	- bin/            : Folder containing the binary file produce when running SIMP.
	- data/ 	  : Folder containing the output data produce by running SIMP.
	- obj/  	  : Folder conatining the object files produce when compiling SIMP.
	- src/ 		  : Folder containing all the source file to run SIMP.
	- ext_force/  	  : Folder containing all the subroutine related to the external forcing in the ice momentum equation.
	- in_out/  	  : Folder containing the subroutines to initialize the run and to store the output data.
	- integration/ 	  : Folder containing the subroutines to integrate the equation in time.
	- int_force/ 	  : Folder containing the subroutines involving the computation of the rheology term in the ice dynamic equation.
	- modules/ 	  : Folder containing the various module to run the programme and the variable declaration.
	- neighbor/ 	  : Folder containing the neighbour search algorithms.
	- phy_quantities/ : Folder containing subroutine to keep track of physical quantities.
	- sph/		  : Folder containing the subroutine specifically related to SPH computation.
	- python/	  : Folder containing all the python related folder and files.
	- plots/	  : Folder containing the output of SIMP run visualization.
	- src/		  : Folder containing source file for SIMP run visualization.


#### Files and subroutines


	- Makefile              : File containing the shell script to run, clean and visualize SIMP.
	- template.txt 		: File use as a template to create new subroutines
	- simp.f90 		: File containing the top layer of the SIMP programme.
	- boundary_force.f90 	: File containing the subroutines to compute normal and tangential force associate to boundaries.
				Contains :
					- subroutine boundary_force
					- subroutine bf_Monaghan
					- subroutine bf_Wang ***
					- subroutine bf_Marquis
	- water_drag.f90 	: File containing the subroutine to compute the water drag force.
				Contains :
					- subroutine water_drag
	- wind_forcing.f90 	: File containing the subroutine to compute the wind drag force.
				Contains :
					- subroutine wind_forcing
	- boundaries.f90 	: File containing the suboutines to initialize the boundaries.
				Contains :
					- subroutine boundaries
					- subroutine boundary_placing_1
					- subroutine boundary_placing_2
					- subroutine boundary_reading ***
					- subroutine boundary_normal ***
	- initialize.f90 	: File containing the suboutines to initialize the ice particles.
				  Contains :
					- subroutine initialize
					- subroutine particle_placing
					- subroutine particle_reading ***
	- output.f90    	: File containing the subroutines to output the data.
				Contains :
					- subroutine output
	- single_step.f90 	: File containing the subroutines to do a single step in time.
				Contains :
					- subroutine single_step
	- time_integration.f90  : Fle containing the subroutines to do the different integration shceme.
				Contains :
					- subroutine time_integration
					- subroutine euler_scheme
					- subroutine predictor_corrector
					- subroutine leapfrog_scheme ***
	- time_step.f90  	: File containting the subroutines for the time step evolution.
				Contains : 
					- subroutine time_step
					- subroutine VP_time_step
					- subroutine Monaghan_time_step
	- bulk_viscosity.f90 	: File containing the subroutines for the bulk visocsity computations of different rheology.
				Contains :
					- subroutine bulk_viscosity
					- subroutine elliptical_bulk_viscosity
					- subroutine MC_bulk_viscosity ***
	- shear_viscosity.f90 	: File containing the subroutines for the shear visocsity computations of different rheology.
				Contains :
					- subroutine shear_viscosity
					- subroutine elliptical_shear_viscosity
					- subroutine MC_bulk_viscosity ***
	- eos.f90 		: File containing the subroutines to compute the pressure term if the rheology equation.
				Contains :
					- subroutine eos
					- subroutine hibler_eos
					- subroutine gutfraind_eos
					- subroutine kreyscher_eos
	- internal_stress.f90   : File containing the subroutine to compute the whole rheology term.
				Contains :
					- subroutine internal_stress
					- subroutine stress_1
					- subroutine stress_2
	- stress_tensor.f90     : File containing the subroutine to compute each component of the stress tensor.
				Contains : 
					- subroutine stress_tensor
	- boundary_domain.f90   : File to initialize the boundary domain shape.
				Contains :
					- module boundary_domain
	- data_variables.f90 	: File to initialize the physical variables and parameter values.
				Contains :
					- module data_variables
	- domain.f90 		: File to initialize the ice domain shape.
				Contains :
					- module domain
	- globals.f90 		: File to chose the simulations properties and subroutine to uses.
				Contains :
					- module globals
	- interfaces.f90 	: File creating the interface to ensure the protability of SIMP.
				Contains :
					- module interfaces
	- mpi_variables.f90     : File to initialize the variables related to OpenMP parallelization.
				Contains :
					- module mpi_variables
	- tree_module.f90       : File containting all the subroutine to make the tree NNPS algorithm.
				Contains :
					- module tree_module
					- subroutine build_tree
					- subroutine p_interaction_tree
					- subroutine q_interaction_tree
					- subroutine interact
					- subroutine tree_loop
					- subroutine delete_tree
					- subroutine boundary_tree_loop
					- subroutine bp_interaction_tree
					- subroutine b_interact
	- NNPS.f90              : File containing the various subroutines to compute the nearest neighbour search.
				Contains : 
					- subroutine NNPS
					- subroutine direct_find
					- subroutine b_tree
					- subroutine bucket_search
					- subroutine interaction_stats
	- density.f90 		: File containing the subroutines to compute the density according various definition .
				Contains : 
					- subroutine density
					- subroutine sum_density ***
					- subroutine con_density ***
					- subroutine diagnostic_density
	- mass_evolution.f90    : File containing the subroutines to compute the particle mass according various definition .
				Contains : 
					- subroutine mass_evolution
					- subroutine mass_update_from_area
	- sound_speed.f90       : File to compute the speed of sound in the ice.
				Contains : 
					- subroutine sound_speed
					- subroutine c_yang
	- boundary_kernel.f90   : File to compute the boundary particle kernel interaction.
				Contains : 
					- subroutine boundary_kernel
					- subroutine b_Quadratic ***
					- subroutine b_Gaussian
					- subroutine b_quintic_spline
					- subroutine b_Wendland_C2 ***
					- subroutine b_Wendland_C4
					- subroutine b_Wendland_C6
	- kernel.f90            : File to compute the ice particle kernel interaction.
				Contains : 
					- subroutine kernel
					- subroutine Gaussian
					- subroutine quintic_spline
					- subroutine Quadratic ***
					- subroutine Wendland_C2 ***
					- subroutine Wendland_C4
					- subroutine Wendland_C6
					- subroutine Marquis ***
					- subroutine Marquis2 ***
	- regularization.f90    : File to compute stabilization technic possible in SPH.
				Contains :
					- subroutine regularization
					- subroutine XSPH ***
					- subroutine art_visc_HK ***
					- subroutine art_visc_Mon ***
	- smoothing_length.f90  : File containing the subroutine to compute the smoothing length evolution.
				Contains :
					- subroutine smoothing_length
					- subroutine Benz_hsml ***
					- subroutine Marquis_hsml
	- vector_divergence.f90 : File containing the subroutine to compute a vector divergence in SPH.
				Contains :
					- subroutine vector_divergence
	- vector_gradient.f90   : File containing the subroutine to compute a vector gradient in SPH.
				Contains :
					- subroutine vector_gradient
	- data_visualization.py : File use to load the data and call the desired plots.
	- plotting_functions.py : File containing the visualization functions.


## Contributing

Everyone can contribute, but to do so please contact Prof. Bruno Tremblay at bruno.tremblay@mcgill.ca.


## Versioning

For now no versionning has been implemented.

## Authors

  - **Oreste Marquis** -  Coded and developped the whole program between 2021 and 2022.


