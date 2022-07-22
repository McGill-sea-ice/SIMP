!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 Module for parameters and constants used                 !
!                   in the entire SPH software packages.                   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module globals
    implicit none

    !!!                       Restriction parameters                      !!!
    !dim : Dimension of the problem 

    integer :: dim
    parameter ( dim = 2)
    
    !!!             Set maximums before hand to manage memory             !!!
    !maxn           : Maximum number of particles
    !maxn           : Maximum number of boundary particles
    !max_interation : Maximum number of interaction pair
    !max_ind_int    : Maximum number of interaction for a single particle

    integer  :: maxn, maxnb, max_interaction, max_ind_int
    parameter( maxn = 80000,                                                &
               maxnb = 5000,                                               &
               max_ind_int = 200,                                          &
               max_interaction = 40 * maxn )

    !resolution : Desired initial spacing between particle               [m]

    DOUBLE PRECISION :: resolution
    parameter( resolution = 7500.)

    !integration_time : Desired time for the simulation                  [s]

    DOUBLE PRECISION :: integration_time = 672. * 3600.

    !!!                    Control parameters for input                  !!!
    !   spin_up : .TRUE. : Start run from input files
    !   path  : Directory path to particle information files

    CHARACTER(len=255) :: path
    logical :: spin_up
    parameter (spin_up =.FALSE.)

    !!!                    Control parameters for output                 !!!
    !   int_stat : .TRUE. : Print statistics about SPH particle interactions
    !                       including virtual particle information.
    !   print_step: Print Timestep (On Screen)                   [iteration]
    !   save_step : Save Timestep (To Disk File)                         [s]

    logical :: int_stat
    parameter (int_stat =.TRUE.)

    integer  :: print_step
    REAL(KIND=8) :: save_step
    parameter ( print_step = 2000 , &
                save_step = 1 * 3600 )

    !!!                        Integration scheme                         !!!
    !   int_algorithm = 1 : Euler scheme
    !                   2 : Predictor-Corrector (Hosseini, 2019)
    !                   3 : Leap-Frog (Monaghan, 2005)

    integer :: int_algorithm
    parameter(int_algorithm = 2)


    !!!        Nearest neighboring particle searching (nnps) method       !!!
    !   nnps_algorithm = 1 : Direct find O(N^2)
    !                    2 : Bucket search O(N + N^2/K + K) 
    !                    3 : B-Tree algorithm 

    integer :: nnps_algorithm
    parameter(nnps_algorithm = 2)
    
    !!!             Smoothing length evolution (sle) algorithm            !!!
    !   sle = 0 : Keep unchanged
    !         1 : dh/dt = (-1/dim)*(h/rho)*(drho/dt) (Benz 1990)
    !         2 : h = sqrt(m/rho)*hsml_multiplicity (Marquis 2022)

    integer :: sle
    parameter(sle = 2)

    !!!                    Timestep evolution criteria                   !!!
    !   tse = 0 : Keep unchanged
    !         1 : SPH usual (Monaghan 2005)
    !         2 : VP material
    
    integer :: tse
    parameter(tse = 2)

    !!!                       Rheology implementation                     !!!
    !   rheology = 1 : VP elliptical yield curve and normal flow rule
    !              2 : VP Mohr-Coulomb yield curve and shear flow rule
    
    integer :: rheology
    parameter(rheology = 1)

    !!!                 Artificial viscosity implementation               !!!
    !   art_visc_alg = 1 : Herquist and Katz (1989)
    !                = 2 : Monaghan (1989),

    integer :: art_visc_alg
    parameter(art_visc_alg = 1)

    !!!                       Boundary initial setup                      !!!
    !   b_setup = 0 : Read from input file
    !           = 1 : One layer
    !           = 2 : Multiple layers

    integer :: b_setup
    parameter(b_setup = 1)

    !!!                     Boundary force algorithm                      !!!
    !   b_setup = 1 : (Monaghan 2009)
    !           = 2 : (Wang 2014)
    !           = 3 : Artificial stopping (Marquis 2022)

    integer :: bfa
    parameter(bfa = 1)

    !!!                      Equation of state algorithm                  !!!
    !   eos_alg = 1 : P = P_star * h * exp( -C * (1 -A)) (Hibler 1979)
    !           = 2 : P = g(rho - rho_0) / (rho_max - rho) (Gutfraind 1997)
    !           = 3 : P = P' * D / (D') (Kreyscher 2000)

    integer :: eos_alg
    parameter(eos_alg = 3)

    !!!                        Particle approximation                     !!!
    !   pa_sph = 1 : (p(i)/rho(i)**2+p(j)/rho(j)**2
    !            2 : (p(i)+p(j))/(rho(i)*rho(j))

    integer :: pa_sph
    parameter(pa_sph = 1)
    
    !!!                     Smoothing kernel function                     !!!
    !   skf = 1 : Quadratic (Johnson 1996)
    !       = 2 : Gaussian (Gingold and Monaghan 1981)
    !       = 3 : Quintic spline (Morris 1997)
    !       = 4 : Wendland_C2 (Dehnen, 2012)
    !       = 5 : Wendland_C4 (Dehnen, 2012)
    !       = 6 : Wendland_C6 (Dehnen, 2012)
    !       = 7 : Step kernel (Marquis, 2022)
    !       = 8 : Triangle kernel (Marquis, 2022)

    integer :: skf
    parameter(skf = 6)

    !!!                         Density definition                        !!!
    !   den_alg = 1 : Normalized density (Randles and Libersky, 1996)
    !           = 2 : Continuity density (Monaghan 2012) 
    !           = 3 : Diagnostic (Staroszyck 2017)

    integer :: den_alg
    parameter(den_alg = 3)

    
    !!!                  Switches for different scenarios                 !!!
    !   viscosity           : .TRUE. : Consider viscosity,
    !                         .FALSE. : No viscosity.
    !   pressure_force      : .TRUE. : Consider pressure between particles,
    !                         .FALSE. : No pressure between particles.
    !   wind_force          : .TRUE. : Consider wind forcing
    !                         .FALSE. : No wind forcing
    !   water_force         : .TRUE. : Consider water drag
    !                         .FALSE. : No water drag
    !   artificial_visc     : .TRUE. : Add artificial viscosity
    !                         .FALSE. : No artificial viscosity
    !   boundary_pressure   : .TRUE. : Add boundary normal force
    !                         .FALSE. : No boundary normal force
    !   boundary_friction   : .TRUE. : Add boundary tangential force
    !                         .FALSE. : No boundary tangential force
    !   xsph_regularization : .TRUE. : Add xsph regularization
    !                         .FALSE. : No xsph regularization
 
    logical :: viscosity, wind_force, water_force, pressure_force,         &
               artificial_visc, boundary_pressure, boundary_friction,      &
               xsph_regularization

    parameter (pressure_force    = .TRUE.,                                 &
               viscosity         = .TRUE.,                                 &
               wind_force        = .TRUE.,                                 &
               water_force       = .TRUE.,                                 &
               boundary_pressure = .TRUE.,                                 &
               boundary_friction = .FALSE.,                                 &
               artificial_visc   = .FALSE.,                                &
               xsph_regularization = .FALSE.)

    !!!                         Systeme variables                         !!!
    !   str : dummy variable for getcwd function
    !   cwd : current working directory
    !   output_folder : name of the folder where the output should be 
                
    CHARACTER(len=255) :: str, output_folder
    CHARACTER(len=:), allocatable :: cwd

    !!! Mathematics constants !!!
    DOUBLE PRECISION :: pi 
    parameter (pi = 3.1415926535897932D0)

    contains
        !Get working directory
        subroutine sys_variables
            implicit none 

            CALL get_environment_variable("PWD", str)
            cwd = trim(str)
        end subroutine sys_variables

end module globals



