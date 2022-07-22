!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Module for defining particles variables                 !
!                          used in the whole program                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module data_variables

    !!!                         Module to include:                       !!!
    use globals

    !!!                       Variables definitions:                     !!!
    implicit none

    !!!                       Particles variables                        !!!
    !   ntotal : Total number of particle to use             
    !   x      : Coordinates of all particles                            [m]
    !   vx     : Velocities of all particles                           [m/s]
    !   mass   : Mass of particles                                      [kg]  
    !   rho    : Densities of particles                             [kg/m^2]
    !   p      : Pressure of particules                                 [Pa]
    !   hsml   : Smoothing lengths of particles                          [m]
    !   A      : Sea ice concentration of particles
    !   h      : Mean ice thickness of particles                         [m]
    !   c      : Sound speed                                           [m/s]

    INTEGER          :: ntotal
    DOUBLE PRECISION :: x(dim, maxn),                                      &
                        vx(dim, maxn),                                     &
                        mass(maxn),                                        &
                        rho(maxn),                                         &
                        p(maxn),                                           &
                        hsml(maxn),                                        &
                        A(maxn),                                           &
                        h(maxn),                                           &
                        c(maxn)
    
    !!!                   Boundary particles variables                   !!!
    !   ntotalb : Total number of boundary particles to use       
    !   xb      : Coordinates of all boundary particles particles        [m]
    !   massb   : Mass of boundary particles                            [kg]
    !   b_normal: Normal unit vector to the boundary

    INTEGER          :: ntotalb
    DOUBLE PRECISION :: xb(dim, maxnb),                                    &
                        massb(maxnb),                                      &
                        b_normal(dim,maxnb)

    !!!                       Integration variables                      !!!
    !   dvxdt     : dvx/dt, speed change rate                        [m/s^2]
    !   dxdt      : dx/dt, position change rate                        [m/s]
    !   intdvxdt  : dvx/dt of internal force                         [m/s^2]
    !   extdvxdt  : dvx/dt of external force                         [m/s^2]
    !   dadt      : Concetration rate of change                       [s^-1]
    !   dhdt      : Thickness rate of change                           [m/s]
    !   drhodt    : Density rate of change                        [kg/m^2/s]
    !   itimestep : Current iteration of the simulation
    !   dt        : Time step of the integration                         [s]
    !   current_time : Current time in the integration                   [s]

    INTEGER          :: itimestep
    DOUBLE PRECISION :: dvxdt(dim,maxn),                                   &
                        dxdt(dim,maxn),                                    &
                        dAdt(maxn),                                        &
                        dhdt(maxn),                                        &
                        drhodt(maxn),                                      &
                        dt                                                
    REAL(kind = 8)   :: current_time

    !!!                       Interaction variables                      !!!
    !   p_i      : List of first particle partner of interaction pair
    !   p_j      : List of second particle partner of interaction pair
    !   nip      : number of interacting particle
    !   r_ij     : distance between the interacting particles            [m]
    !   hsml_ij  : mean smoothing length between particle i and j        [m]
    !   bp_i     : List of boundary particle for boundary force interaction
    !   bp_j     : List of particle for boundary force interaction
    !   bnip     : number of interacting particle for boundary force 
    !   br_ij    : distance between particles for boundary force          [m]
    !   max_hsml : Caping of smoothing length

    INTEGER          :: nip,                                               &
                        p_i(max_interaction),                              &
                        p_j(max_interaction),                              &
                        bnip,                                              &
                        bp_i(max_interaction),                             &
                        bp_j(max_interaction)
    DOUBLE PRECISION :: r_ij(max_interaction),                             &
                        br_ij(max_interaction),                            &
                        hsml_ij(max_interaction),                          &
                        max_hsml

    !!!                     Kernel related variables                     !!!
    !   w        : Kernel value for particle i interacting with j    [1/m^2]
    !   dwdx     : Gradient of kernel with respect to particle i     [1/m^3]
    !   scale_h  : Support domain scaling factor
    !   kernel_0 : Value of the kernel at the given particle location
    !   hsml_multiplicity : Average number of particle in support domain 
    !                       in any direction

    INTEGER          :: scale_h,                                           &
                        hsml_multiplicity 
    DOUBLE PRECISION :: w(max_interaction),                                &
                        dwdx(dim,max_interaction),                         &
                        kernel_zero(maxn)
    parameter         ( hsml_multiplicity = 3)

    !!!                     Stress related variables                     !!!
    !   epsilon_11 : XX componante of the strain rate                 [s^-1]
    !   epsilon_12 : XY and YX componante of the strain rate          [s^-1]
    !   epsilon_22 : YY componante of the strain rate                 [s^-1]
    !   zeta       : Bulk viscosity                                   [Pa s]
    !   eta        : Shear viscosity                                  [Pa s]
    !   sigma_11 : XX componante of the stress tensor                [N/m^2]
    !   sigma_12 : XY and YX componante of the stress tensor         [N/M^2]
    !   sigma_22 : YY componante of the stress tensor                [N/M^2]

    DOUBLE PRECISION :: epsilon_11(maxn),                                  &
                        epsilon_12(maxn),                                  &
                        epsilon_22(maxn),                                  &
                        zeta(maxn),                                        &
                        eta(maxn),                                         &
                        sigma_11(maxn),                                    &
                        sigma_12(maxn),                                    &
                        sigma_22(maxn)                                    

    !!!                         Material parameters                      !!!
    !   Y         : Young modulus                                    [N/m^2]
    !   nu        : Poisson ratio
    !   rho_ice   : Density of ice                                  [kg/m^3]
    !   rho_water : Density of water                                [kg/m^3]
    !   rho_air   : Density of air                                  [kg/m^3]
    !   c_wind    : Dimensionless wind drag coefficient 
    !   c_water   : Dimensionless water drag coefficient 
    !   alpha     : Proportionality constant for area vs smoothing length
    !   e         : Ellipse aspect ratio in VP rheology
    !   c_eos     : Empirical constant for ice strength exponential decay
    !   P_star    : Empirical constant for ice strength in compresison  [Pa]
    !   T_star    : Empirical constant for ice strength in tension      [Pa]
    !   S_star    : Empirical constant for ice strength in shear        [Pa]
    !   zeta_max  : Maximum bulk viscosity                            [Pa s]
    !   eta_max   : Maximum shear viscosity                           [Pa s]
    !   delta_min : Viscosity caping                                   [1/s]
    !   alpha_pi, beta_pi : Artificial viscosity coefficient          
    !   kappa_n   : Boundary normal force parameter
    !   kappa_t   : Boundary tangential force parameter
    !   b_length  : Boundary particle length of interaction              [m]
    !   mu        : Slope of the Mohr-Coulomb limbs
    !   k_t       : Tensile strenght factor
    !   phi       : Angle of the Mohr-Coulomb limbs                    [deg]
    !   chi       : XSPH parameter
    !   xi        : Friction coefficient of boundaries

    DOUBLE PRECISION :: Y, nu, c_eos, rho_ice, rho_air, P_star, e, alpha,  &
                        rho_water, c_wind, c_water ,eta_max ,zeta_max,     &
                        alpha_pi, beta_pi, b_length, kappa_n, kappa_t, mu, &
                        T_star, S_star, k_t , phi , chi, xi, delta_min
    parameter         ( Y         = 1.D9,                                  &
                        nu        = 0.3,                                   &
                        phi       = 45,                                    &
                        c_eos     = 20.0,                                  &
                        P_star    = 2.75D4,                                &
                        T_star    = 0,                                     &
                        S_star    = (P_star+T_star)/(2*2),                &
                        rho_air   = 1.3,                                   &
                        rho_water = 1026.,                                 &
                        rho_ice   = 900.,                                  &
                        c_wind    = 1.2D-3,                                &
                        c_water   = 5.5D-3,                                &
                        alpha     = 1. / ( 2 * hsml_multiplicity) ,        &
                        delta_min = 2.D-9,                                &
                        alpha_pi  = 1.,                                    &
                        beta_pi   = 2.,                                    &
                        chi       = 0.3,                                   &
                        b_length  = 8000.,                                 &
                        xi        = 0.5225,                                &
                        kappa_n   = 5.4D12,                                  &
                        kappa_t   = 5.D-8,                                 &
                        k_t       = T_star / P_star,                       &
                        e         = 0.5 * (P_star + T_star) / S_star,      &
                        mu        = TAN(phi * pi / 180.))

end module data_variables