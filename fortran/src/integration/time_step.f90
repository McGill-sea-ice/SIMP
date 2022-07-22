!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       Subroutine for computing the time step evolution with respect      !
!                          to stability criteria                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine time_step

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Call the desired eos definition
    IF ( tse .EQ. 0 ) THEN
        RETURN

    ELSEIF ( tse .EQ. 1 ) THEN
        call Monaghan_time_step

    ELSEIF ( tse .EQ. 2 ) THEN
        call VP_time_step

    ENDIF

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Subroutine for computing the time step evolution with respect     !
    !            to usual SPH stability criteria (Monaghan 2005)           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Monaghan_time_step

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    !!!              Variables use in the current routine:               !!!
    !   dt : Simulation time step                                      [out]
    !   hsml : Smoothing lengths of particles                           [in]
    !   dvxdt : dvx/dt, force per unit of mass                          [in]
    !   c   : Sound speed                                               [in]
    !   wdvx     : Various products between kernel and speed            [in]
    !   alpha_pi, beta_pi : Artificial viscosity coefficient            [in]
    !   p_i     : List of first particle partner of interaction pair    [in]
    !   p_j     : List of second particle partner of interaction pair   [in]
    !   nip     : number of interacting particle                        [in]
    !   mass   : Mass of particles                                      [in]  
    !   rho    : Densities of particles                                 [in]

    implicit none

    !!!                        Local variables:                          !!!
    !   dt_force : Timestep following the acceleration condition         [s]
    !   dt_CFL : Timestep following the CFL condition                    [s]
    !   fF : Tunning constant for the force criteria lower than 1
    !   acc : acceleration vector                                    [m/s^2]
    !   k : loop index
    !   i : particle i current index
    !   j : particle j current index
    !   div_v : Particle divergence                                   [s^-1]

    INTEGER          :: i, j, k
    DOUBLE PRECISION :: dt_force, dt_CFL, fF, acc(ntotal), div_v(ntotal)
    parameter ( fF = 0.25)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    div_v = 0.
    dt_force = 0.
    dt_CFL = 0.
    acc = 0.

    !Get the force condition
    acc = (dvxdt(1,:ntotal) ** 2 + dvxdt(2,:ntotal) ** 2 ) ** 0.5
    dt_force = fF * MINVAL( hsml(:ntotal) / acc) ** 0.5

    !Get the CFL condition
    IF (artificial_visc) THEN

        !Calculate velocity divergence
        DO k = 1, nip 

            ! Get index
            i = p_i(k)
            j = p_j(k)

            !Get relative velocity
            dvx(:)  = vx(:, i) - vx(:,j)

            !Get Kernel product with particle velocity difference
            wdvx(1) = dvx(1) * dwdx(1,k)
            wdvx(2) = dvx(2) * dwdx(1,k)
            wdvx(3) = dvx(1) * dwdx(2,k)
            wdvx(4) = dvx(2) * dwdx(2,k)

            !Velocity divergence
            div_v(i) = div_v(i) - mass(j) * (wdvx(1) + wdvx(4)) / rho(i)
            div_v(j) = div_v(j) - mass(i) * (wdvx(1) + wdvx(4)) / rho(j)

        ENDDO

        dt_CFL = MINVAL( 0.3 * hsml(:ntotal) / ( c(:ntotal) + hsml(:ntotal)&
                         * ABS(div_v) + 1.2 * ( alpha_pi * c(:ntotal) +    &
                         beta_pi * hsml(:ntotal) * ABS(div_v)) ) )

    ELSE
        dt_CFL = MINVAL( hsml(:ntotal) / c(:ntotal)) 

    ENDIF

    !Get new time step length
    dt = MIN( dt_force, dt_CFL)

    end subroutine Monaghan_time_step

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Subroutine for computing the time step evolution with respect     !
    !                        to VP stability criteria                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine VP_time_step

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    !   dt : Simulation time step                                      [out]
    !   hsml : Smoothing lengths of particles                           [in]
    !   P_star    : Empirical constant for ice strength in compresison  [in]
    !   rho_ice   : Density of ice                                      [in]
    !   zeta_max  : Maximum bulk viscosity                              [in]
    !   h      : Mean ice thickness of particles                        [in]

    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Get new time step length
    dt = MINVAL( rho_ice * hsml(:ntotal) ** 2 * delta_min / (P_star *(1+k_t)))

    end subroutine VP_time_step

end subroutine time_step