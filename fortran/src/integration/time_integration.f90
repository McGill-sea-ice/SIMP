!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for fowarding in time the solution              !
!                          and particle information.                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine time_integration

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    !   integration_time : Desired time for the simulation              [in]
    !   dt               : time step used in the time integration       [in]
    !   itimestep        : current iteration of the simulation         [out]
    !   current_time : Current time in the integration                 [out]

    implicit none    

    !!!                        Local variables:                          !!!
    !   last_save : Last integration time the data was save to file

    REAL(kind=8)        :: last_save

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    last_save = - save_step
    itimestep = 0

    !Print message of the time integration start
    WRITE(*,*) '**********************************************************'
    WRITE(*,*) '                    Starting integration'

    !Loop for the time integration
    DO WHILE (current_time .LT. integration_time)

        !Print information if requested
        IF (mod(itimestep,print_step) .EQ. 0) THEN
            WRITE(*,'(A27,I10)') ' Current iteration number =', itimestep
            WRITE(*,'(A19,f10.4)') ' Current time [h] =', current_time/3600
            WRITE(*,'(A31,f14.6)') ' Current time step length [s] =', dt
        ENDIF


        IF (int_algorithm .EQ. 1) THEN
            call euler_scheme

        ELSEIF (int_algorithm .EQ. 2) THEN
            call predictor_corrector

        ELSEIF (int_algorithm .EQ. 3) THEN
            call leapfrog_scheme

        ENDIF


        !Compile the data if desired
        IF ( current_time .GT. save_step + last_save ) THEN
            call output
            last_save = current_time

        ENDIF

        !Update itimestep
        itimestep = itimestep + 1

        !Update current time
        current_time = current_time + dt

    ENDDO

    !Print message of the time integration end
    WRITE(*,*) '                     Ending integration'
    WRITE(*,*) '**********************************************************'

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme for the time integration with the            !
    !                               Euler scheme                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine euler_scheme

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   x                : coordinates of particles                 [in/out]
    !   vx               : velocities of particles                  [in/out]
    !   A                : Sea ice concentration of particles       [in/out]
    !   h                : Mean ice thickness of particles          [in/out]
    !   rho              : Densities of particles                   [in/out]
    !   ntotal           : total particle number                        [in]
    !   dt               : time step used in the time integration       [in]
    !   dvxdt            : dvx/dt, force per unit of mass               [in]
    !   dadt             : Concetration rate of change                  [in]
    !   dhdt             : Thickness rate of change                     [in]
    !   drhodt           : Density rate of change                       [in]

    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Get single time step integration quantities
    call single_step

    !Do Euler technique for time step integration
    x(:,:ntotal)  = x(:,:ntotal) + dxdt(:,:ntotal) * dt
    vx(:,:ntotal) = vx(:,:ntotal) + dvxdt(:,:ntotal) * dt
    h(:ntotal)    = h(:ntotal) + dhdt(:ntotal) * dt
    A(:ntotal)    = A(:ntotal) + dAdt(:ntotal) * dt

    !Correct concentration
    WHERE (A(:ntotal) .GT. 1)
        A = 1
    ENDWHERE

    !Compute density evolution 
    call density

    !Compute the change in density only if necessary
    IF (den_alg .EQ. 2) THEN
        rho(:ntotal)    = rho(:ntotal) + drhodt(:ntotal) * dt

    ENDIF

    !Update masses
    call mass_evolution

    !Update sound speed
    call sound_speed

    !Update hsml 
    call smoothing_length

    !Update time step
    call time_step

    end subroutine euler_scheme

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme for the time integration with a              !
    !              predictro-corrector scheme (Hosseini,2019)              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine predictor_corrector

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   x                : coordinates of particles                 [in/out]
    !   vx               : velocities of particles                  [in/out]
    !   h                : Mean ice thickness of particles          [in/out]
    !   rho              : Densities of particles                   [in/out]
    !   ntotal           : total particle number                        [in]
    !   dt               : time step used in the time integration       [in]
    !   dvxdt            : dvx/dt, force per unit of mass               [in]
    !   dhdt             : Thickness rate of change                     [in]
    !   drhodt           : Density rate of change                       [in]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   x_12 : coordinates of particles at half time step                [m]
    !   vx_12 : velocities of particles at half time step              [m/s]
    !   h_12 : Mean ice thickness of particles at half time step         [m]
    !   A_12: Mean ice concentration  of particle at half a time step
    !   rho_12: Density of particle at half a time step             [kg/m^2]
    !   x_1 : coordinates of particles at current time step              [m]
    !   vx_1 : velocities of particles at  current time step           [m/s]
    !   h_1 : Mean ice thickness of particles at  current time step      [m]
    !   A_1: Mean ice concentration of particle at current time step  
    !   rho_1: Density of particle at current time step             [kg/m^2]

    DOUBLE PRECISION :: x_12(dim,ntotal), vx_12(dim,ntotal), rho_1(ntotal),&
                        h_12(ntotal), rho_12(ntotal), x_1(dim,ntotal),     &
                        vx_1(dim,ntotal), h_1(ntotal), A_1(ntotal),        &
                        A_12(ntotal)
                        

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Monitor current quantities
    x_1   = x(:,:ntotal)
    vx_1  = vx(:,:ntotal) 
    h_1   = h(:ntotal) 
    A_1   = A(:ntotal) 
    rho_1 = rho(:ntotal)

    !Get single time step integration at current time
    call single_step

    !Predicts the evolution of parameters at half time step
    x(:,:ntotal)  = x(:,:ntotal) + dxdt(:,:ntotal) * dt/ 2
    vx(:,:ntotal) = vx(:,:ntotal) + dvxdt(:,:ntotal) * dt/ 2
    h(:ntotal)    = h(:ntotal) + dhdt(:ntotal) * dt/ 2
    A(:ntotal)    = A(:ntotal) + dAdt(:ntotal) * dt/ 2

    !Correct concentration
    WHERE (A(:ntotal) .GT. 1)
        A = 1
    ENDWHERE

    !Compute density evolution 
    call density

    !Compute the change in density only if necessary
    IF (den_alg .EQ. 2) THEN
        rho(:ntotal)    = rho(:ntotal) + drhodt(:ntotal) * dt/ 2

    ENDIF
    
    !Update masses
    call mass_evolution

    !Update sound speed
    call sound_speed

    !Get single time step at half the time step later
    call single_step

    !Get corrected half time step quantities
    x_12   = x_1 + dxdt(:,:ntotal) * dt/ 2
    vx_12  = vx_1 + dvxdt(:,:ntotal) * dt/ 2
    h_12   = h_1 + dhdt(:ntotal) * dt/ 2
    A_12   = A_1 + dAdt(:ntotal) * dt/ 2

    !Compute density evolution 
    call density

    !Compute the change in density only if necessary
    IF (den_alg .EQ. 2) THEN
        rho_12   = rho_1 + drhodt(:ntotal) * dt/ 2

    ENDIF

    !Get final quantities at the next time step
    x(:,:ntotal)  = 2 * x_12 - x_1
    vx(:,:ntotal) = 2 * vx_12 - vx_1
    h(:ntotal)    = 2 * h_12 - h_1
    A(:ntotal)    = 2 * A_12 - A_1

    !Correct concentration
    WHERE (A(:ntotal) .GT. 1)
        A = 1
    ENDWHERE

    !Compute density evolution 
    call density

    !Compute the change in density only if necessary
    IF (den_alg .EQ. 2) THEN
        rho(:ntotal) = 2 * rho_12 - rho_1

    ENDIF

    !Update masses
    call mass_evolution

    !Update sound speed
    call sound_speed

    !Update hsml 
    call smoothing_length
    
    !Update time step
    call time_step

    end subroutine predictor_corrector

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme for the time integration with the            !
    !                             Leap-Frog scheme                         !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine leapfrog_scheme

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   x                : coordinates of particles                 [in/out]
    !   vx               : velocities of particles                  [in/out]
    !   h                : Mean ice thickness of particles          [in/out]
    !   rho              : Densities of particles                   [in/out]
    !   ntotal           : total particle number                        [in]
    !   dt               : time step used in the time integration       [in]
    !   dvxdt            : dvx/dt, force per unit of mass               [in]
    !   dhdt             : Thickness rate of change                     [in]
    !   drhodt           : Density rate of change                       [in]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   x_1 : coordinates of particles at current time step              [m]
    !   vx_1 : velocities of particles at  current time step           [m/s]
    !   h_1 : Mean ice thickness of particles at  current time step      [m]
    !   rho_1: Density of particle at current time step             [kg/m^2]
    !   A_1: Mean ice concentration of particle at current time step
    
    DOUBLE PRECISION :: x_1(dim,ntotal), vx_1(dim,ntotal), h_1(ntotal),    &
                        rho_1(ntotal), A_1(ntotal)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Monitor current quantities
    x_1   = x(:,:ntotal)
    vx_1  = vx(:,:ntotal)
    h_1   = h(:ntotal) 
    rho_1 = rho(:ntotal) 
    A_1   = A(:ntotal)

    !Get single time step at current time 
    call single_step

    !Do Euler technique for half time step integration
    x(:,:ntotal)  = x(:,:ntotal) + dxdt(:,:ntotal) * dt/ 2
    vx(:,:ntotal) = vx(:,:ntotal) + dvxdt(:,:ntotal) * dt/ 2
    h(:ntotal)    = h(:ntotal) + dhdt(:ntotal) * dt / 2
    A(:ntotal)    = A(:ntotal) + dAdt(:ntotal) * dt / 2

    !Correct concentration
    WHERE (A(:ntotal) .GT. 1)
        A = 1
    ENDWHERE

    !Compute density evolution 
    call density

    !Compute the change in density only if necessary
    IF (den_alg .EQ. 2) THEN
        rho(:ntotal)  = rho(:ntotal) + drhodt(:ntotal) * dt / 2
    ENDIF

    !Update masses
    call mass_evolution

    !Update sound speed
    call sound_speed

    !Get single time step at half the time step later
    call single_step

    !Get final quantities at the next time step
    vx(:,:ntotal) = vx_1 + dvxdt(:,:ntotal) * dt 
    h(:ntotal)    = h_1 + dhdt(:ntotal) * dt 
    x(:,:ntotal)  = x_1 + (vx_1 + dxdt(:,:ntotal)) * dt / 2
    A(:ntotal)    = A_1 + dAdt(:ntotal) * dt 

    !Correct concentration
    WHERE (A(:ntotal) .GT. 1)
        A = 1
    ENDWHERE

    !Compute density evolution 
    call density

    !Compute the change in density only if necessary
    IF (den_alg .EQ. 2) THEN
        rho(:ntotal)  = rho_1 + drhodt(:ntotal) * dt 
    ENDIF

    !Update masses
    call mass_evolution

    !Update sound speed
    call sound_speed

    !Update hsml 
    call smoothing_length
    
    !Update time step
    call time_step

    end subroutine leapfrog_scheme

end subroutine time_integration