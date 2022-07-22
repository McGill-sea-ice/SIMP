!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for loading or generating initial               !
!                            particle information.                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine initialize

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    implicit none

    !!!                        Local variables:                          !!!

    !!!                Call modules necessary properties:                !!!
    CALL sys_variables

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialize data variables to zero to avoid bug
    ntotal       = 0
    x            = 0.
    vx           = 0.
    mass         = 0.
    rho          = 0.
    p            = 0.
    hsml         = 0.
    A            = 0.
    h            = 0.
    c            = 0.
    ntotalb      = 0
    b_normal     = 0.
    xb           = 0.
    massb        = 0.
    itimestep    = 0
    dvxdt        = 0.
    dxdt         = 0.
    dAdt         = 0.
    dhdt         = 0.
    drhodt       = 0.
    dt           = 0.
    current_time = 0.
    nip          = 0
    p_i          = 0
    p_j          = 0
    bnip         = 0
    bp_i         = 0
    bp_j         = 0
    r_ij         = 0.
    br_ij        = 0.
    hsml_ij      = 0.
    max_hsml     = 0.
    scale_h      = 0
    w            = 0.
    dwdx         = 0.
    kernel_zero  = 0.
    epsilon_11   = 0.
    epsilon_12   = 0.
    epsilon_22   = 0.
    zeta         = 0.
    eta          = 0.
    sigma_11     = 0.
    sigma_12     = 0.
    sigma_22     = 0.

    !Start run from scratch
    IF (.NOT. spin_up) THEN

        !Get initial properties of the particles
        CALL particle_placing

    !Start run from previous integration
    ELSEIF (spin_up) THEN
        CALL particle_reading

    ENDIF
    
    !!!                File creation for compilling data:                !!!
    open(1,file= cwd // trim(output_folder) // "particle_momentum.csv",    &
                                action='write',                            &
                                status='replace')
    open(2,file= cwd // trim(output_folder) // "particle_state.csv",       &
                                action='write',                            &
                                status='replace')
    open(3,file= cwd // trim(output_folder) // "particle_other.csv",       &
                                action='write',                            &
                                status='replace')

    CLOSE(1)
    CLOSE(2)
    CLOSE(3)
    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !              Subroutine for generating initial particles              !
    !                       position from a given domain                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine particle_placing

    !!!                         Module to include:                       !!!
    use globals
    use domain
    use data_variables

    !!!              Variables use in the current routine:               !!!
    !   x       : Coordinates of all particles                         [out]
    !   ntotal  : Total number of particle to use                      [out]
    !   hsml    : Smoothing lengths of particles                       [out]
    !   scale_h : Support domain scaling factor                        [out]
    !   vx     : Velocities of all particles                           [out]
    !   mass   : Mass of particles                                     [out]  
    !   rho    : Densities of particles                                [out]
    !   p      : Pressure of particules                                [out]
    !   A      : Sea ice concentration of particles                    [out]
    !   h      : Mean ice thickness of particles                       [out]
    !   c      : Sound speed                                           [out]
    !   dt        : Time step of the integration                       [out]

    implicit none

    !!!                        Local variables:                          !!!
    !   nppg : number of particle per grid cell
    !   res_ratio : resolution ratio between particles and grid
    !   index : particle index
    !   m, n, k, l : indices for do loops
    INTEGER :: nppg, res_ratio, m , n , k, l, index

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    index = 0
    nppg = 0
    res_ratio = 0

    !Get the number of particle per grid cell
    res_ratio = NINT (map_resolution / resolution)
    nppg = res_ratio ** 2
    ntotal = COUNT (map .EQ. 1) * nppg

    !Verify if resolution is allowed
    IF (nppg .EQ. 0) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) '    Resolution desired too coarse by comparison to    '
        WRITE (*,*) '              the domain map resolution.              '
        WRITE (*,*) '******************************************************'

        STOP
    ENDIF
        
    IF (ntotal .GT. maxn) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) ' Resolution desired too high, the number of particles '
        WRITE (*,*) '    generated is higher than the maximum allowed.     '
        WRITE (*,*) '******************************************************'

        STOP
    ENDIF

    !Generate particles for the given domain
    DO n = 0, y_numcell - 1
        DO m = 0, x_numcell - 1

            IF( map(m+1,n+1) .EQ. 1 ) THEN

                !Generate the particle in a square cell subdivided according
                !to the resolution:
                DO k = 0, res_ratio - 1     
                    DO l = 0, res_ratio - 1

                        !Compile the coordinates of the particles:
                        index = index + 1
                        x(1,index) = x_mingeom + (-1./2. + m + (k + 1./2.) &
                                    / res_ratio) * map_resolution
                        x(2,index) = y_maxgeom + ( 1./2. - n - (l + 1./2.) &
                                    / res_ratio) * map_resolution
                    ENDDO
                ENDDO
            ENDIF

        ENDDO
    ENDDO

    !Set initial quantities
    hsml(:ntotal)  = hsml_multiplicity * map_resolution / res_ratio
    max_hsml       = MAXVAL(hsml(:ntotal)) * 10.
    A(:ntotal)     = 1
    h(:ntotal)     = 1
    vx(1,:ntotal)  = 0.
    vx(2,:ntotal)  = 0.
    rho(:ntotal)   = rho_ice * h(:ntotal)
    mass(:ntotal)  = rho(:ntotal) * (hsml(:ntotal) / hsml_multiplicity) ** 2

    call time_step
    call sound_speed

    !Set boundaries only if necessary
    IF (boundary_pressure .OR. boundary_friction) THEN
        CALL boundaries
    ENDIF

    !Defining scale factor according to kernel
    IF ((skf .EQ. 4) .OR. (skf .EQ. 5) .OR. (skf .EQ. 6)) THEN
        scale_h = 1
    
    ELSEIF ( (skf .EQ. 1)) THEN
        scale_h = 2

    ELSEIF ((skf .EQ. 3) .OR. (skf .EQ. 2) .OR. (skf .EQ. 7) .OR.          &
        (skf .EQ. 8)) THEN
        scale_h = 3
    ENDIF

    !Print the resulting particles simulation parameter
    WRITE (*,*) '**********************************************************'
    WRITE (*,*) '          ', ntotal, ' particles generated.                 '
    WRITE (*,'(12x,A35,f8.2)') 'Actual resolution of particles is: ',       &
                                    map_resolution/res_ratio/1000
    WRITE (*,'(15x,A31,f10.3)') 'Corresponding smoothing length:',          &
                                MINVAL(hsml(:ntotal))/1000
    WRITE (*,*) '**********************************************************'

    end subroutine particle_placing

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !              Subroutine for generating initial particles              !
    !                   infromations from a given spin up                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine particle_reading

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    !   x       : Coordinates of all particles                         [out]
    !   ntotal  : Total number of particle to use                      [out]
    !   hsml    : Smoothing lengths of particles                       [out]
    !   scale_h : Support domain scaling factor                        [out]
    !   vx     : Velocities of all particles                           [out]
    !   mass   : Mass of particles                                     [out]  
    !   rho    : Densities of particles                                [out]
    !   p      : Pressure of particules                                [out]
    !   A      : Sea ice concentration of particles                    [out]
    !   h      : Mean ice thickness of particles                       [out]
    !   c      : Sound speed                                           [out]
    !   dt        : Time step of the integration                       [out]

    implicit none

    !!!                        Local variables:                          !!!
    !   i : particle index
    !   n, d : loop index

    INTEGER :: n, i

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Get the path of the folder containing spin up information
    WRITE(*,*) "Enter the path to the folder containing input:"
    read (*,*) path

    !Open corresponding file
    open(1,file= trim(path) // "particle_momentum.csv")
    open(2,file= trim(path) // "particle_state.csv")
    open(3,file= trim(path) // "particle_other.csv")

    !Get global informations
    read(3,*) ntotal, ntotalb, current_time

    !Verify if enough memory is allocated for those input files
    IF (ntotal .GT. maxn) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) '  The number of particles read from the input file    '
        WRITE (*,*) '        is higher than the maximum allowed.           '
        WRITE (*,*) '******************************************************'

        STOP
    ENDIF

    !Convert current_time in secondes
    current_time = current_time * 3600

    !Update integration time to fit with the current time
    integration_time = integration_time + current_time

    WRITE (*,*) '**********************************************************'
    WRITE (*,*) '           Loading initial particle information           '
    WRITE (*,'(15x,A29,I5)') 'Number of particle detected:',ntotal
    WRITE (*,'(7x,A37,I5)') 'Number of boundary particle detected:',ntotalb
    WRITE (*,'(21x,A23,e15.8)') 'Current simulation time', current_time
    WRITE (*,*) '**********************************************************'

    !Reading the particle informations
    do n = 1, ntotal
        read(1,*) i, x(1, i), x(2, i), vx(1, i), vx(2,i)

        read(2,*) i, mass(i), rho(i), p(i), hsml(i), A(i), h(i),       &
                     epsilon_11(i), epsilon_12(i), epsilon_22(i), zeta(i),&
                     eta(i), sigma_11(i), sigma_12(i), sigma_22(i)
    enddo

    !Convert position to meter and speed to m/s
    x(:,:ntotal) = x(:,:ntotal) * 1000
    vx(:, :ntotal) = vx(:, :ntotal) / 100

    !Defining scale factor according to kernel
    IF (skf .EQ. 1) then
        scale_h = 2

    ELSEIF (skf .EQ. 2) then
        scale_h = 3

    ELSEIF (skf .EQ. 3) then
        scale_h = 3

    ENDIF

    CLOSE(1)
    CLOSE(2)
    CLOSE(3)

    end subroutine particle_reading

end subroutine initialize

