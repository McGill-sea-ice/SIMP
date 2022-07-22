!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for loading or generating boundary              !
!                             particles information                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine boundaries

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    !   massb   : Mass of boudary particles                            [out]
    !   ntotalb : Total number of boundary particle to use              [in]
    !   mass   : Mass of particles                                      [in]

    implicit none

    !!!                        Local variables:                          !!!
    !   i : Particle index
    !   d : Dimension
    INTEGER :: i, d

    !!!                Call modules necessary properties:                !!!
    CALL sys_variables

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!                       Formating the output:                      !!!
    1001 format (1x, I5, 3(1x, A1, 1x, e14.8))

    !Get initial properties of the boundary particles
    IF (b_setup .EQ. 0) THEN
        CALL boundary_reading

        !Compute boundary unit normal vector if necessary
        IF (bfa .EQ. 2) THEN
            call boundary_normal
        
        ENDIF

    ELSEIF (b_setup .EQ. 1) THEN
        CALL boundary_placing_1

        !Compute boundary unit normal vector if necessary
        IF (bfa .EQ. 2) THEN
            call boundary_normal
        
        ENDIF
    
        massb(:ntotalb) = SUM(mass(:ntotal))/ntotal

        !!!                File opening for compilling data:                 !!!  
        open(1,file= cwd // trim(output_folder) // "boundary_particles.csv",   &
                                   action='write',                             &
                                   status='replace')
    
        !!!                        Writing the output:                       !!!
        do i = 1, ntotalb
            write(1,1001) i, (',', xb(d, i)/1000 , d = 1, dim),',', massb(i)
        enddo
    
        CLOSE(1)

    ELSEIF (b_setup .EQ. 2) THEN
        CALL boundary_placing_2

        !Compute boundary unit normal vector if necessary
        IF (bfa .EQ. 2) THEN
            WRITE (*,*) '******************************************************'
            WRITE (*,*) '     The boundary force cannot be with bfa 2 if       '
            WRITE (*,*) '  the boundary layer has multiple set of particles.   '
            WRITE (*,*) '******************************************************'
    
            STOP
        
        ENDIF

        massb(:ntotalb) = SUM(mass(:ntotal))/ntotal

        !!!                File opening for compilling data:                 !!!  
        open(1,file= cwd // trim(output_folder) // "boundary_particles.csv",   &
                                   action='write',                             &
                                   status='replace')
    
        !!!                        Writing the output:                       !!!
        do i = 1, ntotalb
            write(1,1001) i, (',', xb(d, i)/1000, d = 1, dim), ',', massb(i)
        enddo
    
        CLOSE(1)

    ENDIF

    contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !              Subroutine for generating initial particles              !
    !                   infromations from a given spin up                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine boundary_reading

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
    !   j : particle index
    !   n : loop index
    INTEGER :: n, j

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Open file with boundary particle information
    open(1,file= trim(path) // "boundary_particles.csv")

    !Verify if enough memory is allocated for those input files
    IF (ntotalb .GT. maxnb) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) '     The number of boundary particles read from       '
        WRITE (*,*) '  the input file is higher than the maximum allowed.  '
        WRITE (*,*) '******************************************************'

        STOP
    ENDIF

    !Reading the particle informations
    do n = 1, ntotalb
        read(1,*) j, xb(1, j), xb(2, j), massb(j)

    enddo

    !Convert position to meter 
    xb(:,:ntotal) = xb(:,:ntotal) * 1000

    CLOSE(1)

    end subroutine boundary_reading

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subroutine for generating a single layer of              !
    !           boundary particles position from a given domain            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine boundary_placing_1

    !!!                         Module to include:                       !!!
    use globals
    use boundary_domain
    use data_variables
    
    !!!              Variables use in the current routine:               !!!
    !   xb      : Coordinates of all boundary particles                [out]
    !   ntotalb : Total number of boundary particle to use             [out]
    !   b_length    : Length of interaction of the boundary             [in]
    !   hsml    : smoothing length of the particles                     [in]
    !   bx_maxgeom : Boundary upper limit of allowed x-regime           [in]
    !   bx_mingeom : Boundary lower limit of allowed x-regime           [in]
    !   by_maxgeom : Boundary upper limit of allowed y-regime           [in]
    !   by_mingeom : Boundary lower limit of allowed y-regime           [in]
    !   bmap_resolution : Boundary side length of the map square grid   [in]
    !   bx_numcell : Boundary number of cell in the x direction         [in]
    !   by_numcell : Boundary number of cell in the y direction         [in]
    !   map_factor : Resolution ratio between domain maps (bmap / map)  [in]
    !   bmap : True/False domain map

    implicit none
    
    !!!                        Local variables:                          !!!
    !   bppg : number of boundary particle per grid cell edge
    !   index : particle index
    !   n, m, k : indices for do loops

    INTEGER :: n , m , k, index ,bppg

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    bppg = 0
    index = 0

    !Get number of boundary particle
    bppg = NINT(bmap_resolution/(b_length / 3))
    
    IF (bppg .EQ. 0) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) '   Boundary map resolution too high by comparison to  '
        WRITE (*,*) '          the boundary particle resolution.           '
        WRITE (*,*) '******************************************************'
    
        STOP
    ENDIF

    !Generate boundary particles
    !Pass through all cell:
    DO CONCURRENT (n = 1: bx_numcell, m = 1 : by_numcell ,                  &
                    bmap(n,m) .EQ. 0) 
    
        !Place particles on the left cell boundary:
        IF ((bmap(n+1,m) .EQ. 1) .AND. ((n+1) .LE. bx_numcell)) THEN 
            DO k = 0, bppg - 1

                !Create layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1./2. * map_factor)        &
                            * bmap_resolution 
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor) &
                            - (k + 1./2.)/ bppg) * bmap_resolution

                !Get normal orientation
                b_normal(1,index) = 1.
                b_normal(2,index) = 0.

            ENDDO
        ENDIF
    
        !Place particles on the right cell boundary:
        IF ((bmap(n-1,m) .EQ. 1) .AND. ((n-1) .GE. 1)) THEN 
            DO k = 0, bppg - 1

                !Create layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor))  &
                            * bmap_resolution 
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor) &
                            - (k + 1./2.)/ bppg) * bmap_resolution

                !Get normal orientation
                b_normal(1,index) = -1.
                b_normal(2,index) =  0.

            ENDDO
        ENDIF
    
        !Place particles on the upper cell boundary:
        IF ((bmap(n,m+1) .EQ. 1) .AND. ((m+1) .LE. by_numcell)) THEN 
            DO k = 0, bppg - 1

                !Create 3 layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor) +  &
                            (k + 1./2.) / bppg) * bmap_resolution
                xb(2,index) = by_maxgeom + ( -m + (1./2. * map_factor) )    &
                            * bmap_resolution 

                !Get normal orientation
                b_normal(1,index) = 0.
                b_normal(2,index) = -1.

            ENDDO
        ENDIF
    
        !Place particles on the lower cell boundary:
        IF ((bmap(n,m-1) .EQ. 1) .AND. ((m-1) .GE. 1)) THEN 
            DO k = 0, bppg - 1

                !Create 3 layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor) +  &
                            (k + 1./2.) / bppg) * bmap_resolution
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor)) &
                            * bmap_resolution 

                !Get normal orientation
                b_normal(1,index) = 0.
                b_normal(2,index) = 1.

            ENDDO
        ENDIF

    ENDDO
    
    !Test if enough memory is allocated
    ntotalb = index
    IF (ntotalb .GT. maxnb) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) '      The number of boundary particles generated      '
        WRITE (*,*) '          is higher than the maximum allowed.         '
        WRITE (*,*) '******************************************************'

        STOP
    ENDIF
    
    !Print the resulting boundary simulation parameter
    WRITE (*,*) '**********************************************************'
    WRITE (*,*) '       ', ntotalb, ' boundary particles generated.        '
    WRITE (*,'(12x,A36,f6.2)') 'Actual resolution of boundaries is: ',  &
                                bmap_resolution/bppg / 1000
    WRITE (*,*) '**********************************************************'
    
    end subroutine boundary_placing_1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subroutine for generating multiple layer of              !
    !           boundary particles position from a given domain            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine boundary_placing_2

    !!!                         Module to include:                       !!!
    use globals
    use boundary_domain
    use data_variables
    
    !!!              Variables use in the current routine:               !!!
    !   xb      : Coordinates of all boundary particles                [out]
    !   ntotalb : Total number of boundary particle to use             [out]
    !   b_length    : Boundary particle length of interaction
    !   hsml    : smoothing length of the particles                     [in]
    !   bx_maxgeom : Boundary upper limit of allowed x-regime           [in]
    !   bx_mingeom : Boundary lower limit of allowed x-regime           [in]
    !   by_maxgeom : Boundary upper limit of allowed y-regime           [in]
    !   by_mingeom : Boundary lower limit of allowed y-regime           [in]
    !   bmap_resolution : Boundary side length of the map square grid   [in]
    !   bx_numcell : Boundary number of cell in the x direction         [in]
    !   by_numcell : Boundary number of cell in the y direction         [in]
    !   map_factor : Resolution ratio between domain maps (bmap / map)  [in]
    !   bmap : True/False domain map

    implicit none
    
    !!!                        Local variables:                          !!!
    !   bppg : number of boundary particle per grid cell edge
    !   index : particle index
    !   n, m, k : indices for do loops
    !   spacing : Spacing between the 3 layer of boundary particle

    INTEGER :: n , m , k, l, index ,bppg
    DOUBLE PRECISION :: spacing
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    bppg = 0
    index = 0

    !Get number of boundary particle
    bppg = NINT(bmap_resolution/(b_length / 3))
    
    IF (bppg .EQ. 0) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) '   Boundary map resolution too high by comparison to  '
        WRITE (*,*) '          the boundary particle resolution.           '
        WRITE (*,*) '******************************************************'
    
        STOP
    ENDIF
    
    !Get biggest smoothing length
    spacing = MAXVAL(hsml(:ntotal)) / hsml_multiplicity

    !Generate boundary particles
    !Pass through all cell:
    DO CONCURRENT (n = 1: bx_numcell, m = 1 : by_numcell ,                  &
                   bmap(n,m) .EQ. 0) 
    
        !Place particles on the left cell boundary:
        IF ((bmap(n+1,m) .EQ. 1) .AND. ((n+1) .LE. bx_numcell)) THEN 
            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

                !Create layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1./2. * map_factor)        &
                            * bmap_resolution - l * spacing
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor) &
                            - (k + 1./2.)/ bppg) * bmap_resolution
                ENDDO
            ENDDO
        ENDIF
    
        !Place particles on the right cell boundary:
        IF ((bmap(n-1,m) .EQ. 1) .AND. ((n-1) .GE. 1)) THEN 
            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

                !Create layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor))  &
                            * bmap_resolution + l * spacing
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor) &
                            - (k + 1./2.)/ bppg) * bmap_resolution

                ENDDO
            ENDDO
        ENDIF
    
        !Place particles on the upper cell boundary:
        IF ((bmap(n,m+1) .EQ. 1) .AND. ((m+1) .LE. by_numcell)) THEN 
            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

                !Create 3 layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor) +  &
                            (k + 1./2.) / bppg) * bmap_resolution
                xb(2,index) = by_maxgeom + ( -m + (1./2. * map_factor) )    &
                            * bmap_resolution + l * spacing

                ENDDO
            ENDDO
        ENDIF
    
        !Place particles on the lower cell boundary:
        IF ((bmap(n,m-1) .EQ. 1) .AND. ((m-1) .GE. 1)) THEN 
            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

                !Create 3 layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor) +  &
                            (k + 1./2.) / bppg) * bmap_resolution
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor)) &
                            * bmap_resolution - l * spacing

                ENDDO
            ENDDO
        ENDIF
    
        !Place particles for the lower left cell corner boundary:

        IF ((bmap(n+1,m-1) .EQ. 1) .AND. (bmap(n,m-1) .NE. 1) .AND.         &
            (bmap(n,m+1) .NE. 1) .AND. (bmap(n-1,m) .NE. 1) .AND.           &
            (bmap(n+1,m) .NE. 1) .AND. (((n+1) .LE. bx_numcell) .AND.       &
            ((m-1) .GE. 1))) THEN 

            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

                !Create layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1./2. * map_factor)        &
                            * bmap_resolution - l * spacing
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor) &
                            - (k + 1./2.)/ bppg) * bmap_resolution

                ENDDO
            ENDDO
            
            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

            !Create 3 layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor) +  &
                            (k + 1./2.) / bppg) * bmap_resolution
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor)) &
                            * bmap_resolution - l * spacing

                ENDDO
            ENDDO
        ENDIF

        !Place particles for the upper right cell corner boundary:

        IF ((bmap(n-1,m+1) .EQ. 1) .AND. (bmap(n,m-1) .NE. 1) .AND.         &
            (bmap(n,m+1) .NE. 1) .AND. (bmap(n-1,m) .NE. 1) .AND.           &
            (bmap(n+1,m) .NE. 1) .AND. (((m+1) .LE. by_numcell) .AND.       &
            ((n-1) .GE. 1))) THEN 

            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1
    
                !Create layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor))  &
                            * bmap_resolution + l * spacing
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor) &
                            - (k + 1./2.)/ bppg) * bmap_resolution

                ENDDO
            ENDDO
            
            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

                !Create 3 layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor) +  &
                            (k + 1./2.) / bppg) * bmap_resolution
                xb(2,index) = by_maxgeom + ( -m + (1./2. * map_factor) )    &
                            * bmap_resolution + l * spacing

                ENDDO
            ENDDO
        ENDIF

        !Place particles for the upper left cell corner boundary:

        IF ((bmap(n+1,m+1) .EQ. 1) .AND. (bmap(n,m-1) .NE. 1) .AND.         &
            (bmap(n,m+1) .NE. 1) .AND. (bmap(n-1,m) .NE. 1) .AND.           &
            (bmap(n+1,m) .NE. 1) .AND. (((n+1) .LE. bx_numcell) .AND.       &
            ((m+1) .LE. by_numcell))) THEN 

            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

                !Create layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1./2. * map_factor)        &
                            * bmap_resolution - l * spacing
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor) &
                            - (k + 1./2.)/ bppg) * bmap_resolution

                ENDDO
            ENDDO
            
            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

                !Create 3 layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor) +  &
                            (k + 1./2.) / bppg) * bmap_resolution
                xb(2,index) = by_maxgeom + ( -m + (1./2. * map_factor) )    &
                            * bmap_resolution + l * spacing

                ENDDO
            ENDDO
        ENDIF

        !Place particles for the lower right cell corner boundary:

        IF ((bmap(n-1,m-1) .EQ. 1) .AND. (bmap(n,m-1) .NE. 1) .AND.         &
            (bmap(n,m+1) .NE. 1) .AND. (bmap(n-1,m) .NE. 1) .AND.           &
            (bmap(n+1,m) .NE. 1) .AND. (((m-1) .GE. 1) .AND.                 &
            ((n-1) .GE. 1))) THEN 

            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

                !Create 3 layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor) +  &
                            (k + 1./2.) / bppg) * bmap_resolution
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor)) &
                            * bmap_resolution - l * spacing

                ENDDO
            ENDDO
            
            DO k = 0, bppg - 1
                DO l = 0, hsml_multiplicity - 1

                !Create layer of boundary according to hsml_multiplicity
                index = index + 1
                xb(1,index) = bx_mingeom + (n - 1 - (1./2. * map_factor))  &
                            * bmap_resolution + l * spacing
                xb(2,index) = by_maxgeom + ( -m + 1 + (1./2. * map_factor) &
                            - (k + 1./2.)/ bppg) * bmap_resolution

                ENDDO
            ENDDO
        ENDIF

    ENDDO
    
    !Test if enough memory is allocated
    ntotalb = index
    IF (ntotalb .GT. maxnb) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) '      The number of boundary particles generated      '
        WRITE (*,*) '          is higher than the maximum allowed.         '
        WRITE (*,*) '******************************************************'

        STOP
    ENDIF

    !Print the resulting boundary simulation parameter
    WRITE (*,*) '**********************************************************'
    WRITE (*,*) '       ', ntotalb, ' boundary particles generated.        '
    WRITE (*,'(12x,A36,f6.2)') 'Actual resolution of boundaries is: ',  &
                                bmap_resolution/bppg / 1000
    WRITE (*,*) '**********************************************************'
    
    end subroutine boundary_placing_2


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subroutine for computing the normal unit vector          !
    !               of the boundary at each boundary particle              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine boundary_normal

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    
    !!!              Variables use in the current routine:               !!!
    !   xb      : Coordinates of all boundary particles                [out]
    !   ntotalb : Total number of boundary particle to use             [out]
    !   b_length    : Length of interaction of the boundary             [in]
    !   hsml    : smoothing length of the particles                     [in]
    !   bx_maxgeom : Boundary upper limit of allowed x-regime           [in]
    !   bx_mingeom : Boundary lower limit of allowed x-regime           [in]
    !   by_maxgeom : Boundary upper limit of allowed y-regime           [in]
    !   by_mingeom : Boundary lower limit of allowed y-regime           [in]
    !   bmap_resolution : Boundary side length of the map square grid   [in]
    !   bx_numcell : Boundary number of cell in the x direction         [in]
    !   by_numcell : Boundary number of cell in the y direction         [in]
    !   map_factor : Resolution ratio between domain maps (bmap / map)  [in]
    !   bmap : True/False domain map

    implicit none
    
    !!!                        Local variables:                          !!!
    !   NN    : The two nearest neighbour
    !   matrix: Use for the slope calculation
    !   matrix_inv: Inverse of matrix
    !   i, j : Loop index
    !   d : distance between current particle
    !   d_n : Keep distance of current neighbour particle to compare
    !   tang : Tangente slope and origin values
    !   vector : Vector to multiply the inverse with
    !   n_vector : Normal vector
    !   orientation: Vector to monitor the orientation

    INTEGER          :: NN(2,ntotalb), i, j, orientation(dim,ntotalb),     &
                        mask(3)
    DOUBLE PRECISION :: matrix(2,2), matrix_inv(2,2), d, d_n(2,ntotalb),   &
                        tang(2), vector(2), xb_temp(dim,ntotalb),          &
                        b_n_temp(dim)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    d_n = 10. ** 20
    orientation = INT(b_normal(:,:ntotalb))
    xb_temp = xb(:,:ntotalb) / 100 !Too large number lead to unprecision

    !Find 2 nearest neighbour boundary particles
    DO i = 1, ntotalb
        DO j = 1, ntotalb

            !Skip the current particle
            IF (j .NE. i) THEN
                d = sqrt((xb(1,i) - xb(1,j))**2 + (xb(2,i) - xb(2,j))**2)

                !Update d_n(1,i)
                IF (d .LT. d_n(1,i)) THEN

                    !Copy d_n(1,i) if necessary
                    IF (d_n(1,i) .LT. d_n(2,i)) THEN

                        d_n(2,i) = d_n(1,i)
                        NN(2,i) = NN(1,i)

                    ENDIF

                    d_n(1,i) = d 
                    NN(1,i) = j

                !Update d_n(2,i)
                ELSEIF (d .LT. d_n(2,i)) THEN

                    d_n(2,i) = d
                    NN(2,i) = j

                ENDIF
            ENDIF
        ENDDO
    ENDDO

    !Get the normal slope for each boundary particle
    DO i = 1, ntotalb

        !Check if all the neighbour are on the same x coordinates
        IF( ABS((xb(1,i) + xb(1,NN(1,i)) + xb(1,NN(2,i)))/3. - xb(1,i))    &
            .LT. 1D-20 ) THEN

            !Get the unit normal vector
            b_normal(1,i) = 1.
            b_normal(2,i) = 0.
        
        !Check if all the neighbour are on the same y coordinates
        ELSEIF( ABS((xb(2,i) + xb(2,NN(1,i)) + xb(2,NN(2,i)))/3. - xb(2,i))&
                .LT. 1D-20 ) THEN

            !Get the unit normal vector
            b_normal(1,i) = 0.
            b_normal(2,i) = 1.

        !If none of the above, compute the normal by solving matrix
        ELSE 
            !Get the square matrix and inverse for this set of coordinates
            matrix(1,1) = xb_temp(1,i) ** 2 + xb_temp(1,NN(1,i)) ** 2 +              &
                          xb_temp(1,NN(2,i))** 2
            matrix(1,2) = xb_temp(1,i) + xb_temp(1,NN(1,i)) + xb(1,NN(2,i)) 
            matrix(2,1) = matrix(1,2)
            matrix(2,2) = 3

            matrix_inv(1,1) = matrix(2,2)
            matrix_inv(1,2) = -matrix(2,1)
            matrix_inv(2,1) = -matrix(1,2)
            matrix_inv(2,2) = matrix(1,1)

            matrix_inv = matrix_inv / (matrix(1,1) * matrix(2,2) -         &
                         matrix(1,2) * matrix(2,1))

            !Compute the tangent slope
            vector(1) = xb_temp(2,i) * xb_temp(1,i) + xb_temp(2,NN(1,i)) * xb_temp(1,NN(1,i)) +&
                        xb_temp(2,NN(2,i)) * xb_temp(1,NN(2,i)) 
            vector(2) = xb_temp(2,i) + xb(2,NN(1,i)) + xb(2,NN(2,i)) 

            tang(1) = matrix_inv(1,1) * vector(1) + matrix_inv(1,2) *      &
                      vector(2)
            tang(2) = matrix_inv(2,1) * vector(1) + matrix_inv(2,2) *      &
                      vector(2)

            !Get the unit normal vector
            b_normal(1,i) = sqrt(1 / ( tang(1) ** (-2) + 1))
            b_normal(2,i) = ABS(b_normal(1,i) /tang(1))

            !Get the proper orientation
            !Shift coeffiecient when necessary
            IF( orientation(1,i) .EQ. 0) THEN

                b_n_temp = b_normal(:,i)
                b_normal(1,i) = b_n_temp(2)
                b_normal(2,i) = b_n_temp(1)

            ENDIF
        ENDIF

        !Check the sign on the x axis
        mask(1) = orientation(1,i)
        mask(2) = orientation(1,NN(1,i))
        mask(3) = orientation(1,NN(2,i))
        IF( ANY(mask .LT. 0) ) THEN

            b_normal(1,i) = -b_normal(1,i)

        ENDIF

        !Check the sign on the y axis
        mask(1) = orientation(2,i)
        mask(2) = orientation(2,NN(1,i))
        mask(3) = orientation(2,NN(2,i))
        IF( ANY(mask .LT. 0) ) THEN
                
            b_normal(2,i) = -b_normal(2,i)

        ENDIF
    ENDDO

    end subroutine boundary_normal

end subroutine boundaries

