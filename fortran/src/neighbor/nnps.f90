!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the nearest neighboring           !
!                                particle search                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine NNPS

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Call the desired nnps algorithm
    IF (nnps_algorithm .EQ. 1) THEN
        call direct_find

    ELSEIF (nnps_algorithm .EQ. 2) THEN
        call bucket_search

    ELSEIF (nnps_algorithm .EQ. 3) THEN
        call b_tree

    ENDIF



    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       Subprograme to calculate the distance and the interacting      !
    !                  pair with the direct find algorithm                 !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine direct_find

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   ntotal : Total number of particle to use                        [in]
    !   x      : Coordinates of all particles                           [in]
    !   hsml   : Smoothing lengths of particles                         [in]
    !   p_i  : List of first particle partner of interaction pair      [out]
    !   p_j  : List of second particle partner of interaction pair     [out]
    !   nip  : number of interacting particle                          [out]
    !   r_ij : distance between the interacting particles              [out]
    !   hsml_ij : mean smoothing length between particle i and j       [out]
    !   dx       : Distances in x and y between particles              [out]
    !   dvx      : Velocity difference in x and y between particles    [out]
    !   bdx     : Distances between particles and boundaries           [out]
    !   bp_i : List of boundary particle for boundary interaction      [out]
    !   bp_j : List of particle for boundary force interaction         [out]
    !   bnip : number of interacting particle for boundary force       [out]
    !   br_ij : distance between particles for boundary force          [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   i, j : index for do loops
    !   count_np : number of interaction per particles                 [out]
    !   d : distance of particle
    !   mhsml : mean smoothing length

    INTEGER          :: i, j, count_np(maxn)
    DOUBLE PRECISION :: d, mhsml

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    bnip = 0.
    nip = 0.
    count_np = 0.
    d = 0.
    mhsml = 0.

    !Compute interaction:
    DO i = 1, ntotal - 1
        DO j = i+1, ntotal
            d     = SQRT(SUM((x(:,i) - x(:,j)) ** 2))
            mhsml = (hsml(i)+hsml(j))/2.

            !Verify if particle are in their support domain:
            IF (d .LE. mhsml) THEN

                !Store interacting particle information
                nip          = nip + 1
                p_i(nip)     = i
                p_j(nip)     = j
                r_ij(nip)    = d
                hsml_ij(nip) = mhsml
                count_np(i)  = count_np(i) + 1
                count_np(j)  = count_np(j) + 1

            ENDIF

        ENDDO
    ENDDO

    !If necessary get boundary interaction
    IF (boundary_pressure .OR. boundary_friction) THEN

        !Compute particle/boundary interaction:
        DO CONCURRENT ( i = 1:ntotalb, j = 1:ntotal)
            d = SQRT(SUM((xb(:,i) - x(:,j)) ** 2))

            !Verify if particle are in their support domain:
            IF (d .LE. b_length) THEN

                !Store interacting particle information
                bnip         = bnip + 1
                bp_i(bnip)   = i
                bp_j(bnip)   = j
                br_ij(bnip)  = d
                count_np(j)  = count_np(j) + 1

            ENDIF

        ENDDO
    ENDIF

    !Verify if the maximum number of interaction has been reach
    IF ( (nip .GT. max_interaction) .OR. (bnip .GT. max_interaction) ) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) '    Too much particle interaction please allocate     '
        WRITE (*,*) '      more memory for arrays in globals module.       '
        WRITE (*,*) '******************************************************'

        STOP
    ENDIF


    !Writing the monitored quantities:
    IF (mod(itimestep,print_step).eq.0) THEN
        call interaction_stats(count_np)

    ENDIF

    end subroutine direct_find

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       Subprograme to calculate the distance and the interacting      !
    !           pair with bucket search algorithm (Rhoades,1991)           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bucket_search

    !!!                         Module to include:                       !!!
    use globals
    use domain
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   ntotal : Total number of particle to use                        [in]
    !   x      : Coordinates of all particles                           [in]
    !   hsml   : Smoothing lengths of particles                         [in]
    !   p_i  : List of first particle partner of interaction pair      [out]
    !   p_j  : List of second particle partner of interaction pair     [out]
    !   nip  : number of interacting particle                          [out]
    !   r_ij : distance between the interacting particles              [out]
    !   hsml_ij : mean smoothing length between particle i and j       [out]
    !   dx       : Distances in x and y between particles              [out]
    !   dvx      : Velocity difference in x and y between particles    [out]
    !   bdx     : Distances between particles and boundaries           [out]
    !   bp_i : List of boundary particle for boundary interaction      [out]
    !   bp_j : List of particle for boundary force interaction         [out]
    !   bnip : number of interacting particle for boundary force       [out]
    !   br_ij : distance between particles for boundary force          [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   i, j, k, l, m : index for do loops
    !   ID : Particle interacting index
    !   count_np : number of interaction per particles                 [out]
    !   d : distance of particle
    !   mhsml : mean smoothing length
    !   x_bin : Bin number of particles in x direction
    !   y_bin : Bin number of particles in y direction
    !   bin_size: Length of the square bin
    !   number_bin : Array containing the number of particle
    !   bins : Bin array containing particle ID in those bins
    !   bucket_coor : min and max bucket ID in x and y direction
    !   max_int_pp  : Maximum interaction per processes

    INTEGER          :: i, j, k, l, count_np(maxn), x_bin, y_bin, ID(2),   &
                        bucket_coor(4), max_int_pp, count_np_temp(maxn)
    INTEGER, ALLOCATABLE :: bins(:,:,:), number_bin(:,:)
    DOUBLE PRECISION :: d, mhsml, bin_size

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    bnip = 0
    nip = 0
    count_np = 0

    !Get bin size and bucket coord
    bin_size = 2 * MINVAL(hsml(:ntotal))
    bucket_coor = (/                                                       &
        NINT(MINVAL([xb(1,:ntotalb), x(1,:ntotal)]) / bin_size),           &
        NINT(MAXVAL([xb(1,:ntotalb), x(1,:ntotal)]) / bin_size),           &
        NINT(MINVAL([xb(2,:ntotalb), x(2,:ntotal)]) / bin_size),           &
        NINT(MAXVAL([xb(2,:ntotalb), x(2,:ntotal)]) / bin_size) /)

    ALLOCATE ( number_bin(bucket_coor(1) - 1 : bucket_coor (2) + 1,        &
    bucket_coor(3) - 1 : bucket_coor (4) + 1) )
    number_bin = 0

    ALLOCATE ( bins(1 : max_ind_int,                                       &
               bucket_coor(1) - 1 : bucket_coor (2) + 1,                   &
               bucket_coor(3) - 1 : bucket_coor (4) + 1) )

    !Initialize parallel region
    !$OMP PARALLEL &
    !$OMP PRIVATE (i, j, k, l, d, mhsml, ID, x_bin, y_bin)&
    !$OMP PRIVATE (count_np_temp, max_int_pp)

    process_rank = OMP_GET_THREAD_NUM()
    cluster_size = OMP_GET_NUM_THREADS()

    !Allocate private mpi variables
    max_int_pp = NINT(REAL(max_interaction)/cluster_size)
    ALLOCATE( p_i_temp(max_int_pp),                &
              p_j_temp(max_int_pp),                &
              r_ij_temp(max_int_pp),               &
              hsml_ij_temp(max_int_pp) )

    !Initialize private variables
    nip_temp = 0
    count_np_temp = 0


    !$OMP SINGLE
    !Get the number of particles and their ID in bins
    DO i = 1, ntotal
        x_bin = NINT(x(1,i) / bin_size)
        y_bin = NINT(x(2,i) / bin_size)
        number_bin(x_bin,y_bin) = number_bin(x_bin,y_bin) + 1
        bins(number_bin(x_bin,y_bin), x_bin, y_bin) = i
    ENDDO
    !$OMP END SINGLE

    !Get chunk size:
    chunk = CEILING( REAL((bucket_coor(2) - bucket_coor(1))) / cluster_size)

    !$OMP DO COLLAPSE(2) SCHEDULE(DYNAMIC, chunk) 
    !Get the interaction parameters
    DO l = bucket_coor(3), bucket_coor(4)
        DO k = bucket_coor(1), bucket_coor(2)

            !Loop on the particle i in the cell
            DO i = 1, number_bin(k,l) - 1

                !Get i particle ID
                ID(1) = bins(i, k, l)

                !Compute interaction within the cell
                DO j = i + 1, number_bin(k,l)

                    !Get j particle ID
                    ID(2) = bins(j, k, l)

                    !Compute distance and interaction length
                    d     = SQRT(SUM((x(:,ID(1)) - x(:,ID(2))) ** 2))
                    mhsml = (hsml(ID(1))+ hsml(ID(2)))/2.
        
                    !Verify if particle are in their support domain:
                    IF (d .LE. mhsml) THEN

                        !Store interacting particle information
                        nip_temp               = nip_temp + 1
                        p_i_temp(nip_temp)     = ID(1)
                        p_j_temp(nip_temp)     = ID(2)
                        r_ij_temp(nip_temp)    = d
                        hsml_ij_temp(nip_temp) = mhsml
                        count_np_temp(ID(1))   = count_np_temp(ID(1)) + 1
                        count_np_temp(ID(2))   = count_np_temp(ID(2)) + 1


                    ENDIF
                
                ENDDO
            
            ENDDO

            !Loop on the particle i in the cell
            DO i = 1, number_bin(k,l)

                !Get i particle ID
                ID(1) = bins(i, k, l)

                !Compute interaction with the righter cell
                DO j = 1, number_bin(k + 1, l)

                    !Get j particle ID
                    ID(2) = bins(j, k + 1, l)

                    !Compute distance and interaction length
                    d     = SQRT(SUM((x(:,ID(1)) - x(:,ID(2))) ** 2))
                    mhsml = (hsml(ID(1))+hsml(ID(2)))/2.
        
                    !Verify if particle are in their support domain:
                    IF (d .LE. mhsml) THEN

                        !Store interacting particle information
                        nip_temp          = nip_temp + 1
                        p_i_temp(nip_temp)    = ID(1)
                        p_j_temp(nip_temp)    = ID(2)
                        r_ij_temp(nip_temp)    = d
                        hsml_ij_temp(nip_temp) = mhsml
                        count_np_temp(ID(1))  = count_np_temp(ID(1)) + 1
                        count_np_temp(ID(2))  = count_np_temp(ID(2)) + 1

                    ENDIF
                
                ENDDO

                !Compute interaction with upper cell
                DO j = 1, number_bin(k, l + 1)

                    !Get j particle ID
                    ID(2) = bins(j, k, l + 1)

                    !Compute distance and interaction length
                    d     = SQRT(SUM((x(:,ID(1)) - x(:,ID(2))) ** 2))
                    mhsml = (hsml(ID(1))+hsml(ID(2)))/2.
        
                    !Verify if particle are in their support domain:
                    IF (d .LE. mhsml) THEN

                        !Store interacting particle information
                        nip_temp          = nip_temp + 1
                        p_i_temp(nip_temp)    = ID(1)
                        p_j_temp(nip_temp)    = ID(2)
                        r_ij_temp(nip_temp)    = d
                        hsml_ij_temp(nip_temp) = mhsml
                        count_np_temp(ID(1))  = count_np_temp(ID(1)) + 1
                        count_np_temp(ID(2))  = count_np_temp(ID(2)) + 1

                    ENDIF
                
                ENDDO

                !Compute interaction with upper righter cell
                DO j = 1, number_bin(k + 1, l + 1)

                    !Get j particle ID
                    ID(2) = bins(j, k + 1, l + 1)

                    !Compute distance and interaction length
                    d     = SQRT(SUM((x(:,ID(1)) - x(:,ID(2))) ** 2))
                    mhsml = (hsml(ID(1))+hsml(ID(2)))/2.
        
                    !Verify if particle are in their support domain:
                    IF (d .LE. mhsml) THEN

                        !Store interacting particle information
                        nip_temp          = nip_temp + 1
                        p_i_temp(nip_temp)    = ID(1)
                        p_j_temp(nip_temp)    = ID(2)
                        r_ij_temp(nip_temp)    = d
                        hsml_ij_temp(nip_temp) = mhsml
                        count_np_temp(ID(1))  = count_np_temp(ID(1)) + 1
                        count_np_temp(ID(2))  = count_np_temp(ID(2)) + 1

                    ENDIF

                ENDDO

                !Compute interaction with upper lefter cell
                DO j = 1, number_bin(k - 1, l + 1)

                    !Get j particle ID
                    ID(2) = bins(j, k - 1, l + 1)

                    !Compute distance and interaction length
                    d     = SQRT(SUM((x(:,ID(1)) - x(:,ID(2))) ** 2))
                    mhsml = (hsml(ID(1))+hsml(ID(2)))/2.
        
                    !Verify if particle are in their support domain:
                    IF (d .LE. mhsml) THEN

                        !Store interacting particle information
                        nip_temp           = nip_temp + 1
                        p_i_temp(nip_temp) = ID(1)
                        p_j_temp(nip_temp) = ID(2)
                        r_ij_temp(nip_temp)    = d
                        hsml_ij_temp(nip_temp) = mhsml
                        count_np_temp(ID(1))  = count_np_temp(ID(1)) + 1
                        count_np_temp(ID(2))  = count_np_temp(ID(2)) + 1
                        
                    ENDIF
                
                ENDDO
            ENDDO
        ENDDO

    ENDDO
    !$OMP END DO NOWAIT

    !Verify if the maximum number of interaction has been reach
    IF ( nip_temp .GT. max_int_pp ) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) ' Too much particle interaction in processe :',        &
                      process_rank
        WRITE (*,*) '       Please allocate more memory for arrays.        '
        WRITE (*,*) '******************************************************'

        STOP
    ENDIF

    !Aggregate the information into the globals array
    !$OMP CRITICAL
    p_i(nip + 1 : nip + nip_temp) = p_i_temp(1:nip_temp)
    p_j(nip + 1 : nip + nip_temp) = p_j_temp(1:nip_temp)
    r_ij(nip + 1 : nip + nip_temp) = r_ij_temp(1:nip_temp)
    hsml_ij(nip + 1: nip + nip_temp) = hsml_ij_temp(1:nip_temp)
    count_np  = count_np + count_np_temp
    nip = nip + nip_temp
    !$OMP END CRITICAL

    !Deallocate everything
    DEALLOCATE( p_i_temp, p_j_temp, r_ij_temp, hsml_ij_temp)


    !If necessary get boundary interaction
    IF (boundary_pressure .OR. boundary_friction) THEN

        !Allocate everything for boundary interaction
        ALLOCATE( bp_i_temp(max_int_pp),                &
                  bp_j_temp(max_int_pp),                &
                  br_ij_temp(max_int_pp) )

        !Initialize private variables
        bnip_temp = 0
        count_np_temp = 0

        !$OMP DO SCHEDULE(DYNAMIC)
        !Get the interaction parameters
        DO i = 1, ntotalb

            !Get i particle ID
            ID(1) = i

            !Get current bin of boundary particle
            x_bin = NINT(xb(1,ID(1)) / bin_size)
            y_bin = NINT(xb(2,ID(1)) / bin_size)

            !Loop on the interacting cell 
            DO k = y_bin - 1, y_bin + 1
                DO j = x_bin - 1, x_bin + 1 

                    !Loop on the particle in the cell
                    DO l = 1, number_bin(j,k)

                        !Get j particle ID
                        ID(2) = bins(l, j, k)

                        !Compute distance and interaction length
                        d = SQRT(SUM((xb(:,ID(1)) - x(:,ID(2))) ** 2))

                        !Verify if particle are in their support domain:
                        IF (d .LE. b_length) THEN
            
                            !Store interacting particle information
                            bnip_temp         = bnip_temp + 1
                            bp_i_temp(bnip_temp)   = ID(1)
                            bp_j_temp(bnip_temp)   = ID(2)
                            br_ij_temp(bnip_temp)  = d
                            count_np_temp(ID(2))  = count_np_temp(ID(2)) + 1
                            
                        ENDIF
                    
                    ENDDO
                
                ENDDO
            ENDDO

        ENDDO
        !$OMP END DO NOWAIT

        !$OMP CRITICAL
        bp_i(bnip + 1 : bnip + bnip_temp) = bp_i_temp(1:bnip_temp)
        bp_j(bnip + 1 : bnip + bnip_temp) = bp_j_temp(1:bnip_temp)
        br_ij(bnip + 1 : bnip + bnip_temp) = br_ij_temp(1:bnip_temp)
        count_np  = count_np + count_np_temp
        bnip = bnip + bnip_temp
        !$OMP END CRITICAL
    
        !Deallocate everything
        DEALLOCATE( bp_i_temp, bp_j_temp, br_ij_temp)


    ENDIF

    !$OMP END PARALLEL

    !Writing the monitored quantities:
    IF (mod(itimestep,print_step).eq.0) THEN
        call interaction_stats(count_np)

    ENDIF

    end subroutine bucket_search

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !       Subprograme to calculate the distance and the interacting      !
    !                     pair with the b-tree algorithm                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine b_tree

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use tree_module
        
    !!!              Variables use in the current routine:               !!!
    !   ntotal : Total number of particle to use                        [in]
    !   x      : Coordinates of all particles                           [in]
    !   hsml   : Smoothing lengths of particles                         [in]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   root : First layer of the Tree strucure
    !   domain : Domaine incorporating all particle
    !   max_c_hsml : Maximal current smoothing length value
    !   count_np : number of interaction per particles                 [out]
    !   switch : Switch to avoid counting interaction twice
    !   ID : Particle index
    !   i : Counter

    INTEGER          :: switch(ntotal), ID(ntotal) , i , count_np(maxn)
    DOUBLE PRECISION :: domain(dim,2), max_c_hsml
    TYPE (node), pointer :: root

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialize
    switch = 1
    ID = (/(i, i=1,ntotal, 1)/)
    bnip = 0.
    nip = 0.
    count_np = 0.
    domain = 0.

    !Get domain
    !The +/- hsml is used because domain needs to encompass all particle
    max_c_hsml = maxval(hsml(:ntotal)) + 1
    domain = reshape( (/ MINVAL(x(1, :ntotal)) - max_c_hsml,             &
                         MINVAL(x(2, :ntotal)) - max_c_hsml,             &
                         MAXVAL(x(1, :ntotal)) + max_c_hsml,             &
                         MAXVAL(x(2, :ntotal)) + max_c_hsml /),(/dim,2/) )

    !Set the first node for the tree structure
    allocate(root)
    root%parent_branch => null()

    !Create the tree structure
    call build_tree(ntotal, ID, domain, root)

    !Compute the particle-particle interactions
    call tree_loop(root, switch, count_np)

    !Compute the boundary-particle interactions
    IF (boundary_pressure .OR. boundary_friction) THEN
        call boundary_tree_loop(root,count_np)
    ENDIF

    !Deallocate tree structure
    call delete_tree(root)

    !Verify if the maximum number of interaction has been reach
    IF ( (nip .GT. max_interaction) .OR. (bnip .GT. max_interaction) ) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) '    Too much particle interaction please allocate     '
        WRITE (*,*) '      more memory for arrays in globals module.       '
        WRITE (*,*) '******************************************************'

        STOP
    ENDIF

    !Writing the monitored quantities:
    IF (mod(itimestep,print_step).eq.0) THEN
        call interaction_stats(count_np)

    ENDIF

    end subroutine b_tree

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme to calculate the particle interaction        !
    !                                 statistics                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    subroutine interaction_stats(count_ip)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
            
    !!!              Variables use in the current routine:               !!!
    !   ntotal : Total number of particle to use                        [in]

    implicit none
        
    !!!                        Local variables:                          !!!
    !   i : index for do loops
    !   sumip : Sum of interacting particles                     
    !   maxip : Maximum of single particule interactions            
    !   minip : Minimum of single particule interactions            
    !   noip  : Number of particle with no interaction
    !   count_np : number of interaction per particles                  [in]
    
    INTEGER :: sumip, maxip, minip, noip, count_ip(maxn)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    sumip = 0
    maxip = 0
    minip = 0
    noip  = 0

    !Calculation of quantities to monitor:
    sumip = SUM(count_ip(:ntotal))
    maxip = MAXVAL(count_ip(:ntotal))
    minip = MINVAL(count_ip(:ntotal))
    noip  = COUNT((count_ip(:ntotal) .EQ. 0))


    WRITE(*,*) '**********************************************************'
    WRITE(*,*) '          Statistics interactions per particle:           '
    WRITE(*,*) 'Maximal number interactions:', maxip
    WRITE(*,*) 'Minimal number interactions:', minip
    WRITE(*,*) 'Average number interactions:',real(sumip)/real(ntotal)
    WRITE(*,*) 'Total number of pairs : ', nip
    WRITE(*,*) 'Number of particles with no interactions:', noip
    WRITE(*,*) '**********************************************************'

    IF (maxip .GT. max_ind_int) THEN
        WRITE (*,*) '******************************************************'
        WRITE (*,*) ' Particle are too close to each other, this can lead  '
        WRITE (*,*) '    to memory issue. Please allow a higher maximum    '
        WRITE (*,*) '     number of interaction for a single particle.     '
        WRITE (*,*) '******************************************************'

        STOP
    ENDIF

    end subroutine interaction_stats

end subroutine NNPS