!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the boundary kernel                !
!                        applied to the particles.                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine boundary_kernel

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Call the desired kernel
    IF (skf .EQ. 1) THEN
        call b_Quadratic

    ELSEIF (skf .EQ. 2) THEN
        call b_Gaussian

    ELSEIF (skf .EQ. 3) THEN
        call b_quintic_spline

    ELSEIF (skf .EQ. 4) THEN
        call b_Wendland_C2

    ELSEIF (skf .EQ. 5) THEN
        call b_Wendland_C4

    ELSEIF (skf .EQ. 6) THEN
        call b_Wendland_C6

    ELSEIF (skf .EQ. 7) THEN
        call b_Gaussian

    ELSEIF (skf .EQ. 8) THEN
        call b_Gaussian

    ENDIF

    CONTAINS


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !               Quadratic kernel calculation developped by             !
    !                            (Johnson, 1996)                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine b_Quadratic

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    
    !!!              Variables use in the current routine:               !!!
    !   mass   : Mass of particles                                      [in]
    !   massb   : Mass of boundary particles                            [in]
    !   vx     : Velocities of all particles                            [in]
    !   bp_i : List of boundary particle for boundary force interaction [in]
    !   bp_j    : List of particle for boundary force interaction       [in]
    !   bnip    : number of interacting particle for boundary force     [in]
    !   br_ij   : distance between particles for boundary force         [in]
    !   bdx      : Distances between particles and boundaries           
    !   dvxdt : Acceleration                                        [in/out]
    !   loop_index : Individual index for loop
    !   cp_i : Current particle i ID when computing interactions
    !   cp_j : Current particle j ID when computing interactions
    !   q : Normalize distance for the interaciton
    !   factor : Kernel function factor related to dimension
    !   x      : Coordinates of all particles                           [in]
    !   xb      : Coordinates of all boundary particles particles       [in]
    !   bw : Current boundary interaction kernel                     
    !   bdwdx : Current boundary interaction kernel derivative       
    !   b_length  : Boundary particle length of interaction             [in]
    !   kappa_n   : Boundary normal force parameter                     [in]
    !   kappa_t   : Boundary tangential force parameter                 [in]

    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Set parallel environment
    !$OMP PARALLEL DO SCHEDULE(STATIC) REDUCTION(+:dvxdt)
    !Loop on all boundary interactions
    DO loop_index = 1, bnip

        ! Get index
        cp_i = bp_j(loop_index) !we want i to be associated with the particle
        cp_j = bp_i(loop_index) !and j to be associated with the boundary particle
        int_id = loop_index

        !Compute their relative distance 
        bdx(:) = xb(:,cp_j) - x(:,cp_i)

        !Compute the kernel
        q = br_ij(loop_index) / b_length * scale_h
        factor =  2. / (pi) * (scale_h / b_length) ** 2

        ![0,2]:
        bw = factor * ( 3./16. * q ** 2 - 3./4. * q + 3./4. ) 
        bdwdx(:) = - factor * ( 3./8. * q - 3./4.) / b_length * (bdx(:) /  &
                     br_ij(loop_index)) * scale_h

        !Compute force
        call boundary_force

    ENDDO
    !$OMP END PARALLEL DO

    end subroutine b_Quadratic

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !          Original Gaussian kernel calculation developped by          !
    !         (Gingold & Monaghan 1981) for boundary interactions          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine b_Gaussian

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    
    !!!              Variables use in the current routine:               !!!
    !   mass   : Mass of particles                                      [in]
    !   massb   : Mass of boundary particles                            [in]
    !   vx     : Velocities of all particles                            [in]
    !   bp_i : List of boundary particle for boundary force interaction [in]
    !   bp_j    : List of particle for boundary force interaction       [in]
    !   bnip    : number of interacting particle for boundary force     [in]
    !   br_ij   : distance between particles for boundary force         [in]
    !   bdx      : Distances between particles and boundaries           
    !   dvxdt : Acceleration                                        [in/out]
    !   loop_index : Individual index for loop
    !   cp_i : Current particle i ID when computing interactions
    !   cp_j : Current particle j ID when computing interactions
    !   q : Normalize distance for the interaciton
    !   factor : Kernel function factor related to dimension
    !   x      : Coordinates of all particles                           [in]
    !   xb      : Coordinates of all boundary particles particles       [in]
    !   bw : Current boundary interaction kernel                     
    !   bdwdx : Current boundary interaction kernel derivative       
    !   b_length  : Boundary particle length of interaction             [in]
    !   kappa_n   : Boundary normal force parameter                     [in]
    !   kappa_t   : Boundary tangential force parameter                 [in]

    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Set parallel environment
    !$OMP PARALLEL DO SCHEDULE(STATIC) REDUCTION(+:dvxdt)
    !Gaussian kernel calculation
    DO loop_index = 1, bnip

        ! Get index
        cp_i = bp_j(loop_index) !we want i to be associated with the particle
        cp_j = bp_i(loop_index) !and j to be associated with the boundary particle
        int_id = loop_index
        
        !Compute their relative distance 
        bdx(:) = xb(:,cp_j) - x(:,cp_i)

        !Gaussian kernel calculation
        q = br_ij(loop_index) / b_length * scale_h
        factor = 1. / pi * (scale_h / b_length) ** 2

        ![0,3]:
        bw = factor * exp(-q ** 2)
        bdwdx(:) = - factor * (- 2 * q * exp(-q ** 2)) /         &
        b_length * (bdx(:) / br_ij(loop_index)) &
        * scale_h

        !Compute force
        call boundary_force

    ENDDO
    !$OMP END PARALLEL DO

    end subroutine b_Gaussian

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !            Quintic spline kernel calculation developped by           !
    !               (Morris 1997) for boundary interactions                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine b_quintic_spline

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    
    !!!              Variables use in the current routine:               !!!
    !   mass   : Mass of particles                                      [in]
    !   massb   : Mass of boundary particles                            [in]
    !   vx     : Velocities of all particles                            [in]
    !   bp_i : List of boundary particle for boundary force interaction [in]
    !   bp_j    : List of particle for boundary force interaction       [in]
    !   bnip    : number of interacting particle for boundary force     [in]
    !   br_ij   : distance between particles for boundary force         [in]
    !   bdx      : Distances between particles and boundaries           
    !   dvxdt : Acceleration                                        [in/out]
    !   loop_index : Individual index for loop
    !   cp_i : Current particle i ID when computing interactions
    !   cp_j : Current particle j ID when computing interactions
    !   q : Normalize distance for the interaciton
    !   factor : Kernel function factor related to dimension
    !   x      : Coordinates of all particles                           [in]
    !   xb      : Coordinates of all boundary particles particles       [in]
    !   bw : Current boundary interaction kernel                     
    !   bdwdx : Current boundary interaction kernel derivative       
    !   b_length  : Boundary particle length of interaction             [in]
    !   kappa_n   : Boundary normal force parameter                     [in]
    !   kappa_t   : Boundary tangential force parameter                 [in]

    implicit none

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    !Set parallel environment
    !$OMP PARALLEL DO SCHEDULE(STATIC) REDUCTION(+:dvxdt)
    !Loop on all boundary interaction
    DO loop_index = 1, bnip
        
        ! Get index
        cp_i = bp_j(loop_index) !we want i to be associated with the particle
        cp_j = bp_i(loop_index) !and j to be associated with the boundary particle
        int_id = loop_index

        !Compute their relative distance 
        bdx(:) = xb(:,cp_j) - x(:,cp_i)

        !Quintic kernel calculation
        q = br_ij(loop_index) / b_length * scale_h
        factor = 15. / (478. * pi) * (scale_h / b_length) ** 2

        ![0,1]:
        IF ((q .GE. 0.) .AND. (q .LE. 1.)) THEN
            bw = factor * ( (3 - q) ** 5 - 6 * (2 - q) ** 5 +   &
                    15 * (1 - q) ** 5 )
            bdwdx(:) = - factor * ( -5 * (3 - q) ** 4 + 30 *     &
                        ( 2 - q ) ** 4 - 75 * ( 1 - q ) ** 4 ) /   &
                        b_length * ( bdx(:) / br_ij(loop_index) ) * scale_h
        ENDIF

        !]1,2]:
        IF ((q .GT. 1.) .AND. (q .LE. 2.)) THEN
            bw = factor * ( (3 - q) ** 5 - 6 * (2 - q) ** 5 )
            bdwdx(:) = - factor * ( -5 * (3 - q ) ** 4 + 30 *    &
                    ( 2 - q ) ** 4 ) / b_length * ( bdx(:) / br_ij(loop_index) ) * scale_h
        
        ENDIF

        !]2,3]:
        IF ((q .GT. 2.) .AND. (q .LE. 3.)) THEN
            bw = factor * ( (3 - q) ** 5 )
            bdwdx(:) = - factor * ( -5 * (3 - q ) ** 4 )          &
                        / b_length * (bdx(:) / br_ij(loop_index)) * scale_h
        ENDIF

        !Compute force
        call boundary_force
        
    ENDDO
    !$OMP END PARALLEL DO

    end subroutine b_quintic_spline

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Wendland C4 kernel calculation developped by             !
    !              (Wendland 1995) for boundary interactions               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine b_Wendland_C2

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    
    !!!              Variables use in the current routine:               !!!
    !   mass   : Mass of particles                                      [in]
    !   massb   : Mass of boundary particles                            [in]
    !   vx     : Velocities of all particles                            [in]
    !   bp_i : List of boundary particle for boundary force interaction [in]
    !   bp_j    : List of particle for boundary force interaction       [in]
    !   bnip    : number of interacting particle for boundary force     [in]
    !   br_ij   : distance between particles for boundary force         [in]
    !   bdx      : Distances between particles and boundaries           
    !   dvxdt : Acceleration                                        [in/out]
    !   loop_index : Individual index for loop
    !   cp_i : Current particle i ID when computing interactions
    !   cp_j : Current particle j ID when computing interactions
    !   q : Normalize distance for the interaciton
    !   factor : Kernel function factor related to dimension
    !   x      : Coordinates of all particles                           [in]
    !   xb      : Coordinates of all boundary particles particles       [in]
    !   bw : Current boundary interaction kernel                     
    !   bdwdx : Current boundary interaction kernel derivative       
    !   b_length  : Boundary particle length of interaction             [in]
    !   kappa_n   : Boundary normal force parameter                     [in]
    !   kappa_t   : Boundary tangential force parameter                 [in]

    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Set parallel environment
    !$OMP PARALLEL DO SCHEDULE(STATIC) REDUCTION(+:dvxdt)
    !Loop on all boundary interactions
    DO loop_index = 1, bnip

        ! Get index
        cp_i = bp_j(loop_index) !we want i to be associated with the particle
        cp_j = bp_i(loop_index) !and j to be associated with the boundary particle
        int_id = loop_index

        !Compute their relative distance 
        bdx(:) = xb(:,cp_j) - x(:,cp_i)

        !Compute the kernel
        q = br_ij(loop_index) / b_length * scale_h
        factor = 3. / (2. * pi) * (scale_h / b_length) ** 2

        ![0,2]:
        bw = factor * ( (1. - q) ** 4 * (1 + 4 * q ) )

        bdwdx(:) = - factor * ( 20. * ( q ** 4 - 3 * q ** 3 + 3 * q ** 2 &
                - q) ) / b_length * (bdx(:) / br_ij(loop_index)) * scale_h

        !Compute force
        call boundary_force

    ENDDO
    !$OMP END PARALLEL DO

    end subroutine b_Wendland_C2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Wendland C4 kernel calculation developped by             !
    !              (Wendland 1995) for boundary interactions               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine b_Wendland_C4

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    
    !!!              Variables use in the current routine:               !!!
    !   mass   : Mass of particles                                      [in]
    !   massb   : Mass of boundary particles                            [in]
    !   vx     : Velocities of all particles                            [in]
    !   bp_i : List of boundary particle for boundary force interaction [in]
    !   bp_j    : List of particle for boundary force interaction       [in]
    !   bnip    : number of interacting particle for boundary force     [in]
    !   br_ij   : distance between particles for boundary force         [in]
    !   bdx      : Distances between particles and boundaries           
    !   dvxdt : Acceleration                                        [in/out]
    !   loop_index : Individual index for loop
    !   cp_i : Current particle i ID when computing interactions
    !   cp_j : Current particle j ID when computing interactions
    !   q : Normalize distance for the interaciton
    !   factor : Kernel function factor related to dimension
    !   x      : Coordinates of all particles                           [in]
    !   xb      : Coordinates of all boundary particles particles       [in]
    !   bw : Current boundary interaction kernel                     
    !   bdwdx : Current boundary interaction kernel derivative       
    !   b_length  : Boundary particle length of interaction             [in]
    !   kappa_n   : Boundary normal force parameter                     [in]
    !   kappa_t   : Boundary tangential force parameter                 [in]

    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Set parallel environment
    !$OMP PARALLEL DO SCHEDULE(STATIC) REDUCTION(+:dvxdt)
    !Loop on all boundary interactions
    DO loop_index = 1, bnip

        ! Get index
        cp_i = bp_j(loop_index) !we want i to be associated with the particle
        cp_j = bp_i(loop_index) !and j to be associated with the boundary particle
        int_id = loop_index

        !Compute their relative distance 
        bdx(:) = xb(:,cp_j) - x(:,cp_i)

        !Compute the kernel
        q = br_ij(loop_index) / b_length * scale_h
        factor = 3. / pi * (scale_h / b_length) ** 2

        ![0,2]:
        bw = factor * ((1-q)**6 * (35*q**2 + 18*q + 3))

        bdwdx(:) = - factor * (-56*q * (5*q + 1) * (1 - q)**5) / &
                b_length * (bdx(:) / br_ij(loop_index)) * scale_h

        !Compute force
        call boundary_force

    ENDDO
    !$OMP END PARALLEL DO

    end subroutine b_Wendland_C4

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Wendland C4 kernel calculation developped by             !
    !              (Wendland 1995) for boundary interactions               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine b_Wendland_C6

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    
    !!!              Variables use in the current routine:               !!!
    !   mass   : Mass of particles                                      [in]
    !   massb   : Mass of boundary particles                            [in]
    !   vx     : Velocities of all particles                            [in]
    !   bp_i : List of boundary particle for boundary force interaction [in]
    !   bp_j    : List of particle for boundary force interaction       [in]
    !   bnip    : number of interacting particle for boundary force     [in]
    !   br_ij   : distance between particles for boundary force         [in]
    !   bdx      : Distances between particles and boundaries           
    !   dvxdt : Acceleration                                        [in/out]
    !   loop_index : Individual index for loop
    !   cp_i : Current particle i ID when computing interactions
    !   cp_j : Current particle j ID when computing interactions
    !   q : Normalize distance for the interaciton
    !   factor : Kernel function factor related to dimension
    !   x      : Coordinates of all particles                           [in]
    !   xb      : Coordinates of all boundary particles particles       [in]
    !   bw : Current boundary interaction kernel                     
    !   bdwdx : Current boundary interaction kernel derivative       
    !   b_length  : Boundary particle length of interaction             [in]
    !   kappa_n   : Boundary normal force parameter                     [in]
    !   kappa_t   : Boundary tangential force parameter                 [in]

    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Set parallel environment
    !$OMP PARALLEL DO SCHEDULE(STATIC) REDUCTION(+:dvxdt)
    !Loop on all boundary interactions
    DO loop_index = 1, bnip

        ! Get index
        cp_i = bp_j(loop_index) !we want i to be associated with the particle
        cp_j = bp_i(loop_index) !and j to be associated with the boundary particle
        int_id = loop_index

        !Compute their relative distance 
        bdx(:) = xb(:,cp_j) - x(:,cp_i)

        !Compute the kernel
        q = br_ij(loop_index) / b_length * scale_h
        factor = 78. / (7. * pi) * ( scale_h / hsml_ij(int_id)) ** 2

        ![0,1]:
        bw = factor * (1 - q)**8 * (32 * q**3 + 25 * q**2 + 8 * q + 1)
        bdwdx(:) = - factor * ( -22.* q * (1 - q)**7 * (16 * q**2 + 7.* q &
        + 1)) / b_length * (bdx(:) / br_ij(loop_index)) * scale_h

        !Compute force
        call boundary_force

    ENDDO
    !$OMP END PARALLEL DO

    end subroutine b_Wendland_C6


end subroutine boundary_kernel

