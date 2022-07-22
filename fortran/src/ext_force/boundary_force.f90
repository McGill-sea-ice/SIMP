!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the boundary force                !
!                        applied to the particles.                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine boundary_force

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
    IF (bfa .EQ. 1) THEN
        call bf_Monaghan

    ELSEIF (bfa .EQ. 2) THEN
        call bf_Wang

    ELSEIF (bfa .EQ. 3) THEN
        call bf_Marquis
    ENDIF

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !          Boundary normal and tangential force developped by          !
    !                 (Monaghan 2009) for boundary interactions            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bf_Monaghan

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    
    !!!              Variables use in the current routine:               !!!
    !   mass   : Mass of particles                                      [in]
    !   massb   : Mass of boundary particles                            [in]
    !   vx     : Velocities of all particles                            [in]
    !   br_ij   : distance between particles for boundary force         [in]
    !   bdx      : Distances between particles and boundaries           
    !   dvxdt : Acceleration                                        [in/out]
    !   int_id : Index of the current interaction
    !   cp_i : Current particle i ID when computing interactions
    !   cp_j : Current particle j ID when computing interactions
    !   bw : Current boundary interaction kernel                     
    !   bdwdx : Current boundary interaction kernel derivative       
    !   b_length  : Boundary particle length of interaction             [in]
    !   kappa_n   : Boundary normal force parameter                     [in]
    !   kappa_t   : Boundary tangential force parameter                 [in]

    implicit none

    !!!                        Local variables:                          !!!
    !   F_b : Normal boundary force variable                         [m/s^2]
    !   pi_aj : Artificial viscosity for tangential boundary force    [Pa s]
    !   v_r : Dot product between v and r                            [m^2/s]

    DOUBLE PRECISION :: F_b(dim), v_r, pi_aj

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !   Calculate boundary normal force if requested
    IF (boundary_pressure) THEN

        !Get artificial pressure
        F_b(:) = kappa_n * bdx(:) * bw * 2 * massb(cp_j) / (&
                ( mass(cp_i) + massb(cp_j) ) * br_ij(int_id) ** 2 )

        !Add to external stress component
        dvxdt(:,cp_i) = dvxdt(:,cp_i) + (-1 * F_b)/rho(cp_i)

    ENDIF

    !Calculate boundary tangential force if requested
    IF (boundary_friction) THEN

        !Calculate dot product
        v_r = vx(1,cp_i) * bdx(1) + vx(2,cp_i) * bdx(2)

        !Get artificial viscosity
        pi_aj = kappa_t * c(cp_i) * v_r / br_ij(int_id) / rho(cp_i)

        !Add to external stress component
        dvxdt(:,cp_i) = dvxdt(:,cp_i) + (-1 * massb(cp_j) * pi_aj * &
                         bdwdx(:))/rho(cp_i)


    ENDIF

    end subroutine bf_Monaghan

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !          Boundary normal and tangential force developped by          !
    !                  (Wang 2014) for boundary interactions               !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bf_Wang

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    
    !!!              Variables use in the current routine:               !!!
    !   mass   : Mass of particles                                      [in]
    !   massb   : Mass of boundary particles                            [in]
    !   vx     : Velocities of all particles                            [in]
    !   br_ij   : distance between particles for boundary force         [in]
    !   bdx      : Distances between particles and boundaries           
    !   dvxdt : Acceleration                                        [in/out]
    !   int_id : Index of the current interaction
    !   cp_i : Current particle i ID when computing interactions
    !   cp_j : Current particle j ID when computing interactions
    !   bw : Current boundary interaction kernel                     
    !   bdwdx : Current boundary interaction kernel derivative       
    !   b_length  : Boundary particle length of interaction             [in]
    !   kappa_n   : Boundary normal force parameter                     [in]
    !   kappa_t   : Boundary tangential force parameter                 [in]

    implicit none

    !!!                        Local variables:                          !!!
    !   F : Boundary force variable (normal,tangent)                 [m/s^2]
    !   v : Velocity in normal and tangent direction                 [m^2/s]

    DOUBLE PRECISION :: F(dim), v(dim)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !   Calculate boundary normal force if requested
    IF (boundary_pressure) THEN

        !Get velocity component
        v(1) = (b_normal(1,cp_j) * vx(1,cp_i) + b_normal(2,cp_j) *         &
               vx(2,cp_i))
        v(2) = sqrt(vx(1,cp_i) ** 2 + vx(2,cp_i) ** 2) - v(1)

        !Get repulsive force from boundary
        F(1) = mass(cp_i) * v(1) / dt
        F(2) = mass(cp_i) * v(2) / dt

        !Calculate boundary friction force 
        IF (boundary_friction) THEN

            !Update the tangential velocity if friction
            IF (ABS(v(2)) .GT. ABS(xi * v(1))) THEN

                v(2) = v(2) * xi * ABS(F(1) / F(2))
                F(2) = mass(cp_i) * v(2) / dt

            ENDIF
        ENDIF

        !Add to external stress component
        dvxdt(1,cp_i) = dvxdt(1,cp_i) + (F(1) * b_normal(1,cp_j) + F(2) *   &
                        b_normal(2,cp_j))
        dvxdt(2,cp_i) = dvxdt(2,cp_i) + (F(1) * b_normal(2,cp_j) - F(2) *   &
                        b_normal(1,cp_j))

    ENDIF

    end subroutine bf_Wang

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !          Boundary normal and tangential force developped by          !
    !                 (Marquis 2022) for boundary interactions             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine bf_Marquis

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    
    !!!              Variables use in the current routine:               !!!
    !   vx     : Velocities of all particles                            [in]
    !   br_ij   : distance between particles for boundary force         [in]
    !   bdx      : Distances between particles and boundaries           
    !   dvxdt : Acceleration                                        [in/out]
    !   int_id : Index of the current interaction
    !   cp_i : Current particle i ID when computing interactions
    !   cp_j : Current particle j ID when computing interactions
    !   b_length  : Boundary particle length of interaction             [in]

    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



    !No pressure force from boundary
    IF (boundary_pressure) THEN
        RETURN

    ENDIF

    !Calculate boundary friction
    IF (boundary_friction) THEN
        dvxdt(:,cp_i) = dvxdt(:,cp_i) - vx(:,cp_i) * (1 -  &
                        (br_ij(int_id) / b_length))

    ENDIF

    end subroutine bf_Marquis

end subroutine boundary_force