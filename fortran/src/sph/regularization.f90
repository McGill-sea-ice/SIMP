!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Subroutines to compute particles position regularization          !
!                     to avoid tensile instabilities                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine regularization

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Call the desired regularization algorithm

    !XSPH
    IF (xsph_regularization) THEN
        CALL XSPH

    ENDIF

    !Artificial viscosity
    IF ( art_visc_alg .EQ. 1 ) THEN
        call art_visc_HK

    ELSEIF ( art_visc_alg .EQ. 2 ) THEN
        call art_visc_Mon

    ENDIF


    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !            Subroutine for computing the XSPH correction              !
    !                      developped by (Monaghan 1989)                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine XSPH

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   x      : coordinates of particles                           [in/out]
    !   vx               : velocities of particles                  [in/out]
    !   mass : Mass of particles                                        [in]
    !   rho  : Particle density                                         [in]
    !   w    : Kernel value for particle i                             [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    INTEGER          :: i, j, k
    DOUBLE PRECISION :: xsph_sum(dim,ntotal)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    xsph_sum = 0.

    !Calculate xsph correction
    DO k = 1, nip 

        ! Get index
        i = p_i(k)
        j = p_j(k)

        !Get relative velocity
        dvx(:)  = vx(:, i) - vx(:,j)

        xsph_sum(1,i) = xsph_sum(1,i) + mass(j) * dvx(1) * &
                                w(k) / (0.5 * (rho(i) + rho(j)))
        xsph_sum(2,i) = xsph_sum(2,i) + mass(j) * dvx(2) * &
                                w(k) / (0.5 * (rho(i) + rho(j)))
        xsph_sum(1,j) = xsph_sum(1,j) + mass(i) * dvx(1) * &
                                w(k) / (0.5 * (rho(i) + rho(j)))
        xsph_sum(2,j) = xsph_sum(2,j) + mass(i) * dvx(2) * &
                                w(k) / (0.5 * (rho(i) + rho(j)))

    ENDDO

    dxdt(:,:ntotal) = dxdt(:,:ntotal) - chi * xsph_sum

    end subroutine XSPH

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !           Subroutine for computing the artifical viscosity           !
    !                 proposed by Herquist and Katz(1989)                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine art_visc_HK

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   ntotal : Total number of particles                              [in]
    !   rho : Particles density                                         [in]
    !   hsml : Particle smoothing length                                [in]
    !   dwdx : Kernel gradient                                          [in]
    !   wdvx : Products between kernel and speed (see kernel.f90)       [in]
    !   mass : Particle mass                                            [in]
    !   p_i  : List of first particle partner of interaction pair       [in]
    !   p_j  : List of second particle partner of interaction pair      [in]
    !   nip  : number of interacting particle                           [in]
    !   c : Particle sound speed                                        [in]
    !   alpha_pi, beta_pi : Artificial viscosity parameters             [in]
    !   dvxdt  : dvx/dt of internal force                        [in/out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   Pi_ij : Artificial viscosity                                  [Pa s]
    !   div_v : Particle divergence                                   [s^-1]
    !   k : loop index
    !   i : particle i current index
    !   j : particle j current index
    !   lv : lisibility variable
    !   criteria: criteria to avoid creation of instability

    INTEGER          :: i, j, k
    DOUBLE PRECISION :: div_v(ntotal), pi_ij, lv(ntotal), criteria

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    pi_ij = 0.
    q = 0.
    div_v = 0.
    criteria = 0.

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

    ! No viscosity for diverging flow
    WHERE (div_v .GE. 0.) 
        div_v = 0.
    ENDWHERE

    !Calculate lv
    lv = (alpha_pi * hsml(:ntotal) * c(:ntotal) * ABS(div_v) + beta_pi *    &
        (ABS(div_v) * hsml(:ntotal)) ** 2 ) * rho(:ntotal)
        
    !Calculate artificial viscosity
    DO k = 1, nip 

        ! Get index
        i = p_i(k)
        j = p_j(k)
        
        !Get relative velocity
        dx(:)   = x(:, i) - x(:, j)
        dvx(:)  = vx(:, i) - vx(:,j)

        !Apply criteria from H & K
        criteria = dvx(1) * dx(1) + dvx(2) * dx(2)

        !Get artificial viscosity
        IF (criteria .LE. 0. ) THEN
            pi_ij = lv(i) / rho(i) ** 2 + lv(j) / rho(j) ** 2
        
        ENDIF

        !Particle i artificial viscosity to internal stress x componante
        dvxdt(1,i) = dvxdt(1,i) - mass(j) * pi_ij * dwdx(1,k)

        !Particle j artificial viscosity to internal stress x componante
        dvxdt(1,j) = dvxdt(1,j) + mass(i) * pi_ij * dwdx(1,k)

        !Particle i artificial viscosity to internal stress y componante
        dvxdt(2,i) = dvxdt(2,i) - mass(j) * pi_ij * dwdx(2,k)

        !Particle j artificial viscosity to internal stress y componante
        dvxdt(2,j) = dvxdt(2,j) + mass(i) * pi_ij * dwdx(2,k)


    ENDDO

    end subroutine art_visc_HK

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !           Subroutine for computing the artifical viscosity           !
    !                      proposed by Monaghan (1989)                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine art_visc_Mon

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   rho : Particles density                                         [in]
    !   hsml : Particle smoothing length                                [in]
    !   dwdx : Kernel gradient                                          [in]
    !   mass : Particle mass                                            [in]
    !   p_i  : List of first particle partner of interaction pair       [in]
    !   p_j  : List of second particle partner of interaction pair      [in]
    !   nip  : number of interacting particle                           [in]
    !   c : Particle Sound speed                                        [in]
    !   alpha_pi, beta_pi : Artificial viscosity parameters             [in]
    !   dvx      : Velocity difference in x and y between particles     [in]
    !   dx       : Distances in x and y between particles               [in]
    !   dvxdt  : dvx/dt of internal force                        [in/out]
    !   r_ij    : distance between the interacting particles            [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   Pi_ij : Artificial viscosity                                  [Pa s]
    !   v_r : Dot product between v and r                            [m^2/s]
    !   k : loop index
    !   i : particle i current index
    !   j : particle j current index
    !   iota  : lisibility variable
    !   s_cor : Avoid singularity parameter

    INTEGER          :: i, j, k
    DOUBLE PRECISION :: v_r, pi_ij, iota, s_cor

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    pi_ij = 0.
    iota = 0.
    v_r = 0.
    s_cor = 0.

    !Calculate artificial viscosity
    DO k = 1, nip 

        ! Get index
        i = p_i(k)
        j = p_j(k)
        
        !Get relative velocity
        dx(:)   = x(:, i) - x(:, j)
        dvx(:)  = vx(:, i) - vx(:,j)

        !Calculate dot product
        v_r = dvx(1) * dx(1) + dvx(2) * dx(2)

        !Calculate singularity corrector
        s_cor = 0.01 * hsml_ij(k)

        !Get artificial viscosity
        IF ( v_r .LE. 0. ) THEN
            iota = hsml_ij(k) * v_r / ( r_ij(k) ** 2 +  s_cor ** 2 ) 

            pi_ij = (beta_pi * iota ** 2 - alpha_pi * iota * (c(i) + c(j))/2 ) &
                    / ((rho(i) + rho(j))/2)
        
        ENDIF

        !Particle i artificial viscosity to internal stress x componante
        dvxdt(1,i) = dvxdt(1,i) - mass(j) * pi_ij * dwdx(1,k)

        !Particle j artificial viscosity to internal stress x componante
        dvxdt(1,j) = dvxdt(1,j) + mass(i) * pi_ij * dwdx(1,k)

        !Particle i artificial viscosity to internal stress y componante
        dvxdt(2,i) = dvxdt(2,i) - mass(j) * pi_ij * dwdx(2,k)

        !Particle j artificial viscosity to internal stress y componante
        dvxdt(2,j) = dvxdt(2,j) + mass(i) * pi_ij * dwdx(2,k)


    ENDDO

    end subroutine art_visc_Mon

end subroutine regularization