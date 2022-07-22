!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 Subroutine for computing the density of                  !
!                               each particles                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine density

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Call the desired density definition
    IF (den_alg .EQ. 1) THEN
        call sum_density

    ELSEIF (den_alg .EQ. 2) THEN
        call con_density

    ELSEIF (den_alg .EQ. 3) THEN
        call diagnostic_density

    ENDIF



    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme to calculate the particle density            !
    !               with summation (Randles and Libersky, 1996)            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine sum_density

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   p_i  : List of first particle partner of interaction pair       [in]
    !   p_j  : List of second particle partner of interaction pair      [in]
    !   nip  : number of interacting particle                           [in]
    !   w    : Kernel value for particle i interacting with j           [in]
    !   mass : Mass of particles                                        [in]
    !   rho  : Densities of particles                                  [out]
    !   kernel_0 : Value of the kernel at the given particle location   [in]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   k : loop index
    !   i : particle i current index
    !   j : particle j current index
    !   numerator: sum(mass*kernel)
    !   denominator: sum(mass/density * kernel)

    INTEGER          :: k, i, j
    DOUBLE PRECISION :: numerator(ntotal), denominator(ntotal)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    numerator = 0.
    denominator = 0.

    !Add particle self density to numerator and denominator 
    numerator   = kernel_zero(:ntotal) * mass(:ntotal) 
    denominator = kernel_zero(:ntotal) * mass(:ntotal) / rho(:ntotal)
    
    !Compute the summation
    DO k = 1, nip

        ! Get index
        i = p_i(k)
        j = p_j(k)

        !Firstly calculate the integration of the kernel over the space
        denominator(i) = denominator(i) + mass(j) * w(k) / rho(j)
        denominator(j) = denominator(j) + mass(i) * w(k) / rho(i)

        !Secondly calculate the rho integration over the space
        numerator(i) = numerator(i) + mass(j) * w(k) 
        numerator(j) = numerator(j) + mass(i) * w(k) 
    ENDDO

    ! Thirdly, calculate the normalized rho, rho=sum(rho)/sum(w)
    rho(:ntotal) = numerator /denominator

    end subroutine sum_density

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme to calculate the particle density            !
    !                   with continuity (Monaghan, 2012)                   !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine con_density

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
            
    !!!              Variables use in the current routine:               !!!
    !   p_i  : List of first particle partner of interaction pair       [in]
    !   p_j  : List of second particle partner of interaction pair      [in]
    !   nip  : number of interacting particle                           [in]
    !   mass : Mass of particles                                        [in]
    !   drhodt    : Density rate of change                              [in]
    !   wdvx     : Various products between kernel and speed            [in]

    implicit none
        
    !!!                        Local variables:                          !!!
    !   k : loop index
    !   i : particle i current index
    !   j : particle j current index
    !   v_ij : Lisibility variable

    INTEGER          :: k, i, j
    DOUBLE PRECISION :: v_ij

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    drhodt(:ntotal) = 0.

    !Compute the change in density
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

        !Compute velocity and kernel dot product
        v_ij = wdvx(1) + wdvx(4)

        !Dot product between grad(kernel) and velocity difference
        drhodt(i) = drhodt(i) + mass(j) * v_ij
        drhodt(j) = drhodt(j) + mass(i) * v_ij

    ENDDO

    end subroutine con_density

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme to calculate the particle density            !
    !                   with Staroszyck 2018 definition                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine diagnostic_density

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
            
    !!!              Variables use in the current routine:               !!!
!   A                : Sea ice concentration of particles          [out]
!   h                : Mean ice thickness of particles          [in/out]
!   rho              : Densities of particles                   [in/out]

    implicit none
        
    !!!                        Local variables:                          !!!


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    rho(:ntotal) = h(:ntotal) * rho_ice 

    end subroutine diagnostic_density


end subroutine density