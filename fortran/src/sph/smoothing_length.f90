!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             Subroutine to compute the evolution of the kernel            !
!                               smoothing length                           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine smoothing_length

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
    IF ( sle .EQ. 0 ) THEN
        RETURN

    ELSEIF ( sle .EQ. 1 ) THEN
        call Benz_hsml

    ELSEIF ( sle .EQ. 2 ) THEN
        call Marquis_hsml

    ENDIF


    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme to compute the change in smoothing           !
    !                 length according to inverse Benz 1990                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Benz_hsml

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   rho  : Densities of particles                                   [in]
    !   hsml  : Particle smoothing length                              [out]
    !   wdvx     : Various products between kernel and speed            [in]
    !   drhodt    : Density rate of change                              [in]
    !   p_i     : List of first particle partner of interaction pair    [in]
    !   p_j     : List of second particle partner of interaction pair   [in]
    !   nip     : number of interacting particle                        [in]
    !   dt        : Time step of the integration                        [in]
    !   mass   : Mass of particles                                      [in]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   k : loop index
    !   i : particle i current index
    !   j : particle j current index
    !   v_ij : Lisibility variable
    !   dhsmldt : Rate of change of smoothing length                   [m/s]

    INTEGER          :: k, i, j
    DOUBLE PRECISION :: v_ij, dhsmldt(ntotal)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Compute the change in density if not already done
    IF (den_alg .NE. 2) THEN
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
    ENDIF

    !Compute the change in smoothing length
    dhsmldt = -hsml(:ntotal) * drhodt(:ntotal) / dim / rho(:ntotal)
    hsml(:ntotal) = hsml(:ntotal) + dhsmldt * dt

    !Ensure the smoothing length doesn't become negative
    WHERE (hsml .LT. 0)
        hsml = hsml - dhsmldt * dt

    ENDWHERE

    end subroutine Benz_hsml

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme to compute the change in smoothing           !
    !        developped by myself but by keeping particle mass constant    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Marquis_hsml

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
            
    !!!              Variables use in the current routine:               !!!
    !   hsml  : Particle smoothing length                              [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Compute 
    hsml(:ntotal) = sqrt(mass(:ntotal) / rho(:ntotal)) * hsml_multiplicity
    
    !Cap smoothing length to limit interaction
    WHERE (hsml(:ntotal) .GT. max_hsml)
        hsml(:ntotal) = max_hsml

    ENDWHERE

    
    end subroutine Marquis_hsml

end subroutine smoothing_length