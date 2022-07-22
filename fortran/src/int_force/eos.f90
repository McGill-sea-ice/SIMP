!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the equation of state             !
!                         i.e. update the pressure                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine eos

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
    IF ( eos_alg .EQ. 1 ) THEN
        call hibler_eos

    ELSEIF ( eos_alg .EQ. 2 ) THEN
        call gutfraind_eos

    ELSEIF ( eos_alg .EQ. 3) THEN
        call kreyscher_eos

    ENDIF

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme for the equation of state based on           !
    !                            Hibler 1979 rheology                      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine hibler_eos

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   A      : Sea ice concentration of particles                     [in]
    !   h      : Mean ice thickness of particles                        [in]
    !   p      : Pressure of particules                                [out]

    implicit none
    
    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Pressure calculation
    p(:ntotal) = P_star * h(:ntotal) * exp(- c_eos * (1-A(:ntotal)))

    end subroutine hibler_eos

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme for the equation of state based on           !
    !                 Gutfraind 1997 article on sea ice SPH                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine gutfraind_eos

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
            
    !!!              Variables use in the current routine:               !!!
    !   rho    : Densities of particles                                 [in]
    !   p      : Pressure of particules                                [out]
    !   A      : Sea ice concentration of particles                     [in]
    !   h      : Mean ice thickness of particles                        [in]
    
    implicit none
        
    !!!                        Local variables:                          !!!
    !   rho_0   : Initial ice density                               [kg/m^3]
    !   rho_max : Free parameter                                    [kg/m^3]
    !   gamma   : Free parameter                                        [Pa]
    
    DOUBLE PRECISION :: rho_0, rho_max , gamma(ntotal)
    parameter          (rho_max = 934.,                                    &
                        rho_0   = rho_ice)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    gamma   = 1.9D5 * h(:ntotal) * exp(-20*(1-A(:ntotal)))

    !Pressure calculation
    p(:ntotal) = gamma * (rho(:ntotal) - rho_0) / (rho_max - rho(:ntotal))

    where (p(:ntotal) .LT. 0.)
        p = 0.
    endwhere


    end subroutine gutfraind_eos

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme for the equation of state based on           !
    !                 Kreyscher 2000 article on sea ice SPH                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine kreyscher_eos

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
            
    !!!              Variables use in the current routine:               !!!
    !   rho    : Densities of particles                                 [in]
    !   p      : Pressure of particules                                [out]
    !   A      : Sea ice concentration of particles                     [in]
    !   h      : Mean ice thickness of particles                        [in]
    
    implicit none
        
    !!!                        Local variables:                          !!!
    !   delta = strain rate variable use for lisibility
    !   delta_min = minimum value of delta [1/s]
    !   zero  = value to replace to avoid division by zero

    DOUBLE PRECISION :: delta(ntotal)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    delta = 0. 

    !Compute delta parameter
    delta = ( (epsilon_11(:ntotal) ** 2 + epsilon_22(:ntotal) ** 2) * (1 + &
               e ** (-2)) + 4 * e ** (-2) * epsilon_12(:ntotal) ** 2 + 2 * &
               epsilon_11(:ntotal) * epsilon_22(:ntotal) * (1 - e ** (-2)) &
               ) ** 0.5

    !Get pressure
    p(:ntotal) = zeta(:ntotal) * 2 * delta / (1 + k_t)



    end subroutine kreyscher_eos
    

end subroutine eos