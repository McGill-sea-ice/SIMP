!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the shear viscosity               !
!                            on each partciles                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine shear_viscosity

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
    IF ( rheology .EQ. 1 ) THEN
        call elliptical_shear_viscosity

    ELSEIF ( rheology .EQ. 2 ) THEN
        call MC_shear_viscosity

    ENDIF

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Calculate shear viscosity for each particle according to      !
    !                         an elliptical yield curve                    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine elliptical_shear_viscosity

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
                    
    !!!              Variables use in the current routine:               !!!
    !   zeta = Bulk viscosity                                           [in]
    !   eta  = Shear viscosity                                         [out]
    !   zeta_max = Bulk viscosity                                       [in]
    !   eta_max  = Shear viscosity                                      [in]
        
    implicit none
                
    !!!                        Local variables:                          !!!
            
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    !Compute shear viscosity
    eta(:ntotal) = zeta(:ntotal) / ( e ** 2)
        
            
    end subroutine elliptical_shear_viscosity

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Calculate shear viscosity for each particle according to      !
    !                         Mohr-Coulomb yield curve                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MC_shear_viscosity

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
                    
    !!!              Variables use in the current routine:               !!!
    !   zeta = Bulk viscosity                                           [in]
    !   eta  = Shear viscosity                                         [out]
        
    implicit none
                
    !!!                        Local variables:                          !!!
    !   epsilon_1,2 = Principal strain rate
    !   delta = strain rate variable use for lisibility
    !   zero  = value to replace to avoid division by zero

    DOUBLE PRECISION :: epsilon_1(ntotal), epsilon_2(ntotal), zero,        &
                        epsilon_II(ntotal)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    epsilon_1 = 0. 
    epsilon_2 = 0.
    epsilon_II = 0.
    zero  = 1.d-18

    !Get principal strain rates
    epsilon_1 = 0.5 * ( (epsilon_11(:ntotal) + epsilon_22(:ntotal)) +      &
                sqrt((epsilon_11(:ntotal) - epsilon_22(:ntotal)) ** 2 +    &
                4 * epsilon_12(:ntotal) ** 2) )
    epsilon_2 = 0.5 * ( (epsilon_11(:ntotal) + epsilon_22(:ntotal)) -      &
                sqrt((epsilon_11(:ntotal) - epsilon_22(:ntotal)) ** 2 +    &
                4 * epsilon_12(:ntotal) ** 2) )

    !Compute the difference between principal stress
    epsilon_II = epsilon_1 - epsilon_2

    !Avoid dividing by 0
    WHERE (ABS(epsilon_II) .LT. zero)
        epsilon_II = zero
                
    ENDWHERE
    
    !Compute shear viscosity
    eta(:ntotal) = mu * ( 0.5 * p(:ntotal) * ( 1 + k_t ) - zeta(:ntotal) * &
                   (epsilon_1 + epsilon_2)) / epsilon_II
        
    !Get minimum between eta and eta_max
    WHERE (eta(:ntotal) .GT. eta_max * p(:ntotal))
        eta(:ntotal) = eta_max * p(:ntotal)
        
    ENDWHERE
            
    end subroutine MC_shear_viscosity

end subroutine shear_viscosity