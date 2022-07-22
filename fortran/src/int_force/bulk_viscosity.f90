!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                Subroutine for computing the bulk viscosity               !
!                           for each particles                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine bulk_viscosity

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
        call elliptical_bulk_viscosity

    ELSEIF ( rheology .EQ. 2 ) THEN
        call MC_bulk_viscosity
    ENDIF

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Calculate bulk viscosity for each particle according to       !
    !                        an elliptical yield curve                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine elliptical_bulk_viscosity

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
                
    !!!              Variables use in the current routine:               !!!
    !   epsilon_11 = XX componante of the strain rate                   [in]
    !   epsilon_12 = XY and YX componante of the strain rate            [in]
    !   epsilon_22 = YY componante of the strain rate                   [in]
    !   zeta = Bulk viscosity                                          [out]
    !   p = Pressure                                                    [in]
     
    implicit none
            
    !!!                        Local variables:                          !!!
    !   delta = strain rate variable use for lisibility
    !   zero  = value to replace to avoid division by zero

    DOUBLE PRECISION :: delta(ntotal), zero
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
    !Initialization
    delta = 0. 
    zero  = 1.d-18
    
    !Compute Bulk viscosity
    delta = ( (epsilon_11(:ntotal) ** 2 + epsilon_22(:ntotal) ** 2) * (1 + &
               e ** (-2)) + 4 * e ** (-2) * epsilon_12(:ntotal) ** 2 + 2 * &
               epsilon_11(:ntotal) * epsilon_22(:ntotal) * (1 - e ** (-2)) &
               ) ** 0.5

    !Avoid dividing by 0
    WHERE (ABS(delta) .LT. delta_min)
        delta = delta_min
       
    ENDWHERE

    zeta(:ntotal) = 0.5 * p(:ntotal) * (1 + k_t) / delta(:ntotal)

        
    end subroutine elliptical_bulk_viscosity

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !        Calculate bulk viscosity for each particle according to       !
    !                      Mohr-Coulomb yield curve                        !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine MC_bulk_viscosity

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
                
    !!!              Variables use in the current routine:               !!!
    !   epsilon_11 = XX componante of the strain rate                   [in]
    !   epsilon_12 = XY and YX componante of the strain rate            [in]
    !   epsilon_22 = YY componante of the strain rate                   [in]
    !   zeta = Bulk viscosity                                          [out]
    !   p = Pressure                                                    [in]
        
    implicit none
            
    !!!                        Local variables:                          !!!
    !   epsilon_1,2 = Principal strain rate
    !   epsilon_I   = Linear combination of epsilon 1 and 2
    !   zero  = value to replace to avoid division by zero

    DOUBLE PRECISION :: epsilon_1(ntotal), epsilon_2(ntotal), zero,        &
                        epsilon_I(ntotal)
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    epsilon_1 = 0. 
    epsilon_2 = 0.
    epsilon_I = 0.
    zero  = 1.d-18

    !Get principal strain rates
    epsilon_1 = 0.5 * ( (epsilon_11(:ntotal) + epsilon_22(:ntotal)) +      &
                sqrt((epsilon_11(:ntotal) - epsilon_22(:ntotal)) ** 2 +    &
                4 * epsilon_12(:ntotal) ** 2) )
    epsilon_2 = 0.5 * ( (epsilon_11(:ntotal) + epsilon_22(:ntotal)) -      &
                sqrt((epsilon_11(:ntotal) - epsilon_22(:ntotal)) ** 2 +    &
                4 * epsilon_12(:ntotal) ** 2) )

    !Compute the combination
    epsilon_I = ABS(epsilon_1 + epsilon_2)

    !Avoid dividing by 0
    WHERE (epsilon_I .LT. zero)
        epsilon_I = zero
                        
    ENDWHERE

    !Compute bulk viscosity
    zeta(:ntotal) = 0.5 * p(:ntotal) * (1 + k_t) / epsilon_I
        
    !Get minimum between zeta and zeta_max
    WHERE (zeta(:ntotal) .GT. zeta_max * p(:ntotal))
        zeta(:ntotal) = zeta_max * p(:ntotal)

    ENDWHERE

    end subroutine MC_bulk_viscosity

end subroutine bulk_viscosity
    