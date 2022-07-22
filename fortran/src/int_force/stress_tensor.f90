!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the vertically intgrated          !
!                                 stress tensor                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine stress_tensor

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
                        
    !!!              Variables use in the current routine:               !!!
    !   epsilon_11 = XX componante of the strain rate                   [in]
    !   epsilon_12 = XY and YX componante of the strain rate            [in]
    !   epsilon_22 = YY componante of the strain rate                   [in]
    !   zeta = Bulk viscosity                                           [in]
    !   eta  = Shear viscosity                                          [in]
    !   sigma_11 = XX componante of the stress tensor                  [out]
    !   sigma_12 = XY and YX componante of the stress tensor           [out]
    !   sigma_22 = YY componante of the stress tensor                  [out]
    !   p = Pressure                                                    [in]
    
    implicit none
            
    !!!                        Local variables:                          !!!
                
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            
    !Calculate stress tensor of each particles
    sigma_11(:ntotal) = (zeta(:ntotal) + eta(:ntotal)) * epsilon_11 +      &
                        (zeta(:ntotal) - eta(:ntotal)) * epsilon_22 -      &
                        p(:ntotal) * (1 - k_t) / 2
        
    sigma_12(:ntotal) = 2 * eta(:ntotal) * epsilon_12(:ntotal)
    
    sigma_22(:ntotal) = (zeta(:ntotal) + eta(:ntotal)) * epsilon_22 +      &
                        (zeta(:ntotal) - eta(:ntotal)) * epsilon_11 -      &
                        p(:ntotal) * (1 - k_t) / 2


end subroutine stress_tensor