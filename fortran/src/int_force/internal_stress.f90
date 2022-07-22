!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the internal stress force         !
!                               on each particles                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine internal_stress

    !!!                         Module to include:                       !!!
    use data_variables
    use globals

    !!!              Variables use in the current routine:               !!!
    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Get viscosity if required
    IF (viscosity) THEN

        !Pressure calculation
        p(:ntotal) = P_star * h(:ntotal) * exp(- c_eos * (1-A(:ntotal)))

        !Get Bulk viscosity
        call bulk_viscosity

        !Get Shear viscosity
        call shear_viscosity

    ENDIF

    !Get Pressure if required
    IF (pressure_force) THEN

        call eos

    ENDIF

    !Get stress tensor
    call stress_tensor
    
    !Get the stress force from the desired formulation
    IF (pa_sph .EQ. 1) THEN

        call stress_1

    ELSEIF (pa_sph .EQ. 2) THEN

        call stress_2

    ENDIF

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Calculate pressure force for the interacting pair of particles    !
    !        according to (p(i)/rho(i)**2+p(j)/rho(j)**2) formulation      !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine stress_1

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
            
    !!!              Variables use in the current routine:               !!!
    !   p_i  : List of first particle partner of interaction pair       [in]
    !   p_j  : List of second particle partner of interaction pair      [in]
    !   nip  : number of interacting particle                           [in]
    !   dwdx : Derivative of kernel with respect to x and y             [in]
    !   rho  : Densities of particles                                   [in]
    !   mass : Mass of particles                                        [in]
    !   sigma_11 : XX componante of the stress tensor                   [in]
    !   sigma_12 : XY and YX componante of the stress tensor            [in]
    !   sigma_22 : YY componante of the stress tensor                   [in]
    !   h      : Mean ice thickness of particles                        [in]
    !   dvxdt  : dvx/dt of internal force                              [out]
    !   rho_ice   : Density of ice                                      [in]

    implicit none
        
    !!!                        Local variables:                          !!!
    !   k : loop index
    !   i : particle i current index
    !   j : particle j current index
    !   cdot_j : dot product of particle 2 between Kernel and stress tensor
    !   cdot_i : dot product of particle 1 between Kernel and stress tensor
    
    INTEGER          :: k, i, j
    DOUBLE PRECISION :: cdot_j, cdot_i
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    !Initialization
    cdot_i = 0.
    cdot_j = 0.

    !Calculate the force from internal stress for each pair of particle
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j,cdot_j, cdot_i) &
    !$OMP REDUCTION(+:dvxdt) 
    DO k = 1, nip
    
        ! Get index
        i = p_i(k)
        j = p_j(k)
    
        !Sum the force on particle i for the x componante 
        cdot_j = (sigma_11(j) * dwdx(1,k) + sigma_12(j) * dwdx(2,k)) /     &
                 rho(j) ** 2
        cdot_i = (sigma_11(i) * dwdx(1,k) + sigma_12(i) * dwdx(2,k)) /     &
                 rho(i) ** 2
        dvxdt(1,i) = dvxdt(1,i) + rho(i) * mass(j) / rho_ice / h(i)  &
                        * (cdot_i + cdot_j)
    
        !Sum the force on particle i for the x componante 
        dvxdt(1,j) = dvxdt(1,j) + rho(j) * mass(i) / rho_ice / h(j)  &
                        * (-cdot_i - cdot_j)

        
        !Sum the force on particle i for the y componante 
        cdot_j = (sigma_12(j) * dwdx(1,k) + sigma_22(j) * dwdx(2,k)) /     &
                 rho(j) ** 2
        cdot_i = (sigma_12(i) * dwdx(1,k) + sigma_22(i) * dwdx(2,k)) /     &
                 rho(i) ** 2
        dvxdt(2,i) = dvxdt(2,i) + rho(i) * mass(j) / rho_ice / h(i)  &
                        * (cdot_i + cdot_j)
           
        !Sum the force on particle i for the y componante 
        dvxdt(2,j) = dvxdt(2,j) + rho(j) * mass(i) / rho_ice / h(j)  &
                        * (-cdot_i - cdot_j)
    
    ENDDO
    !$OMP END PARALLEL DO
    
    end subroutine stress_1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Calculate pressure force for the interacting pair of particles    !
    !         according to (p(i)+p(j)/(rho(j)*rho(i)) formulation          !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine stress_2

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
                
    !!!              Variables use in the current routine:               !!!
    !   p_i  : List of first particle partner of interaction pair       [in]
    !   p_j  : List of second particle partner of interaction pair      [in]
    !   nip  : number of interacting particle                           [in]
    !   dwdx : Derivative of kernel with respect to x and y             [in]
    !   rho  : Densities of particles                                   [in]
    !   mass : Mass of particles                                        [in]
    !   sigma_11 : XX componante of the stress tensor                   [in]
    !   sigma_12 : XY and YX componante of the stress tensor            [in]
    !   sigma_22 : YY componante of the stress tensor                   [in]
    !   h      : Mean ice thickness of particles                        [in]
    !   rho_ice   : Density of ice                                      [in]
    !   dvxdt  : dvx/dt of internal force                           [out]
        
    implicit none
            
    !!!                        Local variables:                          !!!
    !   k : loop index
    !   i : particle i current index
    !   j : particle j current index
    !   cdot_j : dot product of particle 2 between Kernel and stress tensor
    !   cdot_i : dot product of particle 1 between Kernel and stress tensor
        
    INTEGER          :: k, i, j
    DOUBLE PRECISION :: cdot_j, cdot_i
        
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        
    !Initialization
    cdot_i = 0.
    cdot_j = 0.
    
    !Calculate the force from internal stress for each pair of particle
    !$OMP PARALLEL DO SCHEDULE(STATIC) PRIVATE(i,j, cdot_j, cdot_i) &
    !$OMP REDUCTION(+:dvxdt) 
    DO k = 1, nip
        
        ! Get index
        i = p_i(k)
        j = p_j(k)
        
        !Sum the force on particle i for the x componante 
        cdot_j = (sigma_11(j) * dwdx(1,k) + sigma_12(j) * dwdx(2,k)) 
        cdot_i = (sigma_11(i) * dwdx(1,k) + sigma_12(i) * dwdx(2,k)) 
        dvxdt(1,i) = dvxdt(1,i) + mass(j) / rho_ice / h(i) / rho(j)  &
                        * (cdot_i + cdot_j)
        
        !Sum the force on particle j for the x componante 
        dvxdt(1,j) = dvxdt(1,j) - mass(i) / rho_ice / h(j) / rho(i)  &
                        * (cdot_i + cdot_j)
    
            
        !Sum the force on particle i for the y componante 
        cdot_j = (sigma_12(j) * dwdx(1,k) + sigma_22(j) * dwdx(2,k))
        cdot_i = (sigma_12(i) * dwdx(1,k) + sigma_22(i) * dwdx(2,k))
        dvxdt(2,i) = dvxdt(2,i) + mass(j) / rho_ice / h(i) / rho(j)  &
                        * (cdot_i + cdot_j)
               
        !Sum the force on particle j for the y componante 
        dvxdt(2,j) = dvxdt(2,j) - mass(i) / rho_ice / h(j) / rho(i)  &
                        * (cdot_i + cdot_j)
        
    ENDDO
    !$OMP END PARALLEL DO

    end subroutine stress_2

end subroutine internal_stress