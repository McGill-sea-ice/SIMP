!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the water drag on                 !
!                              the particles.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine water_drag

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    
    !!!              Variables use in the current routine:               !!!
    !   ntotal : Total number of particle to use                        [in]
    !   vx     : Velocities of all particles                            [in]
    !   h      : Mean ice thickness of particles                        [in]
    !   dvxdt : Acceleration                                        [in/out]
    !   rho_water : Density of water                                    [in]
    !   rho_ice   : Density of ice                                      [in]
    !   c_water   : Dimensionless water drag coefficient                [in]

    implicit none

    !!!                        Local variables:                          !!!
    !   current_field : Ocean current direction for the forcing      [m/s^2]

    DOUBLE PRECISION :: current_field(dim)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    current_field = 0.

    !Only compute if water forcing is activated in globals
    IF (water_force .AND. current_time .GT. 0. * 3600.) THEN
        current_field(1) = 0
        current_field(2) = 0
        
        !Quadratic drag law :
        !Compute water drag for each particles

        dvxdt(1,:) = dvxdt(1,:) + ( c_water * rho_water   &
                        * (current_field(1)-vx(1,:))                   & 
                        * ABS(current_field(1)-vx(1,:)) )              &
                        / (h(:)*rho_ice)

        dvxdt(2,:) = dvxdt(2,:) + ( c_water * rho_water   &
                        * (current_field(2)-vx(2,:))                   & 
                        * ABS(current_field(2)-vx(2,:)) )              &
                        / (h(:)*rho_ice)


    ENDIF
end subroutine water_drag