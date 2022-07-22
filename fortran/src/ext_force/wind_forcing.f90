!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the wind drag on                  !
!                              the particles.                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine wind_forcing

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    
    !!!              Variables use in the current routine:               !!!
    !   ntotal : Total number of particle to use                        [in]
    !   vx     : Velocities of all particles                            [in]
    !   h      : Mean ice thickness of particles                        [in]
    !   dvxdt : Acceleration                                        [in/out]
    !   rho_air : Density of air                                        [in]
    !   rho_ice   : Density of ice                                      [in]
    !   c_wind   : Dimensionless wind   drag coefficient                [in]

    implicit none

    !!!                        Local variables:                          !!!
    !   wind_field : Wind direction for the forcing                  [m/s^2]

    DOUBLE PRECISION :: wind_field(dim)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    wind_field = 0.

    !Only compute if wind forcing is activated in globals
    IF (wind_force .AND. current_time .GT. 0. * 3600.) THEN
        wind_field(1) = 0
        wind_field(2) = -5.3

        !Compute air drag for each particles
        !Quadratic drag law :
        dvxdt(1,:) = dvxdt(1,:) + ( c_wind * rho_air  &
                            * (wind_field(1)-vx(1,:ntotal))             & 
                            * ABS(wind_field(1)-vx(1,:ntotal)))         &
                            / (h(:ntotal)*rho_ice)

        dvxdt(2,:) = dvxdt(2,:) + ( c_wind * rho_air  &
                            * (wind_field(2)-vx(2,:ntotal))                  & 
                            * ABS(wind_field(2)-vx(2,:ntotal)))            &
                            / (h(:ntotal)*rho_ice)

    ENDIF

end subroutine wind_forcing