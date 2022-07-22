!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                 Subroutine for generating the particles                  !
!                             information output.                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    !   mass   : Mass of particles                                      [in]
    !   ntotal : Total number of particle to use                        [in]
    !   x      : Coordinates of all particles                           [in]
    !   vx     : Velocities of all particles                            [in]
    !   rho    : Densities of particles                                 [in]
    !   p      : Pressure of particules                                 [in] 
    !   hsml   : Smoothing lengths of particles                         [in]
    !   A      : Sea ice concentration of particles                     [in]
    !   h      : Mean ice thickness of particles                        [in]
    !   epsilon_11 : XX componante of the strain rate                   [in]
    !   epsilon_12 : XY and YX componante of the strain rate            [in]
    !   epsilon_22 : YY componante of the strain rate                   [in]
    !   zeta       : Bulk viscosity                                     [in]
    !   eta        : Shear viscosity                                    [in]
    !   sigma_11 : XX componante of the stress tensor                   [in]
    !   sigma_12 : XY and YX componante of the stress tensor            [in]
    !   sigma_22 : YY componante of the stress tensor                   [in]
    !   ntotalb : Total number of boundary particles to use             [in]
    !   current_time : Current time in the integration                  [in]

    implicit none

    !!!                        Local variables:                          !!!
    !   i, d : indices for do loops
    INTEGER :: i, d

    !!!                Call modules necessary properties:                !!!
    CALL sys_variables

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!                File opening for compilling data:                 !!!
    open(1,file= cwd // trim(output_folder) // "particle_momentum.csv",    &
                               action='write',                             &
                               status='old',                               &
                               position='append')
    open(2,file= cwd // trim(output_folder) // "particle_state.csv",       &
                               action='write',                             &
                               status='old',                               &
                               position='append')
    open(3,file= cwd // trim(output_folder) // "particle_other.csv",       &
                               action='write',                             &
                               status='old',                               &
                               position='append')

    !!!                       Formating the output:                      !!!
    1001 format (1x, I7, 4(1x, A1, 1x, e15.8E3))
    1002 format (1x, I7, 11(1x, A1, 1x, e15.8E3))
    1003 format (1x, 2(I7, 1x, A1, 1x), 1(e15.8E3))
                           
    !!!                        Wrinting the output:                      !!!
    do i = 1, ntotal
        write(1,1001) i,                                                   &
                      (',', x(d, i) / 1000, d = 1, dim),                   &
                      (',', vx(d, i) * 100 , d = 1, dim)

        write(2,1002) i, ',', mass(i), ',', p(i), ',', hsml(i),&
                      ',', A(i), ',', h(i), ',', epsilon_11(i), ',',        &
                      epsilon_12(i), ',', epsilon_22(i), ',', &
                       sigma_11(i),',', sigma_12(i), ',',       &
                       sigma_22(i)
    enddo

    write(3,1003) ntotal,',', ntotalb,',', current_time/ 3600

    CLOSE(1)
    CLOSE(2)
    CLOSE(3)


end subroutine output