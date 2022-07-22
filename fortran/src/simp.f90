!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                       SPH Ice Modelisation Project                       !
!                     Oreste Marquis, mai 2021 -  2022                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program SIMP

    !!!                         Module to include:                       !!!
    use interfaces
    use globals
    use data_variables
    use mpi_variables

    !!!              Variables use in the current routine:               !!!

    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of program:                        !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Enter the name of output folder for this run
    output_folder = 'up_arch_W6_5.3ms_E2_ND_R4_C_rho_L/run'
    output_folder = trim('/fortran/data/' // TRIM(output_folder) // '/')

    !Create the directory that will contain the run
    CALL sys_variables
    CALL EXECUTE_COMMAND_LINE("mkdir " // trim(cwd) // TRIM(output_folder))

    !!! Start the run !!!
    !Place the particles in the domain
    call initialize

    !Do the integration
    call time_integration 
    
end program SIMP