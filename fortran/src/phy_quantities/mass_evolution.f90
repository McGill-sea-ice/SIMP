!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the particle masses               !
!      evolution according to the various smoothing length definition      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine mass_evolution

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
        call mass_update_from_area

    ELSEIF ( sle .EQ. 1 ) THEN
        call mass_update_from_area

    ELSEIF ( sle .EQ. 2 ) THEN 
        !Particle mass  is constant and the area it describe vary 
        RETURN

    ENDIF

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Subroutine for computing the mass evolution when the area of a    !
    !                 particle is not described by the mass                !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine mass_update_from_area

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    !   mass   : Mass of particles                                     [out]
    !   hsml   : Smoothing lengths of particles                         [in]
    !   A      : Sea ice concentration of particles                     [in]
    !   h      : Mean ice thickness of particles                        [in]

    implicit none

    !!!                        Local variables:                          !!!
    !   area = Area of ice carried by a particle                       [m^2]
    DOUBLE PRECISION :: area(ntotal)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Initialization
    area = 0.

    !Area calculation
    area(:ntotal) = pi * (alpha * hsml(:ntotal)) ** 2

    !Mass evolution calculation
    mass(:ntotal) = rho_ice * A(:ntotal) * h(:ntotal) * area(:ntotal)

    end subroutine mass_update_from_area

end subroutine mass_evolution