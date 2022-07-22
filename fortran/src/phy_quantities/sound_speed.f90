!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                Subroutine for computing the sound speed of               !
!                       the medium at each particle                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine sound_speed

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    
    !!!              Variables use in the current routine:               !!!
    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Call the desired sound speed definition
    IF ( (artificial_visc) .OR. (tse .EQ. 1) .OR. (boundary_friction) ) THEN
        call c_yang

    ELSE 
        RETURN

    ENDIF



    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !             Subprograme to calculate the sound speed according       !
    !                        to (Yang 2020) definition                     !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine c_yang

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
        
    !!!              Variables use in the current routine:               !!!
    !   Y    : Young modulus                                            [in]
    !   nu   : Poisson ratio                                            [in]
    !   c : Particle sound speed                                       [out]
    !   rho : Particle density                                          [in]

    implicit none
    
    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Calculate sound speed
    c(:ntotal) = ( Y * ( 1 - nu ) / (rho(:ntotal) * (1 + nu) *             &
                 (1 - 2 * nu))) ** 0.5

    end subroutine c_yang

end subroutine sound_speed