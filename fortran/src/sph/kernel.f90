!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the various kernel                !
!                             related quantities.                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine kernel

    !!!                         Module to include:                       !!!
    use globals
    use data_variables

    !!!              Variables use in the current routine:               !!!
    implicit none

    !!!                        Local variables:                          !!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Call the desired kernel
    IF (skf .EQ. 1) THEN
        call Quadratic

    ELSEIF (skf .EQ. 2) THEN
        call Gaussian

    ELSEIF (skf .EQ. 3) THEN
        call quintic_spline
       
    ELSEIF (skf .EQ. 4) THEN
        call Wendland_C2

    ELSEIF (skf .EQ. 5) THEN
        call Wendland_C4

    ELSEIF (skf .EQ. 6) THEN
        call Wendland_C6

    ELSEIF (skf .EQ. 7) THEN
        call Marquis

    ELSEIF (skf .EQ. 8) THEN
        call Marquis2
    
    ENDIF

    CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !            Quartic spline kernel calculation developped by           !
    !                               (Liu 2003)                             !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    subroutine quartic

    !!!                         Module to include:                       !!!
!    use globals
!    use data_variables
!    use mpi_variables

    !!!              Variables use in the current routine:               !!!
    !   nip  : Number of interacting particles                          [in]
    !   dx   : Distances in x and y between particles i and j          [out]
    !   w    : Kernel value for particle i                             [out]
    !   dwdx : Derivative of kernel with respect to x and y            [out]
    !   r_ij : distance between the interacting particles               [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]
    !   wdvx     : Various products between kernel and speed           [out]

!    implicit none
    
    !!!                        Local variables:                          !!!
    !   q : Normalize distance 
    !   factor : Kernel function factor related to dimension

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Kernel quantity calculation
!    q = r_ij(int_id) / hsml_ij(int_id) * scale_h
!    factor = 15. / (7. * pi * hsml_ij(int_id) ** 2) 

    ![0,2]:
!    w_temp = factor * (2./3. - 9./8. * q ** 2 + 19./24. * q ** 3&
!            - 5./32. * q ** 4)
!    dwdx_temp(:) = factor * ( - 9./4. * q + 19./8. * q ** 2      &
!                - 5./8. * q ** 3) / hsml_ij(int_id) * (dx(:) / r_ij(int_id))&
!                * scale_h

    !Store kernel values for future use
!    w(int_id) =  w_temp
!    dwdx(:,int_id) = dwdx_temp

    !Get Kernel product with particle velocity difference
!    wdvx(1) = dvx(1) * dwdx_temp(1)
!    wdvx(2) = dvx(2) * dwdx_temp(1)
!    wdvx(3) = dvx(1) * dwdx_temp(2)
!    wdvx(4) = dvx(2) * dwdx_temp(2)

!   end subroutine quartic

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !          Original Gaussian kernel calculation developped by          !
    !                      (Gingold & Monaghan 1981)                       !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Gaussian

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   nip  : Number of interacting particles                          [in]
    !   dx   : Distances in x and y between particles i and j          [out]
    !   w    : Kernel value for particle i                             [out]
    !   dwdx : Derivative of kernel with respect to x and y            [out]
    !   r_ij : distance between the interacting particles               [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]
    !   wdvx     : Various products between kernel and speed           [out]
    !   kernel_0 : Value of the kernel at the given particle location  [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   q : Normalize distance 
    !   factor : Kernel function factor related to dimension

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Kernel quantity calculation
    q = r_ij(int_id) / hsml_ij(int_id) * scale_h
    factor = 1. / (pi) * (scale_h / hsml_ij(int_id)) ** 2

    !Get the kernel at a particle when it doesn't interact for self density
    kernel_zero(cp_i) = factor 

    ![0,3]:
    w_temp = factor * exp(-q ** 2)
    dwdx_temp(:) = factor * (- 2 * q * exp(-q ** 2)) /          &
                hsml_ij(int_id) * (dx(:) / r_ij(int_id)) * scale_h

    !Store kernel values for future use
    w(int_id) =  w_temp
    dwdx(:,int_id) = dwdx_temp

    !Get Kernel product with particle velocity difference
    wdvx(1) = dvx(1) * dwdx_temp(1)
    wdvx(2) = dvx(2) * dwdx_temp(1)
    wdvx(3) = dvx(1) * dwdx_temp(2)
    wdvx(4) = dvx(2) * dwdx_temp(2)

    end subroutine Gaussian


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !            Quintic spline kernel calculation developped by           !
    !                              (Morris 1997)                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine quintic_spline

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   nip  : Number of interacting particles                          [in]
    !   dx   : Distances in x and y between particles i and j          [out]
    !   w    : Kernel value for particle i                             [out]
    !   dwdx : Derivative of kernel with respect to x and y            [out]
    !   r_ij : distance between the interacting particles               [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]
    !   wdvx     : Various products between kernel and speed           [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   q : Normalize distance 
    !   factor : Kernel function factor related to dimension


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Kernel quantity calculation
    q = r_ij(int_id) / hsml_ij(int_id) * scale_h
    factor = 7.e0 / (478.e0 * pi) * (scale_h / hsml_ij(int_id)) ** 2

    !Get the kernel at a particle when it doesn't interact for self density
    kernel_zero(cp_i) = factor * ( (3) ** 5 - 6 * (2) ** 5 + 15)

    ![0,1]:
    IF ((q .GE. 0.) .AND. (q .LE. 1.)) THEN
        w_temp = factor * ( (3 - q) ** 5 - 6 * (2 - q) ** 5 +   &
                15 * (1 - q) ** 5 )
        dwdx_temp(:) = factor * ( -5 * (3 - q ) ** 4 + 30 *         &
                        ( 2 - q ) ** 4 - 75 * ( 1 - q ) ** 4 )    &
                        / hsml_ij(int_id) * ( dx(:) / r_ij(int_id) ) * scale_h
    ENDIF

    !]1,2]:
    if ((q .GT. 1.) .AND. (q .LE. 2.)) THEN
        w_temp = factor * ( (3 - q) ** 5 - 6 * (2 - q) ** 5 )
        dwdx_temp(:) = factor * ( -5 * (3 - q ) ** 4 + 30 *         &
                            ( 2 - q ) ** 4 ) / hsml_ij(int_id) *          &
                            ( dx(:) / r_ij(int_id) ) * scale_h
    ENDIF

    !]2,3]:
    if ((q .GT. 2.) .AND. (q .LE. 3.)) THEN
        w_temp = factor * ( (3 - q) ** 5 )
        dwdx_temp(:) = factor * ( -5 * (3 - q ) ** 4 )               &
                        / hsml_ij(int_id) * (dx(:) / r_ij(int_id)) * scale_h
    ENDIF

    !Store kernel values for future use
    w(int_id) =  w_temp
    dwdx(:,int_id) = dwdx_temp

    !Get Kernel product with particle velocity difference
    wdvx(1) = dvx(1) * dwdx_temp(1)
    wdvx(2) = dvx(2) * dwdx_temp(1)
    wdvx(3) = dvx(1) * dwdx_temp(2)
    wdvx(4) = dvx(2) * dwdx_temp(2)

    end subroutine quintic_spline


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !               Quadratic kernel calculation developped by             !
    !                            (Johnson, 1996)                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Quadratic

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   nip  : Number of interacting particles                          [in]
    !   dx   : Distances in x and y between particles i and j          [out]
    !   w    : Kernel value for particle i                             [out]
    !   dwdx : Derivative of kernel with respect to x and y            [out]
    !   r_ij : distance between the interacting particles               [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]
    !   wdvx     : Various products between kernel and speed           [out]
    !   kernel_0 : Value of the kernel at the given particle location  [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   q : Normalize distance 
    !   factor : Kernel function factor related to dimension

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Kernel quantity calculation
    q = r_ij(int_id) / hsml_ij(int_id) * scale_h
    factor = 2. / (pi) * ( scale_h / hsml_ij(int_id)) ** 2

    !Get the kernel at a particle when it doesn't interact for self density
    kernel_zero(cp_i) = factor *  3./4. 

    ![0,3]:
    w_temp = factor * ( 3./16. * q ** 2 - 3./4. * q + 3./4. ) 
    dwdx_temp(:) = factor * ( 3./8. * q - 3./4.)  /          &
                hsml_ij(int_id) * (dx(:) / r_ij(int_id)) * scale_h

    !Store kernel values for future use
    w(int_id) =  w_temp
    dwdx(:,int_id) = dwdx_temp

    !Get Kernel product with particle velocity difference
    wdvx(1) = dvx(1) * dwdx_temp(1)
    wdvx(2) = dvx(2) * dwdx_temp(1)
    wdvx(3) = dvx(1) * dwdx_temp(2)
    wdvx(4) = dvx(2) * dwdx_temp(2)

    end subroutine Quadratic

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !          Original Cubic kernel calculation developped by            !
    !                             (Morris 1994)                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    subroutine Cubic_spline

    !!!                         Module to include:                       !!!
!    use globals
!    use data_variables
!    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   nip  : Number of interacting particles                          [in]
    !   dx   : Distances in x and y between particles i and j          [out]
    !   w    : Kernel value for particle i                             [out]
    !   dwdx : Derivative of kernel with respect to x and y            [out]
    !   r_ij : distance between the interacting particles               [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]
    !   wdvx     : Various products between kernel and speed           [out]
    !   kernel_0 : Value of the kernel at the given particle location  [out]

!    implicit none
    
    !!!                        Local variables:                          !!!
    !   q : Normalize distance 
    !   factor : Kernel function factor related to dimension

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Kernel quantity calculation
!    q = r_ij(int_id) / hsml_ij(int_id) * scale_h
!    factor = 15.e0 / (7.e0 * pi * hsml_ij(int_id) ** 2) 

    ![0,1]:
!    IF ((q .GE. 0.) .AND. (q .LE. 1.)) THEN
!        w_temp = factor * ( 2./3. - q ** 2 + 0.5 * q ** 3 )
!        dwdx_temp(:) = factor * ( -2 * q + 1.5 * q ** 2 )     &
!                    / hsml_ij(int_id) * ( dx(:) / r_ij(int_id) ) * scale_h
!    ENDIF

    !]1,2]:
!    if ((q .GT. 1.) .AND. (q .LE. 2.)) THEN
!        w_temp = factor * ( 1./6. * (2 - q) ** 3 )
!        dwdx_temp(:) = factor * ( -0.5 * (2 - q) ** 2 ) / hsml_ij(int_id) * &
!                            ( dx(:) / r_ij(int_id) ) * scale_h
!    ENDIF

    !Store kernel values for future use
!    w(int_id) =  w_temp
!    dwdx(:,int_id) = dwdx_temp

    !Get Kernel product with particle velocity difference
!    wdvx(1) = dvx(1) * dwdx_temp(1)
!    wdvx(2) = dvx(2) * dwdx_temp(1)
!    wdvx(3) = dvx(1) * dwdx_temp(2)
!    wdvx(4) = dvx(2) * dwdx_temp(2)

!    end subroutine Cubic_spline

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                     Wendland C2 kernel calculation                   !
    !                            (Wendland 1995)                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Wendland_C2

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   nip  : Number of interacting particles                          [in]
    !   dx   : Distances in x and y between particles i and j          [out]
    !   w    : Kernel value for particle i                             [out]
    !   dwdx : Derivative of kernel with respect to x and y            [out]
    !   r_ij : distance between the interacting particles               [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]
    !   wdvx     : Various products between kernel and speed           [out]
    !   kernel_0 : Value of the kernel at the given particle location  [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   q : Normalize distance 
    !   factor : Kernel function factor related to dimension

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Kernel quantity calculation
    q = r_ij(int_id) / hsml_ij(int_id) * scale_h
    factor = 3. / (2. * pi) * ( scale_h / hsml_ij(int_id)) ** 2

    !Get the kernel at a particle when it doesn't interact for self density
    kernel_zero(cp_i) = factor
    
    ![0,1]:
    w_temp = factor * ( (1. - q) ** 4 * (1 + 4 * q ) )
    dwdx_temp(:) = factor * ( 20. * ( q ** 4 - 3 * q ** 3 + 3 * q ** 2 - &
                q) ) / hsml_ij(int_id) * ( dx(:) / r_ij(int_id) ) * scale_h

    !Store kernel values for future use
    w(int_id) =  w_temp
    dwdx(:,int_id) = dwdx_temp

    !Get Kernel product with particle velocity difference
    wdvx(1) = dvx(1) * dwdx_temp(1)
    wdvx(2) = dvx(2) * dwdx_temp(1)
    wdvx(3) = dvx(1) * dwdx_temp(2)
    wdvx(4) = dvx(2) * dwdx_temp(2)

    end subroutine Wendland_C2

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                     Wendland C4 kernel calculation                   !
    !                            (Wendland 1995)                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Wendland_C4

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   nip  : Number of interacting particles                          [in]
    !   dx   : Distances in x and y between particles i and j          [out]
    !   w    : Kernel value for particle i                             [out]
    !   dwdx : Derivative of kernel with respect to x and y            [out]
    !   r_ij : distance between the interacting particles               [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]
    !   wdvx     : Various products between kernel and speed           [out]
    !   kernel_0 : Value of the kernel at the given particle location  [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   q : Normalize distance 
    !   factor : Kernel function factor related to dimension

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Kernel quantity calculation
    q = r_ij(int_id) / hsml_ij(int_id) * scale_h
    factor = 3. / pi * (scale_h / hsml_ij(int_id)) ** 2

    !Get the kernel at a particle when it doesn't interact for self density
    kernel_zero(cp_i) = factor * 3

    ![0,1]:
    w_temp = factor * ((1-q)**6 * (35*q**2 + 18*q + 3))
    dwdx_temp(:) = factor * (-56*q * (5*q + 1) * (1 - q)**5) / &
                   hsml_ij(int_id) * ( dx(:) / r_ij(int_id) ) * scale_h

    !Store kernel values for future use
    w(int_id) =  w_temp
    dwdx(:,int_id) = dwdx_temp

    !Get Kernel product with particle velocity difference
    wdvx(1) = dvx(1) * dwdx_temp(1)
    wdvx(2) = dvx(2) * dwdx_temp(1)
    wdvx(3) = dvx(1) * dwdx_temp(2)
    wdvx(4) = dvx(2) * dwdx_temp(2)

    end subroutine Wendland_C4

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                     Wendland C6 kernel calculation                   !
    !                            (Wendland 1995)                           !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Wendland_C6

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   nip  : Number of interacting particles                          [in]
    !   dx   : Distances in x and y between particles i and j          [out]
    !   w    : Kernel value for particle i                             [out]
    !   dwdx : Derivative of kernel with respect to x and y            [out]
    !   r_ij : distance between the interacting particles               [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]
    !   wdvx     : Various products between kernel and speed           [out]
    !   kernel_0 : Value of the kernel at the given particle location  [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   q : Normalize distance 
    !   factor : Kernel function factor related to dimension

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Kernel quantity calculation
    q = r_ij(int_id) / hsml_ij(int_id) * scale_h
    factor = 78. / (7. * pi) * ( scale_h / hsml_ij(int_id)) ** 2

    !Get the kernel at a particle when it doesn't interact for self density
    kernel_zero(cp_i) = factor 
    
    ![0,1]:
    w_temp = factor * (1 - q)**8 * (32 * q**3 + 25 * q**2 + 8 * q + 1)
    dwdx_temp(:) = factor * ( -22.* q * (1 - q)**7 * (16 * q**2 + 7.* q &
         + 1)) / hsml_ij(int_id) * ( dx(:) / r_ij(int_id) ) * scale_h

    !Store kernel values for future use
    w(int_id) =  w_temp
    dwdx(:,int_id) = dwdx_temp

    !Get Kernel product with particle velocity difference
    wdvx(1) = dvx(1) * dwdx_temp(1)
    wdvx(2) = dvx(2) * dwdx_temp(1)
    wdvx(3) = dvx(1) * dwdx_temp(2)
    wdvx(4) = dvx(2) * dwdx_temp(2)

    end subroutine Wendland_C6

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                           My kernel attempt                          !
    !                            (Marquis 2022)                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Marquis

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   nip  : Number of interacting particles                          [in]
    !   dx   : Distances in x and y between particles i and j          [out]
    !   w    : Kernel value for particle i                             [out]
    !   dwdx : Derivative of kernel with respect to x and y            [out]
    !   r_ij : distance between the interacting particles               [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]
    !   wdvx     : Various products between kernel and speed           [out]
    !   kernel_0 : Value of the kernel at the given particle location  [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   q : Normalize distance 
    !   factor : Kernel function factor related to dimension

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Kernel quantity calculation
    q = r_ij(int_id) / hsml_ij(int_id) * scale_h
    factor = 1. / ( pi ) * ( scale_h / hsml_ij(int_id)) ** 2

    !Get the kernel at a particle when it doesn't interact for self density
    kernel_zero(cp_i) = factor * (0.5 - 0.5 * tanh(0.))

    ![0,2]:
    w_temp = factor * (0.5 - 0.5 * tanh(300 * (q - 1)))
    dwdx_temp(:) = factor * ( -150 * cosh(300 * (-1 + q)) ** (-2))             &
                   / hsml_ij(int_id) * ( dx(:) / r_ij(int_id) ) * scale_h

    !Store kernel values for future use
    w(int_id) =  w_temp
    dwdx(:,int_id) = dwdx_temp

    !Get Kernel product with particle velocity difference
    wdvx(1) = dvx(1) * dwdx_temp(1)
    wdvx(2) = dvx(2) * dwdx_temp(1)
    wdvx(3) = dvx(1) * dwdx_temp(2)
    wdvx(4) = dvx(2) * dwdx_temp(2)

    end subroutine Marquis

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                      My second kernel attempt                        !
    !                            (Marquis 2022)                            !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine Marquis2

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
        
    !!!              Variables use in the current routine:               !!!
    !   nip  : Number of interacting particles                          [in]
    !   dx   : Distances in x and y between particles i and j          [out]
    !   w    : Kernel value for particle i                             [out]
    !   dwdx : Derivative of kernel with respect to x and y            [out]
    !   r_ij : distance between the interacting particles               [in]
    !   hsml_ij : mean smoothing length between particle i and j        [in]
    !   wdvx     : Various products between kernel and speed           [out]
    !   kernel_0 : Value of the kernel at the given particle location  [out]

    implicit none
    
    !!!                        Local variables:                          !!!
    !   C1, C2 : constante to shape the curve *factor depends on those
    DOUBLE PRECISION :: C1, C2
    PARAMETER (C1 = 10., C2 = 50.)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Kernel quantity calculation
    q = r_ij(int_id) / hsml_ij(int_id) * scale_h
    factor = 1. / (34.1558 * pi * hsml_ij(int_id) ** 2) 

    !Get the kernel at a particle when it doesn't interact for self density
    kernel_zero(cp_i) = factor * C2 * (log(cosh(C1 ))/C1 + 1 + log(2.)/C1)

    ![0,2]:
    w_temp = factor * C2 * (log(cosh(C1 * (1 - q)))/C1 - q + 1 + log(2.)/C1)
    dwdx_temp(:) = factor * ( -C2 * tanh(C1 * (1 - q)) - 1)             &
                    / hsml_ij(int_id) * ( dx(:) / r_ij(int_id) ) * scale_h

    !Store kernel values for future use
    w(int_id) =  w_temp
    dwdx(:,int_id) = dwdx_temp

    !Get Kernel product with particle velocity difference
    wdvx(1) = dvx(1) * dwdx_temp(1)
    wdvx(2) = dvx(2) * dwdx_temp(1)
    wdvx(3) = dvx(1) * dwdx_temp(2)
    wdvx(4) = dvx(2) * dwdx_temp(2)

    end subroutine Marquis2

end subroutine kernel
