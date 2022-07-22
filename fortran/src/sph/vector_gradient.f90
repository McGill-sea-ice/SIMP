!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                Subroutine for computing the vector gradient              !
!                      summation for the next time step                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vector_gradient(grad_v_sum)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables

    !!!              Variables use in the current routine:               !!!
    !   mass : Mass of particles                                        [in]
    !   wdvx     : Various products between kernel and speed            [in]

    implicit none

    !!!                        Local variables:                          !!!

    DOUBLE PRECISION :: grad_v_sum(4, ntotal)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Modify the xx summation componante
    grad_v_sum(1,cp_i) = grad_v_sum(1,cp_i) - mass(cp_j) * wdvx(1) / rho(cp_j)
    grad_v_sum(1,cp_j) = grad_v_sum(1,cp_j) - mass(cp_i) * wdvx(1) / rho(cp_i)

    !Modify the yx summation componante
    grad_v_sum(2,cp_i) = grad_v_sum(2,cp_i) - mass(cp_j) * wdvx(2) / rho(cp_j)
    grad_v_sum(2,cp_j) = grad_v_sum(2,cp_j) - mass(cp_i) * wdvx(2) / rho(cp_i)

    !Modify the xy summation componante
    grad_v_sum(3,cp_i) = grad_v_sum(3,cp_i) - mass(cp_j) * wdvx(3) / rho(cp_j)
    grad_v_sum(3,cp_j) = grad_v_sum(3,cp_j) - mass(cp_i) * wdvx(3) / rho(cp_i)

    !Modify the yy summation componante
    grad_v_sum(4,cp_i) = grad_v_sum(4,cp_i) - mass(cp_j) * wdvx(4) / rho(cp_j)
    grad_v_sum(4,cp_j) = grad_v_sum(4,cp_j) - mass(cp_i) * wdvx(4) / rho(cp_i)

end subroutine vector_gradient