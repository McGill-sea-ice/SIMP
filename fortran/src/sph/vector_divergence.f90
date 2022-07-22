!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the vector divergence             !
!                      summation for the next time step                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine vector_divergence(div_v_sum)

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables

    !!!              Variables use in the current routine:               !!!
    !   mass : Mass of particles                                        [in]
    !   wdvx     : Various products between kernel and speed            [in]

    implicit none

    !!!                        Local variables:                          !!!

    DOUBLE PRECISION :: div_v_sum(ntotal)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Modify the summation
    div_v_sum(cp_i) = div_v_sum(cp_i) - mass(cp_j) * (wdvx(1) + wdvx(4))
    div_v_sum(cp_j) = div_v_sum(cp_j) - mass(cp_i) * (wdvx(1) + wdvx(4))

end subroutine vector_divergence