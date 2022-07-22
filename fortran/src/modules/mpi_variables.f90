!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                    Module for defining message passing                   !
!               interface variables used in the whole program              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mpi_variables

    !!!                         Module to include:                       !!!
    use globals
    use OMP_LIB

    !!!                       Variables definitions:                     !!!
    implicit none

    !!!                       Processes variables                        !!!
    !   process_rank : ID of the process             
    !   cluster_size: Number of threads in the current routine
    !   chunk    : Size of block when doing loops with multiple threads

    INTEGER          :: process_rank,                                      &
                        cluster_size,                                      &
                        chunk

    !$OMP THREADPRIVATE(process_rank)

    !!!                          NNPS variables                           !!!
    !   p_i_temp   : Temporary list of p_i for multiple thread
    !   p_j_temp   : Temporary list of p_j for multiple thread
    !   nip_temp   : Local thread number of interacting particle
    !   r_ij_temp  : Temporary r_ij for multiple thread                   [m]
    !   hsml_ij_temp : Temporary hsml_ij for multiple thread              [m]
    !   bp_i_temp   : Temporary list of bp_i for multiple thread
    !   bp_j_temp   : Temporary list of bp_j for multiple thread
    !   bnip_temp   : Local thread number of interacting boundary particle
    !   br_ij_temp  : Temporary br_ij for multiple thread                 [m]

    INTEGER          :: nip_temp,                                          &
                        bnip_temp                                         

    INTEGER, ALLOCATABLE :: p_i_temp(:),                                   &
                            p_j_temp(:),                                   &
                            bp_i_temp(:),                                  &
                            bp_j_temp(:)                                   

    REAL(kind=8), ALLOCATABLE :: r_ij_temp(:),                             &
                                 hsml_ij_temp(:),                          &
                                 br_ij_temp(:)

    !$OMP THREADPRIVATE(nip_temp,bnip_temp)
    !$OMP THREADPRIVATE(r_ij_temp, hsml_ij_temp, br_ij_temp)
    !$OMP THREADPRIVATE(p_i_temp, p_j_temp, bp_i_temp, bp_j_temp)

    !!!                  Loop integration variables                      !!!
    !   dx       : Distances in x and y between particles                [m]
    !   bdx      : Distances between particles and boundaries            [m]
    !   dvx      : Velocity difference in x and y between particles    [m/s]
    !   loop_index : Individual index for loop
    !   cp_i : Current particle i ID when computing interactions
    !   cp_j : Current particle j ID when computing interactions
    !   int_id : Current interaction ID
    !   dwdx_temp : Current interaction kernel derivative            [1/m^3]
    !   w_temp : Current interaction kernel                          [1/m^2]
    !   bw : Current boundary interaction kernel                     [1/m^2]
    !   bdwdx : Current boundary interaction kernel derivative       [1/m^3]
    !   wdvx : Current external product between dwdx and dvx        [1/ms^2]
    !   q : Normalize distance for the interaciton
    !   factor : Kernel function factor related to dimension

    INTEGER      :: cp_i,                                                  &
                    cp_j,                                                  &
                    loop_index,                                            &
                    int_id
           
    REAL(kind=8) :: dwdx_temp(dim),                                        &
                    w_temp,                                                &
                    bw,                                                    &
                    bdwdx(dim),                                            &
                    bdx(dim),                                              &
                    dx(dim),                                               &
                    dvx(dim),                                              &
                    wdvx(dim ** 2),                                        &
                    q,                                                     &
                    factor

    !$OMP THREADPRIVATE(int_id, cp_i, cp_j)
    !$OMP THREADPRIVATE(w_temp, dwdx_temp, dx, dvx, wdvx, q, factor, bw)
    !$OMP THREADPRIVATE(bdx, bdwdx)

end module mpi_variables