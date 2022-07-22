!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Subroutine for computing the quantities for                !
!                             a single time step.                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine single_step

    !!!                         Module to include:                       !!!
    use globals
    use data_variables
    use mpi_variables
    
    !!!              Variables use in the current routine:               !!!
    !   dvxdt : dvx/dt, force per unit of mass                         [out]
    !   intdvxdt  : dvx/dt of internal force                            [in]
    !   extdvxdt  : dvx/dt of external force                            [in]
    !   dadt      : Concetration rate of change                        [out]
    !   dhdt      : Thickness rate of change                           [out]
    !   drhodt    : Density rate of change                             [out]
    !   p_i  : List of first particle partner of interaction pair       [in]
    !   p_j  : List of second particle partner of interaction pair      [in]
    !   nip  : number of interacting particle                           [in]

    implicit none

    !!!                        Local variables:                          !!!
    !   div_v_sum : divergence of a vector integration
    !   xsph_sum : sum(mass * delta_vx * kernel)

    DOUBLE PRECISION :: div_v_sum(ntotal), grad_v_sum(4,ntotal)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!                      Begining of subroutine:                     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Reset array values
    dxdt(:,:ntotal)  = 0.
    dvxdt(:,:ntotal) = 0.
    dAdt(:ntotal) = 0.
    dhdt(:ntotal) = 0.
    div_v_sum = 0.
    grad_v_sum = 0.
    epsilon_11(:ntotal) = 0.
    epsilon_12(:ntotal) = 0.
    epsilon_22(:ntotal) = 0.

    !Find interacting particles
    call NNPS

    !Set parallel environment
    !$OMP PARALLEL
    chunk = NINT(REAL(nip) / (100 * OMP_GET_NUM_THREADS()))
    
    !Loop on all interactions:
    !$OMP DO SCHEDULE(DYNAMIC, chunk) &
    !$OMP REDUCTION(-:div_v_sum, grad_v_sum) 
    DO loop_index = 1, nip
        !Get the interacting particles and interaction ID
        cp_i = p_i(loop_index)
        cp_j = p_j(loop_index)
        int_id = loop_index

        !Compute their relative distance and velocity
        dx(:)   = x(:, cp_i) - x(:, cp_j)
        dvx(:)  = vx(:, cp_i) - vx(:,cp_j)

        !Compute interaction parameter
        call kernel

        !Compute speed divergence
        call vector_divergence(div_v_sum)

        !Compute speed gradiant
        call vector_gradient(grad_v_sum)

    ENDDO
    !$OMP END DO
    !$OMP END PARALLEL

    !Compute strain rate
    IF (pressure_force .OR. viscosity) THEN
        epsilon_11(:ntotal) = grad_v_sum(1,:)
        epsilon_22(:ntotal) = grad_v_sum(4,:)
        epsilon_12(:ntotal) = 0.5 * (grad_v_sum(2,:) + grad_v_sum(3,:))
    ENDIF

    !Compute internal force
    call internal_stress 

    !Do the  boundary force computation only if necessary
    IF (boundary_pressure .OR. boundary_friction) THEN
        call boundary_kernel
    ENDIF

    !Compute external force
    call wind_forcing
    call water_drag

    !Stability artifial correction if requested
    IF( xsph_regularization .OR. artificial_visc ) THEN
        call regularization
    ENDIF

    !Get prognostic variables evolutions
    dhdt(:ntotal) = -h(:ntotal) * div_v_sum / rho(:ntotal)
    dAdt(:ntotal) = -A(:ntotal) * div_v_sum / rho(:ntotal)
    dxdt(:,:ntotal) = dxdt(:,:ntotal) + vx(:,:ntotal)


end subroutine single_step