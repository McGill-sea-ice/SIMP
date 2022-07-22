!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                Module to define the interfaces of routines               !
!                        Recommended for portability                       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module interfaces
    implicit none

    !Input routine:
    interface
        subroutine initialize
            use globals
            use data_variables
            implicit none
        end subroutine initialize
    end interface

    interface
        subroutine particle_placing
            use globals
            use domain
            use data_variables
            implicit none
            INTEGER :: nppg, res_ratio, m , n , k, l, index
        end subroutine particle_placing
    end interface

    interface
        subroutine particle_reading
            use globals
            use data_variables
            implicit none
            INTEGER :: n, i
        end subroutine particle_reading
    end interface

    !Stabalization routines:
    interface
        subroutine regularization
            use globals
            use data_variables
            implicit none
        end subroutine regularization
    end interface

    interface
        subroutine XSPH
            use globals
            use data_variables
            use mpi_variables
            implicit none
            INTEGER          :: i, j, k
            DOUBLE PRECISION :: xsph_sum(dim,ntotal)        
        end subroutine XSPH
    end interface

    interface
        subroutine art_visc_HK
            use globals
            use data_variables
            use mpi_variables
            implicit none
            INTEGER          :: i, j, k
            DOUBLE PRECISION :: div_v(ntotal), pi_ij, lv(ntotal), criteria
        end subroutine art_visc_HK
    end interface

    interface
        subroutine art_visc_Mon
            use globals
            use data_variables
            use mpi_variables
            implicit none
            INTEGER          :: i, j, k
            DOUBLE PRECISION :: v_r, pi_ij, iota, s_cor
        end subroutine art_visc_Mon
    end interface

    !Boundaries routine:
    interface
        subroutine boundaries
            use globals
            use data_variables
            implicit none
            INTEGER :: i, d
        end subroutine boundaries
    end interface

    interface
        subroutine boundary_placing_1
            use globals
            use boundary_domain
            use data_variables
            implicit none
            INTEGER :: n , m , k, index ,bppg
        end subroutine boundary_placing_1
    end interface

    interface
        subroutine boundary_placing_2
            use globals
            use boundary_domain
            use data_variables
            implicit none
            INTEGER :: n , m , k, l, index ,bppg
            DOUBLE PRECISION :: spacing
        end subroutine boundary_placing_2
    end interface

    interface
        subroutine boundary_reading
            use globals
            use data_variables
            implicit none
            INTEGER :: n , j
        end subroutine boundary_reading
    end interface

    interface
        subroutine boundary_normal
            use globals
            use boundary_domain
            use data_variables
            implicit none
            INTEGER          :: NN(2,ntotalb), i, j,                       &
                                orientation(dim,ntotalb), mask(3)
            DOUBLE PRECISION :: matrix(2,2), matrix_inv(2,2), d,           &
                                d_n(2,ntotalb), tang(2), vector(2),        &
                                xb_temp(dim,ntotalb), b_n_temp(dim)
        end subroutine boundary_normal
    end interface

    interface
        subroutine boundary_force
            use globals
            use data_variables
        end subroutine boundary_force
    end interface

    interface
        subroutine bf_Monaghan
            use globals
            use data_variables
            use mpi_variables
            implicit none
            DOUBLE PRECISION :: F_b(dim), v_r, pi_aj
        end subroutine bf_Monaghan
    end interface

    interface
        subroutine bf_Wang
            use globals
            use data_variables
            use mpi_variables
            implicit none
            DOUBLE PRECISION :: F_b(dim), v(dim)
        end subroutine bf_Wang
    end interface

    interface
        subroutine bf_Marquis
            use globals
            use data_variables
            use mpi_variables
            implicit none
            DOUBLE PRECISION :: F_b(dim), v(dim)
        end subroutine bf_Marquis
    end interface

    interface
        subroutine boundary_kernel
            use globals
            use data_variables
            implicit none
        end subroutine boundary_kernel
    end interface

    interface
        subroutine b_Quadratic
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine b_Quadratic
    end interface

    interface
        subroutine b_Gaussian
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine b_Gaussian
    end interface

    interface
        subroutine b_quintic_spline
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine b_quintic_spline
    end interface

    interface
        subroutine b_Wendland_C2
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine b_Wendland_C2
    end interface

    interface
        subroutine b_Wendland_C4
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine b_Wendland_C4
    end interface

    interface
        subroutine b_Wendland_C6
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine b_Wendland_C6
    end interface

    !Output routine:
    interface
        subroutine output
            use globals
            use data_variables
            implicit none
            INTEGER :: i, d
        end subroutine output
    end interface

    !Single step routine:
    interface
        subroutine single_step
            use globals
            use data_variables
            use mpi_variables
            implicit none
            DOUBLE PRECISION :: div_v_sum(ntotal), grad_v_sum(4,ntotal)
        end subroutine single_step
    end interface

    !Time integration:
    interface
        subroutine time_integration
            use globals
            use data_variables
            implicit none    
            REAL(kind=8)        :: last_save
        end subroutine time_integration
    end interface

    interface
        subroutine euler_scheme
            use globals
            use data_variables
            implicit none
        end subroutine euler_scheme
    end interface

    interface
        subroutine predictor_corrector
            use globals
            use data_variables
            implicit none
            DOUBLE PRECISION :: x_12(dim,ntotal), vx_12(dim,ntotal),        &
                                h_12(ntotal), rho_12(ntotal),               &
                                vx_1(dim,ntotal), A_1(ntotal), h_1(ntotal), &
                                rho_1(ntotal), A_12(ntotal), x_1(dim,ntotal)
        end subroutine predictor_corrector
    end interface

    interface
        subroutine leapfrog_scheme
            use globals
            use data_variables
            implicit none
            DOUBLE PRECISION :: x_1(dim,ntotal), vx_1(dim,ntotal),          &
                                h_1(ntotal), rho_1(ntotal), A_1(ntotal)
        end subroutine leapfrog_scheme
    end interface

    !Kernel:
    interface
        subroutine kernel
            use globals
            use data_variables
            implicit none
        end subroutine kernel
    end interface
    
    interface
        subroutine Gaussian
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine Gaussian
    end interface

    interface
        subroutine quintic_spline
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine quintic_spline
    end interface

    interface
        subroutine Quadratic
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine Quadratic
    end interface

    interface
        subroutine Wendland_C2
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine Wendland_C2
    end interface

    interface
        subroutine Wendland_C4
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine Wendland_C4
    end interface

    interface
        subroutine Wendland_C6
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine Wendland_C6
    end interface

    interface
        subroutine Marquis
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine Marquis
    end interface

    interface
        subroutine Marquis2
            use globals
            use data_variables
            use mpi_variables
            implicit none
        end subroutine Marquis2
    end interface

    !Nearest neighbor search:
    interface
        subroutine NNPS
            use globals
            use data_variables
            implicit none
        end subroutine NNPS
    end interface

    interface
        subroutine direct_find
            use globals
            use data_variables
            INTEGER          :: i, j, count_np(maxn)
            DOUBLE PRECISION :: d, mhsml
        end subroutine direct_find
    end interface

    interface
        subroutine b_tree
            use globals
            use data_variables
            use tree_module
            INTEGER     :: switch(ntotal), ID(ntotal) , i , count_np(maxn)
            DOUBLE PRECISION :: domain(dim,2), max_c_hsml
            TYPE (node), pointer :: root
        end subroutine b_tree
    end interface

    interface
        subroutine bucket_search
            use globals
            use domain
            use data_variables
            use mpi_variables
            implicit none
            INTEGER  :: i, j, k, l, count_np(maxn), x_bin, y_bin,          &
                    ID(2), bucket_coor(4), max_int_pp, count_np_temp(maxn)
            INTEGER, ALLOCATABLE :: bins(:,:,:), number_bin(:,:)
            DOUBLE PRECISION :: d, mhsml, bin_size
        end subroutine bucket_search
    end interface

    interface
        subroutine interaction_stats(count_ip)
            use globals
            use data_variables
            implicit none
            INTEGER :: sumip, maxip, minip, noip, count_ip(maxn)
        end subroutine interaction_stats
    end interface

    !Wind forcing:
    interface
        subroutine wind_forcing
            use globals
            use data_variables
            implicit none
            DOUBLE PRECISION :: wind_field(dim)
        end subroutine wind_forcing
    end interface

    !Water drag:
    interface
        subroutine water_drag
            use globals
            use data_variables
            implicit none
            DOUBLE PRECISION :: current_field(dim)
        end subroutine water_drag
    end interface

    !Internal stresses:
    interface
        subroutine internal_stress
            use data_variables
            use globals
            implicit none
        end subroutine internal_stress
    end interface

    interface
        subroutine stress_tensor
            use data_variables
            use globals
            implicit none
        end subroutine stress_tensor
    end interface

    interface
        subroutine stress_1
            use globals
            use data_variables
            implicit none
            INTEGER          :: k, i, j
            DOUBLE PRECISION :: cdot_j, cdot_i
        end subroutine stress_1
    end interface

    interface
        subroutine stress_2
            use globals
            use data_variables
            implicit none
            INTEGER          :: k, i, j
            DOUBLE PRECISION :: cdot_j, cdot_i
        end subroutine stress_2
    end interface

    !Equation of state
    interface
        subroutine eos
            use globals
            use data_variables
            implicit none
        end subroutine eos
    end interface

    interface
        subroutine hibler_eos
            use globals
            use data_variables
            implicit none
        end subroutine hibler_eos
    end interface

    interface
        subroutine gutfraind_eos
            use globals
            use data_variables
            implicit none
            DOUBLE PRECISION :: rho_0, rho_max , gamma(ntotal)
            parameter          (rho_max = 934.,                             &
                                rho_0   = rho_ice)
        end subroutine gutfraind_eos
    end interface

    interface
    subroutine kreyscher_eos
        use globals
        use data_variables
        implicit none
        DOUBLE PRECISION :: delta(ntotal)
    end subroutine kreyscher_eos
    end interface

    !Rheology

    interface
        subroutine bulk_viscosity
            use globals
            use data_variables
            implicit none
        end subroutine bulk_viscosity
    end interface

    interface
        subroutine elliptical_bulk_viscosity
            use globals
            use data_variables
            implicit none
            DOUBLE PRECISION :: delta(ntotal), zero
        end subroutine elliptical_bulk_viscosity
    end interface

    interface
        subroutine MC_bulk_viscosity
            use globals
            use data_variables
            implicit none
            DOUBLE PRECISION :: epsilon_1(ntotal), epsilon_2(ntotal), zero, &
            epsilon_I(ntotal)
        end subroutine MC_bulk_viscosity
    end interface

    interface
        subroutine shear_viscosity
            use globals
            use data_variables
            implicit none
        end subroutine shear_viscosity
    end interface

    interface
        subroutine elliptical_shear_viscosity
            use globals
            use data_variables
            implicit none
        end subroutine elliptical_shear_viscosity
    end interface

    interface
    subroutine MC_shear_viscosity
        use globals
        use data_variables
        implicit none
        DOUBLE PRECISION :: epsilon_1(ntotal), epsilon_2(ntotal), zero,     &
        epsilon_II(ntotal)
    end subroutine MC_shear_viscosity
    end interface


    !Density:
    interface
       subroutine density
        use globals
        use data_variables
        implicit none
       end subroutine density
    end interface

    interface
        subroutine sum_density
            use globals
            use data_variables
            implicit none
            INTEGER          :: k, i, j
            DOUBLE PRECISION :: numerator(ntotal), denominator(ntotal)
        end subroutine sum_density
    end interface

    interface
        subroutine con_density
            use globals
            use data_variables
            implicit none
            INTEGER          :: k, i, j
            DOUBLE PRECISION :: v_ij
        end subroutine con_density
    end interface

    interface
        subroutine diagnostic_density
            use globals
            use data_variables
            implicit none
        end subroutine diagnostic_density
    end interface

    !Mass:
    interface
        subroutine mass_evolution
            use globals
            use data_variables
            implicit none
        end subroutine mass_evolution
    end interface

    interface
        subroutine mass_update_from_area
            use globals
            use data_variables
            implicit none
            DOUBLE PRECISION :: area(ntotal)
        end subroutine mass_update_from_area
    end interface

    !Smoothing length evolution
    interface
        subroutine smoothing_length
            use globals
            use data_variables
            implicit none
        end subroutine smoothing_length
    end interface

    interface
        subroutine Benz_hsml
            use globals
            use data_variables
            implicit none
            INTEGER          :: k, i, j
            DOUBLE PRECISION :: v_ij, dhsmldt(ntotal)
        end subroutine Benz_hsml
    end interface

    interface
        subroutine Marquis_hsml
            use globals
            use data_variables
            implicit none
            DOUBLE PRECISION :: m_hsml, m_rho
        end subroutine Marquis_hsml
    end interface

    !Sound speed
    interface
        subroutine sound_speed
            use globals
            use data_variables
            implicit none
        end subroutine sound_speed
    end interface

    interface
        subroutine c_yang
            use globals
            use data_variables
            implicit none
        end subroutine c_yang
    end interface

    !Time step length
    interface
        subroutine time_step
            use globals
            use data_variables
            implicit none
        end subroutine time_step
    end interface

    interface
        subroutine VP_time_step
            use globals
            use data_variables
            implicit none
        end subroutine VP_time_step
    end interface

    interface
        subroutine Monaghan_time_step
            use globals
            use data_variables
            use mpi_variables
            implicit none
            INTEGER          :: i, j, k
            DOUBLE PRECISION :: dt_force, dt_CFL, fF, acc(ntotal),          &
                                div_v(ntotal)
            parameter ( fF = 0.25)
        end subroutine Monaghan_time_step
    end interface

    !Vector operator
    interface
        subroutine vector_divergence
            use globals
            use data_variables
            use mpi_variables
            implicit none
            DOUBLE PRECISION :: div_v_sum(ntotal)
        end subroutine vector_divergence
    end interface

    interface
        subroutine vector_gradient
            use globals
            use data_variables
            use mpi_variables
            implicit none
            DOUBLE PRECISION :: grad_v_sum(4, ntotal)
        end subroutine vector_gradient
    end interface

 end module interfaces


