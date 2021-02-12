! ****************************************************************
! All global variables defined here (cf F77 COMMON block).
! ****************************************************************

MODULE constants

  IMPLICIT NONE

  INTEGER, PARAMETER :: num = KIND(1.d0)
  INTEGER, PARAMETER :: dbl = KIND(1.d0)
  INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
  INTEGER, PARAMETER :: i4  = SELECTED_INT_KIND(9)  ! 4-byte 2^31 ~ 10^9
  INTEGER, PARAMETER :: i8  = SELECTED_INT_KIND(18) ! 8-byte 2^63 ~ 10^18
  REAL(num), PARAMETER :: c_tiny = TINY(1.0_num)
  REAL(num), PARAMETER :: c_largest_number = HUGE(1.0_num)
  REAL(num), PARAMETER :: c_maxexponent = MAXEXPONENT(1.0_num)
  REAL(num), PARAMETER :: c_log2 = 0.69314718055994530941723212145817657_num
  REAL(num), PARAMETER :: c_largest_exp = c_maxexponent * c_log2
  REAL(num), PARAMETER :: c_smallest_exp = (MINEXPONENT(1.0_num)-1.0_num) * c_log2

  INTEGER, PARAMETER :: c_ndims = 3

  ! File unit for deck parser diagnostic output
  INTEGER, PARAMETER :: stdout = 6
  INTEGER, PARAMETER :: du = 40
  INTEGER, PARAMETER :: lu = 41
#ifdef NO_IO
  INTEGER, PARAMETER :: io_units(1) = (/ stdout /)
#else
  INTEGER, PARAMETER :: io_units(2) = (/ stdout,du /)
#endif
  INTEGER, PARAMETER :: nio_units = SIZE(io_units)

  ! Produce time history
  LOGICAL, PARAMETER :: write_time_history = .TRUE.

  ! Boundary type codes
  INTEGER, PARAMETER :: c_bc_periodic = 1
  INTEGER, PARAMETER :: c_bc_other = 2
  INTEGER, PARAMETER :: c_bc_simple_laser = 3
  INTEGER, PARAMETER :: c_bc_simple_outflow = 4
  INTEGER, PARAMETER :: c_bc_open = 5
  INTEGER, PARAMETER :: c_bc_dump = 6
  INTEGER, PARAMETER :: c_bc_zero_gradient = 7
  INTEGER, PARAMETER :: c_bc_clamp = 8
  INTEGER, PARAMETER :: c_bc_reflect = 9
  INTEGER, PARAMETER :: c_bc_conduct = 10
  INTEGER, PARAMETER :: c_bc_thermal = 11

  ! Boundary location codes
  INTEGER, PARAMETER :: c_bd_x_min = 1
  INTEGER, PARAMETER :: c_bd_x_max = 2
  INTEGER, PARAMETER :: c_bd_y_min = 3
  INTEGER, PARAMETER :: c_bd_y_max = 4
  INTEGER, PARAMETER :: c_bd_z_min = 5
  INTEGER, PARAMETER :: c_bd_z_max = 6

  ! Error codes
  INTEGER, PARAMETER :: c_err_none = 0
  INTEGER, PARAMETER :: c_err_unknown_block = 2**0
  INTEGER, PARAMETER :: c_err_unknown_element = 2**1
  INTEGER, PARAMETER :: c_err_preset_element = 2**2
  INTEGER, PARAMETER :: c_err_preset_element_use_later = 2**3
  INTEGER, PARAMETER :: c_err_bad_value = 2**4
  INTEGER, PARAMETER :: c_err_missing_elements = 2**5
  INTEGER, PARAMETER :: c_err_terminate = 2**6
  INTEGER, PARAMETER :: c_err_required_element_not_set = 2**7
  INTEGER, PARAMETER :: c_err_pp_options_wrong = 2**8
  INTEGER, PARAMETER :: c_err_bad_array_length = 2**9
  INTEGER, PARAMETER :: c_err_other = 2**10
  INTEGER, PARAMETER :: c_err_warn_bad_value = 2**11
  INTEGER, PARAMETER :: c_err_generic_warning = 2**12
  INTEGER, PARAMETER :: c_err_generic_error = 2**13
  INTEGER, PARAMETER :: c_err_io = 2**14
  INTEGER, PARAMETER :: c_err_setup = 2**15

  INTEGER, PARAMETER :: c_ds_first = 1
  INTEGER, PARAMETER :: c_ds_last = 2

  ! domain codes
  INTEGER, PARAMETER :: c_do_full = 0
  INTEGER, PARAMETER :: c_do_decomposed = 1

  ! Load balance codes
  INTEGER, PARAMETER :: c_lb_x = 1
  INTEGER, PARAMETER :: c_lb_y = 2
  INTEGER, PARAMETER :: c_lb_z = 4
  INTEGER, PARAMETER :: c_lb_all = c_lb_x + c_lb_y + c_lb_z
  INTEGER, PARAMETER :: c_lb_auto = c_lb_all + 1

  ! Taken from http://physics.nist.gov/cuu/Constants (05/07/2012)
  REAL(num), PARAMETER :: pi = 3.141592653589793238462643383279503_num
  REAL(num), PARAMETER :: q0 = 1.602176565e-19_num ! C (+/- 3.5e-27)
  REAL(num), PARAMETER :: m0 = 9.10938291e-31_num ! kg (+/- 4e-38)
  REAL(num), PARAMETER :: c  = 2.99792458e8_num   ! m/s^2 (exact)
  REAL(num), PARAMETER :: kb = 1.3806488e-23_num  ! J/K (+/- 1.3e-29)
  REAL(num), PARAMETER :: mu0 = 4.e-7_num * pi ! N/A^2 (exact)
  ! epsilon0 = 1.0_num / mu0 / c**2 ! F/m (exact)
  REAL(num), PARAMETER :: epsilon0 = 8.854187817620389850536563031710750e-12_num
  REAL(num), PARAMETER :: h_planck = 6.62606957e-34_num ! J s (+/- 2.9e-41)
  REAL(num), PARAMETER :: ev = q0 ! J

  ! direction parameters
  INTEGER, PARAMETER :: c_dir_x = 1
  INTEGER, PARAMETER :: c_dir_y = 2
  INTEGER, PARAMETER :: c_dir_z = 3
  INTEGER, PARAMETER :: c_dir_px = c_ndims + 1
  INTEGER, PARAMETER :: c_dir_py = c_ndims + 2
  INTEGER, PARAMETER :: c_dir_pz = c_ndims + 3
  INTEGER, PARAMETER :: c_dir_en = c_ndims + 4
  INTEGER, PARAMETER :: c_dir_gamma_m1 = c_ndims + 5
  INTEGER, PARAMETER :: c_dir_xy_angle = c_ndims + 6
  INTEGER, PARAMETER :: c_dir_yz_angle = c_ndims + 7
  INTEGER, PARAMETER :: c_dir_zx_angle = c_ndims + 8

  ! Stagger types
  INTEGER, PARAMETER :: c_stagger_centre = 0
  INTEGER, PARAMETER :: c_stagger_ex = 1
  INTEGER, PARAMETER :: c_stagger_ey = 2
  INTEGER, PARAMETER :: c_stagger_ez = 4
  INTEGER, PARAMETER :: c_stagger_bx = c_stagger_ey + c_stagger_ez
  INTEGER, PARAMETER :: c_stagger_by = c_stagger_ex + c_stagger_ez
  INTEGER, PARAMETER :: c_stagger_bz = c_stagger_ex + c_stagger_ey
  INTEGER, PARAMETER :: c_stagger_jx = c_stagger_ex
  INTEGER, PARAMETER :: c_stagger_jy = c_stagger_ey
  INTEGER, PARAMETER :: c_stagger_jz = c_stagger_ez
  INTEGER, PARAMETER :: c_stagger_max = c_stagger_ex + c_stagger_bx

  ! Length of a standard string
  ! TODO Why so many different values here?
  INTEGER, PARAMETER :: string_length = 256
  INTEGER, PARAMETER :: c_max_string_length = 64
  INTEGER, PARAMETER :: c_id_length = 32

END MODULE constants



MODULE shared_data

  USE mpi
  USE constants

  IMPLICIT NONE

  !----------------------------------------------------------------------------
  ! string handling
  !----------------------------------------------------------------------------
  CHARACTER(LEN=string_length) :: blank
  TYPE string_type
    CHARACTER(string_length) :: value
  END TYPE string_type
  CHARACTER(LEN=string_length) :: extended_error_string

  !----------------------------------------------------------------------------
  ! Choice of problem set-up
  !----------------------------------------------------------------------------
  CHARACTER(LEN=c_max_string_length) :: problem

  !----------------------------------------------------------------------------
  ! Particles
  !----------------------------------------------------------------------------

  LOGICAL, PARAMETER :: particles_uniformly_distributed = .TRUE.
  
  ! Time to start the particle push - 0 by default, can be set in the control
  ! block of the deck using 'particle_tstart'.
  REAL(num) :: particle_push_start_time = 0.0_num

  ! The order for the spline interpolation used as a particle representation.
  ! png is the number of ghost cells needed by the particles
  INTEGER, PARAMETER :: sf_min = -2
  INTEGER, PARAMETER :: sf_max =  2
  INTEGER, PARAMETER :: png =  4

  ! Object representing a particle
  ! If you add or remove from this section then you *must* update the
  ! particle pack and unpack routines
  INTEGER, PARAMETER :: work_ndims=4 ! Ideally this would be per-species based.
  !As long as particle_list knows what species it contains, should be able to use a factory-type method.
  !size of particle is needed wherever a particle is allocated.
  TYPE particle
    REAL(num), DIMENSION(3) :: part_p
    REAL(num), DIMENSION(c_ndims) :: part_pos
    REAL(num) :: weight
    REAL(num) :: charge
    REAL(num) :: mass
    REAL(num) :: pvol
    REAL(num), DIMENSION(work_ndims) :: work
    TYPE(particle), POINTER :: next, prev
  END TYPE particle

  ! Object representing a collection of particles
  ! Used internally by the MPI particle transfer code
  TYPE particle_list
    TYPE(particle), POINTER :: head
    TYPE(particle), POINTER :: tail
    INTEGER(i8) :: count
    INTEGER :: id_update
    ! Pointer is safe if the particles in it are all unambiguously linked
    LOGICAL :: safe

    TYPE(particle_list), POINTER :: next, prev
  END TYPE particle_list

  ! Object representing a particle species
  TYPE particle_species
    ! Core properties
    CHARACTER(string_length) :: name
    TYPE(particle_species), POINTER :: next, prev
    INTEGER :: id
    INTEGER :: count_update_step

    REAL(num) :: charge
    REAL(num) :: mass
    REAL(num) :: weight
    INTEGER(i8) :: count
    TYPE(particle_list) :: attached_list
    LOGICAL :: immobile
    LOGICAL :: is_driftkinetic
    LOGICAL :: use_deltaf = .FALSE.
    
    ! Injection of particles
    REAL(num) :: npart_per_cell
    !TYPE(primitive_stack) :: density_function, temperature_function(3)
    !TYPE(primitive_stack) :: drift_function(3)

    ! Thermal boundaries
    REAL(num), DIMENSION(:,:,:), POINTER :: ext_temp_x_min, ext_temp_x_max
    REAL(num), DIMENSION(:,:,:), POINTER :: ext_temp_y_min, ext_temp_y_max
    REAL(num), DIMENSION(:,:,:), POINTER :: ext_temp_z_min, ext_temp_z_max

    ! Initial conditions of a species
    REAL(num), DIMENSION(:,:,:), POINTER :: density
    REAL(num), DIMENSION(:,:,:,:), POINTER :: temp
    REAL(num), DIMENSION(:,:,:,:), POINTER :: drift
    ! Initial conditions for deltaf
    REAL(num), DIMENSION(:,:,:), POINTER :: density_back
    REAL(num), DIMENSION(:,:,:,:), POINTER :: temp_back
    REAL(num), DIMENSION(:,:,:,:), POINTER :: drift_back

    REAL(num) :: density_min
    REAL(num) :: density_max
  END TYPE particle_species

  INTEGER :: deck_state

  REAL(num) :: time_start, time_stop
  INTEGER :: nstep_start, nstep_stop

  !----------------------------------------------------------------------------
  ! Core code
  !----------------------------------------------------------------------------
  INTEGER :: mpireal = MPI_DOUBLE_PRECISION
  INTEGER :: realsize

  ! ng is the number of ghost cells allocated in the arrays
  ! fng is the number of ghost cells needed by the field solver
  ! jng is the number of ghost cells needed by the current arrays
  INTEGER, PARAMETER :: ng = 3
  INTEGER, PARAMETER :: jng =  MAX(ng,png)
  INTEGER :: fng, nx, ny, nz
  INTEGER :: nx_global, ny_global, nz_global
  INTEGER(i8) :: npart_global
  INTEGER :: nsteps, n_species = -1
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: ex, ey, ez, bx, by, bz, jx, jy, jz
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: ex_back, ey_back, ez_back, bx_back, by_back, bz_back
  REAL(num), ALLOCATABLE, DIMENSION(:,:,:) :: jx_d, jy_d, jz_d !Store drift currents separately
  REAL(r4), ALLOCATABLE, DIMENSION(:,:,:) :: r4array

  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ex_x_min, ex_x_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ey_x_min, ey_x_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ez_x_min, ez_x_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bx_x_min, bx_x_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: by_x_min, by_x_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bz_x_min, bz_x_max

  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ex_y_min, ex_y_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ey_y_min, ey_y_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ez_y_min, ez_y_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bx_y_min, bx_y_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: by_y_min, by_y_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bz_y_min, bz_y_max

  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ex_z_min, ex_z_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ey_z_min, ey_z_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: ez_z_min, ez_z_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bx_z_min, bx_z_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: by_z_min, by_z_max
  REAL(num), ALLOCATABLE, DIMENSION(:,:) :: bz_z_min, bz_z_max

  TYPE(particle_species), POINTER :: species_list, io_list, io_list_data

  REAL(num), ALLOCATABLE, DIMENSION(:) :: x, y, z

  INTEGER, PARAMETER :: data_dir_max_length = 64
  CHARACTER(LEN=data_dir_max_length) :: data_dir

  LOGICAL :: use_random_seed = .FALSE.

  REAL(num) :: dt, t_end, time, dt_multiplier, dt_laser, dt_plasma_frequency
  REAL(num) :: cfl
  ! x_min is the left-hand edge of the simulation domain as specified in
  ! the input deck.
  ! x_grid_min is the location of x(1). Since the grid is cell-centred,
  ! this is usually at x_min + dx/2.
  REAL(num) :: length_x, dx, x_grid_min, x_grid_max, x_min, x_max
  REAL(num) :: x_grid_min_local, x_grid_max_local, x_min_local, x_max_local
  REAL(num) :: length_y, dy, y_grid_min, y_grid_max, y_min, y_max
  REAL(num) :: y_grid_min_local, y_grid_max_local, y_min_local, y_max_local
  REAL(num) :: length_z, dz, z_grid_min, z_grid_max, z_min, z_max
  REAL(num) :: z_grid_min_local, z_grid_max_local, z_min_local, z_max_local
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_grid_mins, x_grid_maxs
  REAL(num), DIMENSION(:), ALLOCATABLE :: y_grid_mins, y_grid_maxs
  REAL(num), DIMENSION(:), ALLOCATABLE :: z_grid_mins, z_grid_maxs

  REAL(num) :: total_ohmic_heating = 0.0_num

  LOGICAL :: need_random_state
  LOGICAL :: allow_cpu_reduce
  LOGICAL :: simplify_deck
  INTEGER, DIMENSION(2*c_ndims) :: bc_field, bc_particle
  INTEGER :: step
  LOGICAL :: drift_kinetic_species_exist = .FALSE.
  !----------------------------------------------------------------------------
  ! MPI data
  !----------------------------------------------------------------------------
  INTEGER :: coordinates(c_ndims), neighbour(-1:1, -1:1, -1:1)
  INTEGER :: x_coords, proc_x_min, proc_x_max
  INTEGER :: y_coords, proc_y_min, proc_y_max
  INTEGER :: z_coords, proc_z_min, proc_z_max
  INTEGER :: errcode, comm, tag, rank
  INTEGER :: nproc, nprocx, nprocy, nprocz
  INTEGER :: nprocdir(c_ndims)
  INTEGER :: status(MPI_STATUS_SIZE)
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nx_each_rank, ny_each_rank, nz_each_rank
  INTEGER(i8), ALLOCATABLE, DIMENSION(:) :: npart_each_rank
  LOGICAL :: x_min_boundary, x_max_boundary
  LOGICAL :: y_min_boundary, y_max_boundary
  LOGICAL :: z_min_boundary, z_max_boundary

  !----------------------------------------------------------------------------
  ! domain and loadbalancing
  !----------------------------------------------------------------------------
  LOGICAL :: use_balance
  REAL(num) :: dlb_threshold
  INTEGER(i8), PARAMETER :: npart_per_it = 1000000
  REAL(num), DIMENSION(:), ALLOCATABLE :: x_global, y_global, z_global
  REAL(num), DIMENSION(:), ALLOCATABLE :: xb_global, yb_global, zb_global
  ! The location of the processors
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_x_min, cell_x_max
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_y_min, cell_y_max
  INTEGER, DIMENSION(:), ALLOCATABLE :: cell_z_min, cell_z_max
  INTEGER, DIMENSION(:), ALLOCATABLE :: old_x_max, old_y_max, old_z_max
  INTEGER :: nx_global_min, nx_global_max
  INTEGER :: ny_global_min, ny_global_max
  INTEGER :: nz_global_min, nz_global_max
  INTEGER :: n_global_min(c_ndims), n_global_max(c_ndims)
  INTEGER :: balance_mode
  LOGICAL :: debug_mode

  !----------------------------------------------------------------------------
  ! laser boundaries
  !----------------------------------------------------------------------------
  TYPE laser_block
    ! Boundary to which laser is attached
    INTEGER :: boundary
    ! A unique id number for the laser (not used directly by EPOCH)
    ! Only used if hard coding time profiles
    INTEGER :: id
    REAL(num), DIMENSION(:,:), POINTER :: profile
    REAL(num), DIMENSION(:,:), POINTER :: phase

    LOGICAL :: use_phase_function, use_profile_function

    REAL(num) :: amp, omega, pol_angle, t_start, t_end

    TYPE(laser_block), POINTER :: next
  END TYPE laser_block

  TYPE(laser_block), POINTER :: laser_x_min, laser_x_max
  TYPE(laser_block), POINTER :: laser_y_min, laser_y_max
  TYPE(laser_block), POINTER :: laser_z_min, laser_z_max
  INTEGER :: n_laser_x_min, n_laser_x_max
  INTEGER :: n_laser_y_min, n_laser_y_max
  INTEGER :: n_laser_z_min, n_laser_z_max
  LOGICAL, DIMENSION(2*c_ndims) :: add_laser = .FALSE.

  REAL(num) :: walltime_start, real_walltime_start
  INTEGER :: stdout_frequency

  LOGICAL, DIMENSION(c_dir_x:c_dir_z,0:c_stagger_max) :: stagger
  INTEGER(i8) :: push_per_field = 5

  ! Absorption diagnostic
  REAL(num) :: laser_inject_local = 0.0_num
  REAL(num) :: laser_absorb_local = 0.0_num
  REAL(num) :: laser_injected = 0.0_num
  REAL(num) :: laser_absorbed = 0.0_num

  REAL(num) :: total_particle_energy = 0.0_num
  REAL(num) :: total_field_energy = 0.0_num

END MODULE shared_data
