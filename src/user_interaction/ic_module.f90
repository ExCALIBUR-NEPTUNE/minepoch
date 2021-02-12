MODULE ic_module

  USE shared_data
  USE helper
  USE setup

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: manual_load, custom_problem_setup

CONTAINS

  SUBROUTINE manual_load

  END SUBROUTINE manual_load



  SUBROUTINE custom_problem_setup(deck_state, problem)

    INTEGER, INTENT(IN) :: deck_state
    CHARACTER(LEN=c_max_string_length), INTENT(IN) :: problem

    SELECT CASE (TRIM(problem))
    CASE('two_stream')
      CALL two_stream_setup(deck_state)
    CASE('one_stream')
      CALL one_stream_setup(deck_state)
    CASE('drift_kin_default')
      CALL dk_setup(deck_state)
    CASE default
      PRINT*, 'Unrecognised set-up: ', TRIM(problem)
      CALL MPI_ABORT(MPI_COMM_WORLD, c_err_setup, errcode)
    END SELECT

  END SUBROUTINE custom_problem_setup



  SUBROUTINE two_stream_setup(deck_state)

    INTEGER, INTENT(IN) :: deck_state
    REAL(num), PARAMETER :: v_drift = 0.2_num * c
    REAL(num), PARAMETER :: v_therm = 0.01_num * c
    REAL(num), PARAMETER :: v_pert = 0.1_num * v_therm
    REAL(num), PARAMETER :: n0 = 8e11
    INTEGER, PARAMETER :: ppc = 16
    REAL(num) :: gamma_drift, temp_x, omega
    INTEGER :: ix
    TYPE(particle_species), POINTER :: current_species

    IF (deck_state == c_ds_first) THEN
      ! Set control variables here
      nx_global = 64
      ny_global = 4
      nz_global = 4
      x_min = 0.0_num
      x_max = 2.0_num * pi
      y_min = x_min
      ! Could do (x_max * ny_global) / nx_global, but be wary of compilers
      ! which don't obey precedence implied by parentheses by default
      ! (e.g. Intel)
      y_max = x_max * REAL(ny_global, num) / REAL(nx_global, num)
      z_min = x_min
      z_max = x_max * REAL(nz_global, num) / REAL(nx_global, num)

      ! Plasma frequency
      omega = SQRT(n0 * q0 * q0 / epsilon0 / m0)
      t_end = 30.0_num / omega
      stdout_frequency = 10

      ! Need to set-up species here
      NULLIFY(current_species)
      CALL setup_species(current_species, 'Right')

      ! mass -- MANDATORY
      current_species%mass = 1.0_num * m0

      ! charge -- MANDATORY
      current_species%charge = -1.0_num * q0

      ! npart_per_cell
      current_species%npart_per_cell = ppc

      ! MANDATORY
      NULLIFY(current_species)
      CALL setup_species(current_species, 'Left')

      ! mass -- MANDATORY
      current_species%mass = 1.0_num * m0

      ! charge -- MANDATORY
      current_species%charge = -1.0_num * q0

      ! npart_per_cell
      current_species%npart_per_cell = ppc

      RETURN
    END IF

    ! Calculate gamma_drift
    ! Strictly should be function of x, but vpert << vdrift
    gamma_drift = 1.0_num / SQRT(1.0_num - (v_drift / c)**2)

    ! Calculate (1 DoF) temperature
    temp_x = v_therm**2 * m0 / kb

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Species

    ! MANDATORY (Is the NULLIFY mandatory?)
    NULLIFY(current_species)
    CALL setup_species(current_species, 'Right')

    ! density
    current_species%density = n0

    ! drift_x
    ! Add on perturbation to seed instability
    DO ix = 1-ng, nx+ng
      current_species%drift(ix,:,:,1) = gamma_drift * m0 &
          * (v_drift + v_pert * SIN(3.0_num * x(ix)))
    END DO

    ! temp_x
    current_species%temp(:,:,:,1) = temp_x

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! MANDATORY
    NULLIFY(current_species)
    CALL setup_species(current_species, 'Left')

    ! density
    current_species%density = n0

    ! drift_x
    ! Add on perturbation to seed instability
    DO ix = 1-ng, nx+ng
      current_species%drift(ix,:,:,1) = gamma_drift * m0 &
          * (-v_drift + v_pert * SIN(3.0_num * x(ix)))
    END DO

    ! temp_x
    current_species%temp(:,:,:,1) = temp_x

  END SUBROUTINE two_stream_setup





  SUBROUTINE one_stream_setup(deck_state)

    INTEGER, INTENT(IN) :: deck_state
    REAL(num), PARAMETER :: v_drift = 0.1_num * c
    REAL(num), PARAMETER :: v_therm = 0.01_num * c
    REAL(num), PARAMETER :: v_pert = 10.0_num * v_therm
    REAL(num), PARAMETER :: n0 = 8e11
    INTEGER, PARAMETER :: ppc = 160
    REAL(num) :: gamma_drift, temp_x, omega
    INTEGER :: ix
    TYPE(particle_species), POINTER :: current_species

    ! Plasma frequency
    omega = SQRT(n0 * q0 * q0 / epsilon0 / m0)

    ! Control

    nx_global = 64
    ny_global = 4
    nz_global = 4
    x_min = 0.0_num
    x_max = 2.0_num * pi
    y_min = x_min
    ! Could do (x_max * ny_global) / nx_global, but be wary of compilers
    ! which don't obey precedence implied by parentheses by default
    ! (e.g. Intel)
    y_max = x_max * REAL(ny_global, num) / REAL(nx_global, num)
    z_min = x_min
    z_max = x_max * REAL(nz_global, num) / REAL(nx_global, num)
    t_end = 0.4_num / omega
    stdout_frequency = 10

    ! Calculate gamma_drift
    ! Strictly should be function of x, but vpert << vdrift
    gamma_drift = 1.0_num / SQRT(1.0_num - (v_drift / c)**2)

    ! Calculate (1 DoF) temperature
    temp_x = v_therm**2 * m0 / kb
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Species

    ! MANDATORY
    NULLIFY(current_species)
    CALL setup_species(current_species, 'Right')

    ! mass -- MANDATORY
    current_species%mass = 1.0_num * m0

    ! charge -- MANDATORY
    current_species%charge = -1.0_num * q0

    current_species%use_deltaf = .true.

    ! npart_per_cell
    current_species%npart_per_cell = ppc

    IF (deck_state /= c_ds_first) THEN
      ! density
      current_species%density = n0

      ! drift_x
      ! Add on perturbation to seed instability
      DO ix = 1-ng, nx+ng
        current_species%drift(ix,:,:,1) = gamma_drift * m0 &
            * (v_drift + v_pert * SIN(3.0_num * x(ix)))
        current_species%density(ix,:,:) = n0 * (1.0 + 0.3*SIN(2.0_num * x(ix)) ) 
      END DO

      ! temp_x
      current_species%temp(:,:,:,1) = temp_x
    ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

END SUBROUTINE one_stream_setup








  SUBROUTINE dk_setup(deck_state)

    INTEGER, INTENT(IN) :: deck_state
    REAL(num), PARAMETER :: v_drift = 0.2_num * c
    REAL(num), PARAMETER :: v_therm = 0.3_num * c ! Since nonrelativistic anyway.
    REAL(num), PARAMETER :: v_pert = 0.1_num * v_therm
    ! Density indirectly sets timestep, take to be low (not self-consistent anyway).
    REAL(num), PARAMETER :: n0 = 8e8             
    INTEGER, PARAMETER :: ppc = 1
    REAL(num) :: gamma_drift, temp_x, omega
    TYPE(particle_species), POINTER :: current_species

    ! Plasma frequency
    omega = SQRT(n0 * q0 * q0 / epsilon0 / m0)

    ! Control

    nx_global = 128
    ny_global = 4
    nz_global = 4
    x_min = 0.0_num
    x_max = 2.0_num * pi
    y_min = x_min
    ! Could do (x_max * ny_global) / nx_global, but be wary of compilers
    ! which don't obey precedence implied by parentheses by default
    ! (e.g. Intel)
    y_max = x_max * REAL(ny_global, num) / REAL(nx_global, num)
    z_min = x_min
    z_max = x_max * REAL(nz_global, num) / REAL(nx_global, num)
    t_end = 10*x_max/v_therm ! For drift kinetics, transit frequensy is a sensible timescale
    stdout_frequency = 10

    ! Calculate gamma_drift
    ! Strictly should be function of x, but vpert << vdrift
    gamma_drift = 1.0_num / SQRT(1.0_num - (v_drift / c)**2)

    ! Calculate (1 DoF) temperature
    temp_x = v_therm**2 * m0 / kb
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Species

    ! MANDATORY
    NULLIFY(current_species)
    CALL setup_species(current_species, 'electrons_dk')

    ! mass -- MANDATORY
    current_species%mass = 1.0_num * m0

    ! charge -- MANDATORY
    current_species%charge = -1.0_num * q0

    ! npart_per_cell
    current_species%npart_per_cell = ppc
    current_species%is_driftkinetic = .true.

    IF (deck_state /= c_ds_first) THEN
      ! density
      current_species%density = n0
 
      ! temp_x
      current_species%temp(:,:,:,1) = temp_x
    ENDIF

 

  END SUBROUTINE dk_setup











!  SUBROUTINE example_problem_setup(deck_state)
!
!    INTEGER, INTENT(IN) :: deck_state
!
!    ! Control
!
!    !allow_cpu_reduce = .TRUE.
!    !timer_collect = .TRUE.
!    !nx_global = -1
!    !ny_global = -1
!    !nz_global = -1
!    !x_min =  c_largest_number
!    !x_max = -c_largest_number
!    !y_min =  c_largest_number
!    !y_max = -c_largest_number
!    !z_min =  c_largest_number
!    !z_max = -c_largest_number
!    !nprocx = 0
!    !nprocy = 0
!    !nprocz = 0
!    !npart_global = 0
!    !nsteps = -1
!    !t_end = c_largest_number
!    !dt_multiplier = 0.95_num
!    !dlb_threshold = 1.0_num
!    !use_balance = .FALSE.
!    !CALL set_field_order(2)
!    !stdout_frequency = 10
!    !use_random_seed = .FALSE.
!    !particle_push_start_time = 0.0_num
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    ! Boundaries
!
!    !bc_field(c_bd_x_min) = c_bc_periodic
!    !bc_particle(c_bd_x_min) = c_bc_periodic
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    ! Species
!
!    ! MANDATORY
!    CALL setup_species(current_species, name)
!
!    ! mass -- MANDATORY
!    !current_species%mass = 1.0_num * m0
!
!    ! charge -- MANDATORY
!    !current_species%charge = 1.0_num * q0
!
!    ! frac / fraction
!    !IF (npart_global >= 0) THEN
!    !  current_species%count = 0.1_num * npart_global
!    !ELSE
!    !  current_species%count = 0
!    !ENDIF
!
!    ! npart
!    !current_species%count = 1000
!
!    ! npart_per_cell
!    !current_species%npart_per_cell = 4
!
!    ! immobile
!    !current_species%immobile = .TRUE.
!
!    ! Initial conditions
!    IF (deck_state /= c_ds_first) THEN
!
!      ! density_min / minrho
!      !current_species%density_min = EPSILON(1.0_num)
!
!      ! density_max / maxrho
!      !current_species%density_max = HUGE(1.0_num)
!
!      ! density / rho / mass_density
!      !mult = 1.0_num
!      !mult = 1.0_num / current_species%mass ! mass_density
!
!      !DO iz = -2, nz+3
!      !zz = z(iz)
!      !DO iy = -2, ny+3
!      !yy = y(iy)
!      !DO ix = -2, nx+3
!      !  xx = x(ix)
!      !  current_species%density(ix,iy,iz) = mult * species_density(xx,yy,zz)
!      !ENDDO
!      !ENDDO
!      !ENDDO
!
!      ! drift_x
!      !DO iz = -2, nz+3
!      !zz = z(iz)
!      !DO iy = -2, ny+3
!      !yy = y(iy)
!      !DO ix = -2, nx+3
!      !  xx = x(ix)
!      !  current_species%drift(ix,iy,iz,1) = species_drift_x(xx,yy,zz)
!      !ENDDO
!      !ENDDO
!      !ENDDO
!
!      ! drift_y
!      !DO iz = -2, nz+3
!      !zz = z(iz)
!      !DO iy = -2, ny+3
!      !yy = y(iy)
!      !DO ix = -2, nx+3
!      !  xx = x(ix)
!      !  current_species%drift(ix,iy,iz,2) = species_drift_y(xx,yy,zz)
!      !ENDDO
!      !ENDDO
!      !ENDDO
!
!      ! drift_z
!      !DO iz = -2, nz+3
!      !zz = z(iz)
!      !DO iy = -2, ny+3
!      !yy = y(iy)
!      !DO ix = -2, nx+3
!      !  xx = x(ix)
!      !  current_species%drift(ix,iy,iz,3) = species_drift_z(xx,yy,zz)
!      !ENDDO
!      !ENDDO
!      !ENDDO
!
!      ! temp / temp_k / temp_ev
!      !mult = 1.0_num
!      !mult = ev / kb ! temp_ev
!      !
!      !DO iz = -2, nz+3
!      !zz = z(iz)
!      !DO iy = -2, ny+3
!      !yy = y(iy)
!      !DO ix = -2, nx+3
!      !  xx = x(ix)
!      !  current_species%temp(ix,iy,iz,1) = species_temperature_x(xx,yy,zz)
!      !  current_species%temp(ix,iy,iz,2) = current_species%temp(ix,iy,iz,1)
!      !  current_species%temp(ix,iy,iz,3) = current_species%temp(ix,iy,iz,1)
!      !ENDDO
!      !ENDDO
!      !ENDDO
!
!      ! temp_x / temp_x_k / temp_x_ev
!      !mult = 1.0_num
!      !mult = ev / kb ! temp_x_ev
!      !
!      !DO iz = -2, nz+3
!      !zz = z(iz)
!      !DO iy = -2, ny+3
!      !yy = y(iy)
!      !DO ix = -2, nx+3
!      !  xx = x(ix)
!      !  current_species%temp(ix,iy,iz,1) = species_temperature_x(xx,yy,zz)
!      !ENDDO
!      !ENDDO
!      !ENDDO
!
!      ! temp_y / temp_y_k / temp_y_ev
!      !mult = 1.0_num
!      !mult = ev / kb ! temp_y_ev
!      !
!      !DO iz = -2, nz+3
!      !zz = z(iz)
!      !DO iy = -2, ny+3
!      !yy = y(iy)
!      !DO ix = -2, nx+3
!      !  xx = x(ix)
!      !  current_species%temp(ix,iy,iz,2) = species_temperature_y(xx,yy,zz)
!      !ENDDO
!      !ENDDO
!      !ENDDO
!
!      ! temp_z / temp_z_k / temp_z_ev
!      !mult = 1.0_num
!      !mult = ev / kb ! temp_z_ev
!      !
!      !DO iz = -2, nz+3
!      !zz = z(iz)
!      !DO iy = -2, ny+3
!      !yy = y(iy)
!      !DO ix = -2, nx+3
!      !  xx = x(ix)
!      !  current_species%temp(ix,iy,iz,3) = species_temperature_z(xx,yy,zz)
!      !ENDDO
!      !ENDDO
!      !ENDDO
!    ENDIF
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    ! Lasers
!
!    IF (deck_state == c_ds_first) THEN
!      ! boundary / direction -- MANDATORY
!      !CALL create_laser(c_bd_x_min, working_laser)
!
!      ! amp -- MANDATORY
!      !working_laser%amp = 1.0_num
!
!      ! SI (W/m^2) - irradiance / intensity
!      !working_laser%amp = SQRT(amp / (c*epsilon0/2.0_num))
!
!      ! irradiance_w_cm2 / intensity_w_cm2
!      !working_laser%amp = SQRT(amp / (c*epsilon0/2.0_num)) * 100_num
!
!      ! omega / freq -- MANDATORY
!      !omega = 0.0_num
!      !working_laser%omega = omega
!
!      ! frequency
!      !working_laser%omega = 2.0_num * pi * omega
!
!      ! lambda
!      !working_laser%omega = 2.0_num * pi * c / omega
!    ELSE
!      ! profile
!      !working_laser%profile = 0.0_num
!      !CALL laser_update_profile(working_laser)
!
!      ! phase
!      !working_laser%phase = 0.0_num
!      !CALL laser_update_phase(working_laser)
!
!      !working_laser%t_start = 0.0_num
!
!      !working_laser%t_end = t_end
!
!      ! pol_angle
!      !working_laser%pol_angle = 0.0_num
!
!      ! pol
!      ! Convert from degrees to radians
!      !working_laser%pol_angle = pi * pol_angle / 180.0_num
!
!      ! id
!      !working_laser%id = 0
!    ENDIF
!
!    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    ! Fields
!
!    IF (deck_state /= c_ds_first) THEN
!      !ex = 0.0_num
!      !ey = 0.0_num
!      !ez = 0.0_num
!      !bx = 0.0_num
!      !by = 0.0_num
!      !bz = 0.0_num
!    ENDIF
!
!  END SUBROUTINE example_problem_setup

END MODULE ic_module
