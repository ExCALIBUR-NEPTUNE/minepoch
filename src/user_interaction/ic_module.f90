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



  SUBROUTINE custom_problem_setup(deck_state)

    INTEGER, INTENT(IN) :: deck_state

    CALL two_stream_setup(deck_state)

  END SUBROUTINE custom_problem_setup



  SUBROUTINE two_stream_setup(deck_state)

    INTEGER, INTENT(IN) :: deck_state
    REAL(num), PARAMETER :: drift_p = 2.5e-24_num
    REAL(num), PARAMETER :: temp = 273.0_num
    REAL(num), PARAMETER :: dens_max = 10.0_num
    REAL(num), PARAMETER :: dens_scale = 2.0e9_num
    INTEGER, PARAMETER :: ppc = 4
    REAL(num) :: xc, yc, zc, dens
    INTEGER :: ix, iy, iz
    TYPE(particle_species), POINTER :: current_species

    ! Control

    nx_global = 64
    ny_global = 3
    nz_global = 3
    x_min = 0.0_num
    x_max = 5.0e5_num
    y_min = x_min
    y_max = x_max
    z_min = x_min
    z_max = x_max
    t_end = 0.07_num
    stdout_frequency = 10

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Species

    ! MANDATORY
    NULLIFY(current_species)
    CALL setup_species(current_species, 'Right')

    ! mass -- MANDATORY
    current_species%mass = 1.0_num * m0

    ! charge -- MANDATORY
    current_species%charge = -1.0_num * q0

    ! npart_per_cell
    current_species%npart_per_cell = ppc

    IF (deck_state /= c_ds_first) THEN
      ! density
      IF (particles_uniformly_distributed) THEN
        current_species%density = dens_max
      ELSE
        DO iz = 1-ng, nz+ng
          zc = (z(iz) - 0.5_num*z_max)**2
        
          DO iy = 1-ng, ny+ng
            yc = (y(iy) - 0.5_num*y_max)**2
          
            DO ix = 1-ng, nx+ng
              xc = (x(ix) - 0.75_num*x_max)**2
                        
              dens = dens_max*EXP(-(xc + yc + zc)/dens_scale)
              current_species%density(ix,iy,iz) = dens
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      ! drift_x
      current_species%drift(:,:,:,1) = drift_p

      ! temp_x
      current_species%temp(:,:,:,1) = temp
    ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! MANDATORY
    NULLIFY(current_species)
    CALL setup_species(current_species, 'Left')

    ! mass -- MANDATORY
    current_species%mass = 1.0_num * m0

    ! charge -- MANDATORY
    current_species%charge = -1.0_num * q0

    ! npart_per_cell
    current_species%npart_per_cell = 4

    IF (deck_state /= c_ds_first) THEN
      ! density
      IF (particles_uniformly_distributed) THEN
        current_species%density = dens_max
      ELSE
        DO iz = 1-ng, nz+ng
          zc = (z(iz) - 0.5_num*z_max)**2
        
          DO iy = 1-ng, ny+ng
            yc = (y(iy) - 0.5_num*y_max)**2
          
            DO ix = 1-ng, nx+ng
              xc = (x(ix) - 0.25_num*x_max)**2
                        
              dens = dens_max*EXP(-(xc + yc + zc)/dens_scale)
              current_species%density(ix,iy,iz) = dens
            ENDDO
          ENDDO
        ENDDO
      ENDIF

      ! drift_x
      current_species%drift(:,:,:,1) = -drift_p

      ! temp_x
      current_species%temp(:,:,:,1) = temp
    ENDIF

  END SUBROUTINE two_stream_setup



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
