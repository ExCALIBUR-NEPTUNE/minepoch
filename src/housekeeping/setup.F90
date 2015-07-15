MODULE setup

  USE fields
  USE version_data
  USE laser
  USE timer
  USE strings

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: after_control, minimal_init
  PUBLIC :: setup_species, after_deck_last, set_dt

CONTAINS

  SUBROUTINE minimal_init

    INTEGER, PARAMETER :: r4  = SELECTED_REAL_KIND(r=30)
    INTEGER, PARAMETER :: r8  = SELECTED_REAL_KIND(r=300)
    INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(r=3000)
    INTEGER :: ierr

    IF (num == r4) THEN
      ! Should use MPI_SIZEOF() but this breaks on scalimpi
      realsize = 4
      mpireal = MPI_REAL4
    ELSE IF (num == r8) THEN
      realsize = 8
      mpireal = MPI_REAL8
    ELSE
      IF (rank == 0) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Cannot determine size of real'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      STOP
    ENDIF

    dt_plasma_frequency = 0.0_num
    dt_multiplier = 0.95_num
    stdout_frequency = 0

    npart_global = -1
    use_balance = .FALSE.
    use_random_seed = .FALSE.
    need_random_state = .FALSE.
    nsteps = -1
    t_end = HUGE(1.0_num)

    NULLIFY(laser_x_min)
    NULLIFY(laser_x_max)
    NULLIFY(laser_y_max)
    NULLIFY(laser_y_min)
    NULLIFY(laser_z_max)
    NULLIFY(laser_z_min)

    CALL set_field_order(2)

    CALL timer_init

    ! This array is true if a field component is staggered in the
    ! given direction.
    stagger = .FALSE.
    stagger(c_dir_x,c_stagger_ex) = .TRUE.
    stagger(c_dir_y,c_stagger_ey) = .TRUE.
    stagger(c_dir_z,c_stagger_ez) = .TRUE.

    stagger(c_dir_x:c_dir_z,c_stagger_bx) = .TRUE.
    stagger(c_dir_x:c_dir_z,c_stagger_by) = .TRUE.
    stagger(c_dir_x:c_dir_z,c_stagger_bz) = .TRUE.
    stagger(c_dir_x,c_stagger_bx) = .FALSE.
    stagger(c_dir_y,c_stagger_by) = .FALSE.
    stagger(c_dir_z,c_stagger_bz) = .FALSE.

    ALLOCATE(x(1), y(1), z(1))
    x = 0.0_num
    y = 0.0_num
    z = 0.0_num

  END SUBROUTINE minimal_init



  SUBROUTINE after_control

    CALL setup_grid
    CALL set_initial_values

  END SUBROUTINE after_control



  SUBROUTINE setup_grid

    INTEGER :: iproc, ix, iy, iz
    REAL(num) :: xb_min, yb_min, zb_min

    length_x = x_max - x_min
    dx = length_x / REAL(nx_global, num)
    x_grid_min = x_min
    x_grid_max = x_max

    length_y = y_max - y_min
    dy = length_y / REAL(ny_global, num)
    y_grid_min = y_min
    y_grid_max = y_max

    length_z = z_max - z_min
    dz = length_z / REAL(nz_global, num)
    z_grid_min = z_min
    z_grid_max = z_max

    ! Shift grid to cell centres.
    ! At some point the grid may be redefined to be node centred.

    xb_min = x_grid_min
    yb_min = y_grid_min
    zb_min = z_grid_min
    x_grid_min = x_grid_min + dx / 2.0_num
    x_grid_max = x_grid_max - dx / 2.0_num
    y_grid_min = y_grid_min + dy / 2.0_num
    y_grid_max = y_grid_max - dy / 2.0_num
    z_grid_min = z_grid_min + dz / 2.0_num
    z_grid_max = z_grid_max - dz / 2.0_num

    ! Setup global grid
    DO ix = -2, nx_global + 3
      x_global(ix) = x_grid_min + (ix - 1) * dx
    ENDDO
    DO ix = 1, nx_global + 1
      xb_global(ix) = xb_min + (ix - 1) * dx
    ENDDO
    DO iy = -2, ny_global + 3
      y_global(iy) = y_grid_min + (iy - 1) * dy
    ENDDO
    DO iy = 1, ny_global + 1
      yb_global(iy) = yb_min + (iy - 1) * dy
    ENDDO
    DO iz = -2, nz_global + 3
      z_global(iz) = z_grid_min + (iz - 1) * dz
    ENDDO
    DO iz = 1, nz_global + 1
      zb_global(iz) = zb_min + (iz - 1) * dz
    ENDDO

    DO iproc = 0, nprocx-1
      x_grid_mins(iproc) = x_global(cell_x_min(iproc+1))
      x_grid_maxs(iproc) = x_global(cell_x_max(iproc+1))
    ENDDO
    DO iproc = 0, nprocy-1
      y_grid_mins(iproc) = y_global(cell_y_min(iproc+1))
      y_grid_maxs(iproc) = y_global(cell_y_max(iproc+1))
    ENDDO
    DO iproc = 0, nprocz-1
      z_grid_mins(iproc) = z_global(cell_z_min(iproc+1))
      z_grid_maxs(iproc) = z_global(cell_z_max(iproc+1))
    ENDDO

    x_grid_min_local = x_grid_mins(x_coords)
    x_grid_max_local = x_grid_maxs(x_coords)
    y_grid_min_local = y_grid_mins(y_coords)
    y_grid_max_local = y_grid_maxs(y_coords)
    z_grid_min_local = z_grid_mins(z_coords)
    z_grid_max_local = z_grid_maxs(z_coords)

    x_min_local = x_grid_min_local - 0.5_num * dx
    x_max_local = x_grid_max_local + 0.5_num * dx
    y_min_local = y_grid_min_local - 0.5_num * dy
    y_max_local = y_grid_max_local + 0.5_num * dy
    z_min_local = z_grid_min_local - 0.5_num * dz
    z_max_local = z_grid_max_local + 0.5_num * dz

    ! Setup local grid
    DO ix = -2, nx + 3
      x(ix) = x_global(nx_global_min+ix-1)
    ENDDO
    DO iy = -2, ny + 3
      y(iy) = y_global(ny_global_min+iy-1)
    ENDDO
    DO iz = -2, nz + 3
      z(iz) = z_global(nz_global_min+iz-1)
    ENDDO

  END SUBROUTINE setup_grid



  SUBROUTINE after_deck_last

    CALL setup_field_boundaries

  END SUBROUTINE after_deck_last



  SUBROUTINE setup_species(next_species, name)

    INTEGER :: i
    CHARACTER(LEN=*), INTENT(IN) :: name
    TYPE(particle_species), POINTER :: species, next_species
    TYPE(particle_species), POINTER :: iolist, next_iolist

    ! First check if the species already exists or set "species" to the end of
    ! the list if it doesn't
    next_species => species_list
    next_iolist => io_list_data
    DO i = 1, n_species
      species => next_species
      IF (str_cmp(species%name, name)) THEN
        next_species => species
        RETURN
      ENDIF
      next_species => species%next
      iolist => next_iolist
      next_iolist => iolist%next
    ENDDO

    CALL setup_single_species(next_species)
    ALLOCATE(next_iolist)
    n_species = n_species + 1
    next_species%id = n_species
    next_species%name = TRIM(name)

    IF (n_species == 1) THEN
      species_list => next_species
      io_list_data => next_iolist
      io_list => species_list
    ELSE
      species%next => next_species
      next_species%prev => species
      iolist%next => next_iolist
      next_iolist%prev => iolist
    ENDIF

  END SUBROUTINE setup_species



  SUBROUTINE setup_single_species(species)

    TYPE(particle_species), POINTER, INTENT(OUT) :: species

    ALLOCATE(species)

    species%name = blank
    species%mass = -1.0_num
    species%charge = 0.0_num
    species%weight = 1.0_num
    species%count = -1
    species%id = 0
    species%npart_per_cell = -1
    species%count_update_step = 0
    species%immobile = .FALSE.
    NULLIFY(species%next)
    NULLIFY(species%prev)
    NULLIFY(species%ext_temp_x_min)
    NULLIFY(species%ext_temp_x_max)
    NULLIFY(species%ext_temp_y_min)
    NULLIFY(species%ext_temp_y_max)
    NULLIFY(species%ext_temp_z_min)
    NULLIFY(species%ext_temp_z_max)

    NULLIFY(species%density)
    NULLIFY(species%temp)
    NULLIFY(species%drift)
    species%density_min = EPSILON(1.0_num)
    species%density_max = HUGE(1.0_num)

    NULLIFY(species%attached_list%next)
    NULLIFY(species%attached_list%prev)
    CALL create_empty_partlist(species%attached_list)

  END SUBROUTINE setup_single_species



  SUBROUTINE setup_field_boundaries

    INTEGER :: nx0, nx1, ny0, ny1, nz0, nz1

    ALLOCATE(ex_x_min(-2:ny+3,-2:nz+3), ex_x_max(-2:ny+3,-2:nz+3))
    ALLOCATE(ey_x_min(-2:ny+3,-2:nz+3), ey_x_max(-2:ny+3,-2:nz+3))
    ALLOCATE(ez_x_min(-2:ny+3,-2:nz+3), ez_x_max(-2:ny+3,-2:nz+3))
    ALLOCATE(bx_x_min(-2:ny+3,-2:nz+3), bx_x_max(-2:ny+3,-2:nz+3))
    ALLOCATE(by_x_min(-2:ny+3,-2:nz+3), by_x_max(-2:ny+3,-2:nz+3))
    ALLOCATE(bz_x_min(-2:ny+3,-2:nz+3), bz_x_max(-2:ny+3,-2:nz+3))

    ALLOCATE(ex_y_min(-2:nx+3,-2:nz+3), ex_y_max(-2:nx+3,-2:nz+3))
    ALLOCATE(ey_y_min(-2:nx+3,-2:nz+3), ey_y_max(-2:nx+3,-2:nz+3))
    ALLOCATE(ez_y_min(-2:nx+3,-2:nz+3), ez_y_max(-2:nx+3,-2:nz+3))
    ALLOCATE(bx_y_min(-2:nx+3,-2:nz+3), bx_y_max(-2:nx+3,-2:nz+3))
    ALLOCATE(by_y_min(-2:nx+3,-2:nz+3), by_y_max(-2:nx+3,-2:nz+3))
    ALLOCATE(bz_y_min(-2:nx+3,-2:nz+3), bz_y_max(-2:nx+3,-2:nz+3))

    ALLOCATE(ex_z_min(-2:nx+3,-2:ny+3), ex_z_max(-2:nx+3,-2:ny+3))
    ALLOCATE(ey_z_min(-2:nx+3,-2:ny+3), ey_z_max(-2:nx+3,-2:ny+3))
    ALLOCATE(ez_z_min(-2:nx+3,-2:ny+3), ez_z_max(-2:nx+3,-2:ny+3))
    ALLOCATE(bx_z_min(-2:nx+3,-2:ny+3), bx_z_max(-2:nx+3,-2:ny+3))
    ALLOCATE(by_z_min(-2:nx+3,-2:ny+3), by_z_max(-2:nx+3,-2:ny+3))
    ALLOCATE(bz_z_min(-2:nx+3,-2:ny+3), bz_z_max(-2:nx+3,-2:ny+3))

    nx0 = 1
    nx1 = nx

    ex_x_min = 0.5_num * (ex(nx0,:,:) + ex(nx0-1,:,:))
    ey_x_min = ey(nx0,:,:)
    ez_x_min = ez(nx0,:,:)
    ex_x_max = 0.5_num * (ex(nx1,:,:) + ex(nx1-1,:,:))
    ey_x_max = ey(nx1,:,:)
    ez_x_max = ez(nx1,:,:)

    bx_x_min = bx(nx0,:,:)
    by_x_min = 0.5_num * (by(nx0,:,:) + by(nx0-1,:,:))
    bz_x_min = 0.5_num * (bz(nx0,:,:) + bz(nx0-1,:,:))
    bx_x_max = bx(nx1,:,:)
    by_x_max = 0.5_num * (by(nx1,:,:) + by(nx1-1,:,:))
    bz_x_max = 0.5_num * (bz(nx1,:,:) + bz(nx1-1,:,:))

    ny0 = 1
    ny1 = ny

    ex_y_min = ex(:,ny0,:)
    ey_y_min = 0.5_num * (ey(:,ny0,:) + ey(:,ny0-1,:))
    ez_y_min = ez(:,ny0,:)
    ex_y_max = ex(:,ny1,:)
    ey_y_max = 0.5_num * (ey(:,ny1,:) + ey(:,ny1-1,:))
    ez_y_max = ez(:,ny1,:)

    bx_y_min = 0.5_num * (bx(:,ny0,:) + bx(:,ny0-1,:))
    by_y_min = by(:,ny0,:)
    bz_y_min = 0.5_num * (bz(:,ny0,:) + bz(:,ny0-1,:))
    bx_y_max = 0.5_num * (bx(:,ny1,:) + bx(:,ny1-1,:))
    by_y_max = by(:,ny1,:)
    bz_y_max = 0.5_num * (bz(:,ny1,:) + bz(:,ny1-1,:))

    nz0 = 1
    nz1 = nz

    bx_z_min = 0.5_num * (bx(:,:,nz0) + bx(:,:,nz0-1))

    ex_z_min = ex(:,:,nz0)
    ey_z_min = ey(:,:,nz0)
    ez_z_min = 0.5_num * (ez(:,:,nz0) + ez(:,:,nz0-1))
    ex_z_max = ex(:,:,nz1)
    ey_z_max = ey(:,:,nz1)
    ez_z_max = 0.5_num * (ez(:,:,nz1) + ez(:,:,nz1-1))

    bx_z_min = 0.5_num * (bx(:,:,nz0) + bx(:,:,nz0-1))
    by_z_min = 0.5_num * (by(:,:,nz0) + by(:,:,nz0-1))
    bz_z_min = bz(:,:,nz0)
    bx_z_max = 0.5_num * (bx(:,:,nz1) + bx(:,:,nz1-1))
    by_z_max = 0.5_num * (by(:,:,nz1) + by(:,:,nz1-1))
    bz_z_max = bz(:,:,nz1)

  END SUBROUTINE setup_field_boundaries



  SUBROUTINE set_initial_values

    INTEGER :: seed

    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num

    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num

    jx = 0.0_num
    jy = 0.0_num
    jz = 0.0_num

    ! Set up random number seed
    seed = 7842432
    IF (use_random_seed) CALL SYSTEM_CLOCK(seed)
    seed = seed + rank

    CALL random_init(seed)

  END SUBROUTINE set_initial_values



  SUBROUTINE set_plasma_frequency_dt

    INTEGER :: ispecies, ix, iy, iz
    REAL(num) :: min_dt, omega2, omega, k_max, fac1, fac2
    TYPE(particle_species), POINTER :: species, next_species

    dt_plasma_frequency = 0.0_num

    IF (n_species < 1) RETURN

    min_dt = 1000000.0_num
    k_max = 2.0_num * pi / MIN(dx, dy, dz)

    ! Identify the plasma frequency (Bohm-Gross dispersion relation)
    ! Note that this doesn't get strongly relativistic plasmas right
    next_species => species_list
    DO ispecies = 1, n_species
      species => next_species
      next_species => species%next

      fac1 = q0**2 / species%mass / epsilon0
      fac2 = 3.0_num * k_max**2 * kb / species%mass
      DO iz = 1, nz
      DO iy = 1, ny
      DO ix = 1, nx
        omega2 = fac1 * species%density(ix,iy,iz) &
            + fac2 * MAXVAL(species%temp(ix,iy,iz,:))
        IF (omega2 <= c_tiny) CYCLE
        omega = SQRT(omega2)
        IF (2.0_num * pi / omega < min_dt) min_dt = 2.0_num * pi / omega
      ENDDO ! ix
      ENDDO ! iy
      ENDDO ! iz
    ENDDO

    CALL MPI_ALLREDUCE(min_dt, dt_plasma_frequency, 1, mpireal, MPI_MIN, &
        comm, errcode)
    ! Must resolve plasma frequency
    dt_plasma_frequency = dt_plasma_frequency / 2.0_num

  END SUBROUTINE set_plasma_frequency_dt



  SUBROUTINE set_dt        ! sets CFL limited step

    CALL set_plasma_frequency_dt
    CALL set_laser_dt

    dt = cfl * dx * dy * dz / SQRT((dx*dy)**2 + (dy*dz)**2 + (dz*dx)**2) / c
    IF (dt_plasma_frequency > c_tiny) dt = MIN(dt, dt_plasma_frequency)
    IF (dt_laser > c_tiny) dt = MIN(dt, dt_laser)
    dt = dt_multiplier * dt

  END SUBROUTINE set_dt

END MODULE setup
