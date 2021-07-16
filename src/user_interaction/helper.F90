MODULE helper

  USE boundary
  USE strings
  USE partlist
  USE deltaf_loader2, ONLY: delstaf_load
  USE utilities
  USE particle_loading

  IMPLICIT NONE

CONTAINS

  SUBROUTINE set_thermal_bcs

    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species, next_species

    next_species => species_list
    DO ispecies = 1, n_species
      species => next_species
      next_species => species%next

      ! Set temperature at boundary for thermal bcs.

      IF (bc_particle(c_bd_x_min) == c_bc_thermal) THEN
        species%ext_temp_x_min(-2:ny+3,-2:nz+3,1:3) = &
            species%temp(1,-2:ny+3,-2:nz+3,1:3)
      ENDIF
      IF (bc_particle(c_bd_x_max) == c_bc_thermal) THEN
        species%ext_temp_x_max(-2:ny+3,-2:nz+3,1:3) = &
            species%temp(nx,-2:ny+3,-2:nz+3,1:3)
      ENDIF
      IF (bc_particle(c_bd_y_min) == c_bc_thermal) THEN
        species%ext_temp_y_min(-2:nx+3,-2:nz+3,1:3) = &
            species%temp(-2:nx+3,1,-2:nz+3,1:3)
      ENDIF
      IF (bc_particle(c_bd_y_max) == c_bc_thermal) THEN
        species%ext_temp_y_max(-2:nx+3,-2:nz+3,1:3) = &
            species%temp(-2:nx+3,ny,-2:nz+3,1:3)
      ENDIF
      IF (bc_particle(c_bd_z_min) == c_bc_thermal) THEN
        species%ext_temp_z_min(-2:nx+3,-2:ny+3,1:3) = &
            species%temp(-2:nx+3,-2:ny+3,1,1:3)
      ENDIF
      IF (bc_particle(c_bd_z_max) == c_bc_thermal) THEN
        species%ext_temp_z_max(-2:nx+3,-2:ny+3,1:3) = &
            species%temp(-2:nx+3,-2:ny+3,nz,1:3)
      ENDIF
    ENDDO

  END SUBROUTINE set_thermal_bcs



  SUBROUTINE auto_load

    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species, next_species

    IF (n_species < 1) RETURN

    CALL set_thermal_bcs

    next_species => species_list
    DO ispecies = 1, n_species
      species => next_species
      next_species => species%next

      IF (particles_uniformly_distributed) THEN
        CALL setup_particle_density(species%density, species, &
            species%density_min, species%density_max)
      ELSE
        CALL non_uniform_load_particles(species%density, species, &
            species%density_min, species%density_max)
      ENDIF
      
      CALL setup_particle_temperature(species%temp(:,:,:,1), c_dir_x, species, &
          species%drift(:,:,:,1))
      CALL setup_particle_temperature(species%temp(:,:,:,2), c_dir_y, species, &
          species%drift(:,:,:,2))
      CALL setup_particle_temperature(species%temp(:,:,:,3), c_dir_z, species, &
          species%drift(:,:,:,3))
      CALL delstaf_load(ispecies, species%temp, species%drift)
      IF(species%is_driftkinetic) THEN
         CALL setup_particle_driftkinetic(species)
         drift_kinetic_species_exist = .true.
      END IF

      IF (rank == 0) THEN
        IF (species%count < 0) THEN
          WRITE(*,*) 'No particles specified for species ', &
              '"' // TRIM(species%name) // '"'
          species%count = 0
        ENDIF
      ENDIF
    ENDDO

  END SUBROUTINE auto_load



  SUBROUTINE setup_particle_driftkinetic(part_species)
    TYPE(particle_species), POINTER, INTENT(INOUT) :: part_species
    TYPE(particle_list), POINTER :: partlist
    REAL(num) :: mass, Bnorm, mu, ppll, pperp
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart

    REAL(num), DIMENSION(3)  :: pvec,bvec,evec, bdir
    REAL(num), DIMENSION(3,3) :: btens

     partlist => part_species%attached_list
     current => partlist%head
     ipart = 0
     DO WHILE(ipart < partlist%count)
        mass = current%mass
        CALL get_fields_at_point(current%part_pos,bvec,evec,btens=btens)
        Bnorm = sqrt(dot_product(Bvec,Bvec))
        bdir = Bvec/Bnorm
        pvec = current%part_p
        ppll = dot_product(bdir,pvec)
        pperp = sqrt(max(dot_product(pvec,pvec)-ppll*ppll,0.0_num))
        mu = pperp*pperp/(2*mass*Bnorm)

        !part_p(3) is not used but could track gyroangle if useful.
        current%part_p(1) = ppll/mass
        current%part_p(2) = mu
        !Use this as a label.
        current%part_p(3) = ipart
        current => current%next
        ipart = ipart + 1
     ENDDO
  END SUBROUTINE setup_particle_driftkinetic

  SUBROUTINE allocate_ic

    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species, next_species

    next_species => species_list
    DO ispecies = 1, n_species
      species => next_species
      next_species => species%next

      ALLOCATE(species%density(-2:nx+3,-2:ny+3,-2:nz+3))
      ALLOCATE(species%temp (-2:nx+3,-2:ny+3,-2:nz+3,1:3))
      ALLOCATE(species%drift(-2:nx+3,-2:ny+3,-2:nz+3,1:3))

      species%density = 1.0_num
      species%temp = 0.0_num
      species%drift = 0.0_num
      species%density_min = EPSILON(1.0_num)
      species%density_max = HUGE(1.0_num)
    ENDDO

  END SUBROUTINE allocate_ic



  SUBROUTINE deallocate_ic

    INTEGER :: ispecies
    TYPE(particle_species), POINTER :: species, next_species

    next_species => species_list
    DO ispecies = 1, n_species
      species => next_species
      next_species => species%next

      DEALLOCATE(species%density)
      DEALLOCATE(species%temp)
      DEALLOCATE(species%drift)
    ENDDO

  END SUBROUTINE deallocate_ic


  
  SUBROUTINE non_uniform_load_particles(density, species, density_min, &
      density_max)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(INOUT) :: density
    TYPE(particle_species), POINTER :: species
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current
    REAL(num), INTENT(INOUT) :: density_min, density_max
    INTEGER(i8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(i8) :: npart_per_cell, npart_left
    REAL(num) :: density_total, density_total_global, density_average
    REAL(num) :: npart_per_cell_average
    INTEGER(i8) :: npart_this_proc_new, ipart, npart_this_species
    INTEGER :: ix, iy, iz
    CHARACTER(LEN=15) :: string
    
    
    num_valid_cells_local = 0
    density_total = 0.0_num

    DO iz = -2, nz+3
    DO iy = -2, ny+3
    DO ix = -2, nx+3
      IF (density(ix,iy,iz) > density_max) density(ix,iy,iz) = density_max
    ENDDO ! ix
    ENDDO ! iy
    ENDDO ! iz

    DO iz = 1, nz
    DO iy = 1, ny
    DO ix = 1, nx
      IF (density(ix,iy,iz) >= density_min) THEN
        num_valid_cells_local = num_valid_cells_local + 1
        density_total = density_total + density(ix,iy,iz)
      ENDIF
    ENDDO ! ix
    ENDDO ! iy
    ENDDO ! iz

    CALL MPI_ALLREDUCE(num_valid_cells_local, num_valid_cells_global, 1, &
        MPI_INTEGER8, MPI_SUM, comm, errcode)

    IF (species%npart_per_cell >= 0) THEN
      npart_per_cell_average = FLOOR(species%npart_per_cell, num)
    ELSE
      npart_per_cell_average = REAL(species%count, num) &
          / REAL(num_valid_cells_global, num)
    ENDIF

    IF (npart_per_cell_average <= 0) RETURN

    CALL MPI_ALLREDUCE(density_total, density_total_global, 1, mpireal, &
        MPI_SUM, comm, errcode)
    density_average = density_total_global / REAL(num_valid_cells_global, num)

    npart_this_proc_new = 0
    DO iz = 1, nz
    DO iy = 1, ny
    DO ix = 1, nx
      npart_per_cell = NINT(density(ix, iy, iz) / density_average &
          * npart_per_cell_average)
      npart_this_proc_new = npart_this_proc_new + npart_per_cell
    ENDDO ! ix
    ENDDO ! iy
    ENDDO ! iz

    partlist => species%attached_list
    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, npart_this_proc_new)
    npart_left = npart_this_proc_new
    current => partlist%head
    
    IF (npart_left > 0) THEN
            
      ! Randomly place npart_per_cell particles into each valid cell
      DO iz = 1, nz
      DO iy = 1, ny
      DO ix = 1, nx
        
        npart_per_cell = NINT(density(ix, iy, iz) / density_average &
            * npart_per_cell_average)

        ipart = 0
        DO WHILE(ASSOCIATED(current) .AND. ipart < npart_per_cell)
          ! Even if particles have per particle charge and mass, assume
          ! that initially they all have the same charge and mass (user
          ! can easily over_ride)
          current%charge = species%charge
          current%mass = species%mass
          current%pvol = npart_per_cell

          CALL init_particle_position(current, x(ix), y(iy), z(iz), dx, dy, dz)

          ipart = ipart + 1
          current => current%next

          ! One particle sucessfully placed
          npart_left = npart_left - 1
        ENDDO
     
      ENDDO ! ix
      ENDDO ! iy
      ENDDO ! iz
    
    ENDIF

    
    CALL MPI_ALLREDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)

    species%count = npart_this_species
    species%weight = density_total_global * dx * dy * dz / npart_this_species
    
    IF (rank == 0) THEN
      CALL integer_as_string(npart_this_species, string)
      WRITE(*,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', '"' // TRIM(species%name) // '"'
    ENDIF

  END SUBROUTINE non_uniform_load_particles



  ! This subroutine automatically loads a uniform density of pseudoparticles
  SUBROUTINE load_particles(species, load_list)

    TYPE(particle_species), POINTER :: species
    LOGICAL, DIMENSION(-2:,-2:,-2:), INTENT(IN) :: load_list
    INTEGER(i8), DIMENSION(:), ALLOCATABLE :: valid_cell_list
    TYPE(particle_list), POINTER :: partlist
    TYPE(particle), POINTER :: current, next
    INTEGER(i8) :: ipart, npart_per_cell, num_int, num_total, idx
    INTEGER(i8) :: num_valid_cells_local, num_valid_cells_global
    INTEGER(i8) :: npart_this_species, num_new_particles, npart_left
    INTEGER(i8), ALLOCATABLE :: num_valid_cells_all(:), num_idx(:)
    REAL(num) :: valid_cell_frac, num_real, f0, f1
    REAL(num), ALLOCATABLE :: num_frac(:)
    INTEGER(i8) :: cell_x
    INTEGER(i8) :: cell_y
    INTEGER(i8) :: cell_z
    INTEGER(i8) :: i, ipos
    INTEGER :: ierr, ix, iy, iz
    CHARACTER(LEN=15) :: string
    LOGICAL :: sweep

    npart_this_species = species%count
    IF (npart_this_species <= 0) RETURN

    num_valid_cells_local = COUNT(load_list(1:nx, 1:ny, 1:nz))

    IF (species%npart_per_cell >= 0) THEN
      npart_per_cell = FLOOR(species%npart_per_cell, KIND=i8)
      num_new_particles = &
          FLOOR(species%npart_per_cell * num_valid_cells_local, KIND=i8)
    ELSE
      ALLOCATE(num_valid_cells_all(nproc), num_idx(nproc), num_frac(nproc))

      ! Calculate global number of particles per cell
      CALL MPI_ALLGATHER(num_valid_cells_local, 1, MPI_INTEGER8, &
          num_valid_cells_all, 1, MPI_INTEGER8, comm, errcode)

      num_valid_cells_global = 0
      DO i = 1,nproc
        num_valid_cells_global = num_valid_cells_global + num_valid_cells_all(i)
      ENDDO

      IF (num_valid_cells_global == 0) THEN
        IF (rank == 0) THEN
          WRITE(*,*) '*** ERROR ***'
          WRITE(*,*) 'Intial condition settings mean that there are no cells ' &
              // 'where particles may'
          WRITE(*,*) 'validly be placed for species "' // TRIM(species%name) &
              // '". ', 'Code will now terminate.'
          CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        ENDIF
      ENDIF

      valid_cell_frac = REAL(num_valid_cells_local, num) &
          / REAL(num_valid_cells_global, num)
      num_real = npart_this_species * valid_cell_frac
      num_new_particles = FLOOR(num_real, KIND=i8)

      ! Work out which processors get the remaining fractional numbers
      ! of particles

      ! Get a list of the fractional part on each processor, along with
      ! the total
      num_total = 0
      DO i = 1,nproc
        valid_cell_frac = REAL(num_valid_cells_all(i), num) &
            / REAL(num_valid_cells_global, num)
        num_real = npart_this_species * valid_cell_frac
        num_int = FLOOR(num_real, KIND=i8)
        num_frac(i) = num_real - num_int
        num_idx (i) = i - 1
        num_total = num_total + num_int
      ENDDO
      num_total = npart_this_species - num_total

      IF (num_total > 0) THEN
        ! Sort the list of fractions into decreasing order using bubble sort
        sweep = .TRUE.
        DO WHILE(sweep)
          sweep = .FALSE.
          f0 = num_frac(1)
          DO i = 2,nproc
            f1 = num_frac(i)
            IF (f1 > f0) THEN
              num_frac(i-1) = f1
              num_frac(i) = f0
              f1 = f0
              idx = num_idx(i-1)
              num_idx(i-1) = num_idx(i)
              num_idx(i) = idx
              sweep = .TRUE.
            ENDIF
            f0 = f1
          ENDDO
        ENDDO

        ! Accumulate fractional particles until they have all been accounted
        ! for. If any of them have been assigned to the current processor,
        ! add them and exit the loop.

        DO i = 1,nproc
          IF (num_idx(i) == rank) THEN
            num_new_particles = num_new_particles + 1
            EXIT
          ENDIF
          num_total = num_total - 1
          IF (num_total <= 0) EXIT
        ENDDO
      ENDIF

      DEALLOCATE(num_valid_cells_all, num_idx, num_frac)

      species%npart_per_cell = &
          REAL(npart_this_species,num) / num_valid_cells_global
      npart_per_cell = FLOOR(species%npart_per_cell, KIND=i8)
    ENDIF

    partlist => species%attached_list

    CALL destroy_partlist(partlist)
    CALL create_allocated_partlist(partlist, num_new_particles)

    ! Randomly place npart_per_cell particles into each valid cell
    npart_left = num_new_particles
    current => partlist%head
    IF (npart_per_cell > 0) THEN

      DO iz = 1, nz
      DO iy = 1, ny
      DO ix = 1, nx
        IF (.NOT. load_list(ix, iy, iz)) CYCLE

        ipart = 0
        DO WHILE(ASSOCIATED(current) .AND. ipart < npart_per_cell)
          ! Even if particles have per particle charge and mass, assume
          ! that initially they all have the same charge and mass (user
          ! can easily over_ride)
          current%charge = species%charge
          current%mass = species%mass
          current%pvol = npart_per_cell

          CALL init_particle_position(current, x(ix), y(iy), z(iz), dx, dy, dz)

          ipart = ipart + 1
          current => current%next

          ! One particle sucessfully placed
          npart_left = npart_left - 1
        ENDDO
      ENDDO ! ix
      ENDDO ! iy
      ENDDO ! iz

    ENDIF

    ! When num_new_particles does not equal
    ! npart_per_cell * num_valid_cells_local there will be particles left
    ! over that didn't get placed.
    ! The following loop randomly place remaining particles into valid cells.
    IF (npart_left > 0) THEN
      ALLOCATE(valid_cell_list(num_valid_cells_local))

      ipos = 0
      DO iz = 1, nz
      DO iy = 1, ny
      DO ix = 1, nx
        IF (load_list(ix,iy,iz)) THEN
          ipos = ipos + 1
          valid_cell_list(ipos) = ix - 1 + nx * (iy - 1 + ny * (iz - 1))
        ENDIF
      ENDDO ! ix
      ENDDO ! iy
      ENDDO ! iz

      DO i = 1, npart_left
        ipos = INT(random() * (num_valid_cells_local - 1)) + 1
        ipos = valid_cell_list(ipos)

        cell_z = ipos / (nx * ny)
        ipos = ipos - (nx * ny) * cell_z
        cell_z = cell_z + 1

        cell_y = ipos / nx
        ipos = ipos - nx * cell_y
        cell_y = cell_y + 1

        cell_x = ipos + 1

        CALL init_particle_position(current, x(cell_x), y(cell_y), z(cell_z), dx, dy, dz)

        current => current%next
      ENDDO

      DEALLOCATE(valid_cell_list)
    ENDIF

    ! Remove any unplaced particles from the list. This should never be
    ! called if the above routines worked correctly.
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL remove_particle_from_partlist(partlist, current)
      DEALLOCATE(current)
      current => next
    ENDDO

    CALL MPI_ALLREDUCE(partlist%count, npart_this_species, 1, MPI_INTEGER8, &
        MPI_SUM, comm, errcode)

    species%count = npart_this_species

    IF (rank == 0) THEN
      CALL integer_as_string(npart_this_species, string)
      WRITE(*,*) 'Loaded ', TRIM(ADJUSTL(string)), &
          ' particles of species ', '"' // TRIM(species%name) // '"'
    ENDIF

    CALL particle_bcs

  END SUBROUTINE load_particles



  SUBROUTINE setup_particle_density(density_in, species, density_min, &
      density_max)

    REAL(num), DIMENSION(-2:,-2:,-2:), INTENT(IN) :: density_in
    TYPE(particle_species), POINTER :: species
    REAL(num), INTENT(IN) :: density_min, density_max
    TYPE(particle), POINTER :: current
    INTEGER(i8) :: ipart
    INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: npart_in_cell
    REAL(num) :: wdata
    TYPE(particle_list), POINTER :: partlist
    INTEGER :: ix, iy, iz, i, j, k, isubx, isuby, isubz
    REAL(num), DIMENSION(:,:,:), ALLOCATABLE :: density
    LOGICAL, DIMENSION(:,:,:), ALLOCATABLE :: density_map
    REAL(num), DIMENSION(sf_min:sf_max) :: gx, gy, gz
    REAL(num) :: cell_x_r, cell_frac_x
    REAL(num) :: cell_y_r, cell_frac_y
    REAL(num) :: cell_z_r, cell_frac_z
    INTEGER :: cell_x, cell_y, cell_z
    REAL(num) :: cx2, cy2, cz2
    REAL(num), PARAMETER :: third = 1.0_num / 3.0_num
    REAL(num), PARAMETER :: fac1 = 0.125_num * third
    REAL(num), PARAMETER :: fac2 = 0.5_num * third
    REAL(num), PARAMETER :: fac3 = 7.1875_num * third

    ALLOCATE(density(-2:nx+3,-2:ny+3,-2:nz+3))
    ALLOCATE(density_map(-2:nx+3,-2:ny+3,-2:nz+3))
    density = density_in
    density_map = .FALSE.

    CALL field_bc(density, ng)

    DO iz = -2, nz+3
    DO iy = -2, ny+3
    DO ix = -2, nx+3
      IF (density(ix,iy,iz) > density_max) density(ix,iy,iz) = density_max
      IF (density(ix,iy,iz) >= density_min) THEN
        density_map(ix,iy,iz) = .TRUE.
      ELSE
        density(ix,iy,iz) = 0.0_num
      ENDIF
    ENDDO ! ix
    ENDDO ! iy
    ENDDO ! iz

    ! Uniformly load particles in space
    CALL load_particles(species, density_map)

    ALLOCATE(npart_in_cell(-2:nx+3,-2:ny+3,-2:nz+3))
    npart_in_cell = 0

    partlist => species%attached_list
    ! If using per particle weighing then use the weight function to match the
    ! uniform pseudoparticle density to the real particle density
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
      IF (.NOT. ASSOCIATED(current)) PRINT *, 'Bad Particle'

      cell_x_r = (current%part_pos(1) - x_grid_min_local) / dx
      cell_y_r = (current%part_pos(2) - y_grid_min_local) / dy
      cell_z_r = (current%part_pos(3) - z_grid_min_local) / dz
      cell_x = FLOOR(cell_x_r + 0.5_num)
      cell_y = FLOOR(cell_y_r + 0.5_num)
      cell_z = FLOOR(cell_z_r + 0.5_num)
      cell_frac_x = REAL(cell_x, num) - cell_x_r
      cell_frac_y = REAL(cell_y, num) - cell_y_r
      cell_frac_z = REAL(cell_z, num) - cell_z_r
      cell_x = cell_x + 1
      cell_y = cell_y + 1
      cell_z = cell_z + 1

      cx2 = cell_frac_x**2
      gx(-2) = fac1 * (0.5_num + cell_frac_x)**4
      gx(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_x &
          + cx2 * (1.5_num - cell_frac_x - cx2))
      gx( 0) = 0.25_num * (fac3 + cx2 * (cx2 - 2.5_num))
      gx( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_x &
          + cx2 * (1.5_num + cell_frac_x - cx2))
      gx( 2) = fac1 * (0.5_num - cell_frac_x)**4

      cy2 = cell_frac_y**2
      gy(-2) = fac1 * (0.5_num + cell_frac_y)**4
      gy(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_y &
          + cy2 * (1.5_num - cell_frac_y - cy2))
      gy( 0) = 0.25_num * (fac3 + cy2 * (cy2 - 2.5_num))
      gy( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_y &
          + cy2 * (1.5_num + cell_frac_y - cy2))
      gy( 2) = fac1 * (0.5_num - cell_frac_y)**4

      cz2 = cell_frac_z**2
      gz(-2) = fac1 * (0.5_num + cell_frac_z)**4
      gz(-1) = fac2 * (1.1875_num + 2.75_num * cell_frac_z &
          + cz2 * (1.5_num - cell_frac_z - cz2))
      gz( 0) = 0.25_num * (fac3 + cz2 * (cz2 - 2.5_num))
      gz( 1) = fac2 * (1.1875_num - 2.75_num * cell_frac_z &
          + cz2 * (1.5_num + cell_frac_z - cz2))
      gz( 2) = fac1 * (0.5_num - cell_frac_z)**4

      ! Calculate density at the particle position
      wdata = 0.0_num
      DO isubz = sf_min, sf_max
        i = cell_x
        j = cell_y
        k = cell_z + isubz
        IF (.NOT. density_map(i,j,k)) THEN
          k = cell_z + isubz / 2
          IF (.NOT. density_map(i,j,k)) k = cell_z - isubz / 2
        ENDIF
        DO isuby = sf_min, sf_max
          i = cell_x
          j = cell_y + isuby
          IF (.NOT. density_map(i,j,k)) THEN
            j = cell_y + isuby / 2
            IF (.NOT. density_map(i,j,k)) j = cell_y - isuby / 2
          ENDIF
          DO isubx = sf_min, sf_max
            i = cell_x + isubx
            IF (.NOT. density_map(i,j,k)) THEN
              i = cell_x + isubx / 2
              IF (.NOT. density_map(i,j,k)) i = cell_x - isubx / 2
            ENDIF
            wdata = wdata + gx(isubx) * gy(isuby) * gz(isubz) * density(i,j,k)
          ENDDO ! isubx
        ENDDO ! isuby
      ENDDO ! isubz

      current%weight = wdata
      npart_in_cell(cell_x,cell_y,cell_z) = &
          npart_in_cell(cell_x,cell_y,cell_z) + 1

      current => current%next
      ipart = ipart + 1
    ENDDO
    DEALLOCATE(density_map)
    DEALLOCATE(density)

    wdata = dx * dy * dz

    partlist => species%attached_list
    ! Second loop renormalises particle weights
    current => partlist%head
    ipart = 0
    DO WHILE(ipart < partlist%count)
      cell_x = FLOOR((current%part_pos(1) - x_grid_min_local) / dx + 1.5_num)
      cell_y = FLOOR((current%part_pos(2) - y_grid_min_local) / dy + 1.5_num)
      cell_z = FLOOR((current%part_pos(3) - z_grid_min_local) / dz + 1.5_num)

      current%weight = current%weight * wdata &
          / npart_in_cell(cell_x,cell_y,cell_z)

      current => current%next
      ipart = ipart + 1
    ENDDO

    DEALLOCATE(npart_in_cell)

  END SUBROUTINE setup_particle_density



  FUNCTION sample_dist_function(axis, dist_fn)

    REAL(num), DIMENSION(:), INTENT(IN) :: axis, dist_fn
    REAL(num), DIMENSION(:), ALLOCATABLE :: cdf
    REAL(num) :: position, d_cdf
    INTEGER :: n_points, ipoint, start, endpoint, current
    REAL(num) :: sample_dist_function

    n_points = SIZE(dist_fn)
    ALLOCATE(cdf(n_points))

    cdf(1) = dist_fn(1)
    DO ipoint = 2, n_points
      cdf(ipoint) = cdf(ipoint-1) + dist_fn(ipoint)
    ENDDO

    cdf = cdf / cdf(n_points)

    position = random()
    sample_dist_function = 0.0_num

    start = 1
    endpoint = n_points
    current = (start + endpoint) / 2

    DO current = 1, n_points-1
      IF (cdf(current) <= position .AND. cdf(current+1) >= position) THEN
        d_cdf = cdf(current+1) - cdf(current)
        sample_dist_function = (axis(current) * (position - cdf(current)) &
            + axis(current+1) * (cdf(current+1) - position)) / d_cdf
        EXIT
      ENDIF
    ENDDO

    DEALLOCATE(cdf)

  END FUNCTION sample_dist_function

END MODULE helper
