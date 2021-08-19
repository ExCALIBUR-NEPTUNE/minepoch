MODULE problem_setup_module

  USE ic_module
  USE shared_data
  USE timer
  USE fields
  USE strings

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: problem_setup

CONTAINS

  SUBROUTINE problem_setup(deck_state, problem)

    INTEGER, INTENT(IN) :: deck_state
    CHARACTER(LEN=c_max_string_length), INTENT(IN) :: problem

    IF (deck_state == c_ds_first) THEN
      CALL control_initialise
      CALL boundary_initialise(problem)
      CALL custom_problem_setup(deck_state, problem)
    ELSE
      CALL fields_initialise(problem)
      CALL custom_problem_setup(deck_state, problem)
      CALL laser_finalise
      CALL species_finalise
    ENDIF

  END SUBROUTINE problem_setup



  SUBROUTINE boundary_initialise(problem)

    CHARACTER(LEN=c_max_string_length), INTENT(IN) :: problem
    INTEGER :: i

    ! Could add problem specific boundary conditions here
    ! Default values
    DO i = c_bd_x_min, c_bd_z_max
      bc_field(i) = c_bc_periodic
      bc_particle(i) = c_bc_periodic
    ENDDO

  END SUBROUTINE boundary_initialise



  SUBROUTINE control_initialise

    ! Default values
    allow_cpu_reduce = .TRUE.
    timer_collect = .TRUE.
    nx_global = -1
    ny_global = -1
    nz_global = -1
    x_min =  c_largest_number
    x_max = -c_largest_number
    y_min =  c_largest_number
    y_max = -c_largest_number
    z_min =  c_largest_number
    z_max = -c_largest_number
    nprocx = 0
    nprocy = 0
    nprocz = 0
    npart_global = 0
    nsteps = -1
    t_end = c_largest_number
    dt_multiplier = 0.95_num
    dlb_threshold = 1.0_num
    use_balance = .FALSE.
    CALL set_field_order(2)
    stdout_frequency = 10
    use_random_seed = .FALSE.
    particle_push_start_time = 0.0_num
    n_species = 0
    fixed_fields = .FALSE.

  END SUBROUTINE control_initialise



  SUBROUTINE fields_initialise(problem)

    CHARACTER(LEN=c_max_string_length), INTENT(IN) :: problem
    INTEGER :: i,j,k
    REAL(num) :: xp
    REAL(num), PARAMETER :: em_wave_mag = 1e6_num


    ! Could include problem specific fields here.
    ! Default to zero
    ex = 0.0_num
    ey = 0.0_num
    ez = 0.0_num
    bx = 0.0_num
    by = 0.0_num
    bz = 0.0_num
    SELECT CASE (TRIM(problem))

    CASE ('one_stream')
       IF (.false.) THEN
          bx = 1.0_num
          ey = 1.0_num

          !To do this right, need to get the offsets to the grid positions.
          do i=1-ng,nx+ng
             do j=1-ng,ny+ng
                do k=1-ng,nz+ng
                   !Correct offset for bz, is there a variable
                   ! storing offset per field/direction?
                   xp = x_grid_min_local+(i-0.5_num)*dx
                   bz(i,j,k) = cos(xp)
                end do
             end do
          end do
       END IF

    CASE ('drift_kin_default')
       ex = 0.0_num
       ey = 0.0_num
       ez = 0.0_num
       bx = 0.0_num
       by = 0.0_num
       bz = 1.0_num

    CASE ('em_wave')
      ex = 0.0_num
      ez = 0.0_num
      bx = 0.0_num
      by = 0.0_num

      DO i = 1-ng, nx+ng
        xp = x(i)
        ey(i,:,:) = em_wave_mag * COS(xp)
        xp = xp + 0.5_num * dx
        bz(i,:,:) = em_wave_mag * COS(xp) / c
      END DO

    CASE default

    END SELECT

  END SUBROUTINE fields_initialise



  SUBROUTINE laser_finalise

    INTEGER :: errcode, ierr
    TYPE(laser_block), POINTER :: current
    INTEGER :: error, io, iu

    errcode = c_err_none

    error = 0
    current => laser_x_min
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%amp < 0.0_num) error = IOR(error, 2)
      current => current%next
    ENDDO

    current => laser_x_max
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%amp < 0.0_num) error = IOR(error, 2)
      current => current%next
    ENDDO

    current => laser_y_min
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%amp < 0.0_num) error = IOR(error, 2)
      current => current%next
    ENDDO

    current => laser_y_max
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%amp < 0.0_num) error = IOR(error, 2)
      current => current%next
    ENDDO

    current => laser_z_min
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%amp < 0.0_num) error = IOR(error, 2)
      current => current%next
    ENDDO

    current => laser_z_max
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%amp < 0.0_num) error = IOR(error, 2)
      current => current%next
    ENDDO

    IF (IAND(error,1) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "lambda" or "omega" for every laser.'
        ENDDO
      ENDIF
      errcode = c_err_missing_elements
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF

    IF (IAND(error,2) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define an "amp" or "irradiance" for every laser.'
        ENDDO
      ENDIF
      errcode = c_err_missing_elements
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF

  END SUBROUTINE laser_finalise



  SUBROUTINE species_finalise

    INTEGER :: i, io, iu
    CHARACTER(LEN=8) :: string
    INTEGER :: errcode, ierr
    TYPE(particle_species), POINTER :: species, next_species

    ! Sanity check
    next_species => species_list
    DO i = 1, n_species
      species => next_species
      next_species => species%next

      IF (rank == 0) THEN
        CALL integer_as_string(i, string)
        PRINT*, 'Name of species ', TRIM(ADJUSTL(string)), ' is ', &
            TRIM(species%name)
      ENDIF

      IF (species%mass < 0) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'Particle species "', TRIM(species%name), &
                '" has negative mass'
          ENDDO
        ENDIF
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      ENDIF

      IF (.NOT. (species%charge < HUGE(1.0_num))) THEN
        IF (rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** ERROR ***'
            WRITE(io,*) 'No charge specified for particle species "', &
                TRIM(species%name), '"'
          ENDDO
        ENDIF
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
      ENDIF

      IF (species%npart_per_cell >= 0) THEN
        IF (species%count >= 0 .AND. rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = io_units(iu)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Two forms of npart used for particle species "', &
                TRIM(species%name), '"'
            WRITE(io,*) 'Just using "npart_per_cell".'
          ENDDO
        ENDIF
        species%count = INT(species%npart_per_cell, i8)
      ENDIF
    ENDDO

  END SUBROUTINE species_finalise

END MODULE problem_setup_module
