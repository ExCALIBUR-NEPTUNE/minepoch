MODULE mpi_routines

  USE helper

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: mpi_initialise, mpi_close, mpi_minimal_init, setup_communicator

  REAL(dbl) :: start_time, end_time

CONTAINS

  SUBROUTINE mpi_minimal_init

    CALL MPI_INIT(errcode)
    CALL MPI_COMM_DUP(MPI_COMM_WORLD, comm, errcode)
    CALL MPI_COMM_SIZE(comm, nproc, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
#ifdef MPI_DEBUG
    CALL mpi_set_error_handler
#endif

  END SUBROUTINE mpi_minimal_init



  SUBROUTINE setup_communicator

    INTEGER, PARAMETER :: ndims = 3
    INTEGER :: dims(ndims), idim, old_comm, ierr
    LOGICAL :: periods(ndims), reorder, op, reset
    INTEGER :: test_coords(ndims)
    INTEGER :: ix, iy, iz
    INTEGER :: nxsplit, nysplit, nzsplit
    INTEGER :: area, minarea, nprocyz
    INTEGER :: ranges(3,1), nproc_orig, oldgroup, newgroup
    CHARACTER(LEN=11) :: str

    nproc_orig = nproc

    IF (nx_global < ng .OR. ny_global < ng .OR. nz_global < ng) THEN
      IF (rank == 0) THEN
        CALL integer_as_string(ng, str)
        PRINT*,'*** ERROR ***'
        PRINT*,'Simulation domain is too small.'
        PRINT*,'There must be at least ' // TRIM(str) // &
            ' cells in each direction.'
      ENDIF
      CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
    ENDIF

    reset = .FALSE.
    IF (MAX(nprocx,1) * MAX(nprocy,1) * MAX(nprocz,1) > nproc) THEN
      reset = .TRUE.
    ELSE IF (nprocx * nprocy * nprocz > 0) THEN
      ! Sanity check
      nxsplit = nx_global / nprocx
      nysplit = ny_global / nprocy
      nzsplit = nz_global / nprocz
      IF (nxsplit < ng .OR. nysplit < ng .OR. nzsplit < ng) &
          reset = .TRUE.
    ENDIF

    IF (reset) THEN
      IF (rank == 0) THEN
        PRINT *, 'Unable to use requested processor subdivision. Using ' &
            // 'default division.'
      ENDIF
      nprocx = 0
      nprocy = 0
      nprocz = 0
    ENDIF

    IF (nprocx * nprocy * nprocz == 0) THEN
      DO WHILE (nproc > 1)
        ! Find the processor split which minimizes surface area of
        ! the resulting domain

        minarea = nx_global * ny_global + ny_global * nz_global &
            + nz_global * nx_global

        DO ix = 1, nproc
          nprocyz = nproc / ix
          IF (ix * nprocyz /= nproc) CYCLE

          nxsplit = nx_global / ix
          ! Actual domain must be bigger than the number of ghostcells
          IF (nxsplit < ng) CYCLE

          DO iy = 1, nprocyz
            iz = nprocyz / iy
            IF (iy * iz /= nprocyz) CYCLE

            nysplit = ny_global / iy
            nzsplit = nz_global / iz
            ! Actual domain must be bigger than the number of ghostcells
            IF (nysplit < ng .OR. nzsplit < ng) CYCLE

            area = nxsplit * nysplit + nysplit * nzsplit + nzsplit * nxsplit
            IF (area < minarea) THEN
              nprocx = ix
              nprocy = iy
              nprocz = iz
              minarea = area
            ENDIF
          ENDDO
        ENDDO

        IF (nprocx > 0) EXIT

        ! If we get here then no suitable split could be found. Decrease the
        ! number of processors and try again.

        nproc = nproc - 1
      ENDDO
    ENDIF

    IF (nproc_orig /= nproc) THEN
      IF (.NOT.allow_cpu_reduce) THEN
        IF (rank == 0) THEN
          CALL integer_as_string(nproc, str)
          PRINT*,'*** ERROR ***'
          PRINT*,'Cannot split the domain using the requested number of CPUs.'
          PRINT*,'Try reducing the number of CPUs to ',TRIM(str)
        ENDIF
        CALL MPI_ABORT(MPI_COMM_WORLD, errcode, ierr)
        STOP
      ENDIF
      IF (rank == 0) THEN
        CALL integer_as_string(nproc, str)
        PRINT*,'*** WARNING ***'
        PRINT*,'Cannot split the domain using the requested number of CPUs.'
        PRINT*,'Reducing the number of CPUs to ',TRIM(str)
      ENDIF
      ranges(1,1) = nproc
      ranges(2,1) = nproc_orig - 1
      ranges(3,1) = 1
      old_comm = comm
      CALL MPI_COMM_GROUP(old_comm, oldgroup, errcode)
      CALL MPI_GROUP_RANGE_EXCL(oldgroup, 1, ranges, newgroup, errcode)
      CALL MPI_COMM_CREATE(old_comm, newgroup, comm, errcode)
      IF (comm == MPI_COMM_NULL) THEN
        CALL MPI_FINALIZE(errcode)
        STOP
      ENDIF
      CALL MPI_GROUP_FREE(oldgroup, errcode)
      CALL MPI_GROUP_FREE(newgroup, errcode)
      CALL MPI_COMM_FREE(old_comm, errcode)
    ENDIF

    dims = (/nprocz, nprocy, nprocx/)
    CALL MPI_DIMS_CREATE(nproc, ndims, dims, errcode)

    periods = .FALSE.
    reorder = .TRUE.

    ! Set boundary to be periodic if *any* boundary condition requires it.
    ! Once there are per-species boundary conditions then this will be true
    ! if any of the species are periodic

    IF (bc_field(c_bd_x_min) == c_bc_periodic &
        .OR. bc_particle(c_bd_x_min) == c_bc_periodic) &
            periods(c_ndims) = .TRUE.

    IF (bc_field(c_bd_y_min) == c_bc_periodic &
        .OR. bc_particle(c_bd_y_min) == c_bc_periodic) &
            periods(c_ndims-1) = .TRUE.

    IF (bc_field(c_bd_z_min) == c_bc_periodic &
        .OR. bc_particle(c_bd_z_min) == c_bc_periodic) &
            periods(c_ndims-2) = .TRUE.

    old_comm = comm
    CALL MPI_CART_CREATE(old_comm, ndims, dims, periods, reorder, comm, errcode)
    CALL MPI_COMM_FREE(old_comm, errcode)
    CALL MPI_COMM_RANK(comm, rank, errcode)
    CALL MPI_CART_COORDS(comm, rank, ndims, coordinates, errcode)
    CALL MPI_CART_SHIFT(comm, 2, 1, proc_x_min, proc_x_max, errcode)
    CALL MPI_CART_SHIFT(comm, 1, 1, proc_y_min, proc_y_max, errcode)
    CALL MPI_CART_SHIFT(comm, 0, 1, proc_z_min, proc_z_max, errcode)

    nprocx = dims(3)
    nprocy = dims(2)
    nprocz = dims(1)
    nprocdir = dims

    IF (rank == 0) THEN
      PRINT *, 'Processor subdivision is ', (/nprocx, nprocy, nprocz/)
    ENDIF

    x_coords = coordinates(c_ndims)
    x_min_boundary = .FALSE.
    x_max_boundary = .FALSE.
    IF (x_coords == 0) x_min_boundary = .TRUE.
    IF (x_coords == nprocx - 1) x_max_boundary = .TRUE.

    y_coords = coordinates(c_ndims-1)
    y_min_boundary = .FALSE.
    y_max_boundary = .FALSE.
    IF (y_coords == 0) y_min_boundary = .TRUE.
    IF (y_coords == nprocy - 1) y_max_boundary = .TRUE.

    z_coords = coordinates(c_ndims-2)
    z_min_boundary = .FALSE.
    z_max_boundary = .FALSE.
    IF (z_coords == 0) z_min_boundary = .TRUE.
    IF (z_coords == nprocz - 1) z_max_boundary = .TRUE.

    neighbour = MPI_PROC_NULL
    DO iz = -1, 1
      DO iy = -1, 1
        DO ix = -1, 1
          test_coords = coordinates
          test_coords(1) = test_coords(1)+iz
          test_coords(2) = test_coords(2)+iy
          test_coords(3) = test_coords(3)+ix
          op = .TRUE.
          ! For some stupid reason MPI_CART_RANK returns an error rather than
          ! MPI_PROC_NULL if the coords are out of range.
          DO idim = 1, ndims
            IF ((test_coords(idim) < 0 &
                .OR. test_coords(idim) >= dims(idim)) &
                .AND. .NOT. periods(idim)) op = .FALSE.
          ENDDO
          IF (op) THEN
            CALL MPI_CART_RANK(comm, test_coords, neighbour(ix,iy,iz), errcode)
          ENDIF
        ENDDO
      ENDDO
    ENDDO

  END SUBROUTINE setup_communicator



  SUBROUTINE mpi_initialise

    INTEGER :: ispecies, idim
    INTEGER :: nx0, nxp
    INTEGER :: ny0, nyp
    INTEGER :: nz0, nzp
    TYPE(particle_species), POINTER :: species, next_species

    CALL setup_communicator

    ALLOCATE(npart_each_rank(nproc))
    ALLOCATE(x_grid_mins(0:nprocx-1), x_grid_maxs(0:nprocx-1))
    ALLOCATE(y_grid_mins(0:nprocy-1), y_grid_maxs(0:nprocy-1))
    ALLOCATE(z_grid_mins(0:nprocz-1), z_grid_maxs(0:nprocz-1))
    ALLOCATE(cell_x_min(nprocx), cell_x_max(nprocx))
    ALLOCATE(cell_y_min(nprocy), cell_y_max(nprocy))
    ALLOCATE(cell_z_min(nprocz), cell_z_max(nprocz))

    nx0 = nx_global / nprocx
    ny0 = ny_global / nprocy
    nz0 = nz_global / nprocz

    ! If the number of gridpoints cannot be exactly subdivided then fix
    ! The first nxp processors have nx0 grid points
    ! The remaining processors have nx0+1 grid points
    IF (nx0 * nprocx /= nx_global) THEN
      nxp = (nx0 + 1) * nprocx - nx_global
    ELSE
      nxp = nprocx
    ENDIF

    IF (ny0 * nprocy /= ny_global) THEN
      nyp = (ny0 + 1) * nprocy - ny_global
    ELSE
      nyp = nprocy
    ENDIF

    IF (nz0 * nprocz /= nz_global) THEN
      nzp = (nz0 + 1) * nprocz - nz_global
    ELSE
      nzp = nprocz
    ENDIF

    DO idim = 1, nxp
      cell_x_min(idim) = (idim - 1) * nx0 + 1
      cell_x_max(idim) = idim * nx0
    ENDDO
    DO idim = nxp + 1, nprocx
      cell_x_min(idim) = nxp * nx0 + (idim - nxp - 1) * (nx0 + 1) + 1
      cell_x_max(idim) = nxp * nx0 + (idim - nxp) * (nx0 + 1)
    ENDDO

    DO idim = 1, nyp
      cell_y_min(idim) = (idim - 1) * ny0 + 1
      cell_y_max(idim) = idim * ny0
    ENDDO
    DO idim = nyp + 1, nprocy
      cell_y_min(idim) = nyp * ny0 + (idim - nyp - 1) * (ny0 + 1) + 1
      cell_y_max(idim) = nyp * ny0 + (idim - nyp) * (ny0 + 1)
    ENDDO

    DO idim = 1, nzp
      cell_z_min(idim) = (idim - 1) * nz0 + 1
      cell_z_max(idim) = idim * nz0
    ENDDO
    DO idim = nzp + 1, nprocz
      cell_z_min(idim) = nzp * nz0 + (idim - nzp - 1) * (nz0 + 1) + 1
      cell_z_max(idim) = nzp * nz0 + (idim - nzp) * (nz0 + 1)
    ENDDO

    nx_global_min = cell_x_min(x_coords+1)
    nx_global_max = cell_x_max(x_coords+1)
    n_global_min(1) = nx_global_min
    n_global_max(1) = nx_global_max

    ny_global_min = cell_y_min(y_coords+1)
    ny_global_max = cell_y_max(y_coords+1)
    n_global_min(2) = ny_global_min
    n_global_max(2) = ny_global_max

    nz_global_min = cell_z_min(z_coords+1)
    nz_global_max = cell_z_max(z_coords+1)
    n_global_min(3) = nz_global_min
    n_global_max(3) = nz_global_max

    nx = nx_global_max - nx_global_min + 1
    ny = ny_global_max - ny_global_min + 1
    nz = nz_global_max - nz_global_min + 1

    DEALLOCATE(x, y, z)
    ALLOCATE(x(1-ng:nx+ng), y(1-ng:ny+ng), z(1-ng:nz+ng))
    ALLOCATE(x_global(1-ng:nx_global+ng))
    ALLOCATE(y_global(1-ng:ny_global+ng))
    ALLOCATE(z_global(1-ng:nz_global+ng))
    ALLOCATE(xb_global(nx_global+1))
    ALLOCATE(yb_global(ny_global+1))
    ALLOCATE(zb_global(nz_global+1))
    ALLOCATE(ex(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(ey(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(ez(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(bx(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(by(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(bz(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ! Current may need an extra layer of ghostcells.
    ALLOCATE(jx(1-jng:nx+jng, 1-jng:ny+jng, 1-jng:nz+jng))
    ALLOCATE(jy(1-jng:nx+jng, 1-jng:ny+jng, 1-jng:nz+jng))
    ALLOCATE(jz(1-jng:nx+jng, 1-jng:ny+jng, 1-jng:nz+jng))

    ! For two-step RK.
    ALLOCATE(ex_back(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(ey_back(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(ez_back(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(bx_back(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(by_back(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(bz_back(1-ng:nx+ng, 1-ng:ny+ng, 1-ng:nz+ng))
    ALLOCATE(jx_d(1-jng:nx+jng, 1-jng:ny+jng, 1-jng:nz+jng))
    ALLOCATE(jy_d(1-jng:nx+jng, 1-jng:ny+jng, 1-jng:nz+jng))
    ALLOCATE(jz_d(1-jng:nx+jng, 1-jng:ny+jng, 1-jng:nz+jng))

    ! If needed, set-up 1D index
    IF (.NOT. explicit_pic) THEN
      ! Number of elements in local problem
      local_elements = nx * ny * nz * 6
      ALLOCATE(linear_index(local_elements))
      linear_index = 0
    END IF

    ! Setup the particle lists

    next_species => species_list
    DO ispecies = 1, n_species
      species => next_species
      next_species => species%next

      IF (bc_particle(c_bd_x_min) == c_bc_thermal) THEN
        ALLOCATE(species%ext_temp_x_min(1-ng:ny+ng, 1-ng:nz+ng, 3))
      ENDIF
      IF (bc_particle(c_bd_x_max) == c_bc_thermal) THEN
        ALLOCATE(species%ext_temp_x_max(1-ng:ny+ng, 1-ng:nz+ng, 3))
      ENDIF
      IF (bc_particle(c_bd_y_min) == c_bc_thermal) THEN
        ALLOCATE(species%ext_temp_y_min(1-ng:nx+ng, 1-ng:nz+ng, 3))
      ENDIF
      IF (bc_particle(c_bd_y_max) == c_bc_thermal) THEN
        ALLOCATE(species%ext_temp_y_max(1-ng:nx+ng, 1-ng:nz+ng, 3))
      ENDIF
      IF (bc_particle(c_bd_z_min) == c_bc_thermal) THEN
        ALLOCATE(species%ext_temp_z_min(1-ng:nx+ng, 1-ng:nz+ng, 3))
      ENDIF
      IF (bc_particle(c_bd_z_max) == c_bc_thermal) THEN
        ALLOCATE(species%ext_temp_z_max(1-ng:nx+ng, 1-ng:nz+ng, 3))
      ENDIF
    ENDDO

    CALL allocate_ic

    start_time = MPI_WTIME()

  END SUBROUTINE mpi_initialise



  SUBROUTINE mpi_close

    INTEGER :: seconds, minutes, hours, total

    IF (rank == 0) THEN
      end_time = MPI_WTIME()
      total = INT(end_time - start_time)
      seconds = MOD(total, 60)
      minutes = MOD(total / 60, 60)
      hours = total / 3600
    ENDIF

    CALL MPI_BARRIER(comm, errcode)

  END SUBROUTINE mpi_close



#ifdef MPI_DEBUG
  SUBROUTINE mpi_set_error_handler

    INTEGER :: errhandler

    CALL MPI_COMM_CREATE_ERRHANDLER(mpi_error_handler, errhandler, errcode)
    CALL MPI_COMM_SET_ERRHANDLER(MPI_COMM_WORLD, errhandler, errcode)

  END SUBROUTINE mpi_set_error_handler



  SUBROUTINE mpi_error_handler(comm, error_code)

    INTEGER :: comm, error_code
    REAL :: tmp1, tmp2
    CHARACTER(LEN=29) :: errstring(0:MPI_ERR_LASTCODE)

    errstring(MPI_SUCCESS                  ) = 'MPI_SUCCESS                  '
    errstring(MPI_ERR_BUFFER               ) = 'MPI_ERR_BUFFER               '
    errstring(MPI_ERR_COUNT                ) = 'MPI_ERR_COUNT                '
    errstring(MPI_ERR_TYPE                 ) = 'MPI_ERR_TYPE                 '
    errstring(MPI_ERR_TAG                  ) = 'MPI_ERR_TAG                  '
    errstring(MPI_ERR_COMM                 ) = 'MPI_ERR_COMM                 '
    errstring(MPI_ERR_RANK                 ) = 'MPI_ERR_RANK                 '
    errstring(MPI_ERR_REQUEST              ) = 'MPI_ERR_REQUEST              '
    errstring(MPI_ERR_ROOT                 ) = 'MPI_ERR_ROOT                 '
    errstring(MPI_ERR_GROUP                ) = 'MPI_ERR_GROUP                '
    errstring(MPI_ERR_OP                   ) = 'MPI_ERR_OP                   '
    errstring(MPI_ERR_TOPOLOGY             ) = 'MPI_ERR_TOPOLOGY             '
    errstring(MPI_ERR_DIMS                 ) = 'MPI_ERR_DIMS                 '
    errstring(MPI_ERR_ARG                  ) = 'MPI_ERR_ARG                  '
    errstring(MPI_ERR_UNKNOWN              ) = 'MPI_ERR_UNKNOWN              '
    errstring(MPI_ERR_TRUNCATE             ) = 'MPI_ERR_TRUNCATE             '
    errstring(MPI_ERR_OTHER                ) = 'MPI_ERR_OTHER                '
    errstring(MPI_ERR_INTERN               ) = 'MPI_ERR_INTERN               '
    errstring(MPI_ERR_IN_STATUS            ) = 'MPI_ERR_IN_STATUS            '
    errstring(MPI_ERR_PENDING              ) = 'MPI_ERR_PENDING              '
    errstring(MPI_ERR_ACCESS               ) = 'MPI_ERR_ACCESS               '
    errstring(MPI_ERR_AMODE                ) = 'MPI_ERR_AMODE                '
    errstring(MPI_ERR_ASSERT               ) = 'MPI_ERR_ASSERT               '
    errstring(MPI_ERR_BAD_FILE             ) = 'MPI_ERR_BAD_FILE             '
    errstring(MPI_ERR_BASE                 ) = 'MPI_ERR_BASE                 '
    errstring(MPI_ERR_CONVERSION           ) = 'MPI_ERR_CONVERSION           '
    errstring(MPI_ERR_DISP                 ) = 'MPI_ERR_DISP                 '
    errstring(MPI_ERR_DUP_DATAREP          ) = 'MPI_ERR_DUP_DATAREP          '
    errstring(MPI_ERR_FILE_EXISTS          ) = 'MPI_ERR_FILE_EXISTS          '
    errstring(MPI_ERR_FILE_IN_USE          ) = 'MPI_ERR_FILE_IN_USE          '
    errstring(MPI_ERR_FILE                 ) = 'MPI_ERR_FILE                 '
    errstring(MPI_ERR_INFO_KEY             ) = 'MPI_ERR_INFO_KEY             '
    errstring(MPI_ERR_INFO_NOKEY           ) = 'MPI_ERR_INFO_NOKEY           '
    errstring(MPI_ERR_INFO_VALUE           ) = 'MPI_ERR_INFO_VALUE           '
    errstring(MPI_ERR_INFO                 ) = 'MPI_ERR_INFO                 '
    errstring(MPI_ERR_IO                   ) = 'MPI_ERR_IO                   '
    errstring(MPI_ERR_KEYVAL               ) = 'MPI_ERR_KEYVAL               '
    errstring(MPI_ERR_LOCKTYPE             ) = 'MPI_ERR_LOCKTYPE             '
    errstring(MPI_ERR_NAME                 ) = 'MPI_ERR_NAME                 '
    errstring(MPI_ERR_NO_MEM               ) = 'MPI_ERR_NO_MEM               '
    errstring(MPI_ERR_NOT_SAME             ) = 'MPI_ERR_NOT_SAME             '
    errstring(MPI_ERR_NO_SPACE             ) = 'MPI_ERR_NO_SPACE             '
    errstring(MPI_ERR_NO_SUCH_FILE         ) = 'MPI_ERR_NO_SUCH_FILE         '
    errstring(MPI_ERR_PORT                 ) = 'MPI_ERR_PORT                 '
    errstring(MPI_ERR_QUOTA                ) = 'MPI_ERR_QUOTA                '
    errstring(MPI_ERR_READ_ONLY            ) = 'MPI_ERR_READ_ONLY            '
    errstring(MPI_ERR_RMA_CONFLICT         ) = 'MPI_ERR_RMA_CONFLICT         '
    errstring(MPI_ERR_RMA_SYNC             ) = 'MPI_ERR_RMA_SYNC             '
    errstring(MPI_ERR_SERVICE              ) = 'MPI_ERR_SERVICE              '
    errstring(MPI_ERR_SIZE                 ) = 'MPI_ERR_SIZE                 '
    errstring(MPI_ERR_SPAWN                ) = 'MPI_ERR_SPAWN                '
    errstring(MPI_ERR_UNSUPPORTED_DATAREP  ) = 'MPI_ERR_UNSUPPORTED_DATAREP  '
    errstring(MPI_ERR_UNSUPPORTED_OPERATION) = 'MPI_ERR_UNSUPPORTED_OPERATION'
    errstring(MPI_ERR_WIN                  ) = 'MPI_ERR_WIN                  '
    errstring(MPI_ERR_LASTCODE             ) = 'MPI_ERR_LASTCODE             '

    PRINT*, "Caught MPI error: ", TRIM(errstring(error_code))
    IF (comm == MPI_COMM_WORLD) THEN
      PRINT*, "Communicator MPI_COMM_WORLD"
    ELSE
      PRINT*, "Communicator ", comm, "(Not MPI_COMM_WORLD)"
    ENDIF

    ! Deliberately raise a divide-by-zero error
    tmp1 = 0.0
    tmp2 = 1.0 / tmp1

  END SUBROUTINE mpi_error_handler
#endif

END MODULE mpi_routines
