MODULE deck

  USE mpi
  USE shared_data
  USE timer

  IMPLICIT NONE

CONTAINS

  SUBROUTINE read_deck

    INTEGER :: ierr, i, errcode
    CHARACTER(LEN=1024), ALLOCATABLE, SAVE :: lines(:)
    LOGICAL, SAVE :: first = .TRUE.
    INTEGER, SAVE :: nlines = 0

    NAMELIST/control/ problem, x_min, x_max, y_min, y_max, z_min, z_max, &
         nx_global, ny_global, nz_global, nprocx, nprocy, nprocz, &
         allow_cpu_reduce, timer_collect, use_balance, use_random_seed, &
         npart_global, nsteps, t_end, dt_multiplier, dlb_threshold, &
         stdout_frequency, particle_push_start_time, n_species

    IF (first) THEN
      ! Set the default problem here
      problem = 'two_stream'

      ! Check for input deck existing.
      ! If exists, read into buffer, otherwise print warning
      IF (rank == 0) THEN
        OPEN(UNIT=lu, STATUS='OLD', FILE=TRIM(data_dir) // '/input.deck', &
             IOSTAT=ierr)

        IF (ierr /= 0) THEN
          PRINT *, 'WARNING: Failed to open file: ' // TRIM(data_dir) // '/input.deck'
        ELSE
          ! Count lines in deck
          DO
            READ(lu, FMT=*, IOSTAT=ierr)
            IF (ierr /= 0) EXIT
            nlines = nlines + 1
          END DO

          ! Allocate the buffer
          ALLOCATE(lines(nlines))

          ! Rewind file, read into buffer
          REWIND(lu)
          DO i = 1, nlines
            READ(lu, FMT='(A)', IOSTAT=ierr) lines(i)
            IF (ierr /= 0) THEN
              PRINT*, 'Error processing input file: ' // TRIM(data_dir) // 'input.deck'
              CALL MPI_ABORT(MPI_COMM_WORLD, c_err_io, errcode)
            END IF
          END DO

          CLOSE(lu)
        END IF
      END IF

      ! Broadcast buffer to other procs
      IF (nproc > 1) THEN

        CALL MPI_BCAST(nlines, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, errcode)

        IF (nlines > 0) THEN
          IF (rank > 0) THEN
            ALLOCATE(lines(nlines))
          END IF

          CALL MPI_BCAST(lines, nlines*LEN(lines(1)), MPI_CHARACTER, 0, &
              MPI_COMM_WORLD, errcode)
        END IF
      END IF

      ! No buffer to read
      IF (nlines == 0) THEN
        first = .FALSE.
        RETURN
      END IF

      ! Read deck
      READ(lines, NML=control, IOSTAT=ierr)
      IF (ierr /= 0) THEN
        IF (rank == 0) PRINT *, 'Error reading namelist'
        CALL MPI_ABORT(MPI_COMM_WORLD, c_err_io, errcode)
      END IF

      first = .FALSE.
    ELSE
      ! No buffer to read
      IF (nlines == 0) THEN
        RETURN
      END IF

      ! Read deck
      READ(lines, NML=control, IOSTAT=ierr)
      IF (ierr /= 0) THEN
        IF (rank == 0) PRINT *, 'Error reading namelist'
        CALL MPI_ABORT(MPI_COMM_WORLD, c_err_io, errcode)
      END IF

      DEALLOCATE(lines)
    END IF

  END SUBROUTINE read_deck

END MODULE deck
