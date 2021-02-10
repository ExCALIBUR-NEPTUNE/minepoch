MODULE deck

  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE read_deck

    INTEGER :: nlines, ierr
    INTEGER :: errcode
    LOGICAL, SAVE :: first = .TRUE.

    NAMELIST/control/ problem

    IF (first) THEN
      nlines = 0
      ! Check for input deck existing.
      ! If exists, read into buffer, otherwise abort
      IF (rank == 0) THEN
        OPEN(UNIT=lu, STATUS='OLD', FILE=TRIM(data_dir) // '/input.deck', &
             IOSTAT=ierr)

        IF (ierr /= 0) THEN
          PRINT *, 'Failed to open file: ' // TRIM(data_dir) // '/input.deck'
          CALL MPI_ABORT(MPI_COMM_WORLD, c_err_io, errcode)
        END IF
      END IF

      first = .FALSE.
    ELSE

    END IF

  END SUBROUTINE read_deck

END MODULE deck
